"""Research verification: locally graded thermal-shock benchmark on a perforated plate.

The separate semi-infinite analytical solution is a *local*, early-time
reference at the midpoint of the hot left edge, not an exact solution of the
finite 2D geometry. This module is not used by the production application.

The optional geometric radial grading moves *parent mesh vertices* before
constructing T6 midside nodes; the Q4 geometry is graded separately. No
common element assembly is introduced into the independent Q4 reference.
"""
from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np
from scipy.special import erfc, erfcinv
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

from .geometry import Mesh, make_mesh, upgrade_to_t6
from .q4_reference import QuadMesh, quad_mesh, assemble_quad
from .fem import assemble_thermal


def thermal_diffusivity(mat: Mapping) -> float:
    k, rho, cp = (float(mat[x]) for x in ('conductivity','density','heat_capacity'))
    if min(k, rho, cp) <= 0 or not np.all(np.isfinite([k,rho,cp])):
        raise ValueError('material thermal properties must be positive and finite')
    return k/(rho*cp)


def penetration_length(mat, t: float, fraction: float = .01) -> float:
    """Distance for a semi-infinite temperature step to decay to fraction."""
    if not (t>0 and 0<fraction<1):
        raise ValueError('t>0 and 0<fraction<1 required')
    return 2*np.sqrt(thermal_diffusivity(mat)*t)*erfcinv(fraction)


def halfspace_temperature(depth, t, mat, *, ambient=298., hot=310.):
    """Exact 1-D semi-infinite step solution at distance >=0 from hot face."""
    depth=np.asarray(depth,dtype=float)
    if np.any(depth<0) or not np.all(np.isfinite(depth)) or not t>0:
        raise ValueError('depth must be nonnegative, finite, and t>0')
    return ambient+(hot-ambient)*erfc(depth/(2*np.sqrt(thermal_diffusivity(mat)*t)))


def halfspace_surface_flux(t, mat, *, ambient=298., hot=310.):
    """Heat influx into the solid, W/m^2; positive for a hot surface."""
    return float(mat['conductivity'])*(hot-ambient)/np.sqrt(
        np.pi*thermal_diffusivity(mat)*t)


def radial_fractions(nradial: int, growth: float) -> np.ndarray:
    """Inner-to-outer ring fractions with geometric layers near outer edge.

    growth=1 reproduces equal radial layers. growth>1 places the smallest
    spacing next to the outer boundary, while preserving the geometry.
    """
    if isinstance(nradial,bool) or not isinstance(nradial,int) or nradial < 2:
        raise ValueError('nradial must be integer >=2')
    if not np.isfinite(growth) or growth < 1 or growth > 2:
        raise ValueError('growth must be finite in [1,2]')
    if growth==1:
        return np.linspace(0,1,nradial+1)
    exponent=np.log(growth)
    layers=np.expm1(exponent*np.arange(nradial,-1,-1))/np.expm1(exponent*nradial)
    return 1-layers


def graded_mesh(config: Mapping, family: str, growth: float):
    """Create a locally graded T3/T6 or independent Q4 conforming mesh."""
    geom=dict(config)
    geom.pop('element_type',None)
    if family in ('T3','T6'):
        m=make_mesh(**geom,element_type='T3')
    elif family=='Q4':
        m=quad_mesh(geom)
    else:
        raise ValueError('family must be T3, T6, Q4')
    nr=geom['radial']; count=4*geom['per_side']
    vals=radial_fractions(nr,growth)
    if np.min(np.diff(vals)) < 1e-12:
        raise ValueError("radial grading creates effectively collapsed layers")
    nodes=np.array(m.nodes,copy=True)
    inner=nodes[:count].copy()
    outer=nodes[nr*count:(nr+1)*count].copy()
    for j,t in enumerate(vals):
        nodes[j*count:(j+1)*count]=(1-t)*inner+t*outer
    if family in ('T3','T6'):
        tri=Mesh(nodes,m.elements.copy(),m.edges,m.node_sets)
        return tri if family=='T3' else upgrade_to_t6(tri)
    return QuadMesh(nodes,m.cells.copy(),m.edges,m.node_sets)


@dataclass
class ThermalShockResult:
    family: str
    mesh: object
    temperature: np.ndarray
    times: np.ndarray
    growth: float
    radial_step_m: float
    surface_gradient_K_per_m: float
    max_temperature_K: float
    min_temperature_K: float
    thermal_scheme: str = "consistent"
    min_temperature_all_steps_K: float = float('nan')
    max_temperature_all_steps_K: float = float('nan')
    boundary_energy_input_J: float = float('nan')
    stored_energy_change_J: float = float('nan')
    heat_balance_absolute_error_J: float = float('nan')



def quadratic_surface_gradient(t_surface: float, t_mid: float, t_inner: float, spacing: float) -> float:
    """Exact first derivative at the face for a quadratic 3-point edge trace."""
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError('spacing must be positive')
    return (-3*t_surface+4*t_mid-t_inner)/spacing


def solve_thermal_shock(config: Mapping, family: str, *, growth=1.,
                        final_time=1.6e-4, steps=128,
                        ambient=298., hot=310.,
                        thermal_scheme="consistent") -> ThermalShockResult:
    """Thermal-only transient benchmark, zero flux elsewhere, no FE coupling.

    The face temperature changes instantaneously at t=0+. Initial bulk
    temperature is ambient. Uses backward Euler with a consistent capacity
    matrix and a single sparse factorization. This does NOT solve coupled
    thermoelastic dynamics or provide mechanical stress validation.
    """
    if not isinstance(steps,int) or isinstance(steps,bool) or steps<1:
        raise ValueError('steps must be positive integer')
    if not np.isfinite(final_time) or final_time<=0:
        raise ValueError('final_time must be positive finite')
    g=config['geometry']; mat=config['material']
    if thermal_scheme not in ('consistent', 'monotone'):
        raise ValueError('thermal_scheme must be consistent or monotone')
    mesh=graded_mesh(g,family,growth)
    if thermal_scheme == 'monotone':
        from .monotone_thermal import assemble_monotone_thermal
        model=assemble_monotone_thermal(mesh,mat,family)
        K,C=model.conductivity,model.capacity
    elif family=='Q4':
        K,C,*_=assemble_quad(mesh,mat)
    else:
        K,C=assemble_thermal(mesh,mat)
    n=len(mesh.nodes); fixed={}
    for label,value in [('left',hot),('right',ambient)]:
        for i in mesh.node_sets[label]:
            i=int(i)
            if i in fixed and abs(fixed[i]-value)>1e-11:
                raise ValueError('inconsistent thermal boundary at shared node')
            fixed[i]=value-ambient
    fidx=np.array(sorted(fixed),dtype=int)
    fval=np.array([fixed[int(i)] for i in fidx],dtype=float)
    free=np.setdiff1d(np.arange(n),fidx)
    dt=final_time/steps
    A=(C/dt+K).tocsc()
    linear=splu(A[free,:][:,free]) if len(free) else None
    coupling=A[free,:][:,fidx]
    theta=np.zeros(n)
    theta[fidx]=fval
    first=theta.copy()
    global_min=float(min(theta+ambient))
    global_max=float(max(theta+ambient))
    total_boundary_energy=0.
    for _ in range(steps):
        previous=theta.copy()
        right=(C@previous)/dt
        if len(free):
            theta[free]=linear.solve(right[free]-coupling@fval)
        residual=C@((theta-previous)/dt)+K@theta
        total_boundary_energy+=dt*float(np.sum(residual[fidx]))
        global_min=min(global_min,float(np.min(theta+ambient)))
        global_max=max(global_max,float(np.max(theta+ambient)))
    stored_energy=float(np.sum(C@(theta-first)))
    # Flux at left mid-edge from derivative of the interpolation on the
    # first radial element; T6 has a midpoint for its quadratic profile.
    depths,ids=radial_path(mesh,family,g)
    vals=theta[ids]+ambient
    spacing=float(depths[1]-depths[0])
    if family=='T6':
        # First midside value is on the radial edge; the monotone scheme
        # uses a P1 subsegment, while the consistent scheme uses P2.
        edge_midpoint=profile_midpoint_temperature(mesh,g,theta+ambient)
        if thermal_scheme=='monotone':
            gradient=2*(edge_midpoint-vals[0])/spacing
        else:
            gradient=quadratic_surface_gradient(vals[0],edge_midpoint,vals[1],spacing)
    else:
        gradient=(vals[1]-vals[0])/spacing
    return ThermalShockResult(family,mesh,theta+ambient,
                              np.linspace(0,final_time,steps+1),growth,spacing,
                              float(gradient),float(max(theta+ambient)),
                              float(min(theta+ambient)),thermal_scheme,
                              global_min,global_max,total_boundary_energy,
                              stored_energy,abs(total_boundary_energy-stored_energy))


def radial_path(mesh, family, geom):
    """Hot-left midpoint profile; T6 base vertices retain parent indexing."""
    count=4*geom['per_side']; nradial=geom['radial']; k=3*geom['per_side']+geom['per_side']//2
    ids=np.arange(nradial+1)*count+k
    surface=float(mesh.nodes[ids[-1],0])
    depths=mesh.nodes[ids[::-1],0]-surface
    # temperatures are supplied separately for this profile
    return depths,ids[::-1]


def profile_midpoint_temperature(mesh,geom,temperature):
    """Retrieve temperature at the first radial T6 midside node."""
    _, ids=radial_path(mesh,'T6',geom)
    midpoint=.5*(mesh.nodes[ids[0]]+mesh.nodes[ids[1]])
    distance,nid=cKDTree(mesh.nodes).query(midpoint)
    if distance>1e-10:
        raise ValueError('T6 midside node absent')
    return float(temperature[nid])


def extract_profile(result: ThermalShockResult, geom, distances):
    """Evaluate continuous T3/Q4 linear or T6 quadratic radial FE trace."""
    distances=np.asarray(distances,dtype=float)
    if np.any(distances<0) or not np.all(np.isfinite(distances)):
        raise ValueError('distances must be nonnegative and finite')
    depths,ids=radial_path(result.mesh,result.family,geom)
    field=result.temperature
    value=field[ids]
    if np.any(distances>depths[-1]):
        raise ValueError('distance exceeds radial path')
    j=np.searchsorted(depths,distances,side='right')-1
    j=np.clip(j,0,len(depths)-2)
    s=(distances-depths[j])/(depths[j+1]-depths[j])
    if result.family!='T6':
        return (1-s)*value[j]+s*value[j+1]
    # global midside indices for each radial interval
    mids=.5*(result.mesh.nodes[ids[:-1]]+result.mesh.nodes[ids[1:]])
    err,mid_ids=cKDTree(result.mesh.nodes).query(mids)
    if np.any(err>1e-10):
        raise ValueError('nonconforming quadratic profile')
    midval=field[mid_ids]
    if result.thermal_scheme=='monotone':
        # P1 on two thermal subsegments. Quadratic interpolation would
        # violate the shape-space and maximum-principle interpretation.
        return np.where(s<=.5,
            (1-2*s)*value[j]+2*s*midval[j],
            (2-2*s)*midval[j]+(2*s-1)*value[j+1])
    return (1-s)*(1-2*s)*value[j]+4*s*(1-s)*midval[j]+s*(2*s-1)*value[j+1]


def shock_error(result, config, *, final_time=1.6e-4, depth_max=.00025):
    """RMS and maximum errors over fixed physical depths near mid left edge."""
    distances=np.linspace(0,depth_max,101)
    ref=halfspace_temperature(distances,final_time,config['material'])
    approx=extract_profile(result,config['geometry'],distances)
    delta=approx-ref
    flux_exact=halfspace_surface_flux(final_time,config['material'])
    flux_fem=-float(config['material']['conductivity'])*result.surface_gradient_K_per_m
    return {'family':result.family,'growth':result.growth,'subdivisions':config['geometry']['per_side'],
            'radial':config['geometry']['radial'],'steps':len(result.times)-1,
            'elements':len(result.mesh.cells if result.family=='Q4' else result.mesh.elements),
            'nodes':len(result.mesh.nodes),'first_radial_spacing_m':result.radial_step_m,
            'rms_temperature_error_K':float(np.sqrt(np.mean(delta*delta))),
            'max_temperature_error_K':float(np.max(abs(delta))),
            'surface_flux_FEM_W_per_m2':flux_fem,
            'surface_flux_exact_W_per_m2':flux_exact,
            'surface_flux_relative_error':abs(flux_fem-flux_exact)/abs(flux_exact),
            'undershoot_K':max(0.,298.-result.min_temperature_K),
            'overshoot_K':max(0.,result.max_temperature_K-310.),
            'thermal_scheme':result.thermal_scheme,
            'min_temperature_all_steps_K':result.min_temperature_all_steps_K,
            'max_temperature_all_steps_K':result.max_temperature_all_steps_K,
            'heat_balance_absolute_error_J':result.heat_balance_absolute_error_J}


def assemble_graded_wave_system(config, family, growth, *, thermal_scheme="consistent"):
    """Experimental graded version of the pre-existing coupled-wave harness.

    Returns the existing WaveSystem with matrices built on a *new* graded
    mesh. This is a diagnostic, not a validated production mesh backend.
    """
    from .fem import assemble_elastic, prescribed_displacements
    from .coupled import assemble_inertia_coupling
    from .q4_reference import _edge_load, _fixed
    from .dynamic_diagnostics import WaveSystem

    mesh=graded_mesh(config['geometry'],family,growth)
    mat=config['material']
    if family=='Q4':
        Kt,Ct,Ku,M,G=assemble_quad(mesh,mat)
        traction=_edge_load(mesh,{'top':[0.,2e6]},mat['thickness'],2)
        fixed={}
        for name,component in [('ux',0),('uy',1)]:
            fixed.update(_fixed(mesh,config['analysis']['displacement_bc'][name],component))
    else:
        Kt,Ct=assemble_thermal(mesh,mat)
        Ku,traction=assemble_elastic(
            mesh,mat,np.full(len(mesh.nodes),mat['T_ref']),
            {'top':[0.,2e6]})
        M,G=assemble_inertia_coupling(mesh,mat)
        fixed=prescribed_displacements(mesh,config['analysis']['displacement_bc'])
    if thermal_scheme == 'monotone':
        from .monotone_thermal import assemble_monotone_thermal
        model=assemble_monotone_thermal(mesh,mat,family)
        Kt,Ct=model.conductivity,model.capacity
        if family=='T6':
            # Use the same piecewise-linear thermal interpolation in the
            # mechanical-to-thermal coupling block; preserve H=T0*G.T.
            from .monotone_thermal import assemble_subtriangle_thermal_coupling
            G=assemble_subtriangle_thermal_coupling(mesh,mat)
    elif thermal_scheme != 'consistent':
        raise ValueError('thermal_scheme must be consistent or monotone')
    return WaveSystem(mesh,Kt,Ct,Ku,M,G,traction,fixed,
                      float(mat['T_ref']),family)
