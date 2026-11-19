"""Independent analytical Kirsch benchmark and explicit circular-hole stress recovery.

This research-only module does not alter application geometry or dispatch.
The exact stresses are for an infinite isotropic plane under remote x tension;
finite rectangular domains require the *corresponding* exact exterior data to
constitute an analytical boundary-value verification problem.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .element import shape_and_gradients, strain_displacement
from .fem import assemble_elastic, plane_stress_matrix, solve_dirichlet


@dataclass(frozen=True)
class HoleStressSamples:
    """Raw, one-sided, element-boundary midpoint stresses (no smoothing)."""
    theta_rad: np.ndarray
    radius_m: np.ndarray
    sigma_rr_pa: np.ndarray
    sigma_tt_pa: np.ndarray
    sigma_rt_pa: np.ndarray
    xy_m: np.ndarray


def _coords(xy, radius, centre):
    points = np.asarray(xy, dtype=float)
    c = np.asarray(centre, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 2 or not np.all(np.isfinite(points)):
        raise ValueError('xy must be a finite array of points with final dimension 2')
    if c.shape != (2,) or not np.all(np.isfinite(c)):
        raise ValueError('centre must contain two finite coordinates')
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError('hole radius must be positive and finite')
    relative = points - c
    r = np.linalg.norm(relative, axis=-1)
    if np.any(r < radius * (1 - 1e-12)):
        raise ValueError('Kirsch solution is defined only outside the circular hole')
    return relative, r, np.arctan2(relative[..., 1], relative[..., 0])


def kirsch_polar_stress(xy, radius, sigma_inf, centre=(0.0, 0.0)):
    """Return [sigma_rr,sigma_theta_theta,tau_r_theta] (Pa) for x-directed tension.

    theta is measured CCW from +x. The hole surface is traction-free.
    """
    if not np.isfinite(sigma_inf):
        raise ValueError('remote stress must be finite')
    _, r, th = _coords(xy, radius, centre)
    q = (radius / r)**2
    c = np.cos(2*th)
    ss = np.sin(2*th)
    rr = .5*sigma_inf*((1-q)+(1-4*q+3*q*q)*c)
    tt = .5*sigma_inf*((1+q)-(1+3*q*q)*c)
    rt = -.5*sigma_inf*(1+2*q-3*q*q)*ss
    return np.stack((rr,tt,rt),axis=-1)


def kirsch_cartesian_stress(xy, radius, sigma_inf, centre=(0.0, 0.0)):
    """Return [sigma_xx,sigma_yy,tau_xy] for the infinite plate (Pa)."""
    _, _, th = _coords(xy, radius, centre)
    rr, tt, rt = np.moveaxis(kirsch_polar_stress(xy,radius,sigma_inf,centre),-1,0)
    c = np.cos(th); s = np.sin(th)
    return np.stack((rr*c*c+tt*s*s-2*rt*c*s,
                     rr*s*s+tt*c*c+2*rt*c*s,
                     (rr-tt)*c*s+rt*(c*c-s*s)),axis=-1)


def kirsch_displacement(xy, radius, sigma_inf, E, nu, centre=(0.0, 0.0)):
    """Plane-stress displacement [ux,uy] (m), zero rigid-body translation/rotation.

    Polar displacements are integrated from the Kirsch stresses under isotropic
    plane stress; the far-field solution is ux=sigma_inf*x/E,
    uy=-nu*sigma_inf*y/E, up to a rigid displacement.
    """
    if not np.isfinite(E) or E<=0 or not np.isfinite(nu) or not (-1<nu<.5):
        raise ValueError('expected a positive Young modulus and valid Poisson ratio')
    _,r,th = _coords(xy,radius,centre)
    a2=radius**2
    a4=a2**2
    factor=sigma_inf/(2*E)
    ur=factor*((1-nu)*r+(1+nu)*a2/r +
               ((1+nu)*r+4*a2/r-(1+nu)*a4/r**3)*np.cos(2*th))
    ut=-factor*((1+nu)*r+2*(1-nu)*a2/r+(1+nu)*a4/r**3)*np.sin(2*th)
    return np.stack((ur*np.cos(th)-ut*np.sin(th),
                     ur*np.sin(th)+ut*np.cos(th)),axis=-1)


def hole_edge_stress_samples(mesh, mat, displacement, centre=(0.,0.)):
    """Evaluate raw physical stresses at the midpoint of each hole edge.

    For T6, the physical evaluation point is the actual curved edge's midside
    node. For T3, it is the polygonal chord midpoint. No nodal smoothing,
    averaging, extrapolation, or inferred global maximum is performed.
    """
    u=np.asarray(displacement, dtype=float)
    c=np.asarray(centre,dtype=float)
    if c.shape!=(2,) or not np.all(np.isfinite(c)):
        raise ValueError('invalid hole centre')
    if u.shape!=(2*len(mesh.nodes),) or not np.all(np.isfinite(u)):
        raise ValueError('displacement must contain finite nodal ux,uy values')
    if 'hole' not in mesh.edges:
        raise ValueError('hole boundary missing')
    kind = 'T3' if mesh.elements.shape[1]==3 else 'T6'
    edge_ref=(((0,1),(.5,0.)),((1,2),(.5,.5)),((2,0),(0.,.5)))
    edge_lookup={}
    for i,el in enumerate(mesh.elements):
        for j,((la,lb),rs) in enumerate(edge_ref):
            key=tuple(sorted((int(el[la]),int(el[lb]))))
            if key in edge_lookup:
                edge_lookup[key]=None  # interior edge: not one-sided
            else:
                edge_lookup[key]=(i,j,rs, int(el[3+j]) if kind=='T6' else None)
    D=plane_stress_matrix(mat['E'],mat['nu'])
    samples=[]
    for edge in mesh.edges['hole']:
        key=tuple(sorted((int(edge[0]),int(edge[-1]))))
        record=edge_lookup.get(key)
        if record is None:
            raise ValueError('hole edge must be adjacent to exactly one triangle')
        idx,j,(r,s),mid=record
        if kind=='T6' and mid!=int(edge[1]):
            raise ValueError('hole midside index inconsistent with local element')
        tri=mesh.elements[idx]
        p=mesh.nodes[tri]
        if kind=='T6':
            x=mesh.nodes[edge[1]]
        else:
            x=.5*(mesh.nodes[edge[0]]+mesh.nodes[edge[1]])
        _,grads,_=shape_and_gradients(p,kind,r,s)
        dofs=np.stack((2*tri,2*tri+1),axis=-1).reshape(-1)
        sx,sy,txy=D@(strain_displacement(grads)@u[dofs])
        phi=np.arctan2(x[1]-c[1],x[0]-c[0]);ca=np.cos(phi);sa=np.sin(phi)
        rr=sx*ca*ca+sy*sa*sa+2*txy*ca*sa
        tt=sx*sa*sa+sy*ca*ca-2*txy*ca*sa
        rt=(sy-sx)*ca*sa+txy*(ca*ca-sa*sa)
        samples.append((float(phi%(2*np.pi)),float(np.linalg.norm(x-c)),rr,tt,rt,x))
    samples.sort(key=lambda v:v[0])
    return HoleStressSamples(theta_rad=np.array([v[0] for v in samples]),
                             radius_m=np.array([v[1] for v in samples]),
                             sigma_rr_pa=np.array([v[2] for v in samples]),
                             sigma_tt_pa=np.array([v[3] for v in samples]),
                             sigma_rt_pa=np.array([v[4] for v in samples]),
                             xy_m=np.array([v[5] for v in samples]))


def canonical_first_quadrant_angle(theta_rad):
    """Return a symmetry-equivalent hole angle in [0, 90] degrees.

    For a centred circular hole in a square under symmetric x-directed
    uniaxial loading, reflections about both coordinate axes leave the
    circumferential stress invariant. A sampled peak may therefore be at
    theta, -theta, pi-theta, or pi+theta, with numerically indistinguishable
    values. The canonical angle is a location descriptor only: it does not
    modify the all-edge sampled maximum stress or concentration factor.
    Rounding to 9 decimal places suppresses roundoff-scale angle
    differences between otherwise equivalent sampled locations.
    """
    if not np.isfinite(theta_rad):
        raise ValueError('peak-sample angle must be finite')
    angle_deg=float(np.rad2deg(theta_rad)) % 180.0
    return round(min(angle_deg, 180.0-angle_deg), 9)


def solve_kirsch_outer_dirichlet(mesh, mat, radius, sigma_inf, centre=(0.,0.)):
    """Exact Kirsch displacement at exterior boundary, natural free-hole traction.

    A finite computational subdomain of the infinite-plane exact solution;
    unlike constant remote-traction plate loading, this is an analytical
    boundary-value problem. Curved edge approximates, not exactly matches, the
    circle: a numerical geometry-discretization error remains.
    """
    T=np.full(len(mesh.nodes),mat['T_ref'],dtype=float)
    K,f=assemble_elastic(mesh,mat,T)
    exterior=np.unique(np.concatenate([mesh.node_sets[k] for k in
                                      ('left','right','top','bottom')]))
    exact=kirsch_displacement(mesh.nodes[exterior],radius,sigma_inf,
                              mat['E'],mat['nu'],centre)
    fixed={2*int(i)+d:float(exact[j,d]) for j,i in enumerate(exterior) for d in (0,1)}
    return solve_dirichlet(K,f,fixed)


def solve_finite_plate_uniform_tension(mesh, mat, sigma_inf):
    """Finite *centred* plate: equal and opposite uniform outer x tractions.

    Reference stress is the applied gross-section traction, not net stress.
    Exactly three corner DOFs remove rigid translations and rotation.
    The named rectangle is centred at (0,0); no eccentric hole is inferred.
    """
    if not np.isfinite(sigma_inf) or sigma_inf<=0:
        raise ValueError('remote tensile stress must be positive and finite')
    p=mesh.nodes
    envelope=.5*(np.max(p,axis=0)+np.min(p,axis=0))
    hole_centre=np.mean(p[mesh.node_sets['hole']],axis=0)
    tolerance=1e-10*max(1.,float(np.max(np.ptp(p,axis=0))))
    if np.linalg.norm(hole_centre-envelope)>tolerance:
        raise ValueError('uniform-traction benchmark requires a centred hole')
    # Use the symmetric lower corners. For a centred hole under symmetric
    # tension their relative y displacement is zero, so no deformation is
    # suppressed by the minimal rigid-body constraints.
    bl=np.intersect1d(mesh.node_sets['bottom'],mesh.node_sets['left'])
    br=np.intersect1d(mesh.node_sets['bottom'],mesh.node_sets['right'])
    if len(bl)!=1 or len(br)!=1:
        raise ValueError('expected one node at each lower plate corner')
    T=np.full(len(p),mat['T_ref'],dtype=float)
    K,f=assemble_elastic(mesh,mat,T,
                         {'left':[-sigma_inf,0.], 'right':[sigma_inf,0.]})
    fixed={2*int(bl[0]):0., 2*int(bl[0])+1:0., 2*int(br[0])+1:0.}
    return solve_dirichlet(K,f,fixed)
