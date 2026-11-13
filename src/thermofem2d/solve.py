"""Heat conduction followed by one-way, quasi-static thermoelastic response."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse.linalg import factorized
from .geometry import make_mesh,Mesh
from .fem import (assemble_thermal,assemble_elastic,boundary_load,
                  prescribed_temperature,prescribed_displacements,
                  solve_dirichlet,stresses)
from .element import TRIANGLE_QUADRATURE,shape_and_gradients,element_type_for_mesh


@dataclass
class Solution:
    mesh: Mesh
    temperature: np.ndarray
    displacement: np.ndarray
    sigma: np.ndarray
    von_mises: np.ndarray
    thermal_reactions: np.ndarray
    mechanical_reactions: np.ndarray
    times: np.ndarray
    temperature_history: np.ndarray
    material: dict
    displacement_history: np.ndarray | None = None
    velocity_history: np.ndarray | None = None
    mechanical_acceleration: np.ndarray | None = None
    two_way_feedback: bool = False


def _heat_load(mesh,material,analysis):
    loads=boundary_load(mesh,analysis.get("inward_heat_flux",{}),material["thickness"],1)
    Q=float(analysis.get("volumetric_heat_source",0.0))
    if not np.isfinite(Q):
        raise ValueError("nonfinite volumetric heat source")
    if Q:
        kind=element_type_for_mesh(mesh)
        for tri in mesh.elements:
            p=mesh.nodes[tri[:3]]
            fe=np.zeros(len(tri))
            for r,s,w in TRIANGLE_QUADRATURE:
                N,_,area=shape_and_gradients(p,kind,r,s)
                fe+=material["thickness"]*Q*area*w*N
            loads[tri]+=fe
    return loads


def initial_temperature_field(mesh: Mesh, material: dict, analysis: dict) -> np.ndarray:
    """Uniform or affine initial temperature in kelvin (spatial gradients K/m).

    Configuration format: ``initial_temperature: 310.0`` (uniform), or
    ``initial_temperature: {type: affine, base: 310, gx: 5, gy: -2}``.
    No arbitrary Python expression is evaluated. Initial Dirichlet values
    are imposed by each solver after evaluating this domain-wide field.
    """
    spec = analysis.get("initial_temperature", material["T_ref"])
    if isinstance(spec, dict):
        if set(spec) != {"type", "base", "gx", "gy"} or spec["type"] != "affine":
            raise ValueError("initial temperature mapping requires type='affine', base, gx, gy")
        try:
            base, gx, gy = (float(spec[k]) for k in ("base", "gx", "gy"))
        except (ValueError, TypeError) as err:
            raise ValueError("initial affine temperature values must be numbers") from err
        if not np.all(np.isfinite([base, gx, gy])):
            raise ValueError("initial affine temperatures must be finite")
        values = base + gx*mesh.nodes[:, 0] + gy*mesh.nodes[:, 1]
    else:
        try:
            value = float(spec)
        except (ValueError, TypeError) as err:
            raise ValueError("initial temperature must be a finite scalar") from err
        values = np.full(len(mesh.nodes), value)
    if np.any(~np.isfinite(values)):
        raise ValueError("initial temperature field contains nonfinite values")
    return np.asarray(values, dtype=float)


def solve_model(config:dict)->Solution:
    """Execute a model with steady or backward-Euler transient heat solution."""
    thermoelastic_mode=config["analysis"].get("thermoelastic_mode", "one_way_quasistatic")
    if thermoelastic_mode == "coupled_dynamic":
        from .coupled import solve_coupled_dynamic
        return solve_coupled_dynamic(config)
    if thermoelastic_mode != "one_way_quasistatic":
        raise ValueError("thermoelastic_mode must be one_way_quasistatic or coupled_dynamic")
    mesh=make_mesh(**config["geometry"])
    mat=config["material"];analysis=config["analysis"]
    K,C=assemble_thermal(mesh,mat)
    fth=_heat_load(mesh,mat,analysis)
    fixedT=prescribed_temperature(mesh,analysis["temperature_bc"])
    if not fixedT:
        raise ValueError("temperature problem needs at least one prescribed temperature")
    mode=analysis.get("thermal_mode","steady")
    if mode=="steady":
        temp=solve_dirichlet(K,fth,fixedT)
        history=np.array([temp.copy()]);times=np.array([0.0])
    elif mode=="transient":
        dt=float(analysis["time_step"]);n=int(analysis["steps"])
        if not np.isfinite(dt) or dt<=0 or n<1 or n>10000 or n!=analysis["steps"]:
            raise ValueError("time step must be positive and step count must be 1..10000")
        temp=initial_temperature_field(mesh,mat,analysis)
        if not np.all(np.isfinite(temp)):
            raise ValueError("invalid initial temperature")
        # Dirichlet conditions are applied instantaneously at t=0+; the
        # initial condition is interpreted in the domain interior.
        ids=np.array(sorted(fixedT),dtype=int)
        temp[ids]=[fixedT[k] for k in ids]
        times=np.arange(n+1)*dt
        history=np.empty((n+1,len(temp)));history[0]=temp
        A=C+dt*K
        free=np.setdiff1d(np.arange(len(temp)),ids)
        if free.size:
            factor=factorized(A[free,:][:,free].tocsc())
            Afc=A[free,:][:,ids]
        for i in range(1,n+1):
            rhs=C@temp+dt*fth
            if free.size:
                temp[free]=factor(rhs[free]-Afc@temp[ids])
            history[i]=temp
    else:
        raise ValueError("thermal_mode must be 'steady' or 'transient'")
    stiffness,load=assemble_elastic(mesh,mat,temp,analysis.get("traction_bc",{}))
    fixedU=prescribed_displacements(mesh,analysis["displacement_bc"])
    if not fixedU:
        raise ValueError("mechanical problem needs enough displacement constraints")
    u=solve_dirichlet(stiffness,load,fixedU)
    sigma,vm=stresses(mesh,mat,temp,u)
    if mode == "transient":
        thermal_residual = C@((history[-1]-history[-2])/dt)+K@temp-fth
    else:
        thermal_residual = K@temp-fth
    return Solution(mesh,temp,u,sigma,vm,thermal_residual,
                    stiffness@u-load,times,history,dict(mat))
