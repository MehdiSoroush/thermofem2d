"""Independent-in-time reference and energy diagnostics for the coupled model.

The ODE reference integrates the *same finite-element semi-discretization* as
``coupled.solve_coupled_dynamic``, but does not use its Newmark/BE time stepper.
It is NOT an independent spatial FEM or external Abaqus validation.
"""
from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import splu

from .geometry import make_mesh
from .fem import (assemble_thermal, assemble_elastic, prescribed_displacements,
                  prescribed_temperature)
from .coupled import assemble_inertia_coupling
from .solve import _heat_load, initial_temperature_field


@dataclass
class ODEReference:
    times: np.ndarray
    temperature_history: np.ndarray
    displacement_history: np.ndarray
    velocity_history: np.ndarray
    max_relative_residual: float


def solve_semidiscrete_ode(config: Mapping, times, *, rtol=2e-10, atol=1e-11):
    """High-accuracy method-of-lines reference using adaptive DOP853.

    Only supports the currently implemented constant Dirichlet conditions,
    constant loads, zero initial displacement/velocity and constant initial
    temperature. The 1st-order ODE is assembled *independently of the
    production Newmark/BE step formulation*.

    A strong temporal check, but it shares FEM element matrices and geometry.
    """
    a=config['analysis']; m=config['material']
    if not np.isfinite(m['T_ref']) or m['T_ref'] <= 0:
        raise ValueError('coupled energy/reference requires T_ref > 0 K')
    if a.get('thermoelastic_mode') != 'coupled_dynamic' or a.get('thermal_mode') != 'transient':
        raise ValueError('ODE reference requires coupled_dynamic transient mode')
    if a.get('initial_displacement',0) != 0 or a.get('initial_velocity',0) != 0:
        raise ValueError('only zero initial kinematic fields are supported')
    t=np.asarray(times,dtype=float)
    if t.ndim!=1 or len(t)<2 or t[0]!=0 or np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError('times must start at 0 and increase strictly')
    mesh=make_mesh(**config['geometry'])
    Kt,Ct=assemble_thermal(mesh,m)
    Ku,f=assemble_elastic(mesh,m,np.full(len(mesh.nodes),m['T_ref']),a.get('traction_bc',{}))
    M,G=assemble_inertia_coupling(mesh,m)
    feedback=a.get('thermoelastic_feedback',True)
    if not isinstance(feedback,bool):
        raise ValueError('thermoelastic_feedback must be bool')
    H=(m['T_ref']*G.T) if feedback else (0*G.T)
    q=_heat_load(mesh,m,a)
    fu=prescribed_displacements(mesh,a.get('displacement_bc',{}))
    ft=prescribed_temperature(mesh,a.get('temperature_bc',{}))
    if not fu:
        raise ValueError('mechanical constraints required')
    n=len(mesh.nodes); nu=2*n
    umask=np.ones(nu,dtype=bool)
    umask[list(fu)]=False
    tmask=np.ones(n,dtype=bool)
    tmask[list(ft)]=False
    uf=np.flatnonzero(umask); tf=np.flatnonzero(tmask)
    ufix=np.zeros(nu)
    for i,val in fu.items(): ufix[i]=val
    thetafix=np.zeros(n)
    for i,val in ft.items(): thetafix[i]=val-m['T_ref']
    thetainit=initial_temperature_field(mesh,m,a)-m['T_ref']
    thetainit[~tmask]=thetafix[~tmask]
    uinit=ufix.copy()
    vinit=np.zeros(nu)
    init=np.r_[uinit[uf],vinit[uf],thetainit[tf]]
    mf=lU=sU=None
    if len(uf):
        mf=splu(M[uf,:][:,uf].tocsc())
    if len(tf):
        sT=splu(Ct[tf,:][:,tf].tocsc())

    def full(z):
        u=ufix.copy(); v=np.zeros(nu); th=thetafix.copy()
        a0=len(uf);a1=2*a0
        u[uf]=z[:a0]
        v[uf]=z[a0:a1]
        th[tf]=z[a1:]
        return u,v,th

    def deriv(_time,z):
        u,v,th=full(z)
        d=np.zeros_like(z)
        k=len(uf)
        if k:
            d[:k]=v[uf]
            d[k:2*k]=mf.solve((f+G@th-Ku@u)[uf])
        if len(tf):
            d[2*k:]=sT.solve((q-Kt@th-H@v)[tf])
        return d

    solution=solve_ivp(deriv,(t[0],t[-1]),init,method='DOP853',t_eval=t,
                       rtol=rtol,atol=atol)
    if not solution.success or len(solution.t)!=len(t):
        raise RuntimeError(f'ODE reference failed: {solution.message}')
    uh=[];vh=[];th=[]
    for z in solution.y.T:
        u,v,theta=full(z);uh.append(u);vh.append(v);th.append(theta+m['T_ref'])
    uh=np.array(uh);vh=np.array(vh);th=np.array(th)
    # Check continuous-time semi-discrete residual at the end, independent
    # of the production integrator; denominator is a scale, not a tolerance.
    u,v,theta=full(solution.y[:,-1])
    rhs=deriv(t[-1],solution.y[:,-1])
    k=len(uf)
    a_full=np.zeros(nu);a_full[uf]=rhs[k:2*k]
    thdot=np.zeros(n);thdot[tf]=rhs[2*k:]
    mr=(M@a_full+Ku@u-G@theta-f)[uf]
    tr=(Ct@thdot+Kt@theta+H@v-q)[tf]
    scale=max(1.0,float(np.linalg.norm(f+G@theta)),float(np.linalg.norm(q+H@v)))
    residual=float(np.linalg.norm(np.r_[mr,tr])/scale)
    return ODEReference(t,th,uh,vh,residual)


@dataclass
class EnergyAudit:
    times: np.ndarray
    mechanical_energy: np.ndarray
    thermal_energy: np.ndarray
    conduction_power: np.ndarray
    cumulative_conduction: np.ndarray
    energy_balance_error: np.ndarray


def insulated_energy_audit(config: Mapping, result) -> EnergyAudit:
    """Compute the zero-source, insulated, fixed-support discrete energy budget.

    For this domain, continuous energy satisfies E(t)-E(0)+integral(D)=0,
    with E=0.5*v'Mv+0.5*u'Ku+0.5*theta'Ct*theta/T_ref and
    D=theta'Kt*theta/T_ref >=0. The actual backward-Euler/Newmark split
    has time-discretization error, so this is a *convergence diagnostic*,
    not an assertion of exact stepwise conservation.
    """
    a=config['analysis']; mat=config['material']
    if a.get('thermoelastic_mode')!='coupled_dynamic' or not a.get('thermoelastic_feedback',True):
        raise ValueError('audit requires fully coupled transient dynamics')
    if a.get('temperature_bc',{}) or a.get('inward_heat_flux',{}) or a.get('traction_bc',{}):
        raise ValueError('audit requires insulated thermal boundary and zero mechanical traction')
    if float(a.get('volumetric_heat_source',0.))!=0:
        raise ValueError('audit requires no heat source')
    if any(float(x)!=0 for bc in a.get('displacement_bc',{}).values() for x in bc.values()):
        raise ValueError('audit requires zero-valued displacement supports')
    if result.displacement_history is None or result.velocity_history is None:
        raise ValueError('audit requires coupled displacement and velocity histories')
    mesh=result.mesh
    Ku,_=assemble_elastic(mesh,mat,np.full(len(mesh.nodes),mat['T_ref']))
    M,_=assemble_inertia_coupling(mesh,mat)
    Kt,Ct=assemble_thermal(mesh,mat)
    U=result.displacement_history
    V=result.velocity_history
    theta=result.temperature_history-mat['T_ref']
    e_mech=np.array([.5*(v@(M@v)+u@(Ku@u)) for u,v in zip(U,V)])
    e_thermal=np.array([.5*(x@(Ct@x))/mat['T_ref'] for x in theta])
    d=np.array([(x@(Kt@x))/mat['T_ref'] for x in theta])
    integral=np.r_[0.,np.cumsum(.5*(d[1:]+d[:-1])*np.diff(result.times))]
    balance=e_mech+e_thermal-e_mech[0]-e_thermal[0]+integral
    return EnergyAudit(result.times,e_mech,e_thermal,d,integral,balance)
