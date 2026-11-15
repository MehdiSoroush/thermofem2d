"""Independent post-processing audit of coupled-wave verification runs.

Derives a stepwise discrete quadratic perturbation-energy budget from the
saved displacement, velocity and temperature histories. Does not change the
solver, impose temperature floors, or alter any simulation's dynamics.

Scope: constant-temperature essential boundaries, zero mechanical support
velocities, no distributed heat-source/heat-flux forcing, full reciprocal
thermoelastic feedback, and uniform time steps. See accompanying report.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse.linalg import spsolve


@dataclass
class CoupledEnergyBudget:
    times: np.ndarray
    kinetic_J: np.ndarray
    elastic_J: np.ndarray
    thermal_perturbation_J: np.ndarray
    energy_change_J: np.ndarray
    mechanical_work_J: np.ndarray
    temperature_boundary_work_J: np.ndarray
    conductivity_dissipation_J: np.ndarray
    backward_euler_dissipation_J: np.ndarray
    coupling_time_quadrature_J: np.ndarray
    identity_residual_J: np.ndarray
    relative_identity_residual: float
    max_mechanical_free_residual_N: float
    max_thermal_free_residual_W: float
    min_temperature_K: float
    min_temperature_node: int
    min_temperature_time_s: float
    cooling_node_count: int
    adiabatic_cooling_indicator_K: np.ndarray


def audit_coupled_energy(system, result, *, thermal_case='reference', feedback=True):
    """Audit Newmark-average-acceleration + backward-Euler coupled dynamics.

    The energy is a positive quadratic *perturbation* measure, not the
    calorimetric sensible heat or a statement of thermodynamic equilibrium.
    Its step balance includes the thermal backward-Euler damping and a finite
    coupling quadrature term, neither of which should be silently omitted.
    Boundary thermal work is reconstructed from the *full* thermal residual.
    """
    if not feedback:
        raise ValueError('energy audit requires reciprocal two-way thermal feedback')
    if thermal_case not in ('reference','boundary_step','preheated_equilibrium'):
        raise ValueError('unsupported thermal Dirichlet specification')
    t = np.asarray(result.times, dtype=float)
    if len(t)<2 or not np.all(np.isfinite(t)):
        raise ValueError('at least two finite time points required')
    dts=np.diff(t); dt=float(dts[0]);
    if dt<=0 or not np.allclose(dts,dt,atol=1e-16,rtol=2e-11):
        raise ValueError('constant positive time step required')
    n=len(system.mesh.nodes); nu=2*n; N=len(t)
    U=np.asarray(result.displacement_history,dtype=float).reshape(N,nu)
    V=np.asarray(result.velocity_history,dtype=float).reshape(N,nu)
    T=np.asarray(result.temperature_history,dtype=float)
    loads=np.asarray(result.load_history,dtype=float)
    if T.shape!=(N,n) or loads.shape!=(N,) or not np.all(np.isfinite(U)) or not np.all(np.isfinite(V)) or not np.all(np.isfinite(T)):
        raise ValueError('invalid state history shape or nonfinite entries')
    T0=float(system.reference_temperature)
    if T0<=0: raise ValueError('positive reference temperature required')
    th=T-T0
    fixed_u=np.asarray(sorted(system.fixed_u),dtype=int)
    if any(abs(system.fixed_u[int(i)])>1e-13 for i in fixed_u) or np.max(np.abs(V[:,fixed_u]),initial=0)>1e-11:
        raise ValueError('mechanical support velocities must vanish')
    free_u=np.setdiff1d(np.arange(nu),fixed_u)
    if thermal_case=='reference':
        fixed_t=np.array([],dtype=int)
    else:
        fixed_t=np.unique(np.r_[system.mesh.node_sets['left'],system.mesh.node_sets['right']]).astype(int)
        if np.max(np.abs(T[:,fixed_t]-T[0,fixed_t]),initial=0)>1e-9:
            raise ValueError('thermal Dirichlet values must be time invariant')
    free_t=np.setdiff1d(np.arange(n),fixed_t)
    M,Ku,G,Ct,Kt=(system.M,system.Ku,system.G,system.Ct,system.Kt)
    forces=np.outer(loads,system.traction)
    # Independent kinematic reconstruction. Mechanical equilibrium at t=0
    # uses the applied initial force amplitude in the verification harness.
    A=np.zeros_like(U)
    if free_u.size:
        A[0,free_u]=spsolve(M[free_u,:][:,free_u].tocsc(),
            (forces[0] + G@th[0] - Ku@U[0])[free_u])
    for j in range(1,N):
        pred=U[j-1]+dt*V[j-1]+0.25*dt*dt*A[j-1]
        A[j]=(U[j]-pred)*(4.0/dt**2)
    v_expected=V[:-1]+.5*dt*(A[:-1]+A[1:])
    velocity_error=float(np.max(np.abs(V[1:]-v_expected)))
    if velocity_error>max(1e-9,2e-7*np.max(abs(V))):
        raise ValueError(f'Newmark velocity reconstruction inconsistent: {velocity_error:.4e} m/s')
    kinetic=.5*np.einsum('ni,ni->n',V,(M@V.T).T)
    elastic=.5*np.einsum('ni,ni->n',U,(Ku@U.T).T)
    thermal=.5/T0*np.einsum('ni,ni->n',th,(Ct@th.T).T)
    deltaE=(kinetic+elastic+thermal)-(kinetic[0]+elastic[0]+thermal[0])
    w_mech=np.zeros(N);w_bound=np.zeros(N);d_cond=np.zeros(N);d_be=np.zeros(N);w_coupling=np.zeros(N)
    max_rm=0.;max_rt=0.
    for j in range(1,N):
        v_mid=.5*(V[j]+V[j-1]);theta_mid=.5*(th[j]+th[j-1]);dth=th[j]-th[j-1]
        therm_res=(Ct@dth)/dt+Kt@th[j]+T0*G.T@V[j]
        mech_res=M@A[j]+Ku@U[j]-G@th[j]-forces[j]
        max_rm=max(max_rm,float(np.max(np.abs(mech_res[free_u]),initial=0.)))
        max_rt=max(max_rt,float(np.max(np.abs(therm_res[free_t]),initial=0.)))
        w_mech[j]=w_mech[j-1]+dt*float(v_mid @ (.5*(forces[j]+forces[j-1])))
        w_bound[j]=w_bound[j-1]+dt/T0*float(th[j,fixed_t]@therm_res[fixed_t])
        d_cond[j]=d_cond[j-1]+dt/T0*float(th[j]@(Kt@th[j]))
        d_be[j]=d_be[j-1]+.5/T0*float(dth@(Ct@dth))
        w_coupling[j]=w_coupling[j-1]+dt*float(theta_mid@(G.T@v_mid)-th[j]@(G.T@V[j]))
    # This expression should vanish by summing the exact two discrete
    # equations. Positive d_cond and d_be are separate from boundary input.
    residual=deltaE+d_cond+d_be-w_mech-w_bound-w_coupling
    scale=max(1.,float(np.max(abs(deltaE))),float(np.max(abs(w_mech))),
              float(np.max(abs(w_bound))),float(np.max(abs(d_cond))))
    # Local energy-equivalent adiabatic temperature tendency from mechanical
    # work conjugate G.T*v; useful as a sign check, *not* a continuum result.
    # We report C^{-1} G.T (u_n-u0) because this is the thermal change
    # predicted by the adiabatic C*dtheta=-T0*G.T*du limiting model.
    from scipy.sparse.linalg import spsolve as _solve
    indicator=-T0*_solve(Ct.tocsc(),G.T@(U[-1]-U[0]))
    imin=int(np.argmin(T)); it,node=np.unravel_index(imin,T.shape)
    return CoupledEnergyBudget(t,kinetic,elastic,thermal,deltaE,w_mech,w_bound,
      d_cond,d_be,w_coupling,residual,float(np.max(abs(residual))/scale),
      max_rm,max_rt,float(T[it,node]),int(node),float(t[it]),
      int(np.count_nonzero(T[-1]<T0-1e-9)),np.asarray(indicator))
