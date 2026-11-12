"""Monolithic linear two-way thermoelastic dynamics (T3 or affine T6).

Small-strain plane stress, Fourier heat conduction, constant material data.
Thermal unknown is theta = T - T_ref. The 2D plane-stress heat capacity is
interpreted as the reduced constant-in-plane-strain heat capacity (see docs).

    M * u_ddot + Ku * u - G * theta = f
    Ct * theta_dot + Kt * theta + T_ref * G.T * u_dot = q

Newmark average acceleration (beta=1/4, gamma=1/2) advances mechanics,
backward Euler advances temperature in one coupled sparse linear solve.
"""
from collections.abc import Mapping

import numpy as np
from scipy.sparse import bmat, coo_matrix
from scipy.sparse.linalg import spsolve, factorized

from .element import TRIANGLE_QUADRATURE, element_type_for_mesh, shape_and_gradients, strain_displacement
from .fem import (
    _validate_material, assemble_thermal, assemble_elastic, boundary_load,
    prescribed_displacements, prescribed_temperature, plane_stress_matrix,
    stresses, solve_dirichlet,
)


def assemble_inertia_coupling(mesh, mat: Mapping):
    """Consistent mechanical mass M and reciprocal thermoelastic coupling G.

    G is 2N x N with mechanical force loading f_thermal = G*(T-T_ref).
    The heat-equation velocity block is H = T_ref * G.T, not an independent
    empirically chosen matrix, enforcing a thermodynamic reciprocity identity.
    """
    _validate_material(mat)
    n = len(mesh.nodes)
    kind = element_type_for_mesh(mesh)
    D = plane_stress_matrix(mat["E"], mat["nu"])
    thermal_constitutive = D @ np.array([mat["alpha"], mat["alpha"], 0.0])
    md, gd = [], []
    mi, mj, gi, gj = [], [], [], []
    mass_density = mat["density"] * mat["thickness"]
    for conn in mesh.elements:
        p = mesh.nodes[conn[:3]]
        k = len(conn)
        local_mass = np.zeros((2*k, 2*k))
        local_g = np.zeros((2*k, k))
        for r, s, w in TRIANGLE_QUADRATURE:
            N, dN, area = shape_and_gradients(p, kind, r, s)
            measure = mat["thickness"] * area * w
            scalar_mass = mass_density * area * w * np.outer(N, N)
            local_mass[0::2, 0::2] += scalar_mass
            local_mass[1::2, 1::2] += scalar_mass
            B = strain_displacement(dN)
            local_g += measure * np.outer(B.T @ thermal_constitutive, N)
        mechanical_dofs = np.ravel([[2*int(a), 2*int(a)+1] for a in conn])
        for a, ia in enumerate(mechanical_dofs):
            for b, ib in enumerate(mechanical_dofs):
                mi.append(int(ia)); mj.append(int(ib)); md.append(local_mass[a, b])
            for b, ib in enumerate(conn):
                gi.append(int(ia)); gj.append(int(ib)); gd.append(local_g[a, b])
    M = coo_matrix((md, (mi, mj)), shape=(2*n, 2*n)).tocsr()
    G = coo_matrix((gd, (gi, gj)), shape=(2*n, n)).tocsr()
    return M, G


def _sparse_dirichlet(A, rhs, prescribed):
    """Solve an assembled, nonsymmetric block system with essential constraints."""
    ndof = A.shape[0]
    if not prescribed:
        return spsolve(A.tocsc(), rhs)
    fixed = np.asarray(sorted(prescribed), dtype=int)
    if fixed.min() < 0 or fixed.max() >= ndof:
        raise ValueError("prescribed coupled DOF outside system")
    out = np.zeros(ndof)
    out[fixed] = np.asarray([prescribed[int(i)] for i in fixed])
    free = np.setdiff1d(np.arange(ndof), fixed)
    if len(free):
        out[free] = spsolve(A[free, :][:, free].tocsc(), rhs[free] - A[free, :][:, fixed] @ out[fixed])
    if not np.all(np.isfinite(out)):
        raise ValueError("coupled system is singular or nonfinite; verify constraints")
    return out


def solve_coupled_dynamic(config):
    """Fully coupled small-strain transient problem, returning a Solution.

    Constant mechanical and thermal Dirichlet data are supported. Thermal
    boundaries may be pure Neumann because transient Ct is positive definite.
    No general time-varying boundary functions or temperature-dependent
    materials are implied by this entry point.
    """
    from .geometry import make_mesh
    from .solve import Solution, _heat_load

    mesh = make_mesh(**config["geometry"])
    mat = config["material"]
    a = config["analysis"]
    if a.get("thermal_mode") != "transient":
        raise ValueError("coupled_dynamic requires thermal_mode='transient'")
    dt = float(a["time_step"])
    nstep_raw = a["steps"]
    if isinstance(nstep_raw, bool) or not isinstance(nstep_raw, int) or not (1 <= nstep_raw <= 10000):
        raise ValueError("steps must be an integer from 1 to 10000")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("time_step must be positive and finite")
    nstep = nstep_raw
    for label in ("initial_displacement", "initial_velocity"):
        if label in a and a[label] != 0:
            raise ValueError(f"{label} currently supports zero only")
    initT = float(a.get("initial_temperature", mat["T_ref"]))
    if not np.isfinite(initT):
        raise ValueError("initial temperature must be finite")

    Kt, Ct = assemble_thermal(mesh, mat)
    # Zero-stress reference temperature isolates external forces from thermal loading.
    ref = np.full(len(mesh.nodes), mat["T_ref"], dtype=float)
    Ku, f_mech = assemble_elastic(mesh, mat, ref, a.get("traction_bc", {}))
    M, G = assemble_inertia_coupling(mesh, mat)
    H = mat["T_ref"] * G.T
    heat_load = _heat_load(mesh, mat, a)
    feedback = a.get("thermoelastic_feedback", True)
    if not isinstance(feedback, bool):
        raise ValueError("thermoelastic_feedback must be boolean")
    if not feedback:
        H = H * 0
    fixed_u = prescribed_displacements(mesh, a.get("displacement_bc", {}))
    if not fixed_u:
        raise ValueError("mechanics requires adequate displacement constraints")
    fixedT = prescribed_temperature(mesh, a.get("temperature_bc", {}))
    nu = 2*len(mesh.nodes)
    nt = len(mesh.nodes)
    fixed = dict(fixed_u)
    fixed.update({nu+i: T-mat["T_ref"] for i, T in fixedT.items()})
    u = np.zeros(nu)
    for i, value in fixed_u.items():
        u[i] = value
    v = np.zeros(nu)
    theta = np.full(nt, initT-mat["T_ref"])
    for i, T in fixedT.items():
        theta[i] = T-mat["T_ref"]
    # Initial acceleration is consistent with moment balance on unconstrained
    # DOFs. Prescribed displacement values are constant, so a_fixed=0.
    initial_rhs = f_mech + G@theta - Ku@u
    accel = np.zeros(nu)
    free_u = np.setdiff1d(np.arange(nu), np.asarray(sorted(fixed_u)))
    if len(free_u):
        accel[free_u] = spsolve(M[free_u,:][:,free_u].tocsc(), initial_rhs[free_u])
    if not np.all(np.isfinite(accel)):
        raise ValueError("initial mechanical acceleration is nonfinite")

    # Standard Newmark average acceleration, BE thermal time stepping.
    beta = 0.25
    gamma = 0.5
    mcoef = 1/(beta*dt**2)
    vcoef = gamma/(beta*dt)
    A = bmat([[Ku + mcoef*M, -G],
              [vcoef*H, Kt + Ct/dt]], format="csr")
    fixed_ids = np.asarray(sorted(fixed), dtype=int)
    fixed_values = np.asarray([fixed[int(i)] for i in fixed_ids])
    free_ids = np.setdiff1d(np.arange(nu+nt), fixed_ids)
    if free_ids.size:
        reduced_solver = factorized(A[free_ids,:][:,free_ids].tocsc())
        fixed_block = A[free_ids,:][:,fixed_ids]
    times = np.arange(nstep+1)*dt
    T_history = np.empty((nstep+1, nt))
    u_history = np.empty((nstep+1, nu))
    v_history = np.empty((nstep+1, nu))
    T_history[0] = theta + mat["T_ref"]
    u_history[0] = u
    v_history[0] = v
    for step in range(1, nstep+1):
        u_pred = u + dt*v + dt**2*(0.5-beta)*accel
        v_pred = v + dt*(1-gamma)*accel
        rhs_m = f_mech + mcoef*(M @ u_pred)
        rhs_t = heat_load + (Ct @ theta)/dt - H @ (v_pred - vcoef*u_pred)
        rhs = np.r_[rhs_m, rhs_t]
        combined = np.zeros(nu+nt)
        combined[fixed_ids] = fixed_values
        if len(free_ids):
            combined[free_ids] = reduced_solver(rhs[free_ids] - fixed_block@fixed_values)
        if not np.all(np.isfinite(combined)):
            raise ValueError("coupled system produced nonfinite fields; check constraints")
        next_u, next_theta = combined[:nu], combined[nu:]
        next_accel = mcoef*(next_u-u_pred)
        next_v = v_pred + gamma*dt*next_accel
        u, theta, v, accel = next_u, next_theta, next_v, next_accel
        T_history[step] = theta+mat["T_ref"]
        u_history[step] = u
        v_history[step] = v
    temperature = theta + mat["T_ref"]
    sigma, vm = stresses(mesh, mat, temperature, u)
    # End-of-step residuals; constrained rows yield boundary reactions/fluxes.
    mechanical_res = M@accel + Ku@u-G@theta-f_mech
    thermal_res = Ct@((T_history[-1]-T_history[-2])/dt) + Kt@theta + H@v - heat_load
    return Solution(
        mesh, temperature, u, sigma, vm, thermal_res, mechanical_res,
        times, T_history, dict(mat),
        displacement_history=u_history, velocity_history=v_history,
        mechanical_acceleration=accel, two_way_feedback=feedback,
    )
