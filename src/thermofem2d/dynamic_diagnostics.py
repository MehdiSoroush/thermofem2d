"""Controlled transient wave benchmarks for T6 and independent Q4 elements.

Dedicated v0.7 verification harness. This is not a general transient-load
boundary-condition interface for the public solver. Loads are applied to a
specified top boundary through independently assembled element matrices.
"""
from dataclasses import dataclass
from copy import deepcopy

import numpy as np
from scipy.sparse import bmat
from scipy.sparse.linalg import eigsh, splu
from scipy.spatial import cKDTree

from .geometry import make_mesh
from .fem import assemble_elastic, assemble_thermal, prescribed_displacements
from .coupled import assemble_inertia_coupling
from .q4_reference import quad_mesh, assemble_quad, _edge_load, _fixed


def ramp_amplitude(t: float, rise_time: float, profile: str) -> float:
    """Amplitude of a step or C1-continuous half-cosine ramp (dimensionless)."""
    if profile == 'step':
        return 1.0  # force active at time zero: initial acceleration uses it
    if profile != 'cosine':
        raise ValueError('load profile must be step or cosine')
    if not np.isfinite(rise_time) or rise_time <= 0:
        raise ValueError('rise time must be positive')
    z = float(np.clip(t / rise_time, 0.0, 1.0))
    return 0.5 * (1.0 - np.cos(np.pi*z))


@dataclass
class WaveSystem:
    mesh: object
    Kt: object
    Ct: object
    Ku: object
    M: object
    G: object
    traction: np.ndarray
    fixed_u: dict
    reference_temperature: float
    method: str


def assemble_wave_system(config, method: str) -> WaveSystem:
    """Use production T3/T6 matrices or independently assembled Q4 matrices."""
    g = deepcopy(config['geometry'])
    mat = config['material']
    if method == 'Q4':
        mesh = quad_mesh(g)
        Kt, Ct, Ku, M, G = assemble_quad(mesh, mat)
        traction = _edge_load(mesh, {'top': [0.0, 2e6]}, mat['thickness'], 2)
        fixed = {}
        for label, component in [('ux', 0), ('uy', 1)]:
            fixed.update(_fixed(mesh, config['analysis']['displacement_bc'][label], component))
    elif method in ('T3', 'T6'):
        mesh = make_mesh(**dict(g, element_type=method))
        Kt, Ct = assemble_thermal(mesh, mat)
        Ku, traction = assemble_elastic(
            mesh, mat, np.full(len(mesh.nodes), mat['T_ref']),
            {'top': [0.0, 2e6]},
        )
        M, G = assemble_inertia_coupling(mesh, mat)
        fixed = prescribed_displacements(mesh, config['analysis']['displacement_bc'])
    else:
        raise ValueError('method must be T3, T6 or Q4')
    if not fixed or any(abs(v) > 1e-15 for v in fixed.values()):
        raise ValueError('benchmark requires zero-valued mechanical supports')
    return WaveSystem(mesh, Kt, Ct, Ku, M, G, traction, fixed,
                      float(mat['T_ref']), method)


@dataclass
class WaveResult:
    method: str
    times: np.ndarray
    displacement_history: np.ndarray  # (nstep+1, nnode, 2)
    temperature_history: np.ndarray  # (nstep+1, nnode), K
    velocity_history: np.ndarray  # (nstep+1, nnode, 2)
    mesh: object
    load_history: np.ndarray
    initial_displacement: np.ndarray


def integrate_wave(
    system: WaveSystem, *, steps: int, final_time: float, profile: str,
    rise_time: float = 1.6e-4, feedback: bool = False,
    thermal_case: str = 'reference',
) -> WaveResult:
    """Monolithic Newmark-average-acceleration / backward-Euler time stepping.

    Ramped top traction is a dedicated verification load. The interpolation
    at every step end is exact for the selected prescribed load function.
    No heat fluxes are applied; thermal initial conditions and fixed edge
    temperatures are selected explicitly through ``thermal_case``.
    """
    if not isinstance(steps, int) or isinstance(steps, bool) or steps < 1:
        raise ValueError('steps must be a positive integer')
    if not np.isfinite(final_time) or final_time <= 0:
        raise ValueError('final_time must be positive')
    n = len(system.mesh.nodes)
    nu = 2*n
    dt = final_time / steps
    fixed = np.array(sorted(system.fixed_u), dtype=int)
    if thermal_case not in ('reference', 'boundary_step', 'preheated_equilibrium'):
        raise ValueError('unknown thermal initialization scenario')
    fixed_th = {}
    if thermal_case != 'reference':
        for label, value in [('left', 310.), ('right', system.reference_temperature)]:
            for node in system.mesh.node_sets[label]:
                fixed_th[int(node)] = value - system.reference_temperature
    thermal_ids = np.array(sorted(fixed_th), dtype=int)
    thermal_values = np.array([fixed_th[int(i)] for i in thermal_ids])
    free_u = np.setdiff1d(np.arange(nu), fixed)
    th = np.zeros(n)
    u = np.zeros(nu)
    if len(thermal_ids):
        th[thermal_ids] = thermal_values
    if thermal_case == 'preheated_equilibrium':
        free_th = np.setdiff1d(np.arange(n), thermal_ids)
        if len(free_th):
            th[free_th] = splu(system.Kt[free_th,:][:,free_th].tocsc()).solve(
                -system.Kt[free_th,:][:,thermal_ids] @ thermal_values
            )
        free_mech = np.setdiff1d(np.arange(nu), fixed)
        u[free_mech] = splu(system.Ku[free_mech,:][:,free_mech].tocsc()).solve(
            (system.G @ th)[free_mech]
        )
    initial_displacement = u.copy().reshape(n, 2)
    v = np.zeros(nu)
    t0_load = ramp_amplitude(0., rise_time, profile)
    acceleration = np.zeros(nu)
    if len(free_u):
        acceleration[free_u] = splu(system.M[free_u, :][:, free_u].tocsc()).solve(
            (t0_load * system.traction + system.G @ th - system.Ku @ u)[free_u]
        )
    beta, gamma = .25, .5
    mcoef = 1. / (beta*dt**2)
    vcoef = gamma / (beta*dt)
    H = system.reference_temperature * system.G.T if feedback else 0 * system.G.T
    blocks = bmat([[system.Ku + mcoef*system.M, -system.G],
                   [vcoef*H, system.Kt + system.Ct/dt]], format='csc')
    fixed_combined = np.r_[fixed, nu+thermal_ids]
    fixed_values = np.r_[np.zeros(len(fixed)), thermal_values]
    free = np.setdiff1d(np.arange(nu+n), fixed_combined)
    factor = splu(blocks[free, :][:, free])
    fixed_coupling = blocks[free, :][:, fixed_combined]
    times = np.linspace(0, final_time, steps+1)
    displacements = np.empty((steps+1, n, 2))
    velocities = np.empty_like(displacements)
    temperature = np.empty((steps+1, n))
    load_history = np.array([ramp_amplitude(t, rise_time, profile) for t in times])
    displacements[0] = u.reshape(n, 2)
    velocities[0] = v.reshape(n, 2)
    temperature[0] = th + system.reference_temperature
    for step in range(1, steps+1):
        up = u + dt*v + (0.5-beta)*dt*dt*acceleration
        vp = v + (1.-gamma)*dt*acceleration
        force = system.traction * load_history[step]
        r1 = force + mcoef*(system.M @ up)
        r2 = (system.Ct @ th) / dt - H @ (vp - vcoef*up)
        state = np.zeros(nu+n)
        state[fixed_combined] = fixed_values
        state[free] = factor.solve(np.r_[r1,r2][free] - fixed_coupling @ fixed_values)
        nextu, th = state[:nu], state[nu:]
        acceleration = mcoef * (nextu-up)
        v = vp + gamma*dt*acceleration
        u = nextu
        displacements[step] = u.reshape(n, 2)
        velocities[step] = v.reshape(n, 2)
        temperature[step] = th + system.reference_temperature
    return WaveResult(system.method, times, displacements, temperature, velocities,
                      system.mesh, load_history, initial_displacement)


def matched_node_indices(tri: WaveResult, quad: WaveResult):
    distances, indices = cKDTree(tri.mesh.nodes).query(quad.mesh.nodes)
    if max(distances) > 1e-10:
        raise ValueError('Meshes do not have exactly matched Q4-vertex coordinates')
    return indices


def rms_vector(values):
    values = np.asarray(values)
    return float(np.sqrt(np.mean(np.sum(values**2, axis=-1))))


def compare_wave(tri: WaveResult, quad: WaveResult, static_scale: float | None = None,
                 perturbation: bool = False):
    """Compare matched vertex histories; report absolute and relative scales."""
    if len(tri.times) != len(quad.times) or not np.allclose(tri.times, quad.times,
                                                          rtol=0, atol=1e-13):
        raise ValueError('Unequal time grids')
    ids = matched_node_indices(tri, quad)
    utri = tri.displacement_history[:, ids, :].copy()
    uquad = quad.displacement_history.copy()
    if perturbation:
        utri -= utri[0]
        uquad -= uquad[0]
    numer = np.sqrt(np.mean(np.sum((utri-uquad)**2, axis=-1), axis=-1))
    denom = np.sqrt(np.mean(np.sum(uquad**2, axis=-1), axis=-1))
    return {'times': tri.times,
            'rms_difference_m': numer,
            'q4_rms_displacement_m': denom,
            'relative_difference': np.divide(numer, denom,
                                             out=np.full_like(numer, np.nan), where=denom>1e-18),
            'static_normalized_difference': numer / static_scale if static_scale else None}


def lowest_frequencies(system: WaveSystem, count: int = 6):
    """Low-order undamped mechanical frequencies with benchmark constraints."""
    nu = 2*len(system.mesh.nodes)
    free = np.setdiff1d(np.arange(nu), np.array(sorted(system.fixed_u), dtype=int))
    k = min(count, len(free)-2)
    eigvals = eigsh(system.Ku[free,:][:,free], M=system.M[free,:][:,free],
                    k=k, sigma=0, which='LM', return_eigenvectors=False)
    if np.any(eigvals <= 0):
        raise ValueError('Benchmark elastic operator not positive definite')
    return np.sqrt(np.sort(eigvals)) / (2*np.pi)
