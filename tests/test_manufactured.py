"""Independent, exact-field benchmarks on meshes containing a rectangular void.

These tests compare the FEM solution to prescribed continuum solutions; they
are not comparisons with the unavailable 2016 MATLAB/Abaqus data.
"""
import numpy as np
import pytest

from thermofem2d.fem import (
    assemble_thermal, assemble_elastic, element_heat_flux, plane_stress_matrix,
    solve_dirichlet, stresses,
)
from thermofem2d.geometry import make_mesh
from thermofem2d.solve import _heat_load, solve_model


def _all_boundary_nodes(mesh):
    return np.unique(np.concatenate(list(mesh.node_sets.values())))


def test_linear_harmonic_temperature_exact_on_offset_hole(base_config):
    """Laplace solution T=320 + 4x - 7y is exactly representable by T3."""
    base_config["geometry"].update(hole_cx=0.19, hole_cy=-0.13)
    mesh = make_mesh(**base_config["geometry"])
    material = base_config["material"]
    K, _ = assemble_thermal(mesh, material)
    exact = 320 + 4 * mesh.nodes[:, 0] - 7 * mesh.nodes[:, 1]
    boundary = _all_boundary_nodes(mesh)
    prescribed = {int(i): float(exact[i]) for i in boundary}
    computed = solve_dirichlet(K, np.zeros(len(exact)), prescribed)
    assert np.max(abs(computed - exact)) < 1e-9
    flux = element_heat_flux(mesh, material, computed)
    expected = -material["conductivity"] * np.array([4, -7])
    np.testing.assert_allclose(flux, np.broadcast_to(expected, flux.shape), atol=1e-9, rtol=0)


def test_affine_plane_stress_patch_exact_on_offset_hole(base_config):
    """Constant strain: displacements affine, stress given by constitutive D."""
    mesh = make_mesh(**base_config["geometry"])
    material = base_config["material"]
    nodes = mesh.nodes
    strain = np.array([1.0e-4, -0.7e-4, 0.35e-4])
    ux = strain[0] * nodes[:, 0] + 0.5 * strain[2] * nodes[:, 1]
    uy = strain[1] * nodes[:, 1] + 0.5 * strain[2] * nodes[:, 0]
    exact = np.column_stack([ux, uy]).reshape(-1)
    T = np.full(len(nodes), material["T_ref"])
    K, F = assemble_elastic(mesh, material, T)
    fixed = {2 * int(i) + d: float(exact[2 * i + d])
             for i in _all_boundary_nodes(mesh) for d in (0, 1)}
    computed = solve_dirichlet(K, F, fixed)
    assert np.max(abs(computed - exact)) < 2e-12
    sigma, vm = stresses(mesh, material, T, computed)
    expected_stress = plane_stress_matrix(material["E"], material["nu"]) @ strain
    np.testing.assert_allclose(sigma, np.broadcast_to(expected_stress, sigma.shape), rtol=1e-8, atol=0.1)
    assert np.all(vm > 0)
    free = np.setdiff1d(np.arange(len(exact)), list(fixed))
    assert np.linalg.norm((K @ computed - F)[free], ord=np.inf) < 0.01


def test_quadratic_conduction_manufactured_solution_converges(base_config):
    """T=320+a(x²+y²), Q=-4ka, boundary T exactly specified at nodes."""
    a = 12.0
    material = base_config["material"]
    k = material["conductivity"]
    errors = []
    for resolution in (3, 6, 12):
        mesh = make_mesh(**dict(base_config["geometry"], per_side=resolution,
                                radial=resolution))
        nodes = mesh.nodes
        exact = 320 + a * np.sum(nodes * nodes, axis=1)
        K, _ = assemble_thermal(mesh, material)
        source = _heat_load(mesh, material, {"volumetric_heat_source": -4 * k * a})
        fixed = {int(i): float(exact[i]) for i in _all_boundary_nodes(mesh)}
        computed = solve_dirichlet(K, source, fixed)
        free = np.setdiff1d(np.arange(len(nodes)), list(fixed))
        assert len(free) > 0
        errors.append(np.linalg.norm(computed[free] - exact[free]) / np.sqrt(len(free)))
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < errors[0] / 4


def test_heat_capacity_energy_decay_under_fixed_temperature(base_config):
    """Backward Euler dissipates the deviation thermal energy for zero sources."""
    analysis = base_config["analysis"]
    analysis.update(thermal_mode="transient", initial_temperature=350.0,
                    temperature_bc={name: 298.0 for name in ("left", "right", "top", "bottom", "hole")},
                    volumetric_heat_source=0.0, inward_heat_flux={},
                    time_step=5000.0, steps=6)
    s = solve_model(base_config)
    _, C = assemble_thermal(s.mesh, base_config["material"])
    energies = np.array([0.5 * ((v - 298) @ (C @ (v - 298))) for v in s.temperature_history])
    assert energies[0] > 0
    assert np.all(np.diff(energies) < 0)
    assert abs(energies[-1]) < energies[0]


def test_flux_rejects_nonfinite_or_incompatible_data(base_config):
    mesh = make_mesh(**base_config["geometry"])
    mat = base_config["material"]
    with pytest.raises(ValueError, match="one finite temperature"):
        element_heat_flux(mesh, mat, np.ones(4))
    t = np.full(len(mesh.nodes), 300.0)
    t[3] = np.nan
    with pytest.raises(ValueError, match="one finite temperature"):
        element_heat_flux(mesh, mat, t)


def test_area_matches_exact_perforated_plate_for_offset_geometry():
    mesh = make_mesh(width=2.6, height=1.8, hole_width=.7, hole_height=.3,
                     hole_cx=.31, hole_cy=-.23, per_side=9, radial=7)
    p = mesh.nodes[mesh.elements]
    v1 = p[:, 1] - p[:, 0]
    v2 = p[:, 2] - p[:, 0]
    area = 0.5 * np.sum(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0])
    assert area == pytest.approx(2.6 * 1.8 - .7 * .3, abs=1e-11)
