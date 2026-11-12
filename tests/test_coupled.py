"""Physical and algebraic checks for v0.3 coupled transient thermomechanics."""
import copy
import numpy as np
import pytest
from thermofem2d.geometry import make_mesh
from thermofem2d.fem import assemble_elastic, assemble_thermal, prescribed_displacements
from thermofem2d.coupled import assemble_inertia_coupling
from thermofem2d.solve import solve_model


def coupled(base, element="T3", *, feedback=True):
    case = copy.deepcopy(base)
    case["geometry"].update(per_side=3, radial=3, element_type=element)
    case["analysis"].update(
        thermoelastic_mode="coupled_dynamic", thermal_mode="transient",
        time_step=4e-5, steps=6,
        initial_temperature=298.0, temperature_bc={},
        thermoelastic_feedback=feedback,
    )
    return case


@pytest.mark.parametrize("kind", ["T3", "T6"])
def test_reciprocal_matrices_and_consistent_mass(base_config, kind):
    c = coupled(base_config, kind)
    m = make_mesh(**c["geometry"])
    mat = c["material"]
    M, G = assemble_inertia_coupling(m, mat)
    assert M.shape == (2*len(m.nodes), 2*len(m.nodes))
    assert G.shape == (2*len(m.nodes), len(m.nodes))
    assert np.max(abs((M-M.T).data), initial=0) < 1e-9
    assert np.linalg.eigvalsh(M.toarray()).min() > 0
    # Rigid uniform translations cannot produce thermal force coupling.
    assert np.linalg.norm(G.T @ np.tile([1.0, 0.0], len(m.nodes))) < 1e-6
    assert np.linalg.norm(G.T @ np.tile([0.0, 1.0], len(m.nodes))) < 1e-6
    temp = 298.0 + 11*m.nodes[:, 0] - 5*m.nodes[:, 1]
    _, f = assemble_elastic(m, mat, temp)
    assert np.allclose(f, G@(temp-mat["T_ref"]), rtol=1e-11, atol=1e-5)
    assert abs(G.toarray()).max()>0


@pytest.mark.parametrize("kind", ["T3", "T6"])
def test_equilibrium_at_reference_temperature(base_config, kind):
    c = coupled(base_config, kind)
    c["analysis"]["traction_bc"] = {}
    s = solve_model(c)
    assert np.max(abs(s.temperature-298)) < 1e-11
    assert np.linalg.norm(s.displacement) < 1e-14
    assert np.max(s.von_mises) < 1e-6
    assert s.displacement_history.shape==(7, 2*len(s.mesh.nodes))
    assert s.velocity_history.shape==s.displacement_history.shape
    assert s.two_way_feedback


def test_feedback_from_mechanical_motion_when_thermally_insulated(base_config):
    c = coupled(base_config, "T3")
    c["analysis"]["traction_bc"] = {"top": [0., 1.e7]}
    s = solve_model(c)
    assert np.max(abs(s.temperature-298.))>1e-6
    assert np.max(abs(s.displacement))>1e-9
    no_feedback = coupled(base_config,"T3",feedback=False)
    no_feedback["analysis"]["traction_bc"]={"top":[0.,1.e7]}
    sf=solve_model(no_feedback)
    assert np.max(abs(sf.temperature-298.))<1e-10
    assert np.max(abs(sf.displacement))>1e-9
    # Thermoelastic feedback changes motion compared with one-way dynamics.
    assert np.linalg.norm(s.displacement-sf.displacement)>1e-13


@pytest.mark.parametrize("kind", ["T3", "T6"])
def test_coupled_block_equilibrium_residuals(base_config, kind):
    c=coupled(base_config,kind)
    c["analysis"].update(
        temperature_bc={"left":310., "right":298.},
        traction_bc={"top":[0.,2.e7]},
    )
    # Dirichlet-consistency at corners is automatic here: left and right disjoint.
    s=solve_model(c)
    M,G=assemble_inertia_coupling(s.mesh,c["material"])
    Kt,Ct=assemble_thermal(s.mesh,c["material"])
    Ku,f=assemble_elastic(s.mesh,c["material"],np.full(len(s.mesh.nodes),298.),c["analysis"]["traction_bc"])
    constraints=prescribed_displacements(s.mesh,c["analysis"]["displacement_bc"])
    u_free=np.setdiff1d(np.arange(len(s.displacement)),np.asarray(list(constraints)))
    th=s.temperature-298.
    lhs=M@s.mechanical_acceleration + Ku@s.displacement - G@th - f
    assert np.max(abs(lhs[u_free]))<5e-3  # N, relative scale ~1e7
    temp_last=(s.temperature_history[-1]-s.temperature_history[-2])/c["analysis"]["time_step"]
    from thermofem2d.solve import _heat_load
    rhs=Ct@temp_last+Kt@th+298.*G.T@s.velocity_history[-1]-_heat_load(s.mesh,c["material"],c["analysis"])
    fixedT=np.unique(np.r_[s.mesh.node_sets["left"],s.mesh.node_sets["right"]])
    freeT=np.setdiff1d(np.arange(len(s.mesh.nodes)), fixedT)
    assert np.max(abs(rhs[freeT]))<1e-3


def test_zero_thermal_expansion_decouples_and_bad_modes(base_config):
    c=coupled(base_config)
    c["material"]["alpha"]=0.0
    c["analysis"]["temperature_bc"]={"left":310.,"right":298.}
    s=solve_model(c)
    assert np.max(abs(s.displacement))<1e-12
    assert abs(s.temperature.max()-310.)<1e-9
    bad=coupled(base_config)
    bad["analysis"]["thermal_mode"]="steady"
    with pytest.raises(ValueError, match="transient"):
        solve_model(bad)
    bad=coupled(base_config)
    bad["analysis"]["time_step"]=0
    with pytest.raises(ValueError, match="time_step"):
        solve_model(bad)


def test_positive_temperature_with_unconstrained_thermal_boundary(base_config):
    c=coupled(base_config)
    c["analysis"].update(initial_temperature=310., traction_bc={})
    s=solve_model(c)
    assert np.all(np.isfinite(s.temperature))
    assert np.all(np.isfinite(s.displacement))
