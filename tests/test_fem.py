import numpy as np
import pytest
from thermofem2d.geometry import make_mesh
from thermofem2d.fem import (assemble_thermal,assemble_elastic,plane_stress_matrix,
    prescribed_temperature,prescribed_displacements,solve_dirichlet,stresses)
from thermofem2d.solve import solve_model


def test_thermal_matrices_symmetric_positive(base_config):
    mesh=make_mesh(**base_config["geometry"])
    K,C=assemble_thermal(mesh,base_config["material"])
    assert np.max(abs((K-K.T).data),initial=0)<1e-9
    assert np.max(abs((C-C.T).data),initial=0)<1e-6
    assert np.linalg.norm(K@np.ones(len(mesh.nodes)))<1e-8
    assert C.diagonal().min()>0


def test_uniform_temperature_has_no_gradient(base_config):
    base_config["analysis"]["temperature_bc"]={x:350 for x in ("left","right","bottom","top","hole")}
    res=solve_model(base_config)
    assert np.max(abs(res.temperature-350))<1e-8
    assert np.max(res.von_mises)<100


def test_free_thermal_expansion_patch(base_config):
    base_config["analysis"]["temperature_bc"]={"left":350,"right":350}
    s=solve_model(base_config)
    xy=s.mesh.nodes; a=base_config["material"]["alpha"]*(350-298)
    ux_expected=a*(xy[:,0]+1)
    uy_expected=a*(xy[:,1]+1)
    assert np.max(abs(s.displacement[::2]-ux_expected))<3e-12
    assert np.max(abs(s.displacement[1::2]-uy_expected))<3e-12
    assert max(s.von_mises)<1


def test_linear_steady_temperature_patch(base_config):
    base_config["analysis"]["temperature_bc"]={"left":380,"right":298}
    s=solve_model(base_config)
    assert s.temperature.min()>=297.99999
    assert s.temperature.max()<=380.00001
    mask=np.ones(len(s.mesh.nodes),bool)
    mask[np.unique(np.r_[s.mesh.node_sets["left"],s.mesh.node_sets["right"]])]=False
    assert np.max(abs(s.thermal_reactions[mask]))<1e-8


def test_heat_equilibrium(base_config):
    s=solve_model(base_config)
    assert abs(s.thermal_reactions.sum())<1e-8
    assert np.all(np.isfinite(s.sigma))
    assert np.all(s.von_mises>=0)


def test_thermal_transient_late_time_approaches_steady(base_config):
    steady=solve_model(base_config)
    a=base_config["analysis"]
    a.update(thermal_mode="transient",initial_temperature=298.0,time_step=25000,steps=180)
    trans=solve_model(base_config)
    assert len(trans.times)==181
    assert np.max(abs(trans.temperature-steady.temperature))<.01


def test_conflicting_boundary_temperature(base_config):
    base_config["analysis"]["temperature_bc"]={"left":310,"bottom":330}
    with pytest.raises(ValueError,match="conflicting"):
        solve_model(base_config)


def test_bad_displacement_label(base_config):
    base_config["analysis"]["displacement_bc"]={"uz":{"left":0}}
    with pytest.raises(ValueError,match="unknown displacement"):
        solve_model(base_config)


def test_plane_stress_constitutive_is_positive():
    D=plane_stress_matrix(200e9,.3)
    assert np.all(np.linalg.eigvalsh(D)>0)
    assert D[0,1]==pytest.approx(D[1,0])


def test_invalid_time_step(base_config):
    base_config["analysis"].update(thermal_mode="transient",time_step=0,steps=4)
    with pytest.raises(ValueError,match="time step"):
        solve_model(base_config)
