"""Independent isoparametric Q4 spatial-reference checks (v0.5)."""
import copy
import numpy as np
import pytest
from scipy.spatial import cKDTree

from thermofem2d.q4_reference import (quad_mesh, assemble_quad, q4_patch_tests,
                                      solve_q4_coupled)
from thermofem2d.solve import solve_model


def _small_coupled(base_config, kind="T3", n=3):
    c=copy.deepcopy(base_config)
    c["geometry"].update(per_side=n,radial=n,element_type=kind,
                         hole_cx=0.12,hole_cy=-0.08)
    c["analysis"].update(thermoelastic_mode="coupled_dynamic",
                         thermal_mode="transient",thermoelastic_feedback=True,
                         time_step=2e-5,steps=8,initial_temperature=298.0,
                         temperature_bc={"left":310.,"right":298.},
                         traction_bc={"top":[0.,2.e6]})
    return c


@pytest.mark.parametrize("n",[3,6])
def test_q4_mesh_area_and_jacobian(base_config,n):
    c=_small_coupled(base_config,n=n)
    g=c["geometry"]
    m=quad_mesh(g)
    assert len(m.cells)==4*n*n
    assert len(m.nodes)==4*n*(n+1)
    from thermofem2d.q4_reference import q4_shape
    area=0.
    for con in m.cells:
        p=m.nodes[con]
        for xi in [-1/np.sqrt(3),1/np.sqrt(3)]:
            for eta in [-1/np.sqrt(3),1/np.sqrt(3)]:
                _,grad=q4_shape(xi,eta)
                det=np.linalg.det(p.T@grad)
                assert det>0
                area+=det
    exact=g["width"]*g["height"]-g["hole_width"]*g["hole_height"]
    assert np.isclose(area,exact,rtol=1e-12,atol=1e-12)


@pytest.mark.parametrize("n",[3,5])
def test_q4_affine_patch_exact(base_config,n):
    c=_small_coupled(base_config,n=n)
    errT,errU=q4_patch_tests(quad_mesh(c["geometry"]),c["material"])
    assert errT<2e-9
    assert errU<2e-16


def test_q4_matrices_symmetric_positive_and_reciprocal(base_config):
    c=_small_coupled(base_config)
    mesh=quad_mesh(c["geometry"])
    Kt,Ct,Ku,M,G=assemble_quad(mesh,c["material"])
    n=len(mesh.nodes)
    assert Kt.shape==(n,n) and Ku.shape==(2*n,2*n)
    for mat in (Kt,Ct,Ku,M):
        assert np.max(abs((mat-mat.T).data),initial=0.) < 1e-12*max(1.,np.max(np.abs(mat.data)))
    assert np.linalg.eigvalsh(Ct.toarray()).min()>0
    assert np.linalg.eigvalsh(M.toarray()).min()>0
    assert abs(G.toarray()).max()>0
    rx=np.tile([1.,0.],n);ry=np.tile([0.,1.],n)
    assert np.linalg.norm(G.T@rx)<1e-6
    assert np.linalg.norm(G.T@ry)<1e-6


def test_q4_thermomechanical_feedback_is_active(base_config):
    c=_small_coupled(base_config)
    c["analysis"].update(temperature_bc={},traction_bc={"top":[0.,2.e6]},
                         initial_temperature=298.)
    active=solve_q4_coupled(c)
    disabled=copy.deepcopy(c)
    disabled["analysis"]["thermoelastic_feedback"]=False
    inactive=solve_q4_coupled(disabled)
    assert abs(active.temperature-298.).max()>1e-8
    assert abs(inactive.temperature-298.).max()<1e-9
    assert np.max(abs(active.displacement))>0


@pytest.mark.parametrize("kind",["T3","T6"])
def test_q4_reference_comparison_returns_finite_nonidentical_fields(base_config,kind):
    c=_small_coupled(base_config,kind=kind,n=3)
    prod=solve_model(c)
    ref=solve_q4_coupled(c)
    dist,idx=cKDTree(prod.mesh.nodes).query(ref.mesh.nodes)
    assert dist.max()<1e-12
    td=ref.temperature-prod.temperature[idx]
    ud=ref.displacement.reshape(-1,2)-prod.displacement.reshape(-1,2)[idx]
    assert np.all(np.isfinite(td)) and np.all(np.isfinite(ud))
    assert 1e-8<np.sqrt(np.mean(td**2))<.1
    assert 1e-12<np.sqrt(np.mean(np.sum(ud**2,axis=1)))<1e-4


def test_t3_vs_q4_comparison_reduces_with_refinement(base_config):
    output=[]
    for n in (3,6):
        c=_small_coupled(base_config,n=n)
        prod=solve_model(c);ref=solve_q4_coupled(c)
        dist,idx=cKDTree(prod.mesh.nodes).query(ref.mesh.nodes)
        assert dist.max()<1e-12
        terr=float(np.sqrt(np.mean((prod.temperature[idx]-ref.temperature)**2)))
        uerr=float(np.sqrt(np.mean(np.sum((prod.displacement.reshape(-1,2)[idx]-ref.displacement.reshape(-1,2))**2,axis=1))))
        output.append((terr,uerr))
    assert output[1][0]<output[0][0]
    assert output[1][1]<output[0][1]


def test_q4_rejects_inverted_or_invalid_hole(base_config):
    c=_small_coupled(base_config)
    c["geometry"]["hole_cx"]=1.1
    with pytest.raises(ValueError): quad_mesh(c["geometry"])
