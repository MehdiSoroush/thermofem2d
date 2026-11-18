"""Independent field-level checks for circular-hole isoparametric elements."""
import numpy as np
import pytest

from thermofem2d.circular_geometry import make_circular_mesh, integrated_geometry_area, hole_boundary_length
from thermofem2d.curved_field import circular_field_mesh
from thermofem2d.coupled import assemble_inertia_coupling
from thermofem2d.element import (
    has_curved_elements, integration_rule, shape_and_gradients,
    CURVED_TRIANGLE_QUADRATURE, strain_displacement,
)
from thermofem2d.fem import (
    assemble_thermal, assemble_elastic, boundary_load,
    solve_dirichlet, plane_stress_matrix, stresses, element_heat_flux,
)

CENTRE = np.array([0.12, -0.07])
RADIUS = 0.25


def curved(kind="T6", subdivisions=4):
    geo = make_circular_mesh(width=2., height=2., hole_radius=RADIUS,
                             hole_cx=float(CENTRE[0]), hole_cy=float(CENTRE[1]),
                             per_side=subdivisions, radial=subdivisions,
                             element_type=kind)
    return geo, circular_field_mesh(geo)


def boundary_ids(mesh):
    return np.unique(np.concatenate(tuple(mesh.node_sets.values())))


def log_heat_error(mesh, mat):
    r = np.linalg.norm(mesh.nodes-CENTRE, axis=1)
    exact = 300. + 8.*np.log(r/RADIUS)
    boundary = boundary_ids(mesh)
    k, _ = assemble_thermal(mesh, mat)
    prediction = solve_dirichlet(k, np.zeros(len(mesh.nodes)),
                                 {int(i): float(exact[i]) for i in boundary})
    return float(np.sqrt(np.mean((prediction-exact)**2)))


@pytest.mark.parametrize("kind", ["T3", "T6"])
def test_explicit_geometry_adapter_preserves_data(kind):
    geometry, mesh = curved(kind)
    assert np.array_equal(geometry.coordinates, mesh.nodes)
    assert np.array_equal(geometry.triangles, mesh.elements)
    assert not np.shares_memory(geometry.coordinates, mesh.nodes)
    assert mesh.edges["hole"].shape[1] == (3 if kind == "T6" else 2)
    assert has_curved_elements(mesh) == (kind == "T6")


def test_circular_geometry_cannot_be_assembled_without_explicit_adapter(base_config):
    geometry, _ = curved()
    with pytest.raises((AttributeError, TypeError)):
        assemble_thermal(geometry, base_config["material"])


def test_invalid_midside_geometry_rejected():
    geometry, _ = curved()
    geometry.coordinates[geometry.triangles[0,3]] += np.array([100., 100.])
    with pytest.raises(ValueError, match="Jacobian"):
        circular_field_mesh(geometry)


def test_curved_shape_partition_gradient_consistency():
    _, mesh = curved()
    tri = next(e for e in mesh.elements if np.any(np.isin(e[3:], mesh.node_sets["hole"])))
    p = mesh.nodes[tri]
    for r,s,_ in CURVED_TRIANGLE_QUADRATURE:
        N,g,measure=shape_and_gradients(p,"T6",r,s)
        assert np.isclose(N.sum(),1.,atol=1e-13)
        np.testing.assert_allclose(g.sum(axis=0),0.,atol=1e-12)
        assert measure>0
        # Physical coordinates must have identity gradient.
        np.testing.assert_allclose(p.T@g, np.eye(2),atol=3e-12)


@pytest.mark.parametrize("kind",["T3","T6"])
def test_thermal_heat_capacity_and_stiffness_energy(kind,base_config):
    geom,mesh=curved(kind)
    mat=base_config["material"]
    k,c=assemble_thermal(mesh,mat)
    ones=np.ones(len(mesh.nodes))
    area=integrated_geometry_area(geom)
    np.testing.assert_allclose(ones@c@ones, area*mat["density"]*mat["heat_capacity"]*mat["thickness"], rtol=2e-11)
    assert abs(ones@k@ones) < 5e-7  # floating-point cancellation, dimensional units
    assert np.max(np.abs(k@ones)) < 1e-5
    np.testing.assert_allclose(k.toarray(),k.T.toarray(),rtol=1e-13,atol=1e-8)
    np.testing.assert_allclose(c.toarray(),c.T.toarray(),rtol=1e-13,atol=1e-7)
    assert c.diagonal().min()>0


def test_analytic_circular_log_heat_convergence(base_config):
    mat=base_config["material"]
    errors={kind: [log_heat_error(curved(kind,n)[1],mat) for n in (3,5,8)]
            for kind in ("T3","T6")}
    assert errors["T3"][-1] < .25*errors["T3"][0]
    assert errors["T6"][-1] < .1*errors["T6"][0]
    assert errors["T6"][-1] < .1*errors["T3"][-1]



def test_affine_temperature_flux_is_exact_on_curved_geometry(base_config):
    mat=base_config["material"]
    _,mesh=curved("T6",5)
    T=301.0+3.0*mesh.nodes[:,0]-2.0*mesh.nodes[:,1]
    result=element_heat_flux(mesh,mat,T)
    reference=-mat["conductivity"]*np.array([3.0,-2.0])
    np.testing.assert_allclose(result,np.broadcast_to(reference,result.shape),rtol=1e-12,atol=1e-10)


def test_affine_elastic_patch_and_constant_strain_stress(base_config):
    mat=base_config["material"]
    _,mesh=curved("T6",5)
    coords=mesh.nodes
    A=np.array([[2.3e-5,-1.4e-5],[0.7e-5,1.1e-5]])
    u=(coords@A.T+np.array([0.003,-0.001])).reshape(-1)
    T=np.full(len(coords),mat["T_ref"])
    K,f=assemble_elastic(mesh,mat,T)
    fixed={2*int(i)+d:u[2*i+d] for i in boundary_ids(mesh) for d in (0,1)}
    solved=solve_dirichlet(K,f,fixed)
    assert np.max(abs(solved-u))<2e-11
    expected=plane_stress_matrix(mat["E"],mat["nu"])@np.array([A[0,0],A[1,1],A[0,1]+A[1,0]])
    sigma,_=stresses(mesh,mat,T,solved)
    np.testing.assert_allclose(sigma,np.broadcast_to(expected, sigma.shape),rtol=2e-8,atol=3e-3)


def test_free_uniform_thermal_expansion_is_stress_free(base_config):
    mat=base_config["material"]
    _,mesh=curved("T6",4)
    dt=18.
    T=np.full(len(mesh.nodes),mat["T_ref"]+dt)
    K,f=assemble_elastic(mesh,mat,T)
    expected=(mat["alpha"]*dt*mesh.nodes).reshape(-1)
    fixed={2*int(i)+d:float(expected[2*i+d]) for i in boundary_ids(mesh) for d in (0,1)}
    solved=solve_dirichlet(K,f,fixed)
    assert np.max(abs(solved-expected)) < 2e-12
    sigma,_=stresses(mesh,mat,T,solved)
    assert np.max(abs(sigma)) < .02


def test_curved_edge_scalar_flux_resultant_and_traction(base_config):
    geo,mesh=curved("T6",6)
    length=hole_boundary_length(geo)
    raw=boundary_load(mesh,{"hole":2.4},components=1)
    np.testing.assert_allclose(raw.sum(),2.4*length,rtol=2e-12)
    applied=boundary_load(mesh,{"hole":[3.,-2.]},multiplier=base_config["material"]["thickness"],components=2)
    np.testing.assert_allclose(applied.reshape(-1,2).sum(axis=0),
                               np.array([3.,-2.])*length*base_config["material"]["thickness"],rtol=2e-12)
    assert np.count_nonzero(raw)>len(mesh.edges["hole"])


def test_inertia_mass_and_thermomechanical_work_identity(base_config):
    mat=base_config["material"]
    geom,mesh=curved("T6",4)
    M,G=assemble_inertia_coupling(mesh,mat)
    mass=mat["density"]*mat["thickness"]*integrated_geometry_area(geom)
    n=len(mesh.nodes)
    for j in (0,1):
        a=np.zeros(2*n);a[j::2]=1
        np.testing.assert_allclose(a@M@a,mass,rtol=2e-11)
    assert np.max(np.abs((M-M.T).data))<1e-7 if (M-M.T).nnz else True
    beta=mat["E"]*mat["alpha"]/(1-mat["nu"])
    u=(mesh.nodes*np.array([1.1e-5, -0.4e-5])).reshape(-1)
    theta=np.full(n,6.0)
    analytic=beta*(1.1e-5-0.4e-5)*6.0*mat["thickness"]*integrated_geometry_area(geom)
    np.testing.assert_allclose(u@(G@theta),analytic,rtol=2e-11)
    np.testing.assert_allclose(u@(G@theta),theta@(G.T@u),rtol=1e-12)


def test_linear_element_stays_compatible(base_config):
    mat=base_config["material"]
    _,mesh=curved("T3",4)
    assert integration_rule(mesh) is not CURVED_TRIANGLE_QUADRATURE
    k,_=assemble_thermal(mesh,mat)
    assert k.shape==(len(mesh.nodes),len(mesh.nodes))
