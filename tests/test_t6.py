"""Analytical verification tests for affine six-node triangular elements."""
import numpy as np
import pytest
from thermofem2d.geometry import make_mesh
from thermofem2d.element import (TRIANGLE_QUADRATURE,shape_and_gradients,
                                 strain_displacement,plotting_triangles)
from thermofem2d.fem import (assemble_thermal,boundary_load,solve_dirichlet,
                              element_heat_flux,stresses,assemble_elastic)
from thermofem2d.solve import solve_model,_heat_load


def test_t6_shape_kronecker_partition_and_gradient():
    vertices=np.array([[0.,0.],[2.,0.],[0.,3.]])
    rst=[(0,0),(1,0),(0,1),(.5,0),(.5,.5),(0,.5)]
    for k,(r,s) in enumerate(rst):
        N,dN,area=shape_and_gradients(vertices,"T6",r,s)
        assert N==pytest.approx(np.eye(6)[k],abs=1e-14)
        assert N.sum()==pytest.approx(1,abs=1e-14)
        assert np.linalg.norm(dN.sum(axis=0))<1e-14
        assert area==pytest.approx(3.)


def test_t6_midside_shared_and_boundaries():
    m3=make_mesh(per_side=3,radial=3)
    m6=make_mesh(per_side=3,radial=3,element_type="T6")
    assert m6.elements.shape==(len(m3.elements),6)
    assert len(m6.nodes)>len(m3.nodes)
    assert all(len(edges[0])==3 for edges in m6.edges.values())
    assert len(plotting_triangles(m6))==4*len(m6.elements)
    # exactly one shared midside node per unique topological T3 edge
    e3=m3.elements
    edges={tuple(sorted((int(t[i]),int(t[(i+1)%3])))) for t in e3 for i in range(3)}
    assert len(m6.nodes)-len(m3.nodes)==len(edges)
    assert len(np.unique(m6.elements[:,3:]))==len(edges)


def test_t6_matrices_spd_and_constant_null_mode(base_config):
    mesh=make_mesh(per_side=3,radial=3,element_type="T6")
    K,C=assemble_thermal(mesh,base_config["material"])
    assert np.max(np.abs((K-K.T).data),initial=0)<1e-8
    assert np.max(np.abs((C-C.T).data),initial=0)<1e-5
    assert np.linalg.norm(K@np.ones(len(mesh.nodes)))<1e-8
    assert np.linalg.eigvalsh(C.toarray()).min()>0
    diag=K.diagonal();assert diag.min()>0


def test_t6_edge_load_integrates_exact_constant_traction():
    m=make_mesh(per_side=3,radial=3,element_type="T6")
    loads=boundary_load(m,{"top":(123.,-45.)},components=2)
    assert loads[0::2].sum()==pytest.approx(2*123.)
    assert loads[1::2].sum()==pytest.approx(-2*45.)
    assert np.count_nonzero(loads)>0


def test_t6_manufactured_quadratic_heat_solution(base_config):
    """T=300+x²+2y²; -div(k gradT)=-6k; all boundaries exact Dirichlet."""
    mat=base_config["material"]
    Q=-6.*mat["conductivity"]
    def exact(p):return 300+p[:,0]**2+2*p[:,1]**2
    errors={}
    for kind in ("T3","T6"):
        mesh=make_mesh(per_side=5,radial=4,element_type=kind)
        K,_=assemble_thermal(mesh,mat)
        f=_heat_load(mesh,mat,{"volumetric_heat_source":Q})
        ids=np.unique(np.concatenate(tuple(mesh.node_sets.values())))
        exactT=exact(mesh.nodes)
        fixed={int(i):float(exactT[i]) for i in ids}
        solved=solve_dirichlet(K,f,fixed)
        errors[kind]=np.linalg.norm(solved-exactT)/np.sqrt(len(solved))
    assert errors["T6"]<1e-9
    assert errors["T3"]>1e-3


def test_t6_exact_linear_strain_patch_and_free_expansion(base_config):
    m=make_mesh(per_side=4,radial=3,element_type="T6")
    p=m.nodes
    disp=np.zeros((len(p),2))
    disp[:,0]=1e-3+.003*p[:,0]-.005*p[:,1]
    disp[:,1]=2e-3+.004*p[:,0]+0.002*p[:,1]
    expected=np.array([.003,.002,-.001])
    for tri in m.elements:
        for r,s,w in TRIANGLE_QUADRATURE:
            _,dN,_=shape_and_gradients(m.nodes[tri[:3]],"T6",r,s)
            eps=strain_displacement(dN)@disp[tri].ravel()
            assert eps==pytest.approx(expected,rel=0,abs=5e-15)
    base_config["geometry"].update(per_side=4,radial=3,element_type="T6")
    base_config["analysis"]["temperature_bc"]={"left":350,"right":350}
    solution=solve_model(base_config)
    a=base_config["material"]["alpha"]*(350-298)
    assert np.max(abs(solution.displacement[::2]-a*(solution.mesh.nodes[:,0]+1)))<5e-12
    assert np.max(abs(solution.displacement[1::2]-a*(solution.mesh.nodes[:,1]+1)))<5e-12
    assert max(solution.von_mises)<10


def test_t6_transient_and_flux_linear_patch(base_config):
    g=base_config["geometry"]
    g.update(per_side=3,radial=3,element_type="T6")
    a=base_config["analysis"]
    a["temperature_bc"]={"left":350,"right":350}
    a.update(thermal_mode="transient",initial_temperature=300,time_step=800.,steps=3)
    sol=solve_model(base_config)
    assert len(sol.temperature_history)==4
    assert np.all(np.isfinite(sol.temperature))
    coords=sol.mesh.nodes
    trial=200+7*coords[:,0]-3*coords[:,1]
    q=element_heat_flux(sol.mesh,base_config["material"],trial)
    assert q==pytest.approx(np.tile([-7*base_config["material"]["conductivity"],3*base_config["material"]["conductivity"]],(len(sol.mesh.elements),1)),rel=0,abs=1e-10)


def test_t6_cli_and_plots(tmp_path,base_config):
    from thermofem2d.plot import make_figures
    base_config["geometry"].update(per_side=3,radial=3,element_type="T6")
    sol=solve_model(base_config)
    pics=make_figures(sol,tmp_path)
    assert len(pics)>=7
    assert all(p.stat().st_size>3000 for p in pics)
    assert len(sol.sigma)==len(sol.mesh.elements)


def test_invalid_element_choice():
    with pytest.raises(ValueError,match="element_type"):
        make_mesh(element_type="Q8")
