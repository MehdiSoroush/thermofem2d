"""Independent infinite-plane Kirsch elasticity and finite plate stress tests."""
import numpy as np
import pytest

from thermofem2d.circular_geometry import make_circular_mesh
from thermofem2d.curved_field import circular_field_mesh
from thermofem2d.fem import plane_stress_matrix
from thermofem2d.kirsch import (
    canonical_first_quadrant_angle, hole_edge_stress_samples,
    kirsch_cartesian_stress, kirsch_displacement,
    kirsch_polar_stress, solve_finite_plate_uniform_tension,
    solve_kirsch_outer_dirichlet,
)

R = .25
SIGMA = 10e6


def make(kind='T6', n=4, width=2.):
    return circular_field_mesh(make_circular_mesh(
        width=width, height=width, hole_radius=R,
        per_side=n, radial=n, element_type=kind))


def test_symmetry_equivalent_finite_plate_peak_angles_are_reproducible():
    # A centred square under x tension has reflection symmetry about x and y.
    # Equivalent samples can differ by roundoff between operating systems.
    for nominal in (0., 5.625, 45., 84.375, 90.):
        images=(nominal, -nominal, 180.-nominal,
                180.+nominal, 360.-nominal)
        for theta in images:
            assert canonical_first_quadrant_angle(np.deg2rad(theta)) == nominal
    assert canonical_first_quadrant_angle(np.deg2rad(95.625)) == 84.375
    assert canonical_first_quadrant_angle(np.deg2rad(275.625)) == 84.375
    with pytest.raises(ValueError, match='finite'):
        canonical_first_quadrant_angle(np.nan)


def test_kirsch_hole_traction_and_tangential_extremes():
    angles=np.linspace(0,2*np.pi,65)
    p=np.stack((R*np.cos(angles),R*np.sin(angles)),axis=-1)
    rr,tt,rt=np.moveaxis(kirsch_polar_stress(p,R,SIGMA),-1,0)
    np.testing.assert_allclose(rr,0,atol=2e-8)
    np.testing.assert_allclose(rt,0,atol=2e-8)
    np.testing.assert_allclose(tt,SIGMA*(1-2*np.cos(2*angles)),atol=1e-8)
    np.testing.assert_allclose(tt[16],3*SIGMA,rtol=2e-15)
    np.testing.assert_allclose(tt[0],-SIGMA,rtol=2e-15)


def test_kirsch_cartesian_remote_stress_and_tensor_conversion():
    p=np.array([[60.,0.],[0.,60.],[60.,60.]])
    v=kirsch_cartesian_stress(p,R,SIGMA)
    np.testing.assert_allclose(v,np.broadcast_to([SIGMA,0.,0.],v.shape),rtol=0,atol=600.)
    point=np.array([[.6,.7],[.4,-.8]])
    theta=np.arctan2(point[:,1],point[:,0]);c=np.cos(theta);s=np.sin(theta)
    sig=kirsch_cartesian_stress(point,R,SIGMA)
    polar=kirsch_polar_stress(point,R,SIGMA)
    np.testing.assert_allclose(sig[:,0],polar[:,0]*c*c+polar[:,1]*s*s-2*polar[:,2]*c*s,
                               atol=5e-9)


def test_exact_displacement_strain_matches_independent_stress_formula():
    points=np.array([[.6,.2],[.5,-.4],[-.4,.9]])
    E=210e9; nu=.27;eps=2e-6
    strains=[]
    for p in points:
        gradient=np.column_stack([(kirsch_displacement(p+eps*axis,R,SIGMA,E,nu)-
                                   kirsch_displacement(p-eps*axis,R,SIGMA,E,nu))/(2*eps)
                                  for axis in np.eye(2)])
        strains.append([gradient[0,0],gradient[1,1],gradient[0,1]+gradient[1,0]])
    computed=np.array(strains)@plane_stress_matrix(E,nu).T
    np.testing.assert_allclose(computed,kirsch_cartesian_stress(points,R,SIGMA),
                               rtol=2e-9,atol=.005)


def test_input_validation_and_exclusion_of_void():
    with pytest.raises(ValueError,match='outside'):
        kirsch_polar_stress(np.array([0.,0.]),R,SIGMA)
    with pytest.raises(ValueError,match='outside'):
        kirsch_displacement(np.array([[.1,.1]]),R,SIGMA,200e9,.3)
    with pytest.raises(ValueError,match='Poisson'):
        kirsch_displacement(np.array([[.6,.1]]),R,SIGMA,200e9,.7)
    with pytest.raises(ValueError,match='finite'):
        kirsch_cartesian_stress(np.array([[.5,np.nan]]),R,SIGMA)


def test_boundary_midpoint_extraction_sorts_and_is_one_sided(base_config):
    mat=base_config['material']
    mesh=make('T6',4)
    imposed=(mesh.nodes@np.array([[.001,0.],[-.0005,.0003]]).T).ravel()
    samples=hole_edge_stress_samples(mesh,mat,imposed)
    assert len(samples.theta_rad)==16
    assert np.all(np.diff(samples.theta_rad)>0)
    np.testing.assert_allclose(samples.radius_m,R,atol=1e-14)
    expected=plane_stress_matrix(mat['E'],mat['nu'])@np.array([.001,.0003,-.0005])
    th=samples.theta_rad;c=np.cos(th);s=np.sin(th)
    np.testing.assert_allclose(samples.sigma_tt_pa,
                               expected[0]*s*s+expected[1]*c*c-2*expected[2]*c*s,
                               atol=.003)


def test_kirsch_displacement_boundary_recovery_converges(base_config):
    mat=base_config['material']
    errors=[];peaks=[]
    for n in (4,6,8):
        mesh=make('T6',n)
        u=solve_kirsch_outer_dirichlet(mesh,mat,R,SIGMA)
        exact=kirsch_displacement(mesh.nodes,R,SIGMA,mat['E'],mat['nu']).reshape(-1)
        relative=float(np.linalg.norm(u-exact)/np.linalg.norm(exact))
        s=hole_edge_stress_samples(mesh,mat,u)
        target=1-2*np.cos(2*s.theta_rad)
        rms=float(np.sqrt(np.mean((s.sigma_tt_pa/SIGMA-target)**2)))
        errors.append((relative,rms))
        peaks.append(np.max(s.sigma_tt_pa)/SIGMA)
    assert errors[-1][0] < errors[0][0]
    assert errors[-1][1] < .5*errors[0][1]
    assert abs(peaks[-1]-3) < .06


def test_t3_and_t6_both_improve_kirsch_hole_stress(base_config):
    mat=base_config['material']
    out={}
    for kind in ('T3','T6'):
        levels=[]
        for n in (4,8):
            mesh=make(kind,n)
            u=solve_kirsch_outer_dirichlet(mesh,mat,R,SIGMA)
            s=hole_edge_stress_samples(mesh,mat,u)
            truth=1-2*np.cos(2*s.theta_rad)
            levels.append(np.sqrt(np.mean((s.sigma_tt_pa/SIGMA-truth)**2)))
        out[kind]=levels
        assert levels[1]<levels[0]
    assert out['T6'][-1]<.3*out['T3'][-1]


def test_finite_width_gross_and_net_normalizations(base_config):
    mat=base_config['material']
    mesh=make('T6',8,width=2.)
    u=solve_finite_plate_uniform_tension(mesh,mat,SIGMA)
    s=hole_edge_stress_samples(mesh,mat,u)
    kt_gross=max(s.sigma_tt_pa)/SIGMA
    kt_net=kt_gross*(1-2*R/2.)
    assert 3.4 < kt_gross < 3.8
    assert 2.5 < kt_net < 2.9
    assert kt_net < kt_gross
    # Uniform-traction finite plate is not the exact-Kirsch boundary problem.
    u_exact=solve_kirsch_outer_dirichlet(mesh,mat,R,SIGMA)
    s_exact=hole_edge_stress_samples(mesh,mat,u_exact)
    assert kt_gross - max(s_exact.sigma_tt_pa)/SIGMA>.35


def test_polygonal_t3_stress_samples_on_chord_not_exact_circle(base_config):
    mat=base_config['material']
    mesh=make('T3',4)
    u=solve_kirsch_outer_dirichlet(mesh,mat,R,SIGMA)
    s=hole_edge_stress_samples(mesh,mat,u)
    assert np.all(s.radius_m<R)


def test_uniform_tension_rejects_eccentric_hole(base_config):
    geom=make_circular_mesh(width=2.,height=2.,hole_radius=R,hole_cx=.1,
                            per_side=3,radial=3,element_type='T6')
    with pytest.raises(ValueError,match='centred hole'):
        solve_finite_plate_uniform_tension(circular_field_mesh(geom),
                                           base_config['material'],SIGMA)
