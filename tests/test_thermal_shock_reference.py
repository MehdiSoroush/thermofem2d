"""Numerical verification tests: boundary-layer heat transfer and mesh grading."""
import copy
import numpy as np
import pytest
from scipy.spatial import cKDTree
from thermofem2d.thermal_shock_reference import (
    thermal_diffusivity, penetration_length, halfspace_temperature,
    halfspace_surface_flux, radial_fractions, graded_mesh,
    solve_thermal_shock, shock_error, extract_profile,
    assemble_graded_wave_system, quadratic_surface_gradient,
)
from thermofem2d.dynamic_diagnostics import integrate_wave, compare_wave


def case(base, per_side=4, radial=22):
    cfg=copy.deepcopy(base)
    cfg['geometry'].update(width=2, height=2, hole_width=.24,
        hole_height=.20, hole_cx=.12, hole_cy=-.08,
        per_side=per_side, radial=radial, element_type='T6')
    cfg['material'].update(conductivity=54., density=7830.,
        heat_capacity=452., thickness=.02, T_ref=298.)
    cfg['analysis'].update(displacement_bc={'ux':{'left':0}, 'uy':{'bottom':0}})
    return cfg


def test_semi_infinite_solution_scales(base_config):
    cfg=case(base_config)
    m=cfg['material'];alpha=54/(7830*452)
    assert thermal_diffusivity(m)==pytest.approx(alpha)
    assert penetration_length(m,1.6e-4,.01)==pytest.approx(.000179986051757,rel=1e-8)
    assert halfspace_temperature(0,1.6e-4,m)==pytest.approx(310)
    depth=penetration_length(m,1.6e-4,.01)
    assert halfspace_temperature(depth,1.6e-4,m)==pytest.approx(298.12)
    assert halfspace_surface_flux(1.6e-4,m)==pytest.approx(7399345.33427,rel=1e-6)


def test_quadratic_surface_gradient_exactness():
    # T(d) = 310 - 2000*d + 400000*d^2; derivative at d=0 is -2000 K/m.
    h=2.5e-4
    fun=lambda x:310.-2000.*x+4e5*x*x
    got=quadratic_surface_gradient(fun(0),fun(h/2),fun(h),h)
    assert got == pytest.approx(-2000.,rel=1e-11,abs=1e-9)
    with pytest.raises(ValueError):
        quadratic_surface_gradient(310.,309.,308.,0.)


def test_reference_input_validation(base_config):
    m=case(base_config)['material']
    with pytest.raises(ValueError):halfspace_temperature(-1.,1.,m)
    with pytest.raises(ValueError):halfspace_temperature(0.,0.,m)
    with pytest.raises(ValueError):penetration_length(m,1.,0.)
    with pytest.raises(ValueError):radial_fractions(1,1.5)
    with pytest.raises(ValueError):radial_fractions(6,.9)


def test_geometric_layers_are_monotone_and_locally_refined():
    linear=radial_fractions(12,1.)
    assert np.allclose(linear,np.linspace(0,1,13))
    values=radial_fractions(28,1.45)
    assert values[0]==pytest.approx(0)
    assert values[-1]==pytest.approx(1)
    assert np.all(np.diff(values)>0)
    assert 1-values[-2]<2e-5


@pytest.mark.parametrize('family',['T3','Q4','T6'])
def test_ungraded_geometry_matches_original(base_config,family):
    from thermofem2d.geometry import make_mesh
    from thermofem2d.q4_reference import quad_mesh
    cfg=case(base_config)
    g=cfg['geometry'];orig=(quad_mesh(g) if family=='Q4' else
                             make_mesh(**dict(g,element_type=family)))
    new=graded_mesh(g,family,1.)
    assert np.allclose(orig.nodes,new.nodes,atol=1e-13,rtol=0)


def test_independent_q4_and_t6_graded_vertex_positions_agree(base_config):
    g=case(base_config,6,28)['geometry']
    q=graded_mesh(g,'Q4',1.45);t=graded_mesh(g,'T6',1.45)
    err,_=cKDTree(t.nodes).query(q.nodes)
    assert np.max(err)<1e-12
    assert len(t.nodes)>len(q.nodes)


@pytest.mark.parametrize('family', ['T3','Q4','T6'])
def test_grading_reduces_analytic_thermal_error(base_config,family):
    cfg=case(base_config,4,22)
    uni=solve_thermal_shock(cfg,family,growth=1.,steps=48)
    graded=solve_thermal_shock(cfg,family,growth=1.45,steps=48)
    e0=shock_error(uni,cfg)['rms_temperature_error_K']
    e1=shock_error(graded,cfg)['rms_temperature_error_K']
    assert e0>7.
    assert e1<.25*e0


def test_high_resolution_quadratic_profile_and_flux(base_config):
    cfg=case(base_config,6,28)
    res=solve_thermal_shock(cfg,'T6',growth=1.45,steps=128)
    e=shock_error(res,cfg)
    assert e['first_radial_spacing_m']<2e-5
    assert e['rms_temperature_error_K']<.035
    assert e['surface_flux_relative_error']<.025
    y=extract_profile(res,cfg['geometry'],[0,5e-5,1.8e-4])
    assert y[0]==pytest.approx(310)
    assert np.all(np.diff(y)<0)


def test_local_success_does_not_hide_global_undershoots(base_config):
    cfg=case(base_config,4,22)
    res=solve_thermal_shock(cfg,'Q4',growth=1.45,steps=64)
    e=shock_error(res,cfg)
    assert e['undershoot_K']==pytest.approx(max(0., 298.-res.temperature.min()))
    assert np.isfinite(res.temperature).all()
    assert 298.<=res.max_temperature_K<=310.00000001


def test_transient_coupled_discrepancy_improves_with_thermal_grading(base_config):
    cfg=case(base_config,4,22)
    ratios=[]
    for growth in (1.,1.45):
        q=integrate_wave(assemble_graded_wave_system(cfg,'Q4',growth),steps=32,
            final_time=1.6e-4,profile='step',thermal_case='boundary_step',feedback=True)
        t=integrate_wave(assemble_graded_wave_system(cfg,'T6',growth),steps=32,
            final_time=1.6e-4,profile='step',thermal_case='boundary_step',feedback=True)
        ratios.append(compare_wave(t,q)['relative_difference'][-1])
    assert ratios[0]>.3
    assert ratios[1]<.08
    assert ratios[1]<ratios[0]/4


def test_invalid_capacity_and_time_settings(base_config):
    cfg=case(base_config)
    with pytest.raises(ValueError):solve_thermal_shock(cfg,'Q4',steps=0)
    with pytest.raises(ValueError):graded_mesh(cfg['geometry'],'Q8',1.45)
