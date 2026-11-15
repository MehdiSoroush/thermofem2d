"""Discrete-maximum-principle and coupled-interpolation checks."""
import copy
import numpy as np
import pytest
from thermofem2d.thermal_shock_reference import (
    graded_mesh, solve_thermal_shock, shock_error,
    assemble_graded_wave_system, extract_profile,
)
from thermofem2d.monotone_thermal import (
    assemble_monotone_thermal, assemble_subtriangle_thermal_coupling,
)
from thermofem2d.dynamic_diagnostics import integrate_wave, compare_wave
from thermofem2d.coupled import assemble_inertia_coupling
from thermofem2d.fem import assemble_thermal
from thermofem2d.q4_reference import assemble_quad


def case(base, n=4, radial=22):
    cfg=copy.deepcopy(base)
    cfg['geometry'].update(width=2, height=2, hole_width=.24,
        hole_height=.20, hole_cx=.12, hole_cy=-.08,
        per_side=n, radial=radial, element_type='T6')
    cfg['material'].update(conductivity=54.,density=7830.,
        heat_capacity=452., thickness=.02, T_ref=298.)
    cfg['analysis'].update(displacement_bc={'ux':{'left':0},'uy':{'bottom':0}})
    return cfg


@pytest.mark.parametrize('family',['T3','Q4','T6'])
def test_monotone_operator_has_m_matrix_sign_pattern_and_positive_mass(base_config,family):
    c=case(base_config)
    m=graded_mesh(c['geometry'],family,1.45)
    ref=assemble_monotone_thermal(m,c['material'],family)
    K,C=ref.conductivity,ref.capacity
    assert (C.diagonal()>0).all()
    assert np.all(np.isfinite(C.diagonal()))
    assert np.linalg.norm((K-K.T).data)<1e-7
    off=K.tocoo();mask=off.row!=off.col
    assert np.max(off.data[mask],initial=0)<1e-9
    assert abs(np.asarray(K.sum(axis=1)).ravel()).max()<1e-8
    assert ref.added_graph_diffusion_L1>=0
    if family=='Q4':
        _,original,*_=assemble_quad(m,c['material'])
    else:
        _,original=assemble_thermal(m,c['material'])
    assert ref.capacity_total == pytest.approx(float(original.sum()),rel=1e-10)


@pytest.mark.parametrize('family',['T3','Q4','T6'])
def test_temperature_bounds_hold_across_entire_graded_plate(base_config,family):
    c=case(base_config)
    r=solve_thermal_shock(c,family,growth=1.45,steps=64,thermal_scheme='monotone')
    assert r.min_temperature_all_steps_K>=298-2e-9
    assert r.max_temperature_all_steps_K<=310+2e-9
    assert r.temperature.min()>=298-2e-9
    assert r.temperature.max()<=310+2e-9
    assert r.thermal_scheme=='monotone'
    # Conservation includes the Dirichlet boundary reaction and actual
    # nodal heat capacity (not just the free-node equation residual).
    assert abs(r.boundary_energy_input_J-r.stored_energy_change_J)<2e-7
    assert r.stored_energy_change_J>0.


def test_high_resolution_reference_preserves_both_local_accuracy_and_bounds(base_config):
    c=case(base_config,6,28)
    reference=solve_thermal_shock(c,'T6',growth=1.45,steps=128,
                                  thermal_scheme='monotone')
    error=shock_error(reference,c)
    assert error['undershoot_K']==0
    assert error['overshoot_K']==0
    assert error['rms_temperature_error_K']<.05
    assert error['surface_flux_relative_error']<.025
    assert reference.heat_balance_absolute_error_J<1e-7
    depth=np.linspace(0,0.00025,51)
    profile=extract_profile(reference,c['geometry'],depth)
    assert profile.min()>=298-1e-8
    assert profile.max()<=310+1e-8


def test_consistent_scheme_retains_comparison_and_exhibits_baseline_undershoot(base_config):
    c=case(base_config,4,22)
    r=solve_thermal_shock(c,'Q4',growth=1.45,steps=64,thermal_scheme='consistent')
    assert shock_error(r,c)['undershoot_K']>.1


def test_quadratic_temperature_cannot_be_mass_lumped_by_row_sums(base_config):
    c=case(base_config)
    m=graded_mesh(c['geometry'],'T6',1.45)
    _,C=assemble_thermal(m,c['material'])
    vertex_count=(c['geometry']['radial']+1)*4*c['geometry']['per_side']
    # P2 vertex shape functions have zero domain integrals in affine cells.
    assert np.max(np.abs(np.asarray(C.sum(axis=1)).ravel()[:vertex_count]))<1e-6
    monotone=assemble_monotone_thermal(m,c['material'],'T6')
    assert np.min(monotone.capacity.diagonal()[:vertex_count])>0


def test_subtriangle_coupling_preserves_affine_thermal_loads_and_reciprocity(base_config):
    c=case(base_config)
    m=graded_mesh(c['geometry'],'T6',1.45)
    glow=assemble_subtriangle_thermal_coupling(m,c['material'])
    ghigh=assemble_inertia_coupling(m,c['material'])[1]
    assert glow.shape==ghigh.shape
    for temperature in [np.ones(len(m.nodes)),310+9*m.nodes[:,0]-4*m.nodes[:,1]]:
        new,original=glow@temperature,ghigh@temperature
        # Tiny absolute cancellation errors can occur at constrained rows
        # on highly graded triangles; use a global load-vector norm.
        assert np.linalg.norm(new-original)/np.linalg.norm(original)<5e-9
    for component in (0,1):
        rigid=np.tile([1.0,0.0] if component==0 else [0.0,1.0],len(m.nodes))
        assert np.linalg.norm(glow.T@rigid)<1e-6


def test_no_feedback_coupled_thermal_step_obeys_bounds(base_config):
    c=case(base_config)
    system=assemble_graded_wave_system(c,'T6',1.45,thermal_scheme='monotone')
    r=integrate_wave(system,steps=32,final_time=1.6e-4,profile='step',
        thermal_case='boundary_step',feedback=False)
    assert np.min(r.temperature_history)>=298-1e-8
    assert np.max(r.temperature_history)<=310+1e-8


def test_coupled_reference_agreement_not_artificially_forced(base_config):
    c=case(base_config)
    q=integrate_wave(assemble_graded_wave_system(c,'Q4',1.45,thermal_scheme='monotone'),
        steps=32,final_time=1.6e-4,profile='step',thermal_case='boundary_step',feedback=True)
    t=integrate_wave(assemble_graded_wave_system(c,'T6',1.45,thermal_scheme='monotone'),
        steps=32,final_time=1.6e-4,profile='step',thermal_case='boundary_step',feedback=True)
    relative=compare_wave(t,q)['relative_difference'][-1]
    assert .01<relative<.1
    # With thermoelastic feedback, a strict heat-equation maximum principle
    # does not apply: deformation-rate coupling can create real cooling.
    assert np.min(t.temperature_history)<298


def test_invalid_thermal_option_rejected(base_config):
    c=case(base_config)
    with pytest.raises(ValueError,match='thermal_scheme'):
        solve_thermal_shock(c,'T3',growth=1.45,steps=1,thermal_scheme='ignored')
    with pytest.raises(ValueError,match='P2/P1'):
        assemble_subtriangle_thermal_coupling(graded_mesh(c['geometry'],'T3',1.45),c['material'])
