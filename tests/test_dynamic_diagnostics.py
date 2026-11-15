"""Controlled dynamic-wave benchmark and regression tests."""
from dataclasses import replace
import copy

import numpy as np
import pytest

from thermofem2d.dynamic_diagnostics import (
    assemble_wave_system, compare_wave, integrate_wave,
    lowest_frequencies, ramp_amplitude, rms_vector,
)
from thermofem2d.q4_reference import solve_q4_coupled
from thermofem2d.solve import solve_model


def demo(base_config, n=4):
    case = copy.deepcopy(base_config)
    case['geometry'].update(width=2.,height=2.,hole_width=.24,hole_height=.20,
        hole_cx=.12,hole_cy=-.08,per_side=n,radial=n,element_type='T6')
    case['material'].update(E=2e11,nu=.3,alpha=11.7e-6,conductivity=54.,
         density=7830.,heat_capacity=452.,thickness=.02,T_ref=298.)
    case['analysis'].update(
         thermoelastic_mode='coupled_dynamic',thermal_mode='transient',
         thermoelastic_feedback=True,time_step=1.6e-4/16,steps=16,
         initial_temperature=298.,temperature_bc={'left':310.,'right':298.},
         inward_heat_flux={},volumetric_heat_source=0.,
         displacement_bc={'ux':{'left':0.},'uy':{'bottom':0.}},
         traction_bc={'top':[0.,2e6]})
    return case


def test_ramp_is_c1_continuous_at_start_and_finish():
    rise=1.6e-4
    assert ramp_amplitude(0,rise,'cosine') == 0
    assert ramp_amplitude(rise,rise,'cosine') == 1
    assert ramp_amplitude(2*rise,rise,'cosine') == 1
    assert ramp_amplitude(.5*rise,rise,'cosine') == pytest.approx(.5)
    for t in np.linspace(0,rise,20):
        assert 0 <= ramp_amplitude(t,rise,'cosine') <= 1
    assert ramp_amplitude(0,rise,'step') == 1


@pytest.mark.parametrize('method', ['T3','T6','Q4'])
def test_no_force_at_reference_temperature_stays_at_rest(base_config,method):
    case=demo(base_config)
    s=assemble_wave_system(case,method)
    s=replace(s,traction=np.zeros_like(s.traction))
    out=integrate_wave(s,steps=4,final_time=.00016,
                       profile='cosine',thermal_case='reference',feedback=True)
    assert np.max(abs(out.displacement_history)) < 1e-15
    assert np.max(abs(out.temperature_history-298.)) < 1e-12


@pytest.mark.parametrize('method', ['T6','Q4'])
def test_preheated_undriven_equilibrium_is_stationary(base_config,method):
    s=assemble_wave_system(demo(base_config),method)
    s=replace(s,traction=np.zeros_like(s.traction))
    r=integrate_wave(s,steps=8,final_time=.00016,
                     profile='cosine',thermal_case='preheated_equilibrium',feedback=True)
    assert np.max(abs(r.displacement_history[-1]-r.displacement_history[0])) < 2e-12
    assert np.max(abs(r.temperature_history[-1]-r.temperature_history[0])) < 1e-9


def test_thermal_step_reproduces_production_t6(base_config):
    case=demo(base_config)
    out=solve_model(case)
    ref=integrate_wave(assemble_wave_system(case,'T6'),steps=16,
                       final_time=.00016,profile='step',
                       thermal_case='boundary_step',feedback=True)
    assert np.allclose(out.displacement_history.reshape(ref.displacement_history.shape),
                       ref.displacement_history,rtol=1e-10,atol=1e-13)
    assert np.allclose(out.temperature_history,ref.temperature_history,rtol=1e-12,atol=1e-11)


def test_thermal_step_reproduces_independent_q4(base_config):
    case=demo(base_config)
    out=solve_q4_coupled(case)
    ref=integrate_wave(assemble_wave_system(case,'Q4'),steps=16,
                       final_time=.00016,profile='step',
                       thermal_case='boundary_step',feedback=True)
    assert np.allclose(out.displacement_history.reshape(ref.displacement_history.shape),
                       ref.displacement_history,rtol=1e-10,atol=1e-13)
    assert np.allclose(out.temperature_history,ref.temperature_history,rtol=1e-12,atol=1e-11)


def test_initial_thermal_step_causes_displacement_discrepancy(base_config):
    case=demo(base_config)
    q=assemble_wave_system(case,'Q4')
    t=assemble_wave_system(case,'T6')
    ref=compare_wave(
        integrate_wave(t,steps=32,final_time=.00016,profile='step',
                       thermal_case='reference',feedback=True),
        integrate_wave(q,steps=32,final_time=.00016,profile='step',
                       thermal_case='reference',feedback=True))
    shocked=compare_wave(
        integrate_wave(t,steps=32,final_time=.00016,profile='step',
                       thermal_case='boundary_step',feedback=True),
        integrate_wave(q,steps=32,final_time=.00016,profile='step',
                       thermal_case='boundary_step',feedback=True))
    assert shocked['relative_difference'][-1] > 3*ref['relative_difference'][-1]


def test_smooth_ramp_spatial_improvement(base_config):
    outcomes=[]
    for n in (4,8):
        case=demo(base_config,n)
        q=integrate_wave(assemble_wave_system(case,'Q4'),steps=32,
              final_time=.00016,profile='cosine',rise_time=.00016)
        t=integrate_wave(assemble_wave_system(case,'T6'),steps=32,
              final_time=.00016,profile='cosine',rise_time=.00016)
        outcomes.append(compare_wave(t,q)['relative_difference'][-1])
    assert outcomes[1] < outcomes[0]


def test_first_modes_are_positive_and_approach_each_other(base_config):
    errors=[]
    for n in (4,8):
        case=demo(base_config,n)
        fq=lowest_frequencies(assemble_wave_system(case,'Q4'),3)
        ft=lowest_frequencies(assemble_wave_system(case,'T6'),3)
        assert np.all(np.diff(fq)>0) and np.all(np.diff(ft)>0)
        errors.append(np.max(abs(fq-ft)/ft))
    assert errors[1] < errors[0]


def test_time_refinement_of_smooth_ramp(base_config):
    s=assemble_wave_system(demo(base_config),'T6')
    a=[integrate_wave(s,steps=n,final_time=.00016,
                      profile='cosine',rise_time=.00016).displacement_history[-1]
       for n in (8,16,32)]
    err1=rms_vector(a[1]-a[0]);err2=rms_vector(a[2]-a[1])
    assert 0 < err2 < .5*err1


def test_invalid_wave_options(base_config):
    s=assemble_wave_system(demo(base_config),'Q4')
    with pytest.raises(ValueError):
        integrate_wave(s,steps=0,final_time=.00016,profile='step')
    with pytest.raises(ValueError):
        integrate_wave(s,steps=4,final_time=.00016,profile='step',thermal_case='unknown')
    with pytest.raises(ValueError):
        ramp_amplitude(.1,0,'cosine')
