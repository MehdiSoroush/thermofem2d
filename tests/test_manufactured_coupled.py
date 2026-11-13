"""Continuum-forced, exact-solution checks of coupled transient T3/T6/Q4."""
import numpy as np
import pytest

from thermofem2d.manufactured import (
    continuum_sources, exact_fields, exact_initial_velocity, verify_manufactured_coupled,
)


def test_manufactured_fields_and_forces_are_nontrivial(base_config):
    mat=base_config['material']
    points=np.array([[.2,-.3],[-.5,.4],[0.,0.]])
    f,q=continuum_sources(points,0.0007,mat)
    u,T=exact_fields(points,0.0007,t_ref=mat['T_ref'])
    v=exact_initial_velocity(points)
    assert np.all(np.isfinite(f)) and np.all(np.isfinite(q))
    assert np.linalg.norm(f)>0 and np.linalg.norm(q)>0
    assert np.linalg.norm(u)>0 and np.linalg.norm(v)>0
    assert np.max(np.abs(T-mat['T_ref']))>0
    # Exact continuum thermal source has a nonzero feedback term proportional
    # to x; evaluate it independently from the no-feedback part.
    x,y=points.T
    b=6000.;a=.02;t=.0007
    q_without=mat['density']*mat['heat_capacity']*b*(x*x+2*y*y)-6*mat['conductivity']*b*t
    expected_feedback=mat['T_ref'] * mat['E']*mat['alpha']/(1-mat['nu'])*3*a*x
    np.testing.assert_allclose(q-q_without,expected_feedback,rtol=1e-11,atol=1e-5)


@pytest.mark.parametrize('n',[3,5])
@pytest.mark.parametrize('steps',[4,8])
def test_t6_exact_quadratic_coupled_manufactured(base_config,n,steps):
    r=verify_manufactured_coupled(base_config,'T6',n,steps=steps)
    assert r.temperature_max_error_K < 2e-9
    assert r.displacement_max_error_m < 1e-12
    assert r.temperature_rms_error_K < 1e-9
    assert r.displacement_rms_error_m < 1e-12


@pytest.mark.parametrize('kind',['T3','Q4'])
def test_lower_order_elements_converge_to_continuum_exact_field(base_config,kind):
    coarse=verify_manufactured_coupled(base_config,kind,3,steps=8)
    fine=verify_manufactured_coupled(base_config,kind,6,steps=8)
    assert fine.temperature_rms_error_K < 0.35*coarse.temperature_rms_error_K
    assert fine.displacement_rms_error_m < 0.25*coarse.displacement_rms_error_m
    assert fine.temperature_rms_error_K > 1e-4


@pytest.mark.parametrize('steps',[4,8])
def test_t6_exact_with_nonzero_inertia_and_acceleration(base_config,steps):
    # u ~ t + curvature*t², thus u_ddot is nonzero throughout the body.
    r=verify_manufactured_coupled(base_config,'T6',4,steps=steps,curvature=5000.0)
    assert r.temperature_max_error_K < 2e-9
    assert r.displacement_max_error_m < 1e-12


def test_coupling_feedback_matters_for_continuum_benchmark(base_config):
    correct=verify_manufactured_coupled(base_config,'T6',4,steps=8)
    negative=verify_manufactured_coupled(base_config,'T6',4,steps=8,disable_feedback=True)
    assert correct.temperature_max_error_K < 1e-9
    assert negative.temperature_rms_error_K > 0.001
    assert negative.displacement_rms_error_m > 1e-9


def test_invalid_manufactured_parameters_rejected(base_config):
    with pytest.raises(ValueError):
        verify_manufactured_coupled(base_config,'Q9')
    with pytest.raises(ValueError):
        verify_manufactured_coupled(base_config,'T6',3,steps=0)
    with pytest.raises(ValueError):
        verify_manufactured_coupled(base_config,'T6',2)
    with pytest.raises(ValueError):
        verify_manufactured_coupled(base_config,'T6',3,end_time=-1)
