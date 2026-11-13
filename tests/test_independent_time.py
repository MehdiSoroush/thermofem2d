"""v0.4: independent-in-time method-of-lines and thermal energy checks."""
from copy import deepcopy
import numpy as np
import pytest

from thermofem2d.solve import solve_model, initial_temperature_field
from thermofem2d.geometry import make_mesh
from thermofem2d.verification import solve_semidiscrete_ode, insulated_energy_audit


def benchmark(base, element='T3', *, n=8, thermal_gradient=True):
    c = deepcopy(base)
    c['geometry'].update(element_type=element, per_side=2, radial=2)
    c['analysis'].update(
        thermoelastic_mode='coupled_dynamic', thermal_mode='transient',
        thermoelastic_feedback=True, time_step=1e-4/n, steps=n,
        initial_temperature=({'type':'affine','base':310.,'gx':6.,'gy':-3.}
                             if thermal_gradient else 310.),
        temperature_bc={}, inward_heat_flux={}, volumetric_heat_source=0.,
        traction_bc={})
    return c


def test_initial_affine_field_and_reject_invalid_specifications(base_config):
    mesh=make_mesh(per_side=2,radial=2)
    a={'initial_temperature':{'type':'affine','base':305.,'gx':4.,'gy':-2.}}
    T=initial_temperature_field(mesh,base_config['material'],a)
    assert np.allclose(T,305 + 4*mesh.nodes[:,0] - 2*mesh.nodes[:,1])
    assert not np.allclose(T,T[0]); assert np.all(np.isfinite(T))
    assert np.allclose(initial_temperature_field(mesh,base_config['material'],{}),298.)
    for spec in [float('inf'),{'type':'affine','base':300.,'gx':1.},
                 {'type':'affine','base':300.,'gx':1.,'gy':float('nan')},
                 {'type':'quadratic','base':300.,'gx':1.,'gy':2.}]:
        with pytest.raises(ValueError):
            initial_temperature_field(mesh,base_config['material'],{'initial_temperature':spec})


@pytest.mark.parametrize('element', ['T3','T6'])
def test_scipy_ode_reference_has_small_semidiscrete_residual(base_config,element):
    case=benchmark(base_config,element,n=12)
    solution=solve_model(case)
    ref=solve_semidiscrete_ode(case,solution.times)
    assert ref.max_relative_residual < 1e-9
    assert np.all(np.isfinite(ref.temperature_history))
    assert np.all(np.isfinite(ref.displacement_history))
    assert np.all(np.isfinite(ref.velocity_history))
    assert np.allclose(solution.temperature_history[0],ref.temperature_history[0])
    assert np.allclose(solution.displacement_history[0],ref.displacement_history[0])
    assert np.allclose(solution.velocity_history[0],ref.velocity_history[0])
    # Modest step: within a few percent of the dynamic displacement scale,
    # and below 0.05 K in temperature for this small manufactured scenario.
    err_u=np.linalg.norm(solution.displacement-ref.displacement_history[-1])
    scale=max(np.linalg.norm(ref.displacement_history[-1]),1e-12)
    assert err_u/scale < .1
    assert np.max(abs(solution.temperature-ref.temperature_history[-1])) < .06


def test_time_refinement_reduces_independent_reference_error(base_config):
    errors=[]
    for n in (8,16,32):
        case=benchmark(base_config,n=n)
        calc=solve_model(case)
        ref=solve_semidiscrete_ode(case,calc.times)
        errors.append((np.max(abs(calc.temperature-ref.temperature_history[-1])),
                       np.max(abs(calc.displacement-ref.displacement_history[-1]))))
    assert errors[2][0]<.65*errors[1][0]<.65**2*errors[0][0]
    assert errors[2][1]<.65*errors[1][1]<.65**2*errors[0][1]


def test_energy_balance_converges_with_time_step(base_config):
    errors=[]
    for n in (4,8,16):
        case=benchmark(base_config,n=n)
        calc=solve_model(case)
        audit=insulated_energy_audit(case,calc)
        E=audit.mechanical_energy+audit.thermal_energy
        assert np.all(np.isfinite(E))
        assert np.min(audit.conduction_power)>-1e-8
        assert audit.cumulative_conduction[-1]>0
        assert E[-1] <= E[0]*(1+1e-8)
        errors.append(abs(audit.energy_balance_error[-1])/E[0])
    assert errors[2]<.65*errors[1]<.65**2*errors[0]


def test_energy_audit_rejects_external_work_or_heat(base_config):
    case=benchmark(base_config,n=5)
    calc=solve_model(case)
    for change in ({'traction_bc':{'top':[0.,1000.]}},
                   {'temperature_bc':{'left':298.}},
                   {'volumetric_heat_source':1.0}):
        bad=deepcopy(case)
        bad['analysis'].update(change)
        with pytest.raises(ValueError):
            insulated_energy_audit(bad,calc)


def test_ode_reference_with_thermoelastic_feedback_disabled(base_config):
    case=benchmark(base_config,n=8)
    case['analysis'].update(thermoelastic_feedback=False,initial_temperature=298.,
                            traction_bc={'top':[0.,1e6]})
    calc=solve_model(case)
    ref=solve_semidiscrete_ode(case,calc.times)
    assert np.max(abs(calc.temperature-298))<1e-10
    assert np.max(abs(ref.temperature_history[-1]-298))<1e-10
    assert np.max(abs(calc.displacement-ref.displacement_history[-1]))<2e-5


def test_coupled_rejects_nonphysical_reference_temperature(base_config):
    case=benchmark(base_config,n=4)
    case['material']['T_ref']=0.0
    with pytest.raises(ValueError,match='positive absolute T_ref'):
        solve_model(case)
    with pytest.raises(ValueError,match='T_ref > 0'):
        solve_semidiscrete_ode(case,np.array([0.,1e-4]))
