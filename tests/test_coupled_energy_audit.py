"""Mechanical/thermal energy verification for the coupled-wave harness."""
import copy
from dataclasses import replace
import numpy as np
import pytest
from scipy.sparse.linalg import spsolve
from thermofem2d.dynamic_diagnostics import integrate_wave
from thermofem2d.thermal_shock_reference import assemble_graded_wave_system
from thermofem2d.coupled_energy_audit import audit_coupled_energy


def energy_case(base, *, n=3, radial=8):
    cfg=copy.deepcopy(base)
    cfg['geometry'].update(width=2.0, height=2.0, hole_width=0.24,
        hole_height=0.20, hole_cx=0.12, hole_cy=-0.08,
        per_side=n, radial=radial)
    cfg['material'].update(T_ref=298., conductivity=54.,density=7830.,
        heat_capacity=452., thickness=.02)
    cfg['analysis']['displacement_bc']={'ux':{'left':0.},'uy':{'bottom':0.}}
    return cfg


@pytest.mark.parametrize('family',['T3','Q4','T6'])
@pytest.mark.parametrize('thermal_scheme',['consistent','monotone'])
def test_discrete_coupled_energy_closure_and_positivity(base_config,family,thermal_scheme):
    c=energy_case(base_config)
    m=assemble_graded_wave_system(c,family,1.35,thermal_scheme=thermal_scheme)
    r=integrate_wave(m,steps=24,final_time=1.6e-4,profile='cosine',
        rise_time=8e-5,thermal_case='boundary_step',feedback=True)
    a=audit_coupled_energy(m,r,thermal_case='boundary_step')
    assert a.relative_identity_residual < 1e-8
    assert a.conductivity_dissipation_J[-1]>0
    assert a.backward_euler_dissipation_J[-1]>=0
    assert a.max_mechanical_free_residual_N < 1e-3
    assert a.max_thermal_free_residual_W < 1e-2
    # Boundary reactions cannot be omitted from a mixed Dirichlet energy audit.
    assert abs(a.temperature_boundary_work_J[-1]) > 1e-7


@pytest.mark.parametrize('family',['T3','Q4','T6'])
@pytest.mark.parametrize('thermal_scheme',['consistent','monotone'])
def test_adiabatic_affine_strain_cooling_matches_constitutive_limit(base_config,family,thermal_scheme):
    c=energy_case(base_config)
    m=assemble_graded_wave_system(c,family,1.3,thermal_scheme=thermal_scheme)
    eps=1e-5
    imposed_displacement=eps*m.mesh.nodes.reshape(-1)
    temperature_change=-m.reference_temperature*spsolve(
        m.Ct.tocsc(),m.G.T@imposed_displacement)
    mat=c['material']
    beta=mat['E']*mat['alpha']/(1-mat['nu'])
    analytic=-m.reference_temperature*beta*(2*eps)/(
        mat['density']*mat['heat_capacity'])
    assert analytic<0
    assert np.max(abs(temperature_change-analytic))<2e-10


def test_energy_audit_catches_corrupt_thermal_history(base_config):
    c=energy_case(base_config)
    m=assemble_graded_wave_system(c,'T3',1.35,thermal_scheme='monotone')
    r=integrate_wave(m,steps=20,final_time=1.6e-4,profile='cosine',
        thermal_case='boundary_step',feedback=True)
    ok=audit_coupled_energy(m,r,thermal_case='boundary_step')
    assert ok.relative_identity_residual<1e-8
    altered=r.temperature_history.copy()
    interior=np.setdiff1d(np.arange(len(m.mesh.nodes)),np.unique(
        np.r_[m.mesh.node_sets['left'],m.mesh.node_sets['right']]))
    altered[10,int(interior[0])]+=.04
    corrupt=audit_coupled_energy(m,replace(r,temperature_history=altered),
        thermal_case='boundary_step')
    assert corrupt.relative_identity_residual>1e-5


def test_energy_audit_rejects_nonreciprocal_feedback(base_config):
    c=energy_case(base_config)
    m=assemble_graded_wave_system(c,'Q4',1.35)
    r=integrate_wave(m,steps=8,final_time=1.6e-4,profile='cosine',
        thermal_case='reference',feedback=False)
    with pytest.raises(ValueError,match='two-way'):
        audit_coupled_energy(m,r,thermal_case='reference',feedback=False)


def test_coarse_cooling_not_labelled_as_converged_physical_response(base_config):
    c=energy_case(base_config,n=3,radial=12)
    m=assemble_graded_wave_system(c,'T6',1.45,thermal_scheme='monotone')
    results=[]
    for steps in [64,256]:
        r=integrate_wave(m,steps=steps,final_time=1.6e-4,profile='step',
            thermal_case='boundary_step',feedback=True)
        e=audit_coupled_energy(m,r,thermal_case='boundary_step')
        assert e.relative_identity_residual < 1e-8
        results.append(e)
    coarse,fine=results
    # Current fixed-mesh experiment: temporal refinement reduces the large
    # early-time cooling excursion and mixed-integrator coupling mismatch.
    assert fine.min_temperature_K > coarse.min_temperature_K
    assert abs(fine.coupling_time_quadrature_J[-1]) < abs(coarse.coupling_time_quadrature_J[-1])
    assert fine.backward_euler_dissipation_J[-1] < coarse.backward_euler_dissipation_J[-1]
