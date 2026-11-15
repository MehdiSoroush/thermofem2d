# ThermoFEM-2D v0.7.0 — Experimental Thermal-Shock Verification

**Mehdi Soroush · November 2026**

Version 0.7.0 extends the verification framework for coupled, small-strain plane-stress thermoelastic finite-element analyses of rectangular perforated plates. The established T3/T6 solver, independent Q4 spatial reference, and ordinary application interface are retained. The new methods are implemented as reproducible research benchmarks rather than application-selectable engineering-design solvers.

## New verification studies

- **Dynamic-wave response:** abrupt and smooth mechanical loading, Q4–T6 displacement comparisons, mesh refinement, selected modal frequencies, and controlled thermal initialization.
- **Resolved near-edge thermal shock:** geometrically graded meshes for T3, Q4 and T6, compared with the analytical semi-infinite transient-conduction temperature profile and surface heat flux.
- **Bound-preserving conduction reference:** positive diagonal capacity and graph diffusion for source-free thermal conduction. On six-node triangular meshes, the thermal field is represented using four linear subtriangles; it is **not** a quadratic T6 monotone method.
- **Two-way thermoelastic energy accounting:** postprocessing of mechanical kinetic and elastic energy, quadratic thermal perturbation energy, thermal boundary work, conductive dissipation, backward-Euler damping, and the finite-step coupling correction.

## Key findings and interpretation

The earlier 63–67% early-time Q4–T6 displacement disagreement is strongly associated with an unresolved thermal boundary layer caused by an abrupt thermal boundary step. In tested graded-mesh cases, the remaining relative cross-discretization displacement differences are approximately 3.5–4.6%. This comparison does **not** provide an absolute continuum-error estimate.

The monotone conduction reference eliminates unphysical global temperature undershoots for the tested *uncoupled, source-free conduction* cases with 298–310 K temperature bounds. The original high-order consistent thermal discretization may still exhibit undershoots on highly graded meshes. In two-way coupled thermoelasticity, deformation-rate heat exchange may induce physical cooling or heating; observed short-time minima remain sensitive to the time step and spatial discretization and have **not** been independently validated.

The two-way energy audit closes the **discrete** perturbation-energy ledger. This demonstrates internal consistency with the assembled equations; it is neither a proof of continuum accuracy nor a validation of thermal-shock design stresses.

## Reproduction and QA

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python examples/scripts/check_release_integrity.py
```

Optional higher-cost studies:

```bash
python examples/scripts/validate_dynamic_waves.py
python examples/scripts/verify_thermal_shock_boundary_layer.py
python examples/scripts/verify_monotone_thermal.py
python examples/scripts/verify_coupled_energy_audit.py
```

The scripts write traceable tables and figures to `results/verification/`. The regular automated suite includes smaller representative numerical regressions; the longer time-refinement runs are preserved in the research scripts and datasets.

## Known limitations

- No independently established full-domain, spatially and temporally converged coupled thermal-shock reference.
- High-frequency wave and thermal response remain sensitive to spatial discretization and time integration.
- Stress concentrations at sharp hole corners are not mesh-independent design stresses.
- No experimental validation, nonlinear material model, curved-hole geometry or complete reconstruction of the 2016 MATLAB numerical results is claimed.
- The new graded-mesh and monotone thermal reference methods are verification utilities, **not** general production solver or Streamlit options.

For the evidence matrix, see [verification status](VALIDATION_STATUS.md); for the technical analyses, see [dynamic-wave study](DYNAMIC_WAVE_VALIDATION.md), [thermal boundary-layer study](THERMAL_SHOCK_BOUNDARY_LAYER.md), [monotone thermal reference](POSITIVITY_AND_THERMAL_CONSERVATION.md) and [coupled energy audit](COUPLED_ENERGY_BUDGET.md).

**Intended use:** numerical research, code verification and reproducible computational mechanics studies; not certified thermal-shock or structural-integrity engineering design.
