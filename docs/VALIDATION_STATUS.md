# Verification Evidence and Unresolved Limitations

**ThermoFEM-2D v0.7.0** · Mehdi Soroush · November 2026

This repository contains research software for small-strain linear plane-stress thermoelastic finite-element calculations. The ordinary application supports straight-sided T3/T6 triangular meshes, steady/transient thermal conduction, and one-way or two-way coupling. The independent Q4 reference and v0.7 graded-mesh studies are separate verification modules, **not** interchangeable production solver options.

## Evidence classification

| Numerical component or question | Verification evidence | Limitation |
| --- | --- | --- |
| Thermal and mechanical element assembly | Analytical affine and polynomial patch tests, independent manufactured fields | Limited to implemented materials, affine geometry and documented boundary data |
| Two-way strain-rate coupling and inertia | Continuum-manufactured T3/Q4/T6 solutions with nonzero acceleration, algebraic reciprocity and adiabatic cooling sign tests | Not independent validation for arbitrary transient loads |
| Mechanical dynamic response | Q4–T6 cross-comparisons under abrupt/smooth loading and first-mode frequency sensitivity | High-frequency wave dispersion and arbitrary-geometry convergence remain open |
| Sudden thermal boundary layer | Semi-infinite conduction temperature and flux reference; graded-mesh refinement | Half-space benchmark measures only the near-edge behavior; global 2D coupled field is not analytically verified |
| Heat-conduction temperature bounds | Diagonal thermal capacity and symmetric graph diffusion with recorded global pure-conduction temperature bounds | Monotone reference adds diffusion; T6-node approach uses **piecewise-linear subtriangles**, not quadratic thermal interpolation |
| Coupled energy budget | Exact stepwise discrete quadratic perturbation energy identity and a negative-control test | Identity closure is algebraic consistency, **not** physical accuracy or continuum thermodynamic validation |
| Apparent coupled cooling | Time-step and thermal-scheme sensitivity studies | Physical cooling cannot be distinguished from integration/discretization effects solely from existing outputs |
| Historical MATLAB coursework | Source/provenance audit of the 2016 report | Original numerical results/matrices have not been independently reproduced |

## Interpretation of v0.7 results

The coarse abruptly loaded Q4–T6 displacement discrepancy (about 63–67% in earlier tests) decreases substantially with thermal-layer refinement, to roughly 3.5–4.6% in selected graded-mesh comparisons. This is a **difference between two numerical discretizations**, not an estimate of the error against an independently established continuum solution.

Uncoupled source-free pure conduction admits temperature bounds; reciprocal two-way thermoelasticity includes deformation-rate heat exchange and can generate temperature variations outside the thermal-only bounds. The remaining short-time coupled temperature minima are sensitive to time step, spatial resolution and interpolation choice. **Do not use these transients as validated thermal-shock design loads or material failure predictions.**

## Reproduction

Run the automated test suite first:

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
```

The high-resolution verification studies are available separately:

```bash
python examples/scripts/validate_dynamic_waves.py
python examples/scripts/verify_thermal_shock_boundary_layer.py
python examples/scripts/verify_monotone_thermal.py
python examples/scripts/verify_coupled_energy_audit.py
```

These scripts may perform large transient solves. Numerical data and plots are written to `results/verification/`. Output metadata explicitly labels the studies as experimental and does not identify results as experiments or historical recovered data.

## Further validation

1. Independent full-domain mesh/time convergence under two-way coupling, including temperature, displacement, strain, boundary heat flux and relevant energy terms.
2. External FEM cross-check of the same loading, initial conditions and spatially resolved thermal boundary layer.
3. A consistently time-centered coupled integrator and comparison with the present Newmark/backward-Euler equations.
4. Careful treatment of corner stress singularities and anisotropic mesh-element conditioning before any design applications.

These open items are scientific development objectives, not evidence that the current internally verified numerical routines are known to be invalid.
