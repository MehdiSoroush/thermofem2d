# ThermoFEM-2D

[![Numerical verification](https://github.com/MehdiSoroush/thermofem2d/actions/workflows/ci.yml/badge.svg)](https://github.com/MehdiSoroush/thermofem2d/actions/workflows/ci.yml)

**Parametric finite element analysis of a perforated plate under thermal loading**

> **v0.1.1 — documented numerical prototype, November 2026.** The current software solves steady or transient **heat conduction** and then a **one-way, quasi-static plane-stress thermoelastic problem**. It is **not yet** a fully coupled thermoelastic solver, and it is **not** an exact numerical reproduction of the 2016 coursework results.

![Six-panel computational overview: temperature, displacement, element-centroid von Mises stress, Fourier heat flux, T3 mesh and displayed deformed geometry](results/figures/overview.png)

## Overview

ThermoFEM-2D is a small, auditable Python finite-element research application. The source study was a **summer 1395 (2016) MATLAB coursework project**, concerned with classical coupled thermoelasticity of a square perforated plate using six-node quadratic triangular elements and three nodal unknowns (`u_x`, `u_y`, `theta`). The archival Persian-language manuscript is not distributed here.

**The 2026 v0.1 solver is a new, deliberately narrower baseline**, written in Python and verified using analytical limiting cases. The historical T6 three-field formulation and its missing MATLAB routines are reserved for reconstruction in later releases. See [source audit](docs/SOURCE_AUDIT.md) and [formulation](docs/FORMULATION.md).

## Current capabilities

- Rectangular plate with one **movable and resizable rectangular through-hole** (geometry-fitted mesh; all coordinates in metres).
- Three-node linear triangular **T3** finite elements, with sparse global assembly and correct element orientation.
- **Steady-state** isotropic Fourier heat conduction and **transient** conduction using a consistent heat capacity matrix and implicit backward Euler.
- Prescribed temperatures on named boundaries; constant imposed inward heat flux and uniform volumetric heat sources.
- Linear, small-strain **plane-stress** elasticity with thermal-expansion loads and optional constant mechanical tractions.
- Separately prescribed `ux` and `uy` on named boundaries with incompatible thermal corner values rejected.
- Temperature, displacement magnitude, and **element-centroid** von Mises stress contour plots; **elementwise Fourier heat-flux maps**, deformed-mesh visualization with a stated magnification factor, mesh inspection, CSV + JSON results.
- Optional interactive Streamlit interface (geometry sliders and interactive parameter configuration).
- Reproducible synthetic verification cases, exact manufactured heat and plane-stress benchmarks, a quantitative conduction-convergence example, mesh sensitivity, and GitHub Actions tests.

## Quick start

Requires Python 3.11+.

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install -e ".[dev]"
python -m pytest -q
python -m thermofem2d.cli --config examples/steady_plate.json --out results
python -m thermofem2d.cli --config examples/offset_hole_transient.json --out results-offset
python examples/scripts/mesh_study.py
python examples/scripts/verify_manufactured.py
```

Optional graphical interface:

```bash
python -m pip install -e ".[app]"
streamlit run app.py
```

**The Streamlit UI is an optional convenience, not the verification baseline.** The command-line solver and tests do not depend on Streamlit. The Streamlit UI smoke test passed in GitHub Actions on Ubuntu/Python 3.12. A separate live deployment has not yet been assessed.

## Change the geometry and boundary conditions

All inputs are in JSON files in `examples/`:

```json
{
  "geometry": {"width": 2.0, "height": 2.0,
               "hole_width": 0.24, "hole_height": 0.20,
               "hole_cx": 0.0, "hole_cy": 0.0,
               "per_side": 16, "radial": 14}
}
```

Move the hole by changing `hole_cx` and `hole_cy`, resize it with `hole_width` and `hole_height`, and refine the T3 mesh with `per_side` and `radial`. The hole must remain strictly inside the plate. **No geometrical edits to the FEM kernel are needed.**

Boundary identifiers are `left`, `right`, `top`, `bottom`, `hole`. Prescribed thermal values are given under `temperature_bc`, and mechanical displacement components under `displacement_bc` (`ux` and `uy`). Any unconstrained thermal boundary has zero natural flux; unconstrained mechanical boundaries are traction-free. The sample fixes horizontal displacements on the left edge and vertical displacements on the bottom edge to remove rigid-body modes. These conditions are **a new verification setup**, not a claim about the 2016 physical boundary conditions.

**Units:** m, kg, s, K, Pa, W, J; temperatures are absolute Kelvin and thermal strains use `T - T_ref`. Plots display displacement in mm and stress in MPa for readability. The example geometric scale and temperature loading are illustrative; material coefficients are selected to mirror values shown in the surviving 2016 appendix.

## Visualization and export

The computed 6-panel overview above shows the **actual perforated geometry**. Individual figures are available under `results/figures/`, including `heat_flux.png` and `deformed_mesh.png`. Thermal flux is reported as the per-element vector `q = -k grad(T)` in W/m². The deformed-mesh diagram **magnifies** the displacement field; its amplification factor is shown in the plot title. Original and displaced meshes are overlaid. Raw element heat flux vectors and centroids are exported to `results/data/element_heat_flux.csv`. Stress contours are constant per linear element, without nodal smoothing.

## Verification status

The automated tests include:

- Mesh orientation, conforming inner boundary, displaced-hole geometry, rejection of invalid geometries.
- Symmetry and constant-field nullspace of heat conduction matrices.
- Isothermal problem: zero temperature gradients and numerically negligible thermal stress under compatible supports.
- Uniform temperature rise: exact free thermal expansion, with zero stress to numerical tolerance (plane-stress patch test).
- Thermal temperature bounds, free-DOF residuals, conservation of global conductive reactions.
- Transient backward Euler convergence to the steady solution and monotonic decay of the thermal energy norm under fixed-temperature homogeneous conditions.
- Exact affine-temperature field and constant-strain plane-stress patch test on a perforated mesh.
- Manufactured quadratic heat equation with explicitly prescribed boundary temperatures and constant volumetric source: **RMS free-node error decreases by 37.9×** across six mesh levels (`results/verification/`).
- Conflicting BC and invalid timestep errors; numerical postprocessing smoke tests.

**Not yet verified:** independent ABAQUS/scikit-fem numerical comparisons; actual 2016 displacement or temperature results; curved-hole geometry; quadratic elements; two-way thermoelastic coupling; dynamic mechanical inertia; contact or temperature-dependent constitutive response. Stress peaks near corners and constraints must be studied through mesh convergence; isolated element-centroid maxima are not certified as mesh-independent physical values. See [verification and limitations](docs/VERIFICATION.md).

## Repository layout

```text
src/thermofem2d/geometry.py   Parametric geometry-fitted T3 mesh
src/thermofem2d/fem.py        Physical element matrices and boundary assembly
src/thermofem2d/solve.py      Thermal and mechanical solvers
src/thermofem2d/plot.py       Scientific post-processing and six-panel dashboard
src/thermofem2d/cli.py        Reproducible command-line runner
app.py                        Optional Streamlit user interface
examples/                     Baseline and offset-hole models
results/                      Computed baseline CSV, plots, verification data
tests/                        Automated unit and physical-limit checks
docs/                         Original-source audit, derivations, verification
```

## Roadmap

- **0.2:** Establish source-specific historical baseline, investigate exact model dimensionality, loading and numerical benchmark references.
- **0.3:** Six-node quadratic triangular thermal and mechanical elements; integration and convergence checks.
- **0.5:** Independent geometry/mesh backend supporting multiple hole shapes and at least one quadrilateral family; richer BC and visualization.
- **1.0:** Source-traceable coupled thermoelastic system, dynamic mechanics, comprehensive benchmarks, numerical diagnostics, reproducible scientific reports.

The scope of these releases is contingent on mathematical verification, not merely adding interface controls.

## Provenance, citation, and license

Historical coursework: **2016 (1395)**. New Python software and numerical tests: **2026**. The new solver is not a literal translation of the original MATLAB appendix. Cite software results as *ThermoFEM-2D v0.1.1 (2026)* with the specific input configuration and numerical verification status. MIT License.