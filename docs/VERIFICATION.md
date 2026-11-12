# Numerical verification and known limitations — 0.1.1

## Current verification cases

1. **Element orientation:** all T3 elements have positive area after orientation enforcement; boundary edges follow actual outer and inner faces.
2. **Mesh geometry:** hole size and eccentricity can vary without editing the assembly code; cells remain outside the hole.
3. **Conductivity nullspace:** a constant temperature vector yields numerical zero `K*1`, consistent with ∇T = 0.
4. **Constant temperature:** uniform thermal boundary data produce uniform temperature to numerical roundoff.
5. **Exact uniform-expansion patch:** for `T=350 K`, `T_ref=298 K`, the exact stress-free solution under roller constraints is `ux=alpha*(T-T_ref)*(x+width/2)`, `uy=alpha*(T-T_ref)*(y+height/2)`. Displacement and stress are compared to that reference.
6. **Heat balance:** imposed thermal reaction fluxes sum to zero for steady conduction without sources.
7. **Transient limit:** backward-Euler heat conduction relaxes toward the independently computed steady-state solution under time-invariant boundaries.
8. **Invalid cases:** conflicting temperature data, degenerate holes and invalid timestep are rejected.

Run `python -m pytest -q` to reproduce automated verification. A separate `python examples/scripts/mesh_study.py` compares coarse-to-fine maximum displacement and **unconverged maximum centroid stress** in the reference model. The displacement varies slightly across the tested meshes; centroid stress peaks increase with refinement. This study is a diagnostic, not an error bound or proof of stress convergence.

## Baseline example

- Plate: 2 m × 2 m, thickness 1 m (numerical benchmark thickness; not a practical design assumption).
- Hole: 0.24 m × 0.20 m, centered.
- Steel-like coefficients follow the original coursework appendix, with `T_ref = 298 K`.
- Left edge at 380 K; right edge at 298 K; top/bottom/hole insulated.
- Mechanical rollers: left edge `ux=0`, bottom edge `uy=0`.
- Plane-stress quasi-static deformation, no imposed mechanical traction.

**All numerical results and pictures are newly generated demonstration outputs.** They are not the historical 2016 report's results. The largest element-centroid von Mises stress is not a converged structural design limit, particularly where mixed constraints and corners create local gradients.

## Validation gaps

- Source-specific six-node coupled system and original MATLAB results are not recovered.
- Independent reference solver comparison (e.g. separate FE code or analytically accessible nontrivial problem) remains outstanding.
- Time-step and mesh-convergence studies are required for the nonuniform demonstration case.
- A posteriori error estimates and high-order or curved elements are not yet implemented.
- Streamlit interface requires separate end-to-end runtime verification.
- Contact, geometric/material nonlinearities, convection, thermal radiation, prescribed time-dependent BC functions, arbitrary holes and dynamic structural inertia are not available.

## Additional manufactured-solution verification (v0.1.1)

These are **analytical continuum checks**, not a claim of independent Abaqus benchmark reproduction. The tests use the same custom T3 assembler, while the reference fields and differential equations are specified analytically.

### Exact harmonic temperature field

The function `T(x,y) = 320 + 4*x - 7*y` K satisfies `div(k grad T) = 0` for constant conductivity. We prescribe this exact temperature on **all** outer and hole-boundary nodes, including the interior hole boundary (which therefore is *not insulated* in this benchmark). T3 reproduces the linear field to numerical precision, and the elementwise heat flux must be `q = -k*[4,-7]`.
This checks both the 2D stiffness assembly and Fourier-flux recovery.

### Exact constant-strain mechanical patch

With `T = T_ref`, we prescribe an affine displacement field for a nonzero constant strain vector `epsilon = [1e-4,-0.7e-4,0.35e-4]`. The solver reproduces the nodal displacements and element stresses against `D*epsilon` to numerical precision. The test includes a perforated plate.

### Quadratic manufactured conduction benchmark

With `T(x,y) = 320 + 12*(x² + y²)` K and constant isotropic conductivity,
`Q = -k*Laplacian(T) = -2592 W/m³`. The exact temperature is imposed on all
outer and hole-boundary nodes. For this field, T3 is **not exact**, so a
systematic reduction in numerical error is expected with refinement.

| Per-side / radial divisions | Nodes | Elements | RMS free-node error (K) |
| ---: | ---: | ---: | ---: |
| 3 | 48 | 72 | 0.554117 |
| 4 | 80 | 128 | 0.287153 |
| 6 | 168 | 288 | 0.116378 |
| 8 | 288 | 512 | 0.062435 |
| 12 | 624 | 1152 | 0.026540 |
| 16 | 1088 | 2048 | 0.014637 |

The observed error reduction is approximately **37.9×** from the first to last tested mesh. This is an error-norm convergence demonstration for one manufactured conduction problem, **not** a proof of convergence of local von Mises stress maxima in the original plate loading case.

The values and plot are regenerated with:

```bash
python examples/scripts/verify_manufactured.py
```

Outputs: `results/verification/quadratic_heat_convergence.csv`,
`quadratic_heat_convergence.png`, and `benchmark_summary.json`.

### Transient diffusion stability

For fixed uniform boundary temperature, zero source, and a warmer initial
interior, the discrete deviation energy `0.5*(T-T_B)^T*C*(T-T_B)` decreases
monotonically with each backward-Euler step in the tested case. This validates
one expected property of implicit diffusion time stepping; it does not
establish temporal accuracy in the full nonuniform loading scenario.

### GUI verification status

The core numerical tests and CLI/postprocessing smoke tests pass in the
present local environment. The optional `tests/test_streamlit_app.py` is
skipped when `streamlit` is unavailable; the GitHub Actions UI smoke test has passed on Ubuntu with Python 3.12. This establishes scripted interface initialization and one numerical solve, but not full browser/deployment validation.
