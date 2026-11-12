# Change history

## v0.1.1 — 12 November 2026

- Added analytically specified harmonic thermal and affine constant-strain patch tests on perforated meshes.
- Added quadratic manufactured heat-conduction benchmark with reproducible six-level convergence data and plot.
- Added transient diffusion energy-decay verification, CLI/export and plotting tests (27 passing numerical/post-processing tests).
- Added Fourier heat-flux recovery and raw element CSV, magnified deformed-mesh visualization, and a six-panel results overview.
- Added optional Streamlit end-to-end smoke test as a dedicated future CI job; not yet run in the local environment.
- Preserved transparent scope: T3, one-way heat-to-mechanics coupling, quasi-static plane stress; historical T6 two-way coupling not yet recovered.

## v0.1.0 — 12 November 2026

- Established a standalone Python package and optional Streamlit interface.
- Implemented a geometry-conforming triangular mesh for a rectangular plate with a movable/resizable rectangular opening.
- Added T3 steady/transient thermal FEM, one-way plane-stress mechanical FEM and thermal-expansion loading.
- Added numerical tests, synthetic verification cases, CSV export and computed scientific plots.
- Documented the 2016 source and differences between the historical coupled T6 formulation and this narrower 2026 baseline.