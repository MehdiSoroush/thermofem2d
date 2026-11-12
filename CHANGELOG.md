# Change history

## v0.3.0 — 12 November 2026

- Added a monolithic **two-way coupled transient plane-stress** solver for T3 and T6.
- Assembled a consistent mechanical mass and reciprocal thermal/mechanical matrices; implemented a combined Newmark average-acceleration and backward-Euler time integrator.
- Added a reproducible coupled T6 case, transient field histories, heat and mechanical-response plots, and a Streamlit mode selector.
- Added dedicated coupling/residual tests and documentation of the plane-stress reduction and historical-source distinctions.
- Preserved the sequential thermal/quasi-static baseline; no historical MATLAB or Abaqus reproduction is claimed.

## v0.2.0 — 12 November 2026

- Added conforming quadratic T6 elements with globally shared midside nodes,
  exact quadratic-edge constant loads and degree-four triangular quadrature.
- Generalized thermal conductivity/capacity, plane-stress stiffness,
  thermal-expansion loading, heat flux and centroid stress recovery to T3/T6.
- Added a T6 selector in the Streamlit interface and CLI JSON configurations.
- Added source-accurate T6 formulation notes and T3-versus-T6 numerical comparison.
- Added nine analytical/numerical T6 tests (36 passed locally, one optional GUI test skipped).
- Retained **one-way** coupling; historical 2016 fully coupled three-field T6 not reconstructed.


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
