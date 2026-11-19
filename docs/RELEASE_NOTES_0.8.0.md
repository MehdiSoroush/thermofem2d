# ThermoFEM-2D v0.8.0 — Circular-hole and curved-element verification

**Mehdi Soroush · November 2026 · Release-candidate documentation**

Version 0.8.0 extends the reproducible **research verification** capability of ThermoFEM-2D to circular-hole geometry, quadratically curved isoparametric T6 element assembly, and static analytical stress-concentration benchmarks. It does **not** make circular holes selectable in the ordinary JSON/CLI or Streamlit application. The application-supported rectangular-hole solver remains available with straight-sided T3 and T6 elements.

## Additions

- **Circular-hole meshes:** an explicitly separate geometry utility for centred and eccentric circular openings, polygonal T3 hole edges and curved quadratic T6 edges. Area, perimeter, sampled Jacobian positivity and conforming connectivity are tested. The quadratic boundary approximates, but does not exactly represent, a circle.
- **Curved T6 field operators:** an explicit geometry-to-assembly adapter and numerical integration using the mapped Jacobian for conductivity, heat capacity, plane-stress elasticity, structural mass, reciprocal thermoelastic coupling and curved-edge loads. Established affine-element paths remain in place.
- **Analytical verification:** harmonic conduction on an eccentric-hole domain, constant-strain/thermal-expansion mechanical patches, mass and coupling reciprocity, and integrated boundary load checks.
- **Static Kirsch benchmark:** analytical plane-stress displacement and stress reference for an infinite medium with a circular hole under remote x-directed tension. The computational subdomain is assigned **exact analytical exterior displacements** and the hole has natural zero traction. One-sided, unsmoothed T3 and curved T6 hole-edge stress samples are compared with the analytical field.
- **Separate finite-width traction study:** square plates loaded through uniform opposing exterior tractions. Gross- and net-section stress factors are separately defined; their values are **not** the infinite-plane analytical result for the Dirichlet benchmark.

## Representative verification results

In the analytical-exterior Kirsch benchmark, the exact infinite-plane hole-edge maximum is `K_t = 3` (gross nominal stress). At 12 subdivisions, the reported **sampled** curved T6 peak is approximately `3.0226`; sampled circumferential-stress RMS error is approximately `0.02270` of the applied stress. These are sampled discretization results, **not** a rigorous continuous maximum or general bound on local stress error. At this mesh, hole-edge radial stress RMS is approximately `0.0812` of the applied stress; weak enforcement of natural traction does not make pointwise stress exactly traction-free.

For the **distinct finite-plate traction** problem, T6 at 12 subdivisions yields a sampled gross stress factor of about `3.6090` at width 2 m (hole diameter/width = 0.25), and about `3.1986` at width 4 m (hole diameter/width = 0.125). Different widths also yield different physical mesh spacings at fixed subdivisions; they are not identical local-resolution tests. The source data and definitions are in [Kirsch stress verification](KIRSCH_STRESS_CONCENTRATION.md).

Curved-element harmonic-conduction, mechanical patch, consistent mass and reciprocal coupling checks are described in [curved field assembly](CURVED_FIELD_ASSEMBLY.md). The new static checks supplement—but do not eliminate—the unresolved dynamic thermal-shock issues documented in the [validation status](VALIDATION_STATUS.md).

## Reproducibility and quality checks

From a clean Python 3.11+ environment at the repository root:

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python examples/scripts/verify_curved_fields.py
python examples/scripts/verify_kirsch.py
python examples/scripts/check_release_integrity.py
```

A separate Windows / CPython 3.14.10 **scratch-checkout regeneration audit** at merged PR #10 commit `ba29e63abb6d84fa267029ea3391ba9a348f4d32` removed 77 target outputs from a disposable checkout and regenerated them using 12 commands. All commands completed; 31 structured outputs agreed within stated tolerances, 42 PNG files were pixel-identical, and four more PNGs were accepted after visual plot comparison despite differing pixels. This is **reproducibility evidence**, not an independent spatial-FEM or experimental validation. See [the documented audit](REPRODUCIBILITY_AUDIT_0.8.0.md).

The release must still pass CI and release-metadata checks **at its final candidate commit**. A subsequent automated results-publishing commit changed three generated files after the audited PR #10 merge snapshot.

## Scope limitations

- Circular-hole geometry and curved T6 calculations remain dedicated verification utilities; they are not supported by the general application dispatcher.
- The Kirsch analytical-exterior case and finite uniform-traction plate are different boundary-value problems. Reported `K_t` values are **largest discrete samples**, not certified continuous or mesh-independent maxima.
- No matched external-FEM solver comparison has been completed for these exact curved-hole models, and no fracture, yielding, fatigue or design-allowable predictions are established.
- Mass/coupling assembly consistency does **not** validate full two-way transient curved-hole thermoelasticity.
- Early-time coupled thermal shock and cooling predictions retain substantial spatial/time-integration sensitivity; previous v0.7 qualifications still apply.
- The original 2016 MATLAB coursework matrices and numerical results have not been independently reproduced. All graphs and tables in this release are new 2026 computations; no new experiments were conducted.

**Intended use:** transparent scientific computing, code verification and reproducible research examples. Results must not be treated as certified structural-integrity or thermal-shock design values.
