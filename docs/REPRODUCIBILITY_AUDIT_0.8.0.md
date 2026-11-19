# v0.8 numerical-output reproducibility audit

**Status:** Accepted with documented non-pixel-identical figures; this is not an external physics-validation certificate.  
**Repository:** `MehdiSoroush/thermofem2d`  
**Audited commit:** [merged PR #10, `ba29e63abb6d84fa267029ea3391ba9a348f4d32`](https://github.com/MehdiSoroush/thermofem2d/commit/ba29e63abb6d84fa267029ea3391ba9a348f4d32)  
**Historical result comparison:** `6d4dcb000c321858ae6a4fb22582041896c7098f` → `91a1566a63172b37767a7412458748dd786c8641`  
**Evidence basis:** recorded local audit terminal logs, generated JSON reports and a four-figure visual comparison reviewed during release preparation. This is a documented reproducibility check, not a third-party independent test.

## Method

An audit utility cloned a disposable checkout at the **exact merged PR #10 commit**, removed 77 tracked numerical and image outputs **only inside the scratch copy**, and regenerated the data using the 12 commands from the repository's results-generation workflow. The 77 files comprised 73 outputs modified in the previous automated results commit and four Kirsch verification outputs introduced by PR #10. The user's working checkout, GitHub `main` and installed virtual environments were not modified by the audit.

**Local comparison environment:** Windows; CPython 3.14.10 (`science314`); NumPy 2.5.4; SciPy 1.18.1; Matplotlib 3.11.2; Pillow 13.0.0. The published results were produced in the repository's GitHub Actions workflow under a potentially different Python and plotting environment. Exact image bytes were not assumed to be portable.

## Results

| Check | Observed result | Disposition |
| --- | ---: | --- |
| Clean regeneration commands | 12 / 12 succeeded | Pass |
| Target outputs present after regeneration | 77 / 77 | Pass |
| CSV/JSON structured output comparisons | 31 / 31 passed under declared tolerances | Pass |
| PNG images identical at pixel level | 42 / 46 | Pass |
| PNG images nonidentical but visually reviewed | 4 / 46 | Qualified acceptance |

The v2 strict comparison used `rtol=1e-8` and general `atol=1e-12`. Only heat-flux components were given a separate near-zero absolute tolerance `2e-11 W/m²`. The original strict comparison had flagged two heat-flux CSV files because reference values of order `1e-11 W/m²` regenerated as zero. The largest absolute differences recorded in those files were below `9e-10 W/m²`, on fields with representative nonzero fluxes of order `10^3 W/m²`. Allowing a near-zero tolerance is a reproducibility convention, **not** an accuracy bound for the continuum model.

## Four graphical differences

The following files were regenerated but were **not pixel-identical**:

1. `results/verification/coupled_manufactured_displacement.png`
2. `results/verification/coupled_manufactured_temperature.png`
3. `results/verification/kirsch_hole_angular_stress.png`
4. `results/verification/kirsch_stress_convergence.png`

A three-column montage (published, regenerated, enhanced difference) was reviewed for curve agreement, axis values, labels, markers and benchmark trend interpretation. The curves and scientific interpretations were visually consistent. The enhanced differences were concentrated at text, grid/line and marker edges and numerically tiny traces. The precise cause of the pixel differences was **not established**, and they are not described as pixel-identical.

Near-machine-precision manufactured-field errors are numerical precision checks and must not be extrapolated to general continuum accuracy.

## Boundary of this conclusion

The audit supports **reproducibility of the selected outputs at the audited commit within documented tolerances**. It does **not** demonstrate agreement against a matched external FEM solver, certify continuous stress maxima, prove curved-hole dynamic convergence, or validate thermal-shock material failure predictions. Original 2016 MATLAB model results were not recovered.

**Final release-head qualification:** The repository's `main` advanced after PR #10 through an automated results-publishing commit that modified three generated verification outputs. This audit was pinned to `ba29e63`; final release CI and integrity checks must therefore be rerun at the reviewed release-candidate commit. The historical three-file change should be independently inspected if byte-level provenance of the final tag is required.

The four image exceptions are documented rather than silently waived. No new scientific features or solver-code changes are implied by this audit record.
