# Source audit: 2016 coursework versus the modern implementation

## Historical source available

Original Persian report: *Finite Element Methods coursework project: Thermoelasticity of a square plate with a square opening, MATLAB*, summer 1395 (2016), 39 PDF pages. It has an incomplete MATLAB appendix. The document is not included in the public repository.

## What is directly supported by the surviving report

- The title and figure identify a square plate with a square/rectangular central opening, a thermal boundary-value problem and mechanical surface tractions.
- The mathematical exposition proposes a **classical coupled thermoelastic** formulation. Its governing equations include dynamic momentum balance and thermal strain-rate feedback into the heat equation.
- The original element is a **quadratic triangle with six nodes and three nodal fields**: two mechanical displacements and temperature.
- The MATLAB appendix identifies illustrative steel properties: `E = 200e9 Pa`, `nu = 0.3`, `alpha = 11.7e-6 K^-1`, `rho = 7830 kg/m^3`, `k = 54 W/(m K)`, `cv = 452 J/(kg K)`, `T0 = 298 K`.
- The appendix shows an initial geometry of square outer bounds x,y∈[-1,1] and half-hole dimensions `a = 0.12`, `b = 0.10`, and a code fragment assembling thermal/mechanical matrices.
- The surviving code ends during boundary-condition specification, before final load assembly, initial conditions, time marching, numerical results and complete post-processing are shown.

## Source ambiguities requiring audit

1. Element interpolation, ordering, Jacobians and quadrature in original code require mathematical re-derivation and matrix-level checks.
2. Original material use of 3D Lamé constants must be reconciled with a consistent 2D **plane-strain/plane-stress** reduction.
3. Original dimensional interpretation of unit-thickness and coordinates must be checked before reproducing any physical numerical output.
4. Historical tractions and thermal boundary conditions are not sufficiently specified by the surviving final script fragment for an exact run.
5. The thermomechanical feedback term in the temperature equation cannot be considered reproduced by a sequential thermal-then-mechanical solver.

## Explicit scope of v0.1.1

A clean T3 mesh, Fourier conductivity/capacity matrices, steady/transient heat equation, and one-way quasi-static plane-stress thermal loading have been implemented. The current example boundary conditions and heating history are **new demonstration inputs** chosen to produce a testable baseline; they are not recovered 2016 cases. No claim of experimental or Abaqus validation is made.

## Source material handling

The original Persian report and its incomplete MATLAB appendix remain private. Public files independently describe the mathematics, provenance, and computed results without reproducing the archival manuscript.
