# Quadratic triangular T6 elements — derivation and implementation

**Release:** 0.2.0 (2026). **Scope:** a **new one-way thermal-to-mechanical**
finite-element implementation. This is **not** an exact reimplementation of the
2016 report's fully coupled 18-component (`u_x`, `u_y`, `T` at six nodes)
transient thermoelastic system.

## Geometry, connectivity and polynomial interpolation

Each T6 triangle has three counter-clockwise vertices followed by shared
midside nodes in this order:

`[v0, v1, v2, m01, m12, m20]`.

Midside coordinates are arithmetic midpoints of the corner coordinates. Thus
this release uses a **straight-sided affine geometry mapping**; it does not
represent curved circular-hole boundaries with true quadratic isoparametric
geometry. Element centroids and areas are determined from the first three
vertices, while the T6 shape functions interpolate the fields at all six
nodes.

In reference coordinates `(r,s)`, define barycentric coordinates
`L1 = 1-r-s`, `L2 = r`, `L3 = s`. The six Lagrange shape functions are

```
N1 = L1 (2 L1 - 1)     N4 = 4 L1 L2
N2 = L2 (2 L2 - 1)     N5 = 4 L2 L3
N3 = L3 (2 L3 - 1)     N6 = 4 L3 L1
```

They satisfy `N_i(node_j)=δ_ij`, `sum N_i=1`, and `sum grad N_i=0`. For each
vertex function `grad[N_i]=(4 L_i-1) grad[L_i]`, while the midside gradients
follow the product rule. A constant affine element Jacobian transforms
reference gradients to the physical `(x,y)` coordinates.

## Conductivity and heat capacity

With homogeneous, isotropic conductivity `k`, thickness `t`, density `ρ`, and
specific heat `c`, the consistent local matrices are

```
K_T(e) = t ∫_Ωe k (∇N)^T (∇N) dΩ
C(e)   = t ∫_Ωe ρ c N^T N dΩ
f_Q(e) = t ∫_Ωe N^T Q dΩ
```

The T6 conductivity integrand is polynomial degree two and the **consistent
capacity** integrand has degree four. Therefore a symmetric **six-point
Dunavant degree-four** triangle rule is used for both. Quadrature weights are
normalized to sum to one; the integral is the triangle area multiplied by the
weighted sum. This also integrates constant volumetric sources correctly:
for T6, the three vertex nodal loads vanish and each midside node receives
`Q t A/3` (in contrast with T3, where each vertex receives `Q t A/3`).

The time discretization remains backward Euler:
`(C + Δt K_T) T^(n+1) = C T^n + Δt (f_Q + f_flux)`.

## Mechanical stiffness and temperature loading

For small-strain **plane stress**, at every integration point:

```
ε = B u,  B = [∂N_i/∂x and ∂N_i/∂y contributions]
ε_thermal = α (N T_e - T_ref) [1, 1, 0]^T
K_u(e) = t ∫_Ωe B^T D B dΩ
f_thermal(e) = t ∫_Ωe B^T D ε_thermal dΩ
```

T6 `B` varies linearly over the element and the temperature field varies
quadratically. The six-point degree-four rule integrates these contributions
for constant material coefficients on straight triangles. The mechanical field
has **12 local displacement DOFs**, assembled separately from the **six thermal
DOFs**. The current solver uses sequential one-way coupling and does not
assemble a monolithic 18×18 three-field thermoelastic tangent.

## Boundary integration

T3 edges use `[end0,end1]` and exact constant-edge load weights `[1/2,1/2]`.
T6 edges use `[end0,mid,end1]` and exact constant-edge load weights
`[1/6, 2/3, 1/6]` multiplied by physical edge length. Both prescribed
thermal and mechanical boundary conditions include midside nodes.

## Postprocessing

For T6, stress and Fourier heat flux are evaluated at the **parent-element
centroid** using quadratic interpolation and physical gradients. Stress/flux
remain spatially varying *within* an element; the reported values are point
samples, not elementwise constant physical fields. To render nodal temperature
and displacement on standard triangular plotting routines, each T6 element is
subdivided into four **visualization triangles**, using the existing six nodes;
this does not change the FEM solution. The single parent-centroid stress value
is repeated across the four subtriangles for a **piecewise-flat diagnostic**,
not claimed to be a resolved local stress distribution.

## Reproduce

```bash
python -m pip install -e '.[dev]'
python -m pytest -q tests/test_t6.py
python -m thermofem2d.cli --config examples/steady_plate_t6.json --out results-t6
python examples/scripts/compare_elements.py
```

The T6 test suite covers Kronecker interpolation, boundary load integration,
capacity definiteness, constant null modes, affine strain patch tests, uniform
stress-free expansion, the exact quadratic manufactured temperature field,
transient simulation and plotting. No claim of numerical reproduction of the
historical MATLAB outputs is made.

## Limitations

- One rectangular opening, straight-sided geometry, homogeneous linear material.
- No true curved T6 boundary elements or automatic mesh adaptation.
- No dynamic mechanical inertia or thermoelastic strain-rate feedback.
- Nonuniform field centroid stress extrema are **not** certified converged.
- The element implementation has analytical and cross-element tests, not yet an
  independent FEM-package comparison for the full perforated-plate load case.

## Independent public formulation references

- [DefElement: quadratic Lagrange triangle](https://defelement.org/elements/examples/triangle-lagrange-equispaced-2.html) describes the same six-node polynomial space (with a different edge-node numbering order).
- [J. Burkardt: six-point degree-four triangle integration](https://people.sc.fsu.edu/~jburkardt/classes/math1800_2023/geometry2/geometry2.pdf), Section “Triangle rule of higher degree”, gives the barycentric quadrature points and normalized weights used here.

These references document standard element mathematics. The implementation
and numerical results in this repository are new Python work.
