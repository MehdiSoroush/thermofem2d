# Mathematical formulation — ThermoFEM-2D v0.2.0

## Domains and approximation

The computational domain is a rectangle Ω from which an axis-aligned rectangular hole has been removed. The mesh consists of conforming three-node linear triangles (T3) or **straight-sided six-node quadratic triangles (T6)**, with interpolation `N_i(x,y)` and oriented area `A > 0`. Every T6 edge midpoint is shared between adjoining triangles. The hole is empty; no elements are generated inside it. `t` is plate thickness. Fields are in SI units.

## Thermal problem

Fourier heat conduction with optional volumetric source `Q`:

`rho * c * dT/dt - div(k grad T) = Q`.

The finite element semidiscrete equation is:

`C * dT/dt + K_T * T = f_Q + f_q`, where for a linear triangle (T3):

- `K_T(e) = t * k * A * (grad N)^T (grad N)`;
- `C(e) = rho*c*t*A/12 * [[2,1,1],[1,2,1],[1,1,2]]`;
- `f_Q(e) = t * Q * A/3 * [1,1,1]` for constant heat generation;
- `f_q` is an integrated **inward** surface heat flux (positive adds heat).

Steady state solves `K_T T = f`. Transient backward Euler solves
`(C + Δt K_T) T^(n+1) = C T^n + Δt f^(n+1)` with strong Dirichlet elimination. Natural heat flux defaults to zero on unlabeled/unspecified surfaces.

## Mechanical response (one-way only)

Under small-strain plane stress:

`eps = [∂ux/∂x, ∂uy/∂y, ∂ux/∂y + ∂uy/∂x]^T`,
`eps_th = alpha * (T - T_ref) * [1,1,0]^T`,
`D = E/(1-nu²) * [[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]]`,
`σ = D(eps - eps_th)`.

The mechanical equilibrium is `K_u u = f_traction + f_thermal`, where:

- `K_u(e) = t*A*B^T*D*B`;
- `f_thermal(e) = t*A*B^T*D*alpha*(mean(T_node)-T_ref)*[1,1,0]^T`.

For T3, the gradient matrix B is constant per element. For T6, the element matrices and thermal loading are integrated over each triangle using a six-point degree-four Gauss rule, with a spatially varying B and quadratic temperature interpolation. See [T6 formulation](T6_ELEMENTS.md). Linear temperature interpolation gives exact integration of uniform-material thermal loading using the element mean temperature. Standard consistent line-edge traction loads are included for constant tractions.

The mechanical step is **quasi-static** and does not contain mass `M * ü` or thermal strain-rate feedback `T_ref*beta*tr(εdot)`. Therefore the method is *not two-way coupled thermoelasticity* and does not reproduce the original classical dynamic formulation.

## Derived outputs

Each element's stress is evaluated at its centroid from the interpolated strain and temperature (the simple nodal mean is exact for T3 only). Plane-stress von Mises equivalent stress:

`σ_vm = sqrt(σxx² - σxx σyy + σyy² + 3 τxy²)`.

The stress plot uses raw per-element constant stresses; it does not smooth/extrapolate them onto mesh nodes. Peaks near corners, sharp-hole vertices, or mixed boundary constraints require mesh convergence and physical interpretation.

## Boundary naming and solvability

Outer edges are `left`, `right`, `top`, `bottom`; the inner boundary is `hole`. Dirichlet values at shared corner nodes must agree. The mechanical system requires independent displacement constraints sufficient to remove all rigid body modes. A single fixed node would not suffice in general. The given example's left `ux=0` and bottom `uy=0` constraints eliminate rigid translations and rotation but also represent idealized roller conditions, not a realistic universal mounting model.


## Postprocessing extensions (v0.1.1)

For T3 elements with constant temperature gradients, the heat-flux density
is `q_e = -k * (grad N)^T * T_e` [W/m²]. For T6 the same expression uses the varying quadratic gradients evaluated at the parent-element centroid. No nodal smoothing or extrapolation is applied. It is exported at
geometric element centroids in `element_heat_flux.csv`.

The displayed deformed mesh uses `x_display = x + scale * u`, with the scalar
magnification factor explicitly written on the figure. It must not be read as
a full-size physical deformation unless `scale = 1`. This does not change
any recovered stresses or other solver results.
