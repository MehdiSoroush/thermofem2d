# Two-way transient thermoelasticity — v0.3.0

## Scope and provenance

This is a **new Python implementation** of classical linear, two-way thermoelasticity
for a plate with an opening. It is informed by the 2016 coursework's classical coupled
balance equations (manuscript equations 21–22 and its weak-form development through
41), but it is **not** a byte-for-byte recovery of the incomplete 2016 MATLAB code or
an independent reproduction of historical plots. The 2016 source uses a 3D isotropic
Lamé representation; this version uses a *reduced plane-stress* constitutive model.
The source report is not included in this public repository.

Available discretizations: conforming T3 and affine T6 triangles, identical shape
functions for temperature and displacement but separate nodal unknown blocks.
Mechanical inertia, thermal capacity, conductivity, and both thermoelastic coupling
directions are assembled consistently with six-point degree-four quadrature.

## Strong form and assumptions

For `theta = T - T_ref` with fixed absolute reference temperature `T_ref > 0`:

```text
  div(sigma) + b = rho * u_ddot
  rho * c * theta_dot - div(k * grad(theta))
      + T_ref * (D * eps_alpha) : eps_dot = Q
  sigma = D * (eps(u) - eps_alpha * theta)
  eps_alpha = alpha * [1, 1, 0]^T  (engineering strain notation)
```

Here `D` is the **plane-stress** isotropic stiffness matrix, `k` is scalar
conductivity, and `rho*c` is the **effective heat capacity at fixed in-plane
strain** for the reduced plane-stress model. This is a modelling choice: a
3D constant-strain heat capacity or a catalog isobaric heat capacity is not
automatically interchangeable with it. A fully consistent 3D-to-plane-stress
reduction would require eliminating the out-of-plane strain in both mechanical
and thermodynamic free-energy terms; users supplying 3D material coefficients
must make the corresponding capacity conversion separately. The software does
not currently perform this conversion.

The linearization evaluates the strain-rate thermal coupling at `T_ref` rather
than the instantaneous `T`. It is most defensible when `|T-T_ref|/T_ref` is
small. Small displacement/strain, constant isotropic coefficients, no plasticity,
no convection/radiation, no damping beyond numerical thermal diffusion, and
Fourier heat transfer are assumed. Plane strain and generalized thermoelastic
heat-flux relaxation are **not supported** in this release.

## Coupled finite-element equations

Define, with thickness `h`, strain matrix `B`, scalar shape-function row `N`,
material density `rho`, and plane-stress expansion vector `a`:

```text
  Ku = int(Omega) h * B.T * D * B dA
  M  = int(Omega) h * rho * Nu.T * Nu dA
  Kt = int(Omega) h * k * grad(N).T * grad(N) dA
  Ct = int(Omega) h * rho * c * N.T * N dA
  G  = int(Omega) h * B.T * D * a * N dA
  H  = T_ref * G.T
```

`Nu` is the two-component vector shape-function matrix. `G` has size
`(2*num_nodes, num_nodes)` and `H` its reciprocal transpose, scaled by
reference temperature. The semidiscrete equations are:

```text
    M*u_ddot + Ku*u - G*theta = f_mech
    Ct*theta_dot + Kt*theta + H*u_dot = f_heat
```

The sign is fixed by the mechanical thermal-expansion stress, and the
transpose relation is checked in the test suite. For an insulated, unforced
system, multiplication by velocities and `theta/T_ref` gives the formal
energy balance

```text
  d/dt [ 0.5*u_dot.T*M*u_dot + 0.5*u.T*Ku*u
       + 0.5*theta.T*Ct*theta / T_ref ]
    = -theta.T*Kt*theta / T_ref
```

under homogeneous compatible essential constraints. This makes the
thermodynamic significance of the feedback sign explicit.

## Time discretization

The mechanical block uses **Newmark average acceleration** with
`beta = 1/4, gamma = 1/2` (unconditionally stable for uncoupled linear
elastodynamics). Temperature uses **backward Euler**. The two fields at each
new step are solved **together** rather than sequentially:

```text
A = [ Ku + M/(beta*dt^2)        -G            ]
    [ H*gamma/(beta*dt)          Kt + Ct/dt   ]
```

The right-hand side uses standard displacement/velocity Newmark predictors
and the previous temperature. Constant prescribed mechanical and thermal
boundary conditions are eliminated from the monolithic block; its sparse
factorization is reused for every step. With `thermoelastic_feedback=false`,
`H=0` while the mechanical thermal-expansion term `-G*theta` remains, providing
a **dynamic, one-way** comparison. The uncoupled heat problem is still transient.

## Initial and boundary conditions

- Uniform interior `initial_temperature` (Kelvin) is supported. Prescribed
  temperature values are enforced at boundary nodes at `t=0+`.
- Mechanical initial displacements and velocities are **zero** by default.
  Nonzero `initial_displacement` and `initial_velocity` fields are not implemented.
- Prescribed displacements and temperatures are **constant in time**.
- Thermal Neumann boundaries default to zero inward flux. A **pure Neumann**
  transient thermal problem is supported because the capacity matrix is nonsingular.
- Prescribed tractions, inward fluxes, and volumetric sources are constant
  throughout the solve. Time-varying boundary functions and loads are future work.
- The mechanical constraints must eliminate rigid-body modes; this is not
  automatically guaranteed by the presence of a single Dirichlet condition.

## Reproducible example

```bash
python -m pip install -e ".[dev]"
python -m pytest -q tests/test_coupled.py
python -m thermofem2d.cli --config examples/coupled_plate_t6.json --out results-coupled
```

The example has an offset rectangular hole, 240 **T6 parent elements**, 528
nodes, transient temperatures fixed on left/right boundaries, roller mechanical
constraints, a fixed top traction, and step-wise simultaneous temperature and
displacement integration. The saved CSVs include nodal fields and full
mechanical displacement/velocity histories. The plots contain temperature,
displacement, centroid stress, Fourier flux, magnified deformed mesh, and
transient mechanical and thermal histories.

**This is a synthetic calculation**, not a reconstructed experimental record.
Stress peaks around opening corners and constrained edges are sensitive to
mesh refinement, and temperature-displacement histories should not be used
for machine design without independent numerical comparison and uncertainty
assessment.

## Verification status and omissions

Automated tests cover reciprocal element matrices, the consistent mechanical
mass, rigid-translation invariance, compatibility between `G*theta` and the
independent thermal-load assembler, stress-free reference equilibrium,
alpha-zero decoupling, mechanical-to-thermal feedback, residuals at free DOFs,
positive finite transient fields, and both T3/T6 discretizations.

Outstanding tasks include external solver comparison (Abaqus, scikit-fem, or
another independently derived model), quantified transient mesh/time-step
convergence, generalized time-varying conditions, curved-edge holes, plane
strain, and a complete re-derivation against the original 2016 matrices and
results. The current numerical tests **do not prove industrial validation**.


## Version 0.4 verification extension

Version 0.4 introduces a separately implemented DOP853 method-of-lines
reference time solver and a linearized thermoelastic energy-budget audit.
See [Independent Verification](INDEPENDENT_VERIFICATION.md). These tools
compare time integration and thermodynamic balances on the **same FE mesh**;
they do not constitute independent spatial FEM or external validation.

The initial temperature may be a uniform absolute temperature in kelvin or
an affine field specified as `{"type":"affine", "base":310, "gx":6, "gy":-3}`
(K/m gradients). Constant mechanical initial displacement and velocity remain
zero in this release. The coupled reference temperature must be positive K.
