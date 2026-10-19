# Experiment parameters and unit mapping

`experiment.json` defines the physical inputs and numerical settings. `parameters.py` derives the solver parameters and records their mapping in the run's `parameters.json`. Values below are rounded for display; conversion uses the input values without intermediate rounding.

## Physical experiment

The body occupies 0 ≤ X ≤ L and 0 ≤ Y ≤ H in its reference configuration. It is initially undeformed and at rest. The left edge is fixed, the top and bottom edges are traction-free, and the right edge carries nominal traction per reference area:

$$
\mathbf T(t)=A\sin(2\pi f t)(1,0).
$$

The model is compressible logarithmic Neo-Hookean elasticity under plane strain, with Fzz = 1. Density is volumetric and the elastic moduli and traction are three-dimensional stress quantities. The two-dimensional calculation is interpreted per unit out-of-plane thickness. The reference density is an illustrative input. There is no body force, explicit damping, plasticity, or damage.

The constitutive law returns Cauchy stress. On the moving right edge, the solver converts nominal traction to current traction by dividing by the boundary stretch, `ds/ds0`, where `ds` and `ds0` are the current and reference segment lengths. The wave speeds below are small-strain values in the reference state.

| Quantity | Physical value |
| --- | ---: |
| Reference length L | 0.06 m |
| Reference height H | 0.018 m |
| Reference density ρ₀ | 7850 kg/m³ |
| First Lamé parameter λ | 115 GPa |
| Shear modulus μ | 77 GPa |
| Young's modulus E = μ(3λ + 2μ)/(λ + μ) | 200.119791667 GPa |
| Poisson ratio ν = λ/[2(λ + μ)] | 0.299479166667 |
| Bulk modulus K = λ + 2μ/3 | 166.333333333 GPa |
| Longitudinal speed cₚ = √[(λ + 2μ)/ρ₀] | 5853.846250421 m/s |
| Shear speed cₛ = √(μ/ρ₀) | 3131.919091779 m/s |
| Traction amplitude A | 1 MPa |
| Loading frequency f = cₚ/(2L) | 48782.052086841 Hz |
| End time 1/f = 2L/cₚ | 20.499342632 µs |
| Physical cell spacing Δx = L/60 | 1 mm |
| Physical time step Δt | 0.085413927632 µs |

## Physical to normalized solver units

Choose the longitudinal Courant number Cₚ = cₚ Δt/Δx = 0.5. It sets the physical streaming speed c = Δx/Δt = cₚ/Cₚ. The normalization scales are

$$
\ell_* = L,\qquad
v_* = c,\qquad
t_* = \frac{L}{c},\qquad
\rho_* = \rho_0,\qquad
\sigma_* = \rho_0c^2 = \frac{\lambda+2\mu}{C_p^2}.
$$

| Scale | Value |
| --- | ---: |
| Length ℓ* | 0.06 m |
| Velocity v* | 11707.692500842 m/s |
| Time t* | 5.124835657896 µs |
| Density ρ* | 7850 kg/m³ |
| Stress σ* | 1.076 × 10¹² Pa |

The suffix `code` denotes the normalized quantities used by the solver:

$$
x_{\rm code}=\frac{x}{\ell_*},\quad
t_{\rm code}=\frac{t}{t_*},\quad
v_{\rm code}=\frac{v}{v_*},\quad
\rho_{\rm code}=\frac{\rho}{\rho_*},\quad
\sigma_{\rm code}=\frac{\sigma}{\sigma_*}.
$$

The same stress scale applies to λ, μ, and nominal traction. The deformation gradient F and its determinant J are dimensionless and unchanged by this conversion.

| Solver quantity | Mapping | Value |
| --- | --- | ---: |
| `length` | L/ℓ* | 1 |
| `height` | H/ℓ* | 0.3 |
| `rho0` | ρ₀/ρ* | 1 |
| `lam` | λ/σ* | 0.10687732342007435 |
| `mu` | μ/σ* | 0.07156133828996282 |
| `amplitude` | A/σ* | 9.29368029739777 × 10⁻⁷ |
| Loading frequency | f t* | 0.25 |
| `omega` | 2π f t* | π/2 |
| `end_time` | (1/f)/t* | 4 |
| `h` | Δx/ℓ* | 1/60 |
| `dt` | Δt/t* | 1/60 |

The solver keeps `dt = h = 1/60`, giving unit streaming speed in normalized coordinates. One loading period takes 240 steps.

## Normalized solver to lattice units

Lattice units count cell spacings and time steps:

$$
x_{\rm lattice}=\frac{x_{\rm code}}{h}=60x_{\rm code},\qquad
t_{\rm lattice}=\frac{t_{\rm code}}{\mathrm{dt}}=60t_{\rm code}.
$$

Thus the reference dimensions are 60 × 18 cells and the lattice spacing and time step are both 1. Because `dt/h = 1`, velocity has the same numerical value in lattice and normalized solver units. Density, stress, moduli, and traction retain the normalization given above. The longitudinal and shear speeds in lattice units are Cₚ = 0.5 and Cₛ = 0.2675095106532903. The loading frequency is 1/240 cycles per step, with angular frequency π/120 radians per step.

The D2Q9 equilibrium constant `CS2 = 1/3` is distinct from the elastic shear-wave speed Cₛ.

| Numerical setting | Value |
| --- | ---: |
| Initial body resolution | 60 × 18 cell centres |
| Background padding | 8 cells on each side |
| Background array size | 76 × 34 |
| Passive material sampling grid | 61 × 19 points |
| `tau_f`, `tau_g` | 0.85, 0.85 |
| Recovery tolerance | 2 × 10⁻¹¹ |
| Maximum recovery iterations | 40 |
| Recovery relaxation | 1 |

`tau_f` and `tau_g` are dimensionless BGK relaxation parameters measured relative to a lattice step. They are numerical settings, not prescribed physical viscosity or material relaxation times. The passive material samples provide output coordinates; they do not set the computational resolution.

## Convert saved fields to physical units

The solver saves normalized quantities. For a material point, form displacement from its saved coordinates:

$$
\mathbf u_{\rm code}=\texttt{positions}-\texttt{reference\_points},\qquad
\mathbf u_{\rm SI}=0.06\,\mathbf u_{\rm code}\;\mathrm{m}.
$$

Displacement in lattice cells is instead

$$
\mathbf u_{\rm lattice}=\frac{\mathbf u_{\rm code}}{h},\qquad
\mathbf u_{\rm SI}=\Delta x\,\mathbf u_{\rm lattice}.
$$

Multiplying the saved normalized displacement directly by Δx would give an incorrect conversion.

| Saved or derived quantity | Physical conversion |
| --- | --- |
| Coordinates and displacement | Multiply by 0.06 m |
| `time` | Multiply by 5.124835657896368 × 10⁻⁶ s |
| `velocity` | Multiply by 11707.692500841807 m/s |
| `density` or `rho` | Multiply by 7850 kg/m³ |
| Stress and nominal traction | Multiply by 1.076 × 10¹² Pa |
| `F`, J | No conversion |

The displacement plots use the reference coordinates for the spatial axes and show the global x and y components and their Euclidean magnitude. The tip history refers to the midpoint of the right reference edge, (L, H/2).

| Plot time in solver units | Physical time (µs) |
| ---: | ---: |
| 0.5 | 2.562417829 |
| 1 | 5.124835658 |
| 1.5 | 7.687253487 |
| 2 | 10.249671316 |
| 3 | 15.374506974 |
| 4 | 20.499342632 |
