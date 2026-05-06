# Turbulence Modelling in OpenFOAM

OpenFOAM implements turbulence modelling through a flexible, runtime-selectable framework organised under the `turbulenceModel` abstraction layer. Models are split across three paradigms — **RANS**, **LES/DES**, and **laminar** — and are compatible with both incompressible (`incompressibleTurbulenceModel`) and compressible (`compressibleTurbulenceModel`) solver families. The library also supports buoyant and multiphase extensions. What follows is a systematic treatment of all major model families.

---

## 1. Architectural Overview

OpenFOAM's turbulence library lives in `$FOAM_SRC/TurbulenceModels/`. The class hierarchy is:

```
turbulenceModel
├── incompressible::turbulenceModel
│   ├── RASModel
│   ├── LESModel
│   └── laminar
└── compressible::turbulenceModel
    ├── RASModel
    ├── LESModel
    └── laminar
```

Model selection is purely runtime-driven via `constant/turbulenceProperties`:

```c++
simulationType  RAS;   // or LES, laminar

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
```

Transport coefficients are stored in `constant/transportProperties` (incompressible) or embedded in thermophysical models (compressible). Wall functions are boundary conditions applied at patch level, completely decoupled from the volume model — a significant architectural strength.

---

## 2. RANS / URANS Models

RANS models solve ensemble-averaged Navier–Stokes equations and introduce closure via modelled Reynolds stresses. All RANS models in OpenFOAM are also valid for URANS (unsteady) when run with a time-marching solver such as `pimpleFoam`.

### 2.1 Linear Eddy-Viscosity Models (LEVM)

These models invoke the Boussinesq hypothesis: $\tau_{ij} = 2\mu_t S_{ij} - \frac{2}{3}\rho k \delta_{ij}$.

#### 2.1.1 One-Equation Models

**Spalart–Allmaras (`SpalartAllmaras`)**
- Solves a single transport equation for the modified eddy viscosity $\tilde{\nu}$.
- Production driven by vorticity $\Omega$; destruction by a wall-proximity function $f_{w}$.
- Intended for attached or mildly separated aerodynamic flows; notoriously poor for flows with massive separation or strong adverse pressure gradients.
- **Wall treatment:** Direct integration to the wall (`y⁺ ≈ 1`) or low-Re formulation. No explicit wall function required.
- **Compressible:** Available as `Spalart-Allmaras` in `compressible::RASModels`.

#### 2.1.2 Two-Equation Models

**Standard k–ε (`kEpsilon`)**
- The workhorse of industry CFD. Solves transport equations for turbulent kinetic energy $k$ and dissipation $\varepsilon$.
- Model constants: $C_\mu = 0.09$, $C_{\varepsilon1} = 1.44$, $C_{\varepsilon2} = 1.92$, $\sigma_k = 1.0$, $\sigma_\varepsilon = 1.3$.
- Performs well in free shear flows; over-predicts $\mu_t$ in adverse pressure gradient flows; does not naturally respect rotational effects.
- **Wall treatment:** Requires wall functions (`kqRWallFunction`, `epsilonWallFunction`) for high-Re meshes ($30 < y^+ < 300$). Low-Re formulations (e.g., Launder–Sharma) available separately.

**RNG k–ε (`RNGkEpsilon`)**
- Derived from Renormalization Group theory (Yakhot & Orszag, 1986). Adds a strain-rate–dependent term $R$ to the $\varepsilon$ equation.
- Superior to standard k–ε in swirling flows and flows with streamline curvature.
- Constants differ: $C_\varepsilon2^* = 1.68$ (effective), $\sigma_k = 0.7194$, $\sigma_\varepsilon = 0.7194$.

**Realizable k–ε (`realizableKE`)**
- Ensures physical realizability constraints: $\overline{u_i^2} \geq 0$ and Schwarz inequality for $\overline{u_i u_j}$.
- $C_\mu$ becomes a function of mean strain and rotation: $C_\mu = f(U^*, \Omega^*)$.
- Better performance in round jets, recirculating flows, and strong rotation compared to standard k–ε.

**Standard k–ω (`kOmega`)**
- Transport equations for $k$ and specific dissipation rate $\omega = \varepsilon / (C_\mu k)$.
- Inherently integrates to the wall without wall functions — good near-wall behaviour.
- Highly sensitive to freestream $\omega$ boundary conditions in external flows.

**k–ω SST (`kOmegaSST`)** *(the de facto industrial standard)*
- Menter's (1994) Shear Stress Transport model blends k–ε in the freestream with k–ω near the wall via blending functions $F_1$ and $F_2$.
- Includes a stress limiter: $\mu_t = \rho k / \max(\omega, S F_2 / a_1)$, preventing the overproduction of $\mu_t$ in adverse pressure gradient regions.
- Outperforms both parent models across a wide range of flows: separated flows, turbomachinery, aerodynamics.
- **Variants in OpenFOAM:**
  - `kOmegaSSTSAS` — Scale-Adaptive Simulation extension (see §2.3)
  - `kOmegaSSTLM` — coupled with the $\gamma$–$Re_\theta$ transition model

**Launder–Sharma k–ε (`LaunderSharmaKE`)**
- Low-Reynolds number formulation with damping functions $f_\mu$, $f_1$, $f_2$ to resolve the viscous sublayer.
- Requires fine near-wall meshing ($y^+ < 1$); computationally expensive but avoids wall function assumptions.

**v2-f model (`v2f`) — via external library or community port**
- Solves for $\overline{v^2}$ (wall-normal stress) and an elliptic relaxation function $f$.
- Physically correct near-wall anisotropy without resolving the buffer layer explicitly.
- Highly accurate for heat transfer and flows dominated by wall-normal turbulence suppression.

**k–kL–ω (Three-equation, `kkLOmega`)**
- Separates turbulent kinetic energy into a large-scale component $k_L$ and a small-scale component $k_T$, alongside $\omega$.
- Designed to handle bypass and natural transition.

### 2.2 Non-Linear / Explicit Algebraic Reynolds Stress Models (EARSM)

These models extend the Boussinesq hypothesis to include quadratic or cubic terms in $S_{ij}$ and $\Omega_{ij}$, partially recovering anisotropic effects without solving full RSM transport equations.

**Non-linear k–ε (`NonlinearKEShih`)**
- Adds quadratic and cubic invariants of strain and rotation tensors to the stress–strain relation.
- Captures secondary flows in ducts, rotation, and streamline curvature effects missed by linear EVMs.

### 2.3 Reynolds Stress Models (RSM / Second-Moment Closure)

RSMs abandon the eddy-viscosity concept entirely and solve transport equations for all six independent components of the Reynolds stress tensor $\overline{u_i u_j}$, plus a scalar equation ($\varepsilon$ or $\omega$). Pressure–strain redistribution is the critical closure term.

**Launder–Gibson RSM (`LRR`)**
- Closure: linear pressure–strain model (Launder, Reece, Rodi, 1975): $\Pi_{ij} = -C_1 \varepsilon b_{ij} + C_2 k S_{ij} + ...$
- Accurately captures anisotropy-driven effects: secondary flows, swirl, curvature.
- Computationally expensive (7 extra transport equations); numerical stiffness a known challenge.

**Speziale–Sarkar–Gatski RSM (`SSG`)**
- Quadratic pressure–strain model (1991); superior to LRR for strongly anisotropic flows.
- Model constants calibrated for a wider range of canonical flows.

**ω-based RSM (`realizableKECp`, community variants)**
- Community-contributed RSMs using $\omega$ as the length-scale equation for better near-wall behaviour.

### 2.4 Scale-Adaptive Simulation (SAS)

**k–ω SST–SAS (`kOmegaSSTSAS`)**
- Menter & Egorov's (2010) formulation introduces the von Kármán length scale $L_{vK}$ as a sensor.
- In unsteady regions (large $L_{vK}$), the model reduces dissipation and allows the solution to develop resolved LES-like content without a fixed filter width.
- Acts as URANS in stable attached regions, transitions to LES-like behaviour in separated wakes.
- No explicit grid-dependency in the model term — a key advantage over DES.

---

## 3. LES Models

In LES, the velocity field is decomposed into resolved (filtered) and sub-grid scale (SGS) components: $\mathbf{u} = \bar{\mathbf{u}} + \mathbf{u}'$. SGS models provide the unclosed SGS stress tensor $\tau_{ij}^{sgs}$.

The filter width $\Delta$ in OpenFOAM defaults to the cube root of cell volume: $\Delta = V_{cell}^{1/3}$, but can be overridden.

### 3.1 Smagorinsky-Type Models

**Smagorinsky (`Smagorinsky`)**
- The classical (1963) algebraic SGS model: $\mu_{sgs} = \rho (C_s \Delta)^2 |\bar{S}|$.
- $C_s \approx 0.1$–$0.2$ depending on flow.
- Over-dissipative near walls; does not vanish correctly in laminar regions or at walls. Requires damping (e.g., Van Driest) near walls.

**Dynamic Smagorinsky (`dynamicKEqn` / `dynamicLagrangian`)**
- Germano's dynamic procedure (1991) computes $C_s$ locally and instantaneously using a test filter ($\hat{\Delta} = 2\Delta$).
- Satisfies correct near-wall asymptotic behaviour and adapts to local flow physics.
- OpenFOAM implements both Eulerian (cell-averaged) and Lagrangian averaging of the Germano identity to avoid numerical instability from negative $C_s$.

**WALE (`WALE` — Wall-Adapting Local Eddy-viscosity)**
- Nicoud & Ducros (1999): constructs $\mu_{sgs}$ from the square of the velocity gradient tensor's traceless symmetric part.
- Correct $y^3$ near-wall scaling without explicit damping.
- Zero SGS viscosity in laminar flow and solid-body rotation — a significant advantage.

**Sigma Model (`SigmaModel` — community/extended)**
- Based on singular values ($\sigma_1 \geq \sigma_2 \geq \sigma_3$) of the velocity gradient tensor.
- Correct near-wall and transitional behaviour; does not require test filtering.

### 3.2 One-Equation SGS Models

**One-equation eddy viscosity (`oneEqEddy`)**
- Solves a transport equation for SGS kinetic energy $k_{sgs}$: $\mu_{sgs} = C_k \rho \Delta k_{sgs}^{1/2}$.
- Retains history effects of SGS energy; more robust than algebraic models in transitional flows.

**Dynamic one-equation (`dynamicKEqn`)**
- Applies Germano's dynamic procedure to determine $C_k$ locally.
- More accurate than static one-equation; slightly more expensive.

### 3.3 Implicit LES / No-Model

**Implicit LES (`Smagorinsky` with $C_s = 0$ or `laminar` in LES context)**
- Relies entirely on numerical dissipation as a surrogate SGS model.
- Appropriate only with dissipation-controlled high-order schemes; not recommended with standard OpenFOAM finite-volume numerics.

---

## 4. Detached Eddy Simulation (DES) and Hybrid RANS-LES

DES uses RANS near walls and transitions to LES in separated regions by comparing the RANS length scale to the local grid size.

**DES based on Spalart–Allmaras (`SpalartAllmarasDES`)**
- Original DES97 formulation. The destruction term's length scale $\tilde{d}$ switches: $\tilde{d} = \min(d_w, C_{DES} \Delta)$.
- Suffers from Modelled Stress Depletion (MSD) and Grey Area at the RANS-LES interface.

**DES based on k–ω SST (`kOmegaSSTDES` / `kOmegaDES`)**
- Replaces the turbulent length scale in the $k$ destruction term: $L_t \to \tilde{L}_t = \min(L_t, C_{DES}\Delta)$.
- Inherits SST's better RANS performance over SA.

**Delayed DES (DDES)**
- Introduces a shielding function $f_d$ to protect attached boundary layers from premature LES switching.
- Available as `SpalartAllmarasDDES` and `kOmegaSSTDDES`.

**Improved DDES (IDDES)**
- Combines DDES with a Wall-Modelled LES (WMLES) branch, allowing LES operation inside the boundary layer when the grid is sufficiently refined.
- Available as `SpalartAllmarasIDDES` and `kOmegaSSTIDDES`.
- The preferred DES variant for modern high-fidelity industrial simulations.

---

## 5. Transition Models

**$\gamma$–$Re_\theta$ model (`gammaReTheta`, `kOmegaSSTLM`)**
- Langtry–Menter (2009) four-equation model: $k$, $\omega$, intermittency $\gamma$, and momentum-thickness Reynolds number $\tilde{Re}_\theta$.
- Handles bypass, natural, and separation-induced transition through empirical correlations for $Re_{\theta c}$ and $F_{length}$.
- Fully local — no line integrals required; compatible with parallel decomposition.
- Available as a standalone model or tightly coupled to SST (`kOmegaSSTLM`).

---

## 6. Wall Treatment

Wall functions in OpenFOAM are applied as boundary conditions on patches, completely independently of the volume model. This modularity is a distinguishing feature.

| Wall Function | Variable | Description |
|---|---|---|
| `kqRWallFunction` | $k$ | Standard $k$ wall function |
| `epsilonWallFunction` | $\varepsilon$ | Standard $\varepsilon$ wall function (Launder–Spalding) |
| `omegaWallFunction` | $\omega$ | Blended formulation (Menter, 2001) |
| `nutWallFunction` | $\nu_t$ | Standard log-law-based $\nu_t$ |
| `nutkWallFunction` | $\nu_t$ | Based on $k$, for use with $k$-based models |
| `nutLowReWallFunction` | $\nu_t$ | For low-Re / resolved wall ($y^+ < 1$) |
| `nutUSpaldingWallFunction` | $\nu_t$ | Spalding's unified law valid across all $y^+$ |
| `v2WallFunction` | $\overline{v^2}$ | For v2-f models |

The `nutUSpaldingWallFunction` is particularly powerful for meshes with mixed $y^+$ values, as Spalding's unified law blends viscous sublayer and log-law behaviour continuously.

---

## 7. Compressible Extensions

All models have compressible counterparts accessible via `compressible::RASModels` and `compressible::LESModels`. Key differences:

- **Favre averaging** replaces Reynolds averaging: $\tilde{\phi} = \overline{\rho \phi} / \bar{\rho}$.
- Density fluctuation effects are partially accounted for through Favre-averaged production terms.
- Compressible Smagorinsky and dynamic models account for density weighting in the SGS stress.
- Dilatation dissipation corrections (Sarkar, Zeman) can be patched in for high-Mach LES.

---

## 8. Buoyancy and Thermal Extensions

For buoyant flows (`buoyantPimpleFoam`, `buoyantSimpleFoam`):

- **Turbulent Prandtl number** approach: $q_i = -(\mu_t / Pr_t) \partial \tilde{T} / \partial x_i$ with $Pr_t \approx 0.85$.
- **k–ε with buoyancy production:** $G_b = -g_i (\mu_t / \rho \sigma_T) \partial \rho / \partial x_i$ added to the $k$ and $\varepsilon$ equations.
- Stable stratification can be handled by multiplying $G_b$ by $C_3$ (often set to 0 or 1 depending on flow orientation).

---

## 9. Model Selection Guidance

| Flow Regime | Recommended Model |
|---|---|
| Fully attached, mild separation | k–ω SST |
| Strongly separated / bluff bodies | k–ω SST SAS or IDDES |
| Transition-sensitive (turbomachinery, wind energy) | k–ω SST–LM ($\gamma$–$Re_\theta$) |
| Anisotropic / swirling flows | SSG RSM or Non-linear k–ε |
| External aerodynamics (mid-to-high fidelity) | SA-IDDES or SST-IDDES |
| Wall-bounded, high-accuracy heat transfer | LES with dynamic Smagorinsky or WALE |
| Near-wall resolution impractical (high Re) | IDDES / Wall-Modelled LES |
| Industrial internal flows, general-purpose | k–ω SST with `nutUSpaldingWallFunction` |

---

## 10. Practical Implementation Notes

- **Numerical schemes matter as much as the model.** LES requires bounded central differencing (`filteredLinear`, `limitedLinear`) for convection; RANS tolerates upwind-biased schemes (`linearUpwind`, `limitedLinear`).
- **Time-stepping in LES/DES:** Courant number should be $Co < 0.5$ (ideally near 0.2) for resolved turbulence. `pimpleFoam` with `nOuterCorrectors 1` is effectively equivalent to `pisoFoam`.
- **Mesh quality for LES:** Cell aspect ratios near walls should approach unity in the $y$–$z$ plane; $\Delta x^+ \approx 50$–100, $\Delta z^+ \approx 15$–30 for wall-resolved LES.
- **`turbulenceProperties` vs `RASProperties`:** The legacy `RASProperties`/`LESProperties` dictionaries were consolidated into `turbulenceProperties` from OpenFOAM v3+. Always check the version.
- **`printCoeffs on`:** Invaluable during development — dumps all active model coefficients to `log` at startup so coefficient overrides in the dictionary are confirmed.
- **Custom models** can be compiled as user libraries and linked via `libs ("libmyTurbulenceModels.so")` in `controlDict`, exploiting OpenFOAM's runtime selection table without touching the core source.

---

This framework makes OpenFOAM one of the most extensible turbulence modelling environments available, enabling everything from quick industrial RANS sweeps to rigorous scale-resolving simulations within a single, unified codebase.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 4.6
>
> Prompt: Provide a thorough description of the turbulence models available in OpenFOAM (Open Field Operation and Manipulation). This description is intended for a senior application engineer who specializes in computational fluid dynamics (CFD).
