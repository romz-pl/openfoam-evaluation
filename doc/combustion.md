# Combustion Models in OpenFOAM

OpenFOAM implements combustion modelling through a layered, object-oriented framework that decouples the thermochemistry, turbulence–chemistry interaction, and the flow solver. The following is a systematic treatment of the available model families, their theoretical basis, governing assumptions, and practical applicability.

---

## 1. Framework Architecture

Before discussing individual models, it is important to understand the structural hierarchy. Combustion models in OpenFOAM are housed within the `combustionModels` library and are templated against a `ReactionThermo` type (e.g., `psiReactionThermo`, `rhoReactionThermo`). This means every combustion model is tightly coupled to:

- A **thermodynamic package** (e.g., `janaf`, `polynomial` coefficient tables for Cp, h, s)
- A **transport model** (e.g., `sutherland`, `const`)
- A **mixture model** (e.g., `homogeneousMixture`, `multiComponentMixture`, `singleStepReactingMixture`)
- A **turbulence model** (via `turbulenceModel` pointer within the solver)

The actual call to combustion happens through the `Qdot()` (heat release rate) and `R()` (species source terms) methods, which are overridden in each derived class. Solvers such as `reactingFoam`, `rhoReactingFoam`, `fireFoam`, `XiFoam`, and `PDRFoam` each selects the appropriate model via runtime selection (`combustionModel::New()`).

---

## 2. Infinitely Fast Chemistry / Mixture Fraction Models

### 2.1 Burke–Schumann (Equilibrium) — `infinitelyFastChemistry`

This is the simplest non-trivial combustion model. It assumes that the chemical time scales are infinitely faster than the turbulent mixing time scales (Da → ∞), collapsing the entire thermochemistry onto the **mixture fraction Z** (Bilger definition). The flame structure is assumed to be in chemical equilibrium and the flame sheet exists at Z = Z_st.

- **Species and temperature** are algebraically related to Z via look-up tables derived from equilibrium calculations.
- **No transport equations** for individual species are solved; only one scalar (Z) is convected.
- **Source term**: The energy equation contains the heat release; species are post-processed.
- **Validity**: Diffusion flames at high Damköhler numbers (e.g., gas-turbine pilot flames, furnace combustion at high temperature), where extinction and re-ignition are not important.
- **Implementation note**: Used primarily in `diffusionFlamelet` context or as a limiting case; requires `singleStepReactingMixture`.

### 2.2 `singleStepCombustion`

A global one-step irreversible reaction (e.g., fuel + oxidiser → products) with an Arrhenius rate that is typically set very large to approximate infinitely fast chemistry. This is often the entry-level model in OpenFOAM tutorials (`reactingFoam` with `CH4 + 2O2 → CO2 + 2H2O`).

- Only **one reaction progress variable** is effectively needed.
- Useful for **initial solution development** or coarse feasibility studies.
- Cannot predict CO, NOx, or intermediate species accurately.

---

## 3. Finite-Rate Chemistry Models

### 3.1 Laminar Finite-Rate — `laminarFlameSpeed` / `PaSR` with unity Lewis

When the **laminar** keyword is selected or in conjunction with DNS-type simulations, OpenFOAM solves the full set of species transport equations:

$$\frac{\partial(\rho Y_k)}{\partial t} + \nabla \cdot (\rho \mathbf{u} Y_k) = \nabla \cdot (\rho D_k \nabla Y_k) + \dot{\omega}_k$$

where $\dot{\omega}_k$ is computed directly from Arrhenius expressions for each elementary reaction via **Chemkin-format** mechanism input (`chemkin` or `foamChemistryReader`). This is exact in laminar flows and serves as the foundation for turbulence–chemistry interaction closures.

### 3.2 `PaSR` — Partially Stirred Reactor

The PaSR model (Golovitchev / Chomiak) acknowledges that within a computational cell, mixing and reaction occur on different sub-grid time scales. The cell is conceptually divided into:
- A **reacting fraction** κ where mixing is complete and chemistry proceeds at the full laminar rate.
- A **non-reacting fraction** (1 − κ) representing poorly mixed fluid.

The effective source term is:

$$\dot{\omega}_k^{\text{eff}} = \kappa \cdot \dot{\omega}_k^{\text{lam}}$$

where:

$$\kappa = \frac{\tau_c}{\tau_c + \tau_{mix}}, \quad \tau_{mix} = C_{mix}\sqrt{\frac{\nu}{\varepsilon}}$$

- **τ_c** is the chemical time scale (estimated from the stiffest reaction), **τ_mix** is the sub-grid mixing time scale (Kolmogorov or turbulent), and **C_mix** is a model constant (typically 0.03–0.1).
- **Strength**: Retains the full mechanism without pre-tabulation; captures ignition, extinction, and slow chemistry effects.
- **Weakness**: High computational cost (ODE integration per cell per time step); stiff chemistry requires implicit or semi-implicit solvers (CVODE, SEULEX via `chemistryModel`).
- **Applicable solvers**: `reactingFoam`, `rhoReactingFoam`, `dieselFoam`.

### 3.3 `EDC` — Eddy Dissipation Concept (Magnussen)

The EDC model is one of the most widely used models for turbulent reacting flows with finite-rate chemistry. It is based on the concept that combustion occurs in fine structures (the "fine structure reactor") where turbulent energy cascades down to the Kolmogorov scale and energy dissipation drives mixing.

The **fine structure fraction** (volume of reacting zones):

$$\gamma^* = C_\gamma \left(\frac{\nu \varepsilon}{k^2}\right)^{1/4}$$

The **residence time** within fine structures:

$$\tau^* = C_\tau \left(\frac{\nu}{\varepsilon}\right)^{1/2}$$

The source term is:

$$\dot{\omega}_k = \frac{\rho (\gamma^*)^2}{\tau^* [1-(\gamma^*)^3]} (\tilde{Y}_k^* - Y_k)$$

where $\tilde{Y}_k^*$ is the species mass fraction after reaction within the fine structure (solved by a PSR/PFR sub-model).

- **Model constants**: C_γ = 2.1377, C_τ = 0.4082 (Magnussen, 1981); these are not universal and may require tuning.
- **Chemistry**: Full detailed or reduced mechanisms via Chemkin format, integrated using `chemistryModel` with CVODE.
- **Strength**: Physically motivated; applicable across a wide range of Damköhler numbers; effective for industrial burners, gas turbines, and IC engines.
- **Weakness**: The fine structure constants are empirical; computationally expensive with large mechanisms; sensitive to mesh resolution in the reaction zone.
- **OpenFOAM implementation**: `EDC` class within `combustionModels`; requires k-ε or k-ω SST turbulence model for the ε and k fields.

---

## 4. Eddy Dissipation Model Family

### 4.1 `eddyDissipationModel` (EDM — Magnussen & Hjertager, 1977)

The classic EDM closes the reaction rate by assuming it is limited by the **turbulent mixing rate**, not chemical kinetics. The reaction rate for fuel consumption is:

$$\dot{\omega}_F = -A \frac{\varepsilon}{k} \min\left(Y_F, \frac{Y_O}{s}, B \frac{Y_P}{1+s}\right)$$

where A ≈ 4.0 and B ≈ 0.5 are empirical constants, s is the stoichiometric oxidiser-to-fuel ratio, and Y_P is the product mass fraction.

- **No Arrhenius chemistry**: The model is purely mixing-limited and will predict combustion wherever fuel and oxidiser coexist — it cannot predict ignition delay, extinction, or slow chemistry phenomena.
- **Computationally inexpensive**: Single algebraic expression per cell.
- **Use cases**: Fully turbulent, high-temperature diffusion flames where Da >> 1; initial solution guesses; industrial furnace estimates.
- **Caution**: Significantly overestimates reaction rates in premixed zones and near walls; inappropriate for ignition/extinction studies or slow reactions (e.g., NOx formation).

### 4.2 `eddyDissipationDiffusionModel`

An extension that accounts for both diffusion and reaction limitations, blending the EDM rate with a diffusion-limited rate. Less commonly used in practice.

---

## 5. Flamelet-Based Models

Flamelet models exploit the concept that a turbulent flame, at sufficiently high Damköhler number, consists of an ensemble of locally laminar flame structures (flamelets) embedded in turbulent flow. All thermochemical information is pre-tabulated as a function of a small set of scalars.

### 5.1 Steady Laminar Flamelet Model — `steadyLaminarFlamelet`

Pre-computed laminar **counterflow diffusion flame** solutions are stored as functions of mixture fraction Z and scalar dissipation rate χ (which characterises local strain):

$$\chi = 2D|\nabla Z|^2$$

The flamelet library contains T, Y_k = f(Z, χ). In the CFD solver, only Z and its variance $\widetilde{Z''^2}$ (or g = $\widetilde{Z''^2}/\widetilde{Z}(1-\widetilde{Z})$) are transported. A presumed PDF (typically a β-PDF for Z) is used to integrate the flamelet solutions:

$$\tilde{\phi} = \int_0^1 \phi(Z, \chi_{st}) \tilde{P}(Z; \tilde{Z}, \widetilde{Z''^2}) \, dZ$$

- **Tabulation**: Flamelets are computed offline (e.g., with `flameD`, Cantera, or FlameMaster) and stored in `XiFoam`/`rhoReactingFoam` compatible format.
- **Strength**: Captures detailed chemistry cheaply; can include NOx, CO, soot precursors.
- **Weakness**: Assumes attached flamelets; breaks down near extinction (high χ_st), in strongly premixed regimes, or in recirculation zones with slow chemistry.
- **OpenFOAM solver**: `rhoReactingFoam` with `flameletModel` or `XiFoam` for premixed variant.

### 5.2 Unsteady Flamelet / Representative Interactive Flamelet (RIF)

Not natively available as a standalone model in standard OpenFOAM releases, but implementable via the `flameletModel` infrastructure. Solves the unsteady flamelet equations in mixture fraction space, coupled back to the CFD via χ field. Important for diesel autoignition and staged combustion.

---

## 6. Premixed and Partially Premixed Combustion Models

### 6.1 `XiFoam` — Wrinkling Factor Model (Coherent Flame Model variant)

`XiFoam` is the dedicated premixed/partially-premixed solver in OpenFOAM. It solves for the **regress variable b** (b = 1 in unburnt, b = 0 in burnt) and the **flame surface wrinkling factor Ξ**:

$$\frac{\partial(\rho b)}{\partial t} + \nabla \cdot (\rho \mathbf{u} b) = \nabla \cdot (\rho D \nabla b) - \rho_u S_u \Ξ |\nabla b|$$

The wrinkling factor Ξ accounts for the sub-grid flame surface area increase due to turbulence and is transported via a semi-empirical model (e.g., Weller's model):

$$\Ξ = \Ξ_{eq} = 1 + \frac{u'}{S_L}$$

- **Laminar flame speed S_L**: Supplied as a function of equivalence ratio, pressure, and temperature via `laminarFlameSpeedModels` (Gulders correlation, Metghalchi–Keck, or tabulated).
- **Sub-models**: `Weller`, `RaviPitsch`, `HuijnenVanVeen` wrinkling models; `standardPhiV`, `standardPhi` laminar flame speed correlations.
- **Use cases**: Spark-ignition engines, gas turbine premix burners, lean premixed combustors.
- **PDRFoam**: A variant for **Porosity/Distributed Resistance** modelling — important for explosion hazard assessment in congested geometries (offshore platforms, hydrogen safety).

### 6.2 `FSD` — Flame Surface Density Model

Not shipped as a named model in vanilla OpenFOAM but implementable via the `XiFoam` framework. Transports the flame surface density Σ directly:

$$\Sigma = \frac{dA_{flame}}{dV}$$

and closes the mean reaction rate as $\bar{\dot{\omega}} = \rho_u S_L \Sigma$.

### 6.3 `thickened Flame Model (TFM)`

Available in some distributions (notably foam-extend and specialized forks). The flame is artificially thickened by a factor F to resolve it on the LES mesh, while an efficiency function E compensates for the reduced wrinkling:

$$\dot{\omega}_{TFM} = \frac{E}{F} \dot{\omega}_{lam}$$

Critical for LES of premixed combustion where the flame thickness (O(0.1 mm)) is much smaller than the LES filter width.

---

## 7. Conditional Moment Closure (CMC)

CMC is not in the standard OpenFOAM release but has been implemented by research groups (notably Edinburgh, TU Delft). It transports conditional means of species and temperature $\langle Y_k | Z = \eta \rangle$, i.e., the expected value conditioned on a specific mixture fraction value η. This captures local extinction and re-ignition more accurately than flamelet models.

---

## 8. Tabulated Chemistry Approaches

### 8.1 `TDACChemistryModel` — Tabulation of Dynamic Adaptive Chemistry

TDAC is a chemistry acceleration framework rather than a combustion model per se. It operates on top of any finite-rate model (PaSR, EDC, laminar) and provides two acceleration strategies:

- **ISAT (In-Situ Adaptive Tabulation)**: Stores previously computed ODE solutions in a binary tree; retrieves solutions for similar states via linear interpolation with error control (Pope, 1997). Achieves speed-ups of 10–100× for complex mechanisms.
- **DAC (Dynamic Adaptive Chemistry)**: On-the-fly mechanism reduction — identifies active species for each computational cell at each time step and solves only the relevant subset of the ODE system. Typically reduces the mechanism by 30–70% in stiff regions.

TDAC is particularly powerful for diesel combustion with mechanisms of 50–200 species, where direct integration would be prohibitive.

### 8.2 FGM — Flamelet Generated Manifolds

FGM (van Oijen & de Goey, 2000) is available in several OpenFOAM variants (OpenFOAM+/ESI, foam-extend). A low-dimensional manifold is constructed in composition space using detailed flamelet solutions. The manifold is parameterized by:
- Mixture fraction Z (for partially premixed)
- Progress variable c (for premixed)

All thermochemical quantities are stored in a look-up table as $\phi = \phi(Z, c)$, with a β-PDF integration applied for turbulence–chemistry interaction. This is the industry-preferred method for premixed gas turbine combustion simulations.

---

## 9. Solid/Liquid Fuel Combustion

### 9.1 Coal / Biomass — `coalCombustion`

OpenFOAM's `coalChemistryFoam` handles:
- **Particle drying** (evaporation of moisture)
- **Devolatilization** (single or multiple competing rates — e.g., Kobayashi model)
- **Char combustion** (diffusion/kinetics limited surface reactions — e.g., Field model or Gibb model)

The gas-phase combustion (volatiles + char oxidation products) is handled by EDM or finite-rate models via the Lagrangian particle–gas coupling (`coalParcel`).

### 9.2 Liquid Spray Combustion — `sprayFoam` / `reactingParcelFoam`

`sprayFoam` combines:
- **Lagrangian droplet tracking** with breakup (KH-RT, TAB models), evaporation (Ranz–Marshall, Abramzon–Sirignano), and collision models.
- **Gas-phase combustion** handled via any combustion model compatible with `reactingFoam` (PaSR, EDC, or EDM).
- **Coupling**: Mass, momentum, and energy source terms from droplet evaporation feed the gas-phase combustion.

---

## 10. NOx and Pollutant Post-Processing

While not combustion models in the strict sense, OpenFOAM provides post-processing utilities and runtime function objects for:
- **Thermal NOx**: Zeldovich mechanism (extended) — implemented as a post-processing `thermalnox` function object or included in the mechanism.
- **Prompt NOx**: Fenimore mechanism — requires CH radical, which needs at least a reduced mechanism.
- **Soot**: Two-equation soot model (Moss–Brookes) available in `fireFoam`; moment methods in research distributions.
- **CO**: Accurately requires finite-rate chemistry (PaSR or EDC with a mechanism including CO ↔ CO₂ equilibrium); EDM will always predict zero CO in hot zones.

---

## 11. Selection Guide for Practitioners

| Flow Regime | Damköhler # | Recommended Model | Solver |
|---|---|---|---|
| Laminar diffusion flame | Any | Laminar finite-rate | `reactingFoam` |
| Turbulent diffusion (high Da) | >> 1 | EDM or steady flamelet | `rhoReactingFoam` |
| Turbulent diffusion (moderate Da) | ~ 1–10 | PaSR or EDC | `reactingFoam` |
| Turbulent diffusion (extinction) | < 1 | EDC + detailed chem or CMC | `reactingFoam` |
| Lean premixed (GT combustor) | >> 1 | XiFoam + Weller / FGM | `XiFoam` |
| SI engine | Variable | XiFoam + ignition model | `XiFoam` |
| Diesel HCCI | Low (ignition) | PaSR + TDAC/ISAT | `reactingFoam` |
| Spray combustion | Variable | PaSR/EDC + Lagrangian | `sprayFoam` |
| Fire/pool fire | Low–moderate | EDM or PaSR | `fireFoam` |
| Explosion hazard | Premixed | PDRFoam + Xi | `PDRFoam` |

---

## 12. Key Implementation Considerations

**ODE Solvers for Stiff Chemistry**: OpenFOAM integrates with CVODE (SUNDIALS library) for stiff ODE integration. The `odeChemistryModel` and `EulerImplicit` solvers handle different stiffness levels. For very large mechanisms (>100 species), TDAC/ISAT is effectively mandatory.

**Mechanism Format**: OpenFOAM natively reads Chemkin II/III format mechanisms (`chemkinToFoam` conversion utility) and can also use native `foamChemistryReader` format. Thermodynamic data follow the 7-coefficient NASA polynomial format.

**LES Combustion**: Most models described above are RANS-based. For LES, only PaSR (with τ_mix set to the SGS time scale) and the thickened flame model are physically consistent; EDM and EDC with RANS-calibrated constants are theoretically inconsistent in LES, though widely misused in practice. For DNS, the laminar finite-rate model is the only appropriate choice.

**Radiation Coupling**: For luminous flames and high-temperature applications, the combustion model must be coupled to a radiation model (`P1`, `fvDOM`, `viewFactor`) via the `Qr` field, which feeds back into the energy equation. This is critical for accurate temperature prediction in furnaces and fires.

**Mesh Requirements**: Reaction zones in premixed flames require y⁺-independent resolution down to the flame thickness δ_L; for typical hydrocarbon/air at atmospheric conditions, δ_L ~ 0.3–0.5 mm. This drives mesh sizes in practical LES of premixed combustion, and is one of the primary motivations for the thickened flame and FGM approaches.

---

This framework gives OpenFOAM considerable flexibility in combustion modelling, ranging from simple global-rate models suitable for rapid engineering estimates to detailed finite-rate chemistry with tabulation acceleration suitable for high-fidelity predictive simulations. The selection of an appropriate model requires careful assessment of the combustion regime (premixed vs. non-premixed, Da number), the quantities of interest (temperature field vs. CO/NOx prediction), and the available computational budget.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 4.6
>
> Prompt: Provide a thorough description of the combustion models available in OpenFOAM (Open Field Operation and Manipulation). This description is intended for a senior application engineer specializing in computational fluid dynamics (CFD).
