# Core Solver Categories

OpenFOAM organizes its solver library around the governing equation set and the physical regime being modeled. Here is a thorough technical walkthrough of each major category.

---

## Incompressible flow solvers

These target the incompressible Navier-Stokes equations where density is constant and the pressure-velocity coupling is the central numerical challenge.

**`simpleFoam`** is the workhorse for steady-state turbulent flow. It implements the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm with under-relaxation. Turbulence closure is modular — RANS models (k-ε, k-ω SST, Spalart-Allmaras) or LES/DES plug in via the `turbulenceProperties` dictionary without changing the solver binary. Suitable for internal aerodynamics, HVAC, and pipe networks.

**`pimpleFoam`** is the transient counterpart, merging PISO (Pressure-Implicit with Splitting of Operators) with SIMPLE outer corrections into the PIMPLE loop. The `nOuterCorrectors`, `nCorrectors`, and `nNonOrthogonalCorrectors` control levers let you tune stability vs. computational cost per timestep. It handles LES and DNS directly, and its large-timestep stability (Courant numbers above 1 are feasible) makes it preferred for unsteady RANS and scale-resolving simulations.

**`icoFoam`** is a stripped-down PISO solver for laminar, isothermal flow — primarily useful for validation, teaching, and benchmarking against analytical solutions.

**`potentialFoam`** solves the potential flow (irrotational, inviscid) equations to generate a divergence-free velocity field for initializing viscous solvers, a technique that dramatically accelerates RANS convergence on complex geometries.

---

## Compressible flow solvers

Here density couples to the pressure and energy fields, and the equation of state is explicit.

**`rhoPimpleFoam`** and **`rhoSimpleFoam`** are the compressible analogues of their incompressible counterparts, using the density-based formulation. They solve the full set of compressible Navier-Stokes equations with the ideal gas law (or real-gas options) and an energy transport equation. They cover the subsonic-to-transonic regime effectively.

**`sonicFoam`** targets sonic and supersonic transient flows where strong shocks are present. It employs first-order temporal discretization to maintain stability near discontinuities and is typically paired with Gauss limitedLinear or vanLeer schemes to avoid spurious oscillations across shocks.

**`rhoCentralFoam`** implements a density-based, central-upwind (Kurganov-Tadmor) scheme designed explicitly for high-speed compressible flow. It handles strong shocks, expansion fans, and contact discontinuities with greater fidelity than pressure-based segregated approaches and is the standard choice for hypersonic and ballistic applications.

**`chtMultiRegionFoam`** — although primarily known as a heat transfer solver — belongs here conceptually because it solves the compressible energy equation in coupled fluid and solid regions simultaneously, making it central to conjugate heat transfer analysis in turbomachinery and electronics cooling.

---

## Multiphase flow solvers

These handle flows involving two or more immiscible phases, using distinct interface-tracking or interface-capturing strategies.

**`interFoam`** (and its turbulent extension **`interFoam`** with the `kOmegaSST` model active) implements the Volume of Fluid (VoF) method. The phase-fraction transport equation for α is solved with the interface-compression term (MULES — Multidimensional Universal Limiter for Explicit Solution) to keep the interface sharp without algebraic diffusion. This is the canonical solver for free-surface flows: sloshing, wave propagation, dam break, and coastal engineering.

**`interIsoFoam`** introduces the isoAdvector algorithm for VoF advection, providing a geometrically exact reconstruction of the interface on arbitrary unstructured meshes — a significant accuracy gain over MULES for highly distorted cells.

**`compressibleInterFoam`** extends VoF to compressible two-phase systems, solving individual phase energy and mass equations — essential for cavitation and flash evaporation problems.

**`reactingTwoPhaseEulerFoam`** and **`twoPhaseEulerFoam`** implement the Euler-Euler two-fluid model, where both phases are treated as interpenetrating continua with separate velocity fields and exchange terms for mass, momentum, and energy. This is the appropriate framework for dispersed flows (bubbly flow, fluidized beds, sprays) where individual interface resolution is computationally intractable.

**`driftFluxFoam`** employs the drift-flux mixture model — a simplification between the mixture model and full two-fluid, effective when the relative velocity between phases can be expressed algebraically via a slip relation.

---

## Heat transfer and buoyancy-driven flow

**`buoyantSimpleFoam`** and **`buoyantPimpleFoam`** solve the Boussinesq-approximated or fully compressible buoyant flow equations. They add a gravitational body force term proportional to local density deviation and solve an enthalpy or temperature transport equation. These are the standard solvers for natural convection, building thermal comfort analysis, and fire safety CFD (coupled with radiation models: P1, fvDOM, viewFactor).

**`chtMultiRegionFoam`** (mentioned above) handles solid conduction regions — solved via the Laplace equation — coupled to fluid regions through matched thermal boundary conditions at region interfaces, with automatic interface heat flux conservation.

---

## Combustion and reacting flow solvers

This category is physically the most complex, coupling fluid mechanics to chemical kinetics and multi-species transport.

**`reactingFoam`** is the base solver for chemically reacting flows. It integrates multi-species mass fractions, the energy equation with enthalpy of formation, and chemical kinetics via the OpenFOAM `chemistryModel` layer. Stiff ODE integrators (VODE, Seulex) handle the chemistry sub-step. Suitable for laminar flames and RANS combustion with closure models.

**`fireFoam`** extends `reactingFoam` for fire dynamics: it adds thermal radiation (fvDOM), soot transport, pyrolysis boundary conditions, and buoyancy-driven flow — making it analogous to FDS but with full OpenFOAM meshing flexibility.

**`PDRFoam`** (Porosity/Distributed Resistance) targets gas explosion scenarios in congested industrial geometries, where obstacles below grid resolution are parameterized via a distributed porosity/resistance field rather than being resolved explicitly.

**`XiFoam`** and **`engineFoam`** address premixed and partially premixed turbulent combustion using the flame wrinkling factor Ξ as the combustion progress variable, with particular application to internal combustion engines. `engineFoam` adds dynamic mesh motion (piston/valve motion via prescribed kinematics).

---

## Lagrangian particle and spray solvers

**`reactingParcelFoam`**, **`sprayFoam`**, and related solvers implement the Euler-Lagrange Discrete Phase Model (DPM). A Eulerian carrier phase (gas) is solved on the mesh; parcels of particles/droplets are tracked individually via Newton's law with drag, heat/mass exchange, and breakup models (TAB, KHRT, Reitz-Diwakar). `sprayFoam` includes evaporation, secondary breakup, and droplet-wall interaction specifically tuned for fuel sprays in engines and gas turbine combustors.

---

## Electromagnetics and MHD

**`mhdFoam`** solves the incompressible magnetohydrodynamics equations — the Navier-Stokes system augmented with the Lorentz body force and Faraday induction equation — for conducting fluids in magnetic fields. Applications include liquid metal cooling of fusion blankets and electromagnetic stirring in continuous casting.

---

## Solid mechanics and FSI

**`solidDisplacementFoam`** and **`solidEquilibriumDisplacementFoam`** solve linear elastic solid mechanics using the finite volume method — an unusual but consistent approach within the OpenFOAM paradigm. They compute displacement fields under applied loads, supporting thermal stresses and pressure boundary conditions.

These can be loosely coupled to fluid solvers for Fluid-Structure Interaction (FSI) through explicit or partitioned schemes, though tightly coupled FSI in OpenFOAM typically requires third-party coupling libraries (preCICE, OpenFSI).

---

## DNS and direct LES

While not a separate binary in all cases, OpenFOAM supports Direct Numerical Simulation via `pimpleFoam` or `icoFoam` with turbulence models switched off entirely and fine enough mesh resolution (Kolmogorov-scale resolved). The `channelFoam` solver is purpose-built for DNS of turbulent channel flow, automatically computing and writing turbulence statistics and applying cyclic+driving pressure-gradient boundary conditions.

---

Now here is a structured map of these categories and their principal solvers:

![OpenFOAM — solver taxonomy by physical regime](./solvers.png)

---

## Cross-cutting architectural notes relevant to solver selection

**Turbulence model portability.** The `turbulenceModel` abstraction layer means that RANS, LES, and DES closures are not baked into solver binaries. Any solver that calls `turbulence->correct()` can switch between k-ω SST, dynamic Smagorinsky, or wall-adapting local eddy-viscosity (WALE) without recompilation — a critical design decision that distinguishes OpenFOAM from legacy codes where turbulence was solver-coupled.

**Dynamic mesh.** Several solvers (`pimpleFoam`, `engineFoam`, `interFoam`) support AMI (Arbitrary Mesh Interface) and topoChangerMesh for sliding interfaces and topological remeshing. This is how rotating machinery (MRF steady, sliding mesh transient), piston motion, and free-surface mesh deformation are handled natively.

**fvOptions framework.** Solvers can be augmented at runtime with source terms — momentum porosity, Coriolis rotation, actuator disk models, heat sources — via the `fvOptions` mechanism, decoupling physical modeling from solver code.

**Parallel scalability.** All production solvers use PETSc-style domain decomposition via `decomposePar`. The `GAMG` (Geometric-Algebraic Multi-Grid) preconditioner in the pressure equation is the dominant solver for large incompressible cases and scales to tens of thousands of MPI ranks when properly configured with agglomeration settings matched to mesh topology.

**Solver selection heuristic for a senior practitioner:** the governing equation set (incompressible vs. compressible, single vs. multiphase) determines the category; the temporal character (steady vs. transient) and the Reynolds number regime determine which binary; turbulence closure, radiation, and reactions are then layered on via the runtime dictionary system rather than by changing the solver.


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 4.6
>
> Prompt: Provide a thorough description of the core solver categories available in OpenFOAM (Open Field Operation and Manipulation). This description is intended for a senior application engineer specializing in computational fluid dynamics (CFD).
