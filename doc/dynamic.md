# Dynamic Mesh Capabilities in OpenFOAM

OpenFOAM provides a rich, extensible framework for dynamic mesh operations, implemented primarily through the `dynamicFvMesh` library hierarchy. What follows is a comprehensive technical overview structured around the major capability domains.

---

## 1. Foundational Architecture

The dynamic mesh system is rooted in the abstract base class `dynamicFvMesh`, which inherits from `fvMesh`. All dynamic mesh variants override the `update()` method, which is called once per time step (or sub-step in some formulations) prior to solving the governing equations. The mesh motion velocity **U**_mesh is computed internally and passed to the ALE (Arbitrary Lagrangian-Eulerian) flux correction via the Geometric Conservation Law (GCL) term in the transport equations.

The GCL is enforced through the `fvc::ddt` and `fvc::div` operators acting on `phi` (the face flux), ensuring that a uniform field advected by mesh motion remains exactly conserved — a non-trivial requirement that OpenFOAM satisfies by default through `meshPhi` tracking.

---

## 2. Mesh Motion Solvers (`motionSolver` Hierarchy)

The `motionSolver` framework decouples **how the mesh moves** from **what causes it to move**. Solvers are specified in `constant/dynamicMeshDict`.

### 2.1 Laplacian / Diffusion-Based Solvers

The most widely used solvers for smooth, topology-preserving deformation:

- **`displacementLaplacian`**: Solves a Laplace equation for the point displacement field **d**:

  ∇·(γ ∇**d**) = **0**

  The diffusivity field γ can be set to `uniform`, `inverseDistance`, `inverseVolume`, or `quadratic`, with higher diffusivity near boundaries of interest to concentrate deformation there and preserve cell quality in the far field.

- **`velocityLaplacian`**: The velocity-formulation equivalent, better suited when boundary conditions are expressed as velocities rather than displacements.

- **`displacementComponentLaplacian`**: Applies the Laplacian diffusion independently per spatial component — useful for cases with strongly directional motion constraints.

### 2.2 Radial Basis Function (RBF) Interpolation

Available via `displacementSBRStress` and third-party or extended RBF frameworks. In this approach, control points (typically on moving boundaries) have prescribed displacements which are globally interpolated to all mesh points using a kernel function (e.g., compact support radial basis functions). This is highly accurate for large rigid-body motions of embedded surfaces but scales poorly (O(N²) to O(N³)) for large point counts without truncated support kernels.

### 2.3 Solid-Body Motion Solvers

For prescribed rigid-body kinematics without mesh deformation:

- **`solidBodyMotionFvMesh`**: Applies a pure rigid-body transformation (translation + rotation) to the entire mesh or subzones. Motion functions include:
  - `linearMotion`
  - `rotatingMotion`
  - `oscillatingLinearMotion`
  - `oscillatingRotatingMotion`
  - `SDA` (Ship motion for seakeeping)
  - `tabulated6DoFMotion` (lookup table)

- **`multiSolidBodyMotionFvMesh`**: Extends `solidBodyMotionFvMesh` to allow different zones within a single mesh to undergo independent prescribed rigid motions. Useful in multi-component machinery simulations.

### 2.4 Linear Elasticity / Stress-Based Solvers

- **`displacementSBRStress`** (Stress-Based Relaxation): Solves a pseudo-linear-elasticity equation to compute mesh point displacements, providing better volume preservation under large deformations compared to pure Laplacian smoothing. The governing equation includes a Lamé-type stress tensor:

  ∇·(μ(∇**d** + (∇**d**)ᵀ) + λ tr(∇**d**)**I**) = **0**

  Coefficients μ and λ are chosen heuristically to maximize mesh quality.

---

## 3. Mesh Topology Change: Refinement & Coarsening

This is where OpenFOAM's capabilities become particularly powerful, enabling **run-time adaptive mesh refinement (AMR)**.

### 3.1 `dynamicRefineFvMesh`

The workhorse AMR class. At each refinement interval it:

1. Evaluates a user-specified **refinement field** (any scalar field — typically a gradient, vorticity magnitude, interface indicator, or Mach number).
2. Applies `hexRef8` (isotropic 2:1 refinement in 3D, splitting one hex into 8) or `hexRef4` (2D, 1→4 quads) where the field exceeds `upperRefineLevel`.
3. **Unrefinement (coarsening)**: Collapses previously split cells when the field drops below `lowerRefineLevel`, subject to the constraint that only cells whose entire split family is eligible can be merged back. This constraint is managed by the `refinementHistory` object.
4. Balances the load across MPI ranks via `fvMeshDistribute` after each topology change.

Key parameters in `dynamicMeshDict`:
```
maxRefinement       3;        // max split levels per cell
nBufferLayers       2;        // layers of transitional cells around refined zone
refineInterval      5;        // time steps between refinement passes
dumpLevel           false;
```

The `refinementHistory` object serializes to disk, enabling restart of AMR simulations.

### 3.2 `dynamicMultiMotionSolverFvMesh` + AMR

These can be combined in a single simulation: boundary motion drives Laplacian mesh deformation, while AMR simultaneously adapts the bulk mesh topology. The interaction requires careful ordering of the `update()` calls.

### 3.3 Snappy-based Dynamic Refinement

`snappyHexMesh` geometry-driven refinement is a pre-processing step, but its data structures (`refinementSurfaces`, `refinementRegions`) conceptually underpin the AMR machinery — both rely on `hexRef8` under the hood.

---

## 4. Sliding Meshes and Arbitrary Mesh Interfaces (AMI)

For relative motion between non-conforming mesh regions (rotors, valves, pistons) without remeshing.

### 4.1 `cyclicAMI` Patch Pairs

The Arbitrary Mesh Interface uses weighted interpolation between non-conformal patch pairs. The interpolation weights are computed geometrically by projecting source face polygons onto target faces:

- **`faceAreaWeightAMI`**: Most robust; weights proportional to projected face area overlap.
- **`partialFaceAreaWeightAMI`**: Handles partially overlapping patches.
- **`sweptFaceAreaWeightAMI`**: Swept-area correction for fast-moving interfaces, reducing interpolation error in time-accurate simulations.

For rotating machinery:
- The rotor patch rotates each time step by Δθ = ω·Δt.
- AMI weights are **recomputed every time step** (or cached with incremental updates in newer releases).
- Both `MRF` (frozen-rotor) and transient sliding mesh approaches use AMI as the interpolation backbone.

### 4.2 `cyclicACMI` (Conditionally Coupled AMI)

Extends AMI to handle **partial opening/closing** interfaces — the primary tool for valve simulations. Cells near a valve seat that would become degenerate are instead blocked via a masking function, and flux is conditionally coupled through the ACMI patch:
- When the valve is open: full AMI coupling.
- When closed: the patch behaves as a wall (`blockage` field → 1).
- Partial open states are handled by intermediate masking values, preserving mass conservation.

### 4.3 Overset (Chimera) Mesh — `oversetFvMesh`

Available since OpenFOAM v1806 (ESI) and in the Foundation version via extensions:

- Background and overset component meshes are defined separately.
- The `oversetFvMesh` constructs donor-acceptor relationships between overlapping regions at run time.
- Interpolation stencils (typically tri-linear or radial basis) transfer solution data across the overlap zone.
- Cells in the background mesh that fall inside the overset body are **blanked** (removed from the active domain via the `cellTypes` field).
- Particularly powerful for: underwater vehicle manoeuvring, store separation, multi-body hydrodynamics.

Limitations: the overset approach requires careful attention to overlap zone sizing (generally ≥2 cell layers), and the interpolation introduces a small conservation error that must be monitored.

---

## 5. Topological Mesh Changes: Layer Addition/Removal and Mesh Stitching

### 5.1 `layerAdditionRemoval`

Available via the `topoChangerFvMesh` family. This capability adds or removes entire **cell layers** at a specified face zone when that zone crosses a displacement threshold:

- **Addition**: A new layer of cells is inserted when the gap between moving and static surfaces would otherwise cause cell collapse.
- **Removal**: A layer is deleted when cells become too thin, preventing aspect ratio degradation.

Commonly used in: piston-cylinder simulations, bellows, and peristaltic pump geometries. Requires a structured layer topology in the region of addition/removal.

### 5.2 `slidingInterface`

An older mechanism (predating AMI) that stitches and un-stitches mesh patches along a sliding interface. Largely superseded by AMI/ACMI but still present in the codebase and useful for certain legacy workflows.

---

## 6. Six Degree-of-Freedom (6-DOF) Fluid-Structure Interaction

### 6.1 `sixDoFRigidBodyMotion`

Solves Newton-Euler equations of motion for a rigid body immersed in (or floating on) a fluid:

- **Forces**: Integrated from the pressure and viscous stress fields over the body surface patch(es).
- **Torques**: Computed about the body's centre of mass.
- **Constraints**: Translational and rotational DOFs can be selectively freed or locked. Available constraints include `plane`, `axis`, `line`, `fixedPoint`, and custom user-defined constraints.
- **Restraints**: Spring-damper, linear spring, linearAxialAngularSpring, mooring-line models (for offshore applications).
- **Integration**: Symplectic (leapfrog), Crank-Nicolson, or Newmark-β time integrators. The symplectic integrator is the default — second-order accurate and energy-conservative for undamped systems.

### 6.2 Coupling to Mesh Motion

The computed body displacement/rotation at each time step is passed as the boundary condition to the chosen `motionSolver` (typically `displacementLaplacian` or `solidBodyMotionFvMesh`), which then deforms the surrounding mesh accordingly. This creates a tightly coupled FSI loop within a single OpenFOAM solver.

---

## 7. Mesh Quality Control and Enforcement

Regardless of the motion mechanism, deforming meshes must satisfy minimum quality thresholds. OpenFOAM provides several safeguards:

- **`meshQualityDict` checks**: Maximum non-orthogonality, maximum skewness, minimum face weight, and minimum cell determinant are evaluated after each mesh update. Violations trigger warnings and can optionally abort the run.
- **`pointSmoother`**: Post-motion point smoothing passes to recover quality metrics (laplacian, centroidal, or constrained surface smoothing).
- **Relaxation of mesh motion**: The motion solver can under-relax the computed displacements to prevent single-step quality deterioration.
- **`checkMesh` at runtime**: Dynamic mesh solvers internally call the mesh checking routines; the `writeInterval` parameter for mesh quality can be set separately from field output.

---

## 8. Parallel Decomposition and Dynamic Load Balancing

Dynamic mesh changes — particularly AMR and topology changes — alter cell counts per MPI rank, causing load imbalance. OpenFOAM handles this via:

- **`fvMeshDistribute`**: Redistributes cells among MPI ranks post-topology change using the SCOTCH, METIS, or hierarchical decomposition algorithms.
- **`redistributePar`**: A utility (also callable as a library function at runtime) that can redecompose a running simulation.
- AMR simulations should specify `refinementInterval` large enough to amortize the redistribution cost.

---

## 9. Solver Integration and the ALE Formulation

All dynamic mesh solvers in OpenFOAM use the ALE form of the transport equations. The key modification to the standard conservation law is the **mesh flux** `meshPhi`, defined as:

φ_mesh = ∫_face **U**_mesh · d**S**

The convective flux for any transported quantity becomes:
φ_corrected = φ_fluid − φ_mesh

This correction is applied automatically in the `fvc::ddt` schemes (`Euler`, `backward`, `CrankNicolson`) when the mesh is detected as moving, requiring no user intervention beyond selecting a dynamic mesh class.

---

## 10. Available Dynamic Mesh Solvers (Application-Level)

| Solver | Dynamic Mesh Capability |
|---|---|
| `pimpleFoam` | Full dynamic mesh support (most common for external FSI/AMR) |
| `interFoam` | VOF + dynamic mesh (waves, floating bodies) |
| `rhoPimpleFoam` | Compressible + dynamic mesh |
| `overInterDyMFoam` | Overset + VOF (ship hydrodynamics, wave energy) |
| `mpPimpleFoam` | Multi-physics with dynamic mesh |
| `engineFoam` | Specialised for IC engines (topological layer add/remove) |
| `SRFPimpleFoam` | Single rotating frame variant |

---

## 11. Practical Workflow Considerations for Senior Engineers

- **Laplacian vs. SBRStress**: For moderate deformations (< ~20% of mesh scale), `displacementLaplacian` with `inverseDistance` diffusivity is usually sufficient and fast. For large deformations, switch to `displacementSBRStress`.
- **Time-step control**: Dynamic mesh simulations are sensitive to the Courant number of the mesh motion itself. A Courant number based on mesh velocity (mesh Co) should be monitored alongside the fluid Co.
- **AMI recomputation cost**: At high mesh counts, AMI weight recomputation dominates wall-clock time. Cache strategies and reducing the `matchTolerance` are key tuning levers.
- **Overset hole-cutting**: The blanking algorithm is sensitive to the overlap zone geometry; always verify the `cellTypes` field (0=calculated, 1=interpolated, 2=blanked) in post-processing before interpreting flow results.
- **Restart fidelity**: Both `refinementHistory` (AMR) and `sixDoFRigidBodyMotionState` must be present in the time directory for a fully consistent restart.
- **Testing mesh motion**: Run with a `zeroField` fluid solver first to validate mesh kinematics in isolation before engaging the full CFD solve — a common and time-saving debugging practice.

---

OpenFOAM's dynamic mesh framework is among the most complete open-source offerings available, offering a continuum from simple rigid-body motion to full AMR-coupled FSI with overset capability. Its object-oriented design means that new motion solvers, refinement criteria, and constraint models can be introduced as runtime-selectable classes without modifying the core solver, making it highly amenable to research and advanced industrial workflows alike.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 4.6
>
> Prompt: Provide a thorough description of the dynamic mesh capabilities available in OpenFOAM (Open Field Operation and Manipulation). This description is intended for a senior application engineer specializing in computational fluid dynamics (CFD).
