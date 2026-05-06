# Mesh Infrastructure

## 1. Foundational Philosophy: The Finite Volume Method on Arbitrary Polyhedral Meshes

OpenFOAM's mesh infrastructure is built around a **fully unstructured, polyhedral finite volume** discretization. Unlike solvers that restrict topology to hexahedral or tetrahedral elements, OpenFOAM treats any convex (or mildly non-convex) polyhedron as a valid cell. This generality is not cosmetic — it permeates every layer of the codebase, from storage layout to gradient reconstruction to parallelism. The canonical data structure is the **face-based (face-addressing) mesh**, sometimes called the *fvMesh*, which unifies all mesh operations under a single owner–neighbour topology.

---

## 2. Core Data Structures

### 2.1 `polyMesh` — The Topological Foundation

All mesh classes ultimately inherit from `polyMesh`, which stores the raw combinatorial topology:

| Data Array | Type | Description |
|---|---|---|
| `points_` | `pointField` | $(x,y,z)$ coordinates of all mesh vertices |
| `faces_` | `faceList` | Ordered node-index lists for every face (internal + boundary) |
| `owner_` | `labelList` | Index of the owner cell for every face |
| `neighbour_` | `labelList` | Index of the neighbour cell for every *internal* face |
| `boundary_` | `polyBoundaryMesh` | Ordered list of patch descriptors |

**Ordering convention (critical for performance):** Internal faces come first, indexed `[0, nInternalFaces)`. Boundary patch faces follow contiguously, so patch $p$ spans `[patch.start(), patch.start()+patch.size())`. This contiguous layout enables vectorised boundary loops without pointer chasing.

The face normal orientation convention is:

$$\hat{n}_f \text{ points from owner to neighbour (internal) or outward (boundary)}$$

All flux sign conventions across the entire solver ecosystem flow from this single rule.

### 2.2 `fvMesh` — The Finite Volume Layer

`fvMesh` extends `polyMesh` with precomputed geometric quantities consumed by discretisation schemes:

- **`Sf()`** — face area vectors $\mathbf{S}_f = |A_f|\hat{n}_f$
- **`magSf()`** — face area magnitudes $|A_f|$
- **`V()`** — cell volumes $V_P$
- **`C()`** — cell centroid positions $\mathbf{x}_P$
- **`Cf()`** — face centroid positions $\mathbf{x}_f$
- **`delta()`** — the non-orthogonality correction vector $\mathbf{d}_{PN}$ (owner-to-neighbour centroid vector minus the face-normal component)
- **`weights()`** — linear interpolation weights $w_f$ based on the ratio of distances from face centroid to owner and neighbour centroids

These are stored as `GeometricField` objects, updating automatically on mesh motion via the `movePoints()` / `updateMesh()` call chain.

### 2.3 `primitiveMesh` — Geometry Kernel

Between `polyMesh` and `fvMesh` sits `primitiveMesh`, which computes all geometry from the raw `points_` and `faces_` arrays without knowledge of boundary conditions. It provides:
- Cell-to-face addressing (`cellFaces()`)
- Point-to-cell addressing (`pointCells()`)
- Point-to-face addressing (`pointFaces()`)
- Edge lists and edge addressing (on demand, since edges are rarely needed in FVM)

These addressing lists are computed lazily and cached — a key performance pattern throughout the mesh layer.

---

## 3. Boundary Mesh Architecture

### 3.1 `polyBoundaryMesh` and Patch Hierarchy

Boundary conditions in OpenFOAM are mesh-level concepts, not field-level ones. `polyBoundaryMesh` stores an ordered list of `polyPatch` objects. The patch hierarchy is:

```
polyPatch
└── fvPatch
    ├── basicFvPatch (wall, empty, symmetry, wedge, cyclic)
    ├── coupledFvPatch
    │   ├── processorFvPatch       ← inter-process halo
    │   ├── cyclicFvPatch          ← periodic, possibly with AMI
    │   └── cyclicAMIFvPatch       ← Arbitrary Mesh Interface
    └── mappedPatchBase            ← non-conformal sampling
```

Each patch knows its face range in the global face list, its physical type (for BC dispatch), and carries a `patchField<T>` pointer for each registered field.

### 3.2 Coupled Patches and AMI

**Cyclic patches** require face-to-face pairing; when the two cyclic halves are geometrically conformal (one-to-one face matching), this is trivial. When the surfaces are non-conformal (e.g., rotational sliding interfaces, different mesh refinement on each side), the **Arbitrary Mesh Interface (AMI)** machinery performs a surface intersection to compute interpolation weights. AMI supports:
- `faceAreaWeightAMI` — area-weighted for conservative flux transfer
- `partialFaceAreaWeightAMI` — for partially overlapping interfaces
- `swakAMI` and user extensions via the `AMIMethod` factory

For **Moving Reference Frame (MRF)** and **sliding mesh** (GGI, sliding AMI), the interpolation weights are rebuilt at each time step inside `AMIInterpolation::update()`.

---

## 4. Non-Orthogonality and Mesh Quality Metrics

Because OpenFOAM operates on arbitrary polyhedra, the solver must quantify and correct for mesh defects at runtime. The key quality metric is **non-orthogonality** $\phi$:

$$\phi = \arccos\!\left(\frac{\mathbf{d}_{PN} \cdot \mathbf{S}_f}{|\mathbf{d}_{PN}||\mathbf{S}_f|}\right)$$

When $\phi > 0$, the face normal does not align with the cell-centroid connector, requiring a **non-orthogonality correction** to the Laplacian operator. The face gradient is decomposed as:

$$(\nabla \phi)_f \cdot \mathbf{S}_f = \underbrace{\frac{\phi_N - \phi_P}{|\mathbf{d}_{PN}|}|\mathbf{S}_f|}_{\text{orthogonal (implicit)}} + \underbrace{(\nabla \phi)_f \cdot \mathbf{k}_f}_{\text{correction (explicit)}}$$

where $\mathbf{k}_f$ is the non-orthogonality correction vector. OpenFOAM supports three decomposition strategies selectable in `fvSolution`:

| Strategy | Decomposition of $\mathbf{S}_f$ | Stability |
|---|---|---|
| `corrected` | $\mathbf{S}_f = \overline{\mathbf{d}}_{PN}|\mathbf{S}_f| + \mathbf{k}_f$ | Stable for $\phi < 70°$ |
| `limited` | Clips correction magnitude | Robust for $\phi < 85°$ |
| `uncorrected` | Ignores $\mathbf{k}_f$ | Only for near-orthogonal meshes |
| `orthogonal` | Pure orthogonal (no correction) | Structured Cartesian only |

The `nNonOrthogonalCorrectors` keyword in `fvSolution` controls the number of explicit correction sub-iterations within each pressure-velocity coupling loop.

Other quality metrics computed by `checkMesh` (and accessible programmatically via `polyMeshTools`):

- **Skewness** — distance of face interpolation point from face centroid, normalised by cell-to-cell distance
- **Face flatness** — deviation of face vertices from the best-fit plane (relevant for highly non-planar polyhedral faces)
- **Aspect ratio** — maximum/minimum edge-length ratio within a cell
- **Cell determinant** — measure of cell volume relative to its bounding box

---

## 5. Mesh Motion and Dynamic Meshes

### 5.1 `dynamicFvMesh` Hierarchy

For moving boundary problems (FSI, 6-DOF, free surface, sliding mesh), OpenFOAM provides the `dynamicFvMesh` base class with a virtual `update()` hook called at the beginning of each time step:

```
dynamicFvMesh
├── dynamicMotionSolverFvMesh   ← general motion solver dispatch
│   ├── velocityLaplacian       ← Laplacian smoothing of mesh velocity
│   ├── displacementLaplacian   ← point displacement formulation
│   ├── displacementSBRStress   ← solid-body rotation stress analogy
│   └── RBFMeshMotionSolver     ← Radial Basis Function interpolation
├── solidBodyMotionFvMesh       ← prescribed rigid-body motion (6-DOF)
└── rawTopoChangerFvMesh        ← topological changes (layer addition/removal)
```

### 5.2 Geometric Conservation Law (GCL)

On moving meshes, the Geometric Conservation Law mandates that the discrete mesh motion is consistent with the transported volume:

$$\frac{d V_P}{dt} = \oint_{\partial V_P} \mathbf{u}_g \cdot d\mathbf{S}$$

OpenFOAM enforces GCL through the **swept-volume flux** $\phi_g$ on moving faces, computed in `fvMesh::phi()` and added to the continuity equation. Failure to maintain GCL on rapidly deforming meshes introduces spurious mass sources — a common source of divergence in FSI simulations.

### 5.3 Topological Changes: `polyTopoChanger`

When mesh motion is insufficient (e.g., large displacements, layer addition in engine simulations, dynamic refinement), `polyTopoChanger` performs in-situ topological modifications:

- **`layerAdditionRemoval`** — adds/removes cell layers at a sliding interface (engine stroke simulation)
- **`slidingInterface`** — merges/demerges patches at a sliding plane
- **`refinementHistory`** — tracks parent–child relationships for AMR

After each topological event, all addressing lists, geometric fields, and boundary conditions are rebuilt via `updateMesh(mapPolyMesh)`.

---

## 6. Adaptive Mesh Refinement (AMR)

### 6.1 `hexRef8` — Isotropic Octree Refinement

The standard AMR engine in OpenFOAM is `hexRef8`, which refines hexahedral cells by splitting each cell into 8 children (2×2×2). The algorithm:

1. **Mark cells** for refinement/unrefinement based on a field-based criterion (e.g., gradient magnitude, iso-surface proximity)
2. **Balance** the refinement level across neighbours (maximum 2:1 level jump enforced to avoid hanging nodes)
3. **Split** faces and edges, introducing **hanging nodes** at refinement boundaries
4. **Project** hanging nodes onto their parent face to maintain geometric consistency
5. **Update** `refinementHistory` for consistent unrefinement

The constraint that refinement produces *valid polyhedra* means split faces on coarser neighbours become non-planar pentagons/hexagons — entirely legal in OpenFOAM's polyhedral framework.

### 6.2 `snappyHexMesh` Refinement vs. Runtime AMR

It is important to distinguish between:

- **`snappyHexMesh` (offline)** — surface-based refinement during mesh generation, producing the final polyhedral mesh including castellated cells, snapped surface cells, and boundary layers
- **`dynamicRefineFvMesh` / `dynamicMultiLevelRefine` (runtime)** — field-driven refinement during the simulation, using `hexRef8` under the hood

---

## 7. Parallel Mesh Decomposition and Halo Exchange

### 7.1 Decomposition Strategies

Parallel decomposition is handled by `decomposePar` using methods registered in the `decompositionMethod` factory:

| Method | Algorithm | Best For |
|---|---|---|
| `scotch` / `ptscotch` | Graph partitioning (min-cut) | General unstructured meshes |
| `metis` | Multilevel k-way partitioning | Large meshes, good load balance |
| `simple` | Geometric slicing along a coordinate | Structured/block meshes |
| `hierarchical` | Multi-level geometric slicing | Multi-socket HPC nodes |
| `multiLevel` | Scotch + hierarchical composite | NUMA-aware decomposition |
| `manualDecomp` | User-supplied cell → processor map | Custom topologies |

The decomposition produces `processor*` sub-directories, each containing a subdomain `polyMesh` with:
- Internal faces of the subdomain
- `processorCyclicFvPatch` faces at inter-processor boundaries

### 7.2 Processor Patches and Halo Exchange

Each processor boundary is represented as a **processor patch pair** (patch on owner proc + patch on neighbour proc). Halo exchange is implemented through `Pstream` (built on MPI):

```
Pstream::exchange<List<T>>()   ← non-blocking all-to-all
combineReduce()                ← reduction operations
```

The `processorFvPatch` computes the interpolated face value from the remote processor's ghost cell. Since OpenFOAM stores each face only once (on the owning processor), the halo mechanism is strictly one-sided for flux computation but requires two-sided communication for gradient reconstruction.

### 7.3 `globalMeshData` — Shared Point/Edge Topology

For operations requiring knowledge of the globally shared mesh topology (e.g., point-based smoothing in mesh motion solvers, consistent normals at processor boundaries), `globalMeshData` tracks:

- Shared points (owned by multiple processors)
- Shared edges
- Global point/face/cell numbering via `globalIndex`

---

## 8. Mesh Generation Ecosystem

### 8.1 `blockMesh`

`blockMesh` constructs structured hexahedral meshes from a `blockMeshDict`. Key features:
- Multi-block connectivity with shared edges/faces
- Non-uniform grading via `simpleGrading` or `edgeGrading` (geometric, parabolic, or user-defined expansion ratios)
- **Curved edges** via `arc`, `spline`, `BSpline`, `polyLine` primitives mapped onto block edges
- **Merge patch pairs** for conformal multi-block interfaces that become internal faces
- Wedge topology for axisymmetric 2D problems (single-cell-thick $5°$ sector)

### 8.2 `snappyHexMesh`

The workhorse surface-conforming mesh generator for complex geometry:

**Phase 1 — Castellated mesh:** Background `blockMesh` hex mesh is refined toward the STL/OBJ/triSurface geometry using `hexRef8` (surface proximity + curvature + feature-angle refinement regions). Result is a staircase approximation.

**Phase 2 — Snapping:** Surface vertices are projected onto the target geometry using a constrained Laplacian displacement, with feature-line preservation (sharp edges, corners) controlled by `nFeatureSnapIter` and `explicitFeatureSnap`.

**Phase 3 — Layer addition:** Prismatic boundary layer cells are extruded inward from the snapped surface using `addLayers`, with `expansionRatio` and `finalLayerThickness` controlling the near-wall resolution. Layer collapse detection prevents degenerate cells on concave geometric features.

### 8.3 Third-Party Import: `fluentMeshToFoam`, `gmshToFoam`, `cgnsToFoam`

OpenFOAM maintains conversion utilities for external mesh formats. The conversion pipeline:
1. Read foreign format (CGNS, Fluent `.msh`, Gmsh `.msh v2/v4`)
2. Construct `polyMesh` arrays from the foreign cell/element connectivity
3. Identify and merge duplicate faces to build the owner/neighbour arrays
4. Map zone/group names to OpenFOAM patch names

A notable subtlety: CGNS and Fluent meshes use *element-based* topology (each element owns its faces), requiring a face-deduplication pass to produce OpenFOAM's face-based representation.

---

## 9. `topoSet` and `faceZone`/`cellZone` Framework

Beyond patch boundaries, OpenFOAM uses **zones** to mark subsets of the mesh for localised physics:

- **`cellZone`** — used by MRF zones, porosity models, per-zone material properties, and source terms
- **`faceZone`** — used by AMI interfaces, explicit porosity jump faces (`porousBafflePressure`), phase-change interface tracking
- **`pointZone`** — used by mesh motion solvers for prescribed-displacement regions

`topoSet` operations (`cellToCell`, `faceToFace`, `surfaceToCell`, `cylinderToCell`, etc.) provide a composable set algebra for building zones from geometric primitives or field thresholds — analogous to a Boolean CSG system for mesh subsets.

---

## 10. `checkMesh` — Mesh Validation Pipeline

`checkMesh` runs a systematic quality audit, the results of which directly govern solver stability:

```
Mesh stats:
  points    non-orthogonality: max 45.3, average 8.2
  faces     skewness:          max 0.81, average 0.11
  cells     aspect ratio:      max 12.4, average 2.1
  patches   volume ratio:      max 8.3
```

Critical thresholds for robust RANS/LES simulations:

| Metric | Warning | Fatal |
|---|---|---|
| Non-orthogonality | > 70° | > 85° |
| Skewness | > 4 | > 20 |
| Aspect ratio | > 100 | (problem-dependent) |
| Face flatness | > 0.5 | > 1.0 |

Internally, `checkMesh` calls `primitiveMeshTools::faceOrthogonality()`, `primitiveMeshTools::faceSkewness()`, and related static methods — the same functions used by the solver at runtime to compute correction tensors.

---

## 11. Mesh Representation in Memory: A Systems Perspective

The entire `fvMesh` is held in process memory as flat `Field<T>` arrays (aliased to `std::vector<T>` with alignment guarantees). Key layout considerations:

- **Face ordering** (internal-first, then patches) enables cache-friendly inner product loops in the matrix assembly
- **Lazy demand-driven geometry** (`DemandDrivenData<T>`) avoids recomputing cell volumes during time-stepping unless `movePoints()` is called
- **`IOobject` registration** means every mesh entity participates in the `objectRegistry`, enabling runtime introspection and automated I/O without manual bookkeeping
- **`mapPolyMesh`** carries point/face/cell mapping tables through topological changes, allowing all registered fields to remap conservatively (via `mapFields()`) without explicit field-level intervention

---

## Summary

OpenFOAM's mesh infrastructure is a carefully layered system in which a minimalist face-based topological kernel (`polyMesh`) is progressively augmented with geometry (`primitiveMesh`), FVM discretisation support (`fvMesh`), boundary condition dispatch (`polyBoundaryMesh`), and parallel decomposition — all without any element-type restrictions. Its polyhedral generality, built-in non-orthogonality correction framework, and runtime-extensible topological change machinery make it uniquely suited to industrial CFD problems involving complex geometry, moving parts, and adaptive resolution, at the cost of requiring more careful mesh quality management than topology-restricted solvers.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 4.6
>
> Prompt: Provide a thorough description of the mesh infrastructure in OpenFOAM (Open Field Operation and Manipulation). This description is intended for a senior application engineer specializing in computational fluid dynamics (CFD).
