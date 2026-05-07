# Parallel Execution in OpenFOAM: A Technical Overview

---

## 1. Foundational Architecture

OpenFOAM's parallel framework is built entirely on the **Message Passing Interface (MPI)** standard, typically realized through OpenMPI, MPICH, or vendor-optimized implementations (Intel MPI, Cray MPICH, etc.). The framework is not a thin wrapper — it is deeply woven into the class hierarchy through the `Pstream` abstraction layer, which makes parallelism largely transparent to solver and library developers writing at higher levels of the code.

The central design philosophy is **domain decomposition parallelism (DDP)**: the computational mesh is partitioned into sub-domains, each assigned to exactly one MPI rank (process). All ranks execute the same binary concurrently, following the **Single Program, Multiple Data (SPMD)** model. Each rank operates on its local sub-domain and communicates only at shared inter-processor boundaries.

---

## 2. Domain Decomposition

### 2.1 The `decomposePar` Utility

Before a parallel run, the global mesh and field data must be split into per-processor directories (`processor0/`, `processor1/`, …). This is handled by `decomposePar`, which reads the `system/decomposeParDict` dictionary. The key entries are:

```cpp
numberOfSubdomains  N;
method              scotch; // or hierarchical, simple, metis, manual, multiLevel
```

### 2.2 Decomposition Methods

| Method | Strategy | Typical Use Case |
|---|---|---|
| `simple` | Cuts along user-specified coordinate planes | Structured, Cartesian-like meshes |
| `hierarchical` | Sequential directional cuts (x→y→z) | Structured meshes; controllable aspect ratio |
| `scotch` | Graph-partitioning via PT-Scotch library | General unstructured meshes; minimizes inter-proc faces |
| `metis` | Graph-partitioning via METIS library | Similar to scotch; alternative weighting |
| `kahip` | KaHIP partitioner | High-quality cuts for large partitions |
| `multiLevel` | Hierarchical combination of methods | Clusters/nodes with NUMA topology awareness |
| `manual` | User-supplied `cellToProc` field | Full control; useful for restarts with known good decomposition |

The quality metric that matters most is the **cut-edge count** (inter-processor faces), since every such face requires MPI communication. Scotch/METIS minimize this using multilevel k-way graph partitioning.

### 2.3 Load Imbalance Considerations

Pure geometric decomposition ignores computational cost per cell. For cases with:
- **Combustion / chemistry** (stiff ODE solvers per cell)
- **Lagrangian parcels** (DEM, spray, particle tracking)
- **Adaptive mesh refinement (AMR)**

…geometric balance in cell count does not guarantee temporal load balance. Scotch supports a `processorWeights` list and a `strategy` string to bias partitioning.

---

## 3. Inter-Processor Boundary Handling

### 3.1 Processor Patches

When `decomposePar` cuts the mesh, shared faces between two sub-domains become **processor patches** — a paired boundary condition type. Each inter-processor face appears in *both* the owning processor's patch list (`processorFvPatch`) and the neighbour's, as a mirror image.

In the mesh topology, these are stored as:
```
type  processor;
myProcNo   <ownerRank>;
neighbProcNo  <neighbourRank>;
```

The field values on processor patches are communicated via non-blocking MPI calls (`MPI_Isend` / `MPI_Irecv`) before each evaluation requiring neighbour cell data.

### 3.2 The `Pstream` Class Hierarchy

OpenFOAM's communication infrastructure:

```
Pstream (abstract interface)
├── UPstream           ← low-level MPI rank/size queries
├── IPstream           ← input (receive) stream
├── OPstream           ← output (send) stream
├── PstreamBuffers     ← buffered non-blocking scatter/gather
└── combineReduce      ← collective reduction operations
```

Key global flags: `Pstream::parRun()` returns `true` during parallel execution, enabling conditional logic without code duplication.

### 3.3 Communication Patterns

OpenFOAM uses three primary MPI communication modes:

**Blocking point-to-point** — Used for initialization and setup; avoids deadlock complexity at the cost of serialization.

**Non-blocking (asynchronous)** — The dominant pattern during time-stepping. `initEvaluateCoupled()` posts all `MPI_Isend`/`MPI_Irecv`, then `evaluateCoupled()` calls `MPI_Waitall`. This overlaps communication with computation where possible.

**Collective operations** — `Pstream::gather`, `Pstream::scatter`, `Pstream::allReduce` wrap `MPI_Gather`, `MPI_Bcast`, `MPI_Allreduce`, etc. Used for global residual norms, time-step control (`CoNum`), and convergence monitoring.

---

## 4. Linear Algebra and Parallel Solvers

### 4.1 Matrix Assembly

The `fvMatrix<T>` assembled on each rank is inherently local. Off-diagonal contributions from across processor boundaries are handled implicitly: the processor patch contributions enter the source vector (`b`) via the boundary condition interface mechanism, not as explicit matrix entries. This keeps the local matrix square and banded.

### 4.2 Iterative Solver Parallelism

All iterative solvers in OpenFOAM (GAMG, PCG, PBiCGStab, smoothers) are parallel by design:

- **Dot products / norms** → `MPI_Allreduce` for global sums
- **Matrix-vector products** → local SpMV + halo exchange for processor boundary cells
- **GAMG (Geometric-Algebraic MultiGrid)** — the most sophisticated: it agglomerates cells across processor boundaries at coarser levels, requiring inter-processor communication at every multigrid level. The coarsest level can be solved by a direct method on rank 0 (after gathering) or by further iteration.

### 4.3 Direct Solvers

OpenFOAM does not include a native parallel direct solver. For workflows requiring one (e.g., very stiff systems), external libraries such as **MUMPS**, **SuperLU_DIST**, or **PETSc** are integrated via adapter layers in some distributions (e.g., foam-extend, OpenFOAM-ESI extensions).

---

## 5. Running in Parallel

### 5.1 Invocation

```bash
mpirun -np <N> <solver> -parallel
```

The `-parallel` flag switches `Pstream::parRun()` to `true` and directs I/O to per-processor directories. Under SLURM or PBS:

```bash
srun --ntasks=<N> --ntasks-per-node=<ppn> <solver> -parallel
```

### 5.2 Parallel I/O Strategies

OpenFOAM supports several I/O modes, configurable via `system/controlDict`:

```cpp
// In controlDict:
runTimeModifiable  true;
writeFormat        binary;  // ascii or binary
```

**Distributed I/O (default):** Each rank writes its own `processorN/` subdirectory. Simple, scales well, but produces many small files — problematic on parallel filesystems (Lustre, GPFS) at high rank counts.

**Collated I/O (OpenFOAM v1706+):** All ranks write to a single file per time directory using `fileHandler collated`. A dedicated thread on each node (or rank 0) aggregates data, drastically reducing metadata operations. This is strongly recommended for runs on >256 cores.

```cpp
OptimisationSwitches
{
    fileHandler collated;
    maxThreadFileBufferSize 2e9;
}
```

**masterUncollated:** A hybrid — one file per time step per field, written by the master rank. Good for moderate scale.

---

## 6. Post-Processing and Reconstruction

### 6.1 `reconstructPar`

Recombines per-processor field data into a single case directory for serial post-processing. It reads `constant/polyMesh/cellProcAddressing` (and face/point/boundary addressing files) to map local indices back to global.

Options of note:
```bash
reconstructPar -time '500:'     # only recent times
reconstructPar -fields '(U p)'  # selective field reconstruction
reconstructPar -lagrangianFields # include Lagrangian clouds
```

### 6.2 In-Situ Parallel Post-Processing

For very large cases, reconstruction is impractical. Alternatives:
- **ParaView + pvserver / Catalyst:** OpenFOAM writes `*.foam` dummy files; ParaView's OpenFOAM reader directly reads distributed `processorN/` data in parallel
- **`foamToVTK -parallel`:** Converts in parallel; VTK handles multiblock distributed datasets
- **runTimePostProcessing function objects:** Execute on the fly during the simulation without writing full fields

---

## 7. Special Parallel Considerations

### 7.1 Cyclic AMI in Parallel

**Arbitrary Mesh Interface (AMI)** patches (sliding meshes, GGI interfaces) are among the most communication-intensive features. When an AMI interface is split across multiple processors:
- A **distributed AMI** algorithm triangulates the patch face intersection across ranks
- Requires multiple rounds of `MPI_Alltoallv` to communicate face overlap data
- The `cyclicAMI` weight computation is done once at startup (or each timestep for moving meshes), and the resulting interpolation weights are stored locally per rank

### 7.2 Dynamic Mesh and AMR

`dynamicRefineFvMesh` performs local refinement/coarsening. In parallel:
- Refinement decisions are made locally but must be globally consistent at processor boundaries
- After topology changes, the mesh is **not** re-decomposed; load imbalance accumulates
- **Dynamic load balancing** via `zoltanRenumber` or external libraries (e.g., ParMETIS) can periodically re-partition; this is available in some forks (foam-extend) but not yet fully standard in the ESI/Foundation releases

### 7.3 Lagrangian Particle Tracking

Parcels are owned by the processor whose cell they currently occupy. When a parcel crosses a processor boundary:
1. It is flagged for transfer
2. `Cloud::move()` calls `patchInteractionModel` which invokes `Pstream::exchange` to transfer parcel data
3. The parcel is deleted from the source rank and instantiated on the destination rank

This creates inherent load imbalance when particles cluster (e.g., dense spray zones), since processor boundary crossings also trigger communication overhead.

### 7.4 Overset / Chimera Meshes

Overset meshes (available via `overPimpleFoam`) add another layer: interpolation donor/acceptor pairs may span processor boundaries, requiring explicit communication to transfer field values between overlapping mesh regions in parallel.

---

## 8. Performance Engineering

### 8.1 Scalability Behavior

OpenFOAM's parallel efficiency follows a characteristic profile:

- **Strong scaling** degrades when sub-domains become too small (< ~5,000–10,000 cells/rank for typical incompressible flows), as MPI overhead dominates useful compute
- **Weak scaling** is generally good for structured problems; complex geometries with high surface-to-volume ratios per sub-domain increase communication-to-computation ratio
- GAMG exhibits super-linear coarsening at high rank counts because the coarse-level problem grows with processor count

### 8.2 Profiling Tools

| Tool | What it measures |
|---|---|
| `foamLog` + residual plots | Convergence quality; indirect load imbalance indicator |
| `timer` / `profilingTrigger` function objects | Per-region wallclock instrumentation |
| Scalasca + Score-P | MPI call analysis, communication bottlenecks |
| Intel VTune / Arm MAP | Load imbalance, NUMA effects, memory bandwidth |
| `mpiP` | Lightweight MPI profiling with minimal overhead |

### 8.3 Practical Tuning Recommendations

**Decomposition:** Prefer `scotch` or `metis` for unstructured meshes. Inspect inter-processor face counts with `checkMesh -allTopology` post-decomposition. Aim for < 5–8% of total faces being processor faces.

**MPI transport:** On InfiniBand clusters, ensure OpenMPI is built with `--with-verbs` or use UCX transport. Avoid shared-memory MPI for cross-node communication.

**GAMG coarsening:** Set `nCellsInCoarsestLevel` to approximately `(totalCells / N_ranks)` to prevent the coarsest level from collapsing to a single cell — this degrades parallel efficiency significantly.

**I/O:** Use `collated` format on Lustre/GPFS. Set write interval generously; I/O is often the dominant wall-time consumer at high cell counts.

**Pinning:** Use `--bind-to core` or equivalent to prevent MPI rank migration across NUMA domains.

---

## 9. Emerging and Non-Standard Parallel Paradigms

### 9.1 GPU Acceleration

The ESI OpenFOAM v2212+ includes experimental support for GPU offloading via:
- **AmgX** (NVIDIA): Drop-in replacement for GAMG/PCG on GPU clusters
- **Ginkgo**: Vendor-neutral linear algebra library with CUDA/HIP/SYCL backends
- **PETSc + CUDA**: Through the PETSc interface layer

These target the linear algebra phase only; mesh operations and boundary conditions remain CPU-bound.

### 9.2 Hybrid MPI + OpenMP

Some community forks and research branches explore **MPI + OpenMP** threading (one MPI rank per NUMA domain, multiple threads per rank). The core `fvMesh` data structures are not thread-safe without careful locking, so adoption is limited. The `functionObject` system and some chemistry solvers have thread-safe implementations.

### 9.3 Task-Based Parallelism

Experimental integration with **HPX** and **oneTBB** has been explored in the academic literature for irregular workloads (chemistry integration, particle tracking), but none have been merged into mainline releases.

---

## Summary

OpenFOAM's parallel execution model is mature, robust, and scales well to thousands of cores for mesh-dominated workloads when properly configured. The key levers for a senior CFD engineer are: **decomposition quality** (minimize inter-processor faces, respect load balance), **I/O strategy** (collated format at scale), **linear solver tuning** (GAMG coarsening parameters), and **MPI transport configuration** (fabric-specific tuning). The framework's `Pstream` abstraction means that custom solvers and function objects written at the `fvMesh` level inherit parallel capability almost automatically — but workloads with strong spatial clustering (particles, AMR, combustion fronts) require explicit attention to dynamic load balance, which remains an open engineering challenge in the standard distribution.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 4.6
>
> Prompt: Provide a thorough description of the parallel execution of the code available in OpenFOAM (Open Field Operation and Manipulation). This description is intended for a senior application engineer specializing in computational fluid dynamics (CFD).
