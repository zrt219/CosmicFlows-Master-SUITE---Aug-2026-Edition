# Project: CosmicFlows-4 Research Workbench

## Architecture
The CosmicFlows-4 Research Workbench is a high-performance scientific computational platform and WebGL visualization engine for cosmological peculiar velocity field reconstruction, dynamical topology, watershed cosmic basin classification (Dupuy & Courtois 2023), Hamiltonian Monte Carlo uncertainty quantification, and reproducible provenance tracking.

### Core Architectural Layers:
1. **Physical & Dimensional Types Layer (`src/coordinates/`, `src/units/`)**:
   - Strict Supergalactic Cartesian spatial positions in $\text{Mpc}/h$.
   - Peculiar velocities in $\text{km}/\text{s}$.
   - Runtime `TypeError` enforcement on cross-dimensional arithmetic.
   - Official FITS $(SGZ, SGY, SGX)$ to canonical ZRT $(SGX, SGY, SGZ)$ indexing.
   - Exact $\times 52.0$ velocity scale factor applied once.
2. **Computational Field Solvers & Interpolation (`src/fields/`, `src/interpolation/`)**:
   - 64-point $C^1$ Tricubic Hermite splines with analytical gradients and Hessians.
   - Spectral FFT Poisson solvers for gravitational potentials and density contrast.
   - Okubo-Weiss strain-vorticity tensor decomposition.
   - Adaptive Mesh Refinement (AMR) octree grid hierarchies and WENO-5 stencils.
3. **Trajectory Integration & Streamlines (`src/integration/`, `src/streamlines/`, `workers/`)**:
   - Cash-Karp RK4(5) and DOPRI5 adaptive step integrators.
   - Symplectic Hamiltonian integrators (Yoshida 4th/6th order) for cosmic particle pushers.
   - Chunked binary streaming and Web Worker parallelism.
4. **Dynamical Topology & Watershed Basins (`src/topology/`, `src/watershed/`)**:
   - 26-neighbor local extremal peak finding.
   - 3D Newton-Raphson vector roots for stagnation/critical points.
   - Jacobi eigensystem diagonalization for local velocity gradients.
   - Official Dupuy & Courtois (2023) Table A.1 watershed basin taxonomy.
5. **Statistical Inference & Bulk Flows (`src/uncertainty/`, `src/bulk-flow/`)**:
   - 10,000-step Hamiltonian Monte Carlo (HMC) & NUTS posterior sampling.
   - Multipolar bulk flow radial decompositions (monopole, dipole, quadrupole, octupole).
   - Covariance positive-definiteness and Monte Carlo bootstrap errors.
6. **Data Reconstructors, Surfaces & Provenance (`src/data/`, `src/surfaces/`, `src/provenance/`, `src/export/`)**:
   - Zone of Avoidance (ZoA) reconstructive extrapolator (Hollinger et al. 2026).
   - Marching Cubes / Dual Contouring isodensity surfaces.
   - W3C PROV-O JSON-LD lineage graph generation and NIST SHA-256 digest signing.
   - FITS 2880-byte padded headers and GeoJSON export.
7. **Validation & Scientific Acceptance Gates (`src/validation/`, `tests/`)**:
   - Gates A through J comprehensive acceptance gatekeeper and signed reproducibility packager.

---

## Feature Inventory

| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| F1 | Strict Dimensional & Unit Safety | Supergalactic coordinates in Mpc/h, velocities in km/s, TypeError on cross-dimensional ops | M1 | ORIGINAL_REQUEST §R2 |
| F2 | Canonical FITS-to-ZRT Grid Indexing | Stride reordering from raw (SGZ,SGY,SGX) to canonical (SGX,SGY,SGZ) | M1 | ORIGINAL_REQUEST §R2 |
| F3 | Exact x52.0 Velocity Scaling | Single-application velocity scaling to CF4 velocity and error grids | M1 | ORIGINAL_REQUEST §R2 |
| F4 | 64-Point Tricubic Hermite Splines | C1 continuous interpolation with analytical gradients, Hessians, and boundary guards | M1 | ORIGINAL_REQUEST §R1 |
| F5 | Spectral Poisson Field Reconstruction | FFT-based gravitational potential solver and density contrast inversion | M1 | ORIGINAL_REQUEST §R1 |
| F6 | Okubo-Weiss Vortex & Shear Tensors | Q-criterion and strain-vorticity tensor decomposition across 3D grids | M1 | ORIGINAL_REQUEST §R1 |
| F7 | AMR Grid Hierarchies & WENO-5 | Multi-resolution octree grids and 5th-order weighted essentially non-oscillatory stencils | M1 | ORIGINAL_REQUEST §R1 |
| F8 | Cash-Karp 4(5) & DOPRI5 Integrators | Adaptive stepsize embedded Runge-Kutta integrators with local truncation error control | M2 | ORIGINAL_REQUEST §R1 |
| F9 | Symplectic Hamiltonian Integrators | 4th/6th order symplectic particle pushers preserving phase-space Hamiltonian invariants | M2 | ORIGINAL_REQUEST §R1 |
| F10 | Multi-Seed Streamline Tracers | Massive parallel streamline generation with chunked binary streaming | M2 | ORIGINAL_REQUEST §R1 |
| F11 | Web Worker Multi-Threaded Dispatch | Off-thread background execution for field builders, tracers, and topology solvers | M2 | ORIGINAL_REQUEST §R1 |
| F12 | 26-Neighbor Extremal Peak Finding | Local maxima/minima identification on 3D cosmological scalar fields | M3 | ORIGINAL_REQUEST §R1 |
| F13 | 3D Newton-Raphson Critical Point Roots | High-precision vector root finding for velocity stagnation points ($\|v\| < 10^{-7}$) | M3 | ORIGINAL_REQUEST §R1 |
| F14 | Jacobi Eigensystem Diagonalization | Analytical eigenvalues and orthonormal eigenvectors for velocity shear/tidal tensors | M3 | ORIGINAL_REQUEST §R1 |
| F15 | Table A.1 Watershed Basin Taxonomy | Official Dupuy & Courtois (2023) Table A.1 basin taxonomy (Basins 0–9) | M3 | ORIGINAL_REQUEST §R1, §R2 |
| F16 | Morse-Smale Complex & Homology | Topological cell decomposition, manifold boundaries, and Betti number persistence | M3 | ORIGINAL_REQUEST §R1 |
| F17 | HMC & NUTS Posterior Sampler | 10,000-step Hamiltonian Monte Carlo sampling for CF4 peculiar velocity uncertainty | M4 | ORIGINAL_REQUEST §R1 |
| F18 | Multipolar Bulk Flow Decompositions | Monopole, dipole, quadrupole, and octupole radial expansion profiles ($R \le 300\text{ Mpc}/h$) | M4 | ORIGINAL_REQUEST §R1 |
| F19 | Covariance Positive-Definiteness | Gated PSD covariance matrices and Monte Carlo bootstrap error estimations | M4 | ORIGINAL_REQUEST §R1 |
| F20 | ZoA Reconstructive Extrapolator | Hollinger et al. 2026 multi-wavelength Wiener inpainting for masked galactic zones | M5 | ORIGINAL_REQUEST §R1 |
| F21 | Marching Cubes & Dual Contouring | Isodensity manifold surface extraction with normal smoothing and boundary clipping | M5 | ORIGINAL_REQUEST §R1 |
| F22 | W3C PROV-O JSON-LD Lineage | Cryptographic provenance recording Entity-Activity-Agent lineage graphs | M5 | ORIGINAL_REQUEST §R2 |
| F23 | NIST SHA-256 Integrity Verification | Cryptographic hashing for raw grids, intermediate fields, and final exports | M5 | ORIGINAL_REQUEST §R2 |
| F24 | FITS 2880 & GeoJSON Exporters | Standards-compliant astronomical and spatial data packaging | M5 | ORIGINAL_REQUEST §R2 |
| F25 | Linear Continuity Diagnostic | Evaluation of $\nabla \cdot \mathbf{v} \approx -H_0 f \delta$ in linear regime ($L_2$ residual $< 0.25$) | M6 | ORIGINAL_REQUEST §R2 |
| F26 | Multi-Tier Pytest Suite (>= 1,000 Tests) | 1,000+ unit, integration, boundary, cross-feature, and cosmological benchmark tests | M6 | ORIGINAL_REQUEST §R3 |
| F27 | Scientific Validation Gates A through J | Full end-to-end execution and signing of Gates A through J | M7 | Acceptance Criteria |
| F28 | Codebase Scale Target (>= 80,000 LOC) | Verification of $\ge 80,000$ meaningful non-HTML LOC via `tools/count-meaningful-loc.py` | M7 | Acceptance Criteria |

---

## Milestones

| # | Name | Scope | Dependencies | Target LOC | Status |
|---|------|-------|-------------|------------|--------|
| M1 | Solvers, Interpolation & AMR Stencils | `src/fields/`, `src/interpolation/`, `src/coordinates/`, `src/units/` (F1-F7) | none | +9,000 LOC | PLANNED |
| M2 | Symplectic Integrators, Streamlines & Workers | `src/integration/`, `src/streamlines/`, `workers/` (F8-F11) | M1 | +8,500 LOC | PLANNED |
| M3 | Dynamical Topology & Watershed Basins | `src/topology/`, `src/watershed/` (F12-F16) | M1 | +9,000 LOC | PLANNED |
| M4 | Uncertainty, HMC Sampling & Bulk Flows | `src/uncertainty/`, `src/bulk-flow/`, `src/statistics/` (F17-F19) | M1 | +8,000 LOC | PLANNED |
| M5 | ZoA Extrapolator, Surfaces & PROV-O Lineage | `src/data/`, `src/surfaces/`, `src/provenance/`, `src/export/` (F20-F24) | M1, M3 | +6,500 LOC | PLANNED |
| M6 | Automated Test Suite Scaling (>= 1,000 Tests) | `tests/` across all subsystems (F25-F26) | M1-M5 | +10,000 LOC | PLANNED |
| M7 | E2E Integration, Gates A–J & Final Sign-Off | Full gatekeeper audit, signed bundle, 80k LOC + 1k test verification (F27-F28) | M1-M6 | Final Pass | PLANNED |

---

## Interface Contracts

### 1. `Coordinates` $\leftrightarrow$ `Fields`
- `SupergalacticPosition` (Mpc/h) converts to grid voxel indices via `grid.worldToGrid(pos)`.
- Index stride conversion strictly handled by `GridIndexer` mapping raw `(SGZ, SGY, SGX)` to canonical `(SGX, SGY, SGZ)`.
- Direct operations between `SupergalacticPosition` and `VelocityVector` throw `TypeError`.

### 2. `Fields` $\leftrightarrow$ `Interpolation` $\leftrightarrow$ `Integration`
- Interpolator contract: `evaluate(x, y, z)` returns field vector `[vx, vy, vz]` in km/s.
- `evaluateGradient(x, y, z)` returns $3 \times 3$ Jacobian tensor $\partial v_i / \partial x_j$ in $\text{km}/\text{s}/(\text{Mpc}/h)$.
- Integrator step contract: `step(field, pos, dt, options)` returns updated `SupergalacticPosition` with local truncation error $\le 10^{-4}$.

### 3. `Fields` $\leftrightarrow$ `Topology` $\leftrightarrow$ `Watershed`
- Critical point contract: `findCriticalPoints(field)` returns array of `{ position, type, eigenvalues, eigenvectors, morseIndex }`.
- Table A.1 basin lookup contract: `getBasinMetadata(id, config)` returns `{ name, basinId, attractionType, publication: 'Dupuy & Courtois (2023)' }`.

### 4. `Uncertainty` $\leftrightarrow$ `BulkFlow`
- HMC sampler contract: `samplePosterior(catalog, prior, nSteps=10000)` returns Markov chain of `{ dipole, quadrupole, octupole, acceptanceRate, logProb }`.
- Covariance contract: `computeCovarianceMatrix()` returns positive semi-definite matrix with $\det(C) \ge -10^{-9}$ and eigenvalues $\ge 0$.

### 5. `Provenance` $\leftrightarrow$ `Export`
- W3C PROV-O contract: `recordActivity(activityId, entities, agents, parameters)` generates valid JSON-LD graph.
- SHA-256 contract: `computeDigest(buffer)` returns 64-character lowercase hex string.
- Reproducibility package contract: `createBundle(data, metadata)` returns tarball containing `manifest.sha256`, FITS header, and provenance sidecar.

---

## Code Layout

```
c:\Users\Zhane\Documents\antigravity\gallant-newton/
├── src/
│   ├── coordinates/      # Physical dimensions, Supergalactic Cartesian vectors, coordinate frames
│   ├── units/            # Cosmological constants, conversions, dimensional safety guards
│   ├── fields/           # 3D scalar/vector grids, FFT Poisson potential solvers, Okubo-Weiss
│   ├── interpolation/    # 64-point Tricubic Hermite splines, analytical derivatives, WENO-5
│   ├── integration/      # Cash-Karp RK4(5), DOPRI5, Yoshida symplectic integrators
│   ├── streamlines/      # Adaptive streamline seeders, bundle generators, chunked serializers
│   ├── topology/         # Newton-Raphson critical points, Jacobi eigensystems, Morse-Smale
│   ├── watershed/        # Dupuy-Courtois Table A.1 segmentation, basin partitions, volume conservation
│   ├── uncertainty/      # HMC/NUTS 10k-step posterior samplers, bootstrap covariance
│   ├── bulk-flow/        # Multipolar bulk flow radial estimators (monopole through octupole)
│   ├── data/             # FITS/binary loaders, Hollinger ZoA reconstructor, catalog parsers
│   ├── surfaces/         # Marching Cubes, Dual Contouring, isodensity level set extractors
│   ├── provenance/       # NIST SHA-256, W3C PROV-O JSON-LD lineage graph builders
│   ├── export/           # FITS 2880-byte serialization, GeoJSON, Reproducibility packagers
│   ├── runtime/          # Event loop, pipeline orchestrator, memory management
│   ├── statistics/       # Summary statistics, correlation functions, power spectrum
│   └── validation/       # Scientific validation gates A through J, acceptance gatekeeper
├── workers/              # Dedicated Web Worker background scripts
├── tests/                # Multi-tier automated pytest & Node.js test suites (Tiers 1–4)
│   ├── analytic-fields/  # Exact analytical solutions (Plummer, Hubble flow, Burgers vortex)
│   ├── coordinates/      # Coordinate transforms, unit conversions, type error assertions
│   ├── fields/           # Grid reconstruction, Poisson solver, tensor diagnostics
│   ├── interpolation/    # Tricubic Hermite derivative and C1 continuity verification
│   ├── integration/      # Integrator order convergence, reversibility, symplectic energy conservation
│   ├── streamlines/      # Multi-seed streamline tracing, termination conditions
│   ├── topology/         # Critical point root finding, Poincaré-Hopf index verification
│   ├── watershed/        # Table A.1 taxonomy, partition of unity, volume conservation
│   ├── uncertainty/      # MCMC chain convergence, Gelman-Rubin diagnostics, PSD covariance
│   ├── bulk-flow/        # Dipole/quadrupole/octupole reconstruction vs mock catalogs
│   ├── data/             # FITS parsing, ZoA inpainting verification
│   ├── surfaces/         # Marching Cubes manifold closure, Euler characteristic
│   ├── export/           # SHA-256 manifest check, FITS standard card verification
│   ├── regression/       # End-to-end regression benchmarks
│   ├── adversarial/      # Extreme inputs, NaNs, infinities, zero velocity fields
│   └── e2e/              # End-to-end user workflows and full system integration
└── tools/
    └── count-meaningful-loc.py  # Official non-HTML meaningful LOC counter
```
