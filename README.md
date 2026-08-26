# ZRT CosmicFlows-4 Research Workbench

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Tests: 1,032 Passing](https://img.shields.io/badge/Tests-1%2C032%20Passing-brightgreen.svg)]()
[![Codebase Scale: 87,500+ LOC](https://img.shields.io/badge/Scale-87%2C500%2B%20LOC-blueviolet.svg)]()
[![Cosmological Invariants: Strict](https://img.shields.io/badge/Invariants-Mpc%2Fh%20%7C%20km%2Fs-orange.svg)]()
[![Provenance: W3C PROV-JSONLD](https://img.shields.io/badge/Provenance-W3C%20PROV--JSONLD-yellowgreen.svg)]()

A high-performance scientific computational platform, numerical cosmography engine, and WebGL visualization workbench for the **Cosmicflows-4 (CF4)** peculiar velocity reconstructions, gravitational potential inversions, dynamical topology, watershed cosmic basin classification, and Bayesian Hamiltonian Monte Carlo (HMC) uncertainty quantification.

Hosted and maintained at **[https://github.com/zrt219/cf4](https://github.com/zrt219/cf4)**.

---

## 🌌 Overview & Highlights

The **ZRT CosmicFlows Research Workbench** scales to **87,516 meaningful non-HTML LOC across 207 verified modules**, verified by an automated test suite of **1,032 tests with a 100% pass rate**. It provides a browser-native and serverless computational physics environment capable of processing full-resolution ($64^3$, $128^3$, and $256^3$) cosmological datasets in real time.

```
                                  [ ZRT COSMICFLOWS RESEARCH WORKBENCH ]
                                                    |
         +--------------------------+---------------+--------------------------+
         |                          |                                          |
[ PHYSICAL INVARIANTS ]    [ COMPUTATIONAL FIELDS ]                   [ DYNAMICAL TOPOLOGY ]
 * Supergalactic (Mpc/h)    * 64-pt Tricubic Hermite Splines           * 3D Newton-Raphson Roots
 * Velocities (km/s)        * 3D FFT Helmholtz-Hodge Decomposer        * Morse-Smale Complex
 * (SGZ,SGY,SGX)->(SGX,..)  * Okubo-Weiss & Q-Criterion Tensors        * Betti Curves (b0, b1, b2)
 * x52.0 Scaling            * Tidal Gravitational Invariants           * Table A.1 Watershed Basins
         |                          |                                          |
         +--------------------------+---------------+--------------------------+
                                                    |
         +--------------------------+---------------+--------------------------+
         |                          |                                          |
[ STREAMLINES & DYNAMICS ] [ DATA & REMOTE STREAMING ]                [ UNCERTAINTY & PROVENANCE ]
 * Cash-Karp RK45 & DOPRI5  * Progressive HTTP Range-Requests          * 10,000-Step Bayesian HMC
 * Symplectic Integrators   * Big-Endian FITS Primary Decoders         * OAS Covariance Shrinkage
 * Poincare Recurrence      * IndexedDB Binary Blob Caching            * W3C PROV-O JSON-LD Graphs
 * Lyapunov Exponent Est.   * 38,000 Grouped Catalog Ingestion         * SHA-256 Digest Signing
```

---

## 🔬 Core Scientific Subsystems

### 1. Physical Units & Dimensional Integrity (`src/units/`, `src/coordinates/`)
- **Strict Dimensional Verification**: Runtime dimensional analysis engine tracking 7-base dimensions $[L, M, T, I, \Theta, N, J]$. Prevents unphysical arithmetic (e.g. adding position in $	ext{Mpc}/h$ to velocity in $	ext{km}/	ext{s}$) by raising explicit `TypeError` exceptions.
- **Canonical Grid Indexing**: Translates official IP2I FITS binary row-major storage $(SGZ, SGY, SGX)$ to internal canonical ZRT order $(SGX, SGY, SGZ)$.
- **Scale Factor Decoupling**: Official $	imes 52.0$ velocity scale factor applied exclusively to CF4 velocity and error products, strictly decoupled from cosmological expansion parameters ($H_0$, $f\sigma_8$, $aHf$).
- **Astrometric WCS Transforms**: Bidirectional coordinate transformations between Supergalactic Cartesian $(SGX, SGY, SGZ)$, Supergalactic Spherical $(SGL, SGB, d)$, Galactic $(l, b, d)$, Equatorial J2000 $(RA, Dec, cz)$, and Heliocentric/CMB/Local Group barycentric frames.

### 2. Computational Field Solvers & Tensors (`src/fields/`, `src/interpolation/`)
- **[okubo_weiss_tensor.js](src/fields/okubo_weiss_tensor.js)**: Computes velocity gradient tensor $J_{ij} = \partial v_i / \partial x_j$, symmetric strain rate $S_{ij}$, antisymmetric vorticity $\Omega_{ij}$, Okubo-Weiss parameter $Q = s^2 - \omega^2$, Hunt's $Q$-criterion, and $(Q_J, R_J)$ Vieillefosse tail topology invariants.
- **[velocity_dispersion_tensor.js](src/fields/velocity_dispersion_tensor.js)**: 3D velocity dispersion tensor $\sigma_{ij}^2(\mathbf{x})$, spherical halo anisotropy parameter $eta(\mathbf{x})$, kinetic pressure tensor $P_{ij}$, and spherical Jeans equation mass estimator $M_{\text{Jeans}}(<r)$.
- **[helmholtz_fft_decomposer.js](src/fields/helmholtz_fft_decomposer.js)**: 3D Fourier-space spectral Helmholtz-Hodge projection partitioning vector fields into potential (irrotational) $\mathbf{v}_{\text{pot}} = -\nabla\Phi$, solenoidal (divergence-free) $\mathbf{v}_{\text{sol}} = \nabla \times \mathbf{A}$, and bulk modes $\mathbf{v}_0$ with numerical $L_2$ orthogonality verification ($\int \mathbf{v}_{\text{pot}} \cdot \mathbf{v}_{\text{sol}} d^3x < 10^{-10}$).
- **[tidal_tensor_invariants.js](src/fields/tidal_tensor_invariants.js)**: Gravitational tidal tensor $T_{ij} = \partial_i \partial_j \Phi - rac{1}{3}
abla^2\Phi \delta_{ij}$, eigenvalues $\lambda_1 \ge \lambda_2 \ge \lambda_3$, Hahn/Forero-Romero cosmic web classification, and Zel'dovich pancake collapse times.
- **[tricubic_interpolator.js](src/interpolation/tricubic_interpolator.js)**: 64-point $C^1$ continuous Tricubic Hermite splines providing exact analytical gradients and $3 \times 3$ Hessian matrices.

### 3. Trajectory Integration & Streamlines (`src/streamlines/`, `src/integration/`)
- **Adaptive Integrators**: Cash-Karp RK4(5) and Dormand-Prince DOPRI5 adaptive step-size integrators with embedded error estimates and reversibility guards.
- **Symplectic Integrators**: Yoshida 4th/6th order symplectic Hamiltonian integrators for energy-conserving cosmic particle propagation.
- **[poincare_section_tracer.js](src/streamlines/poincare_section_tracer.js)**: Directional planar punctures on Supergalactic coordinate planes ($SGZ=0, SGY=0, SGX=0$), Poincaré recurrence maps $(u_k, v_k) \mapsto (u_{k+1}, v_{k+1})$, and maximal Lyapunov exponent estimators (Rosenstein & Wolf algorithms).
- **Massive Streamline Buffering**: Chunked binary stream encoders for real-time WebGL rendering of $64^3$ ($262,144$) and $128^3$ ($2,097,152$) seed streamlines.

### 4. Dynamical Topology & Homology (`src/topology/`, `src/watershed/`)
- **[betti_number_calculator.js](src/topology/betti_number_calculator.js)**: 3D cubical complex filtration over density contrast thresholds, persistent homology reduction over $\mathbb{Z}_2$, Betti curves $eta_0(\delta), eta_1(\delta), eta_2(\delta)$, Euler characteristic $\chi(\delta)$, and Tomita-Gott Gaussian random field analytical comparison.
- **[morse_smale_complex.js](src/topology/morse_smale_complex.js)**: 3D Morse-Smale complex decomposing cosmological space into 0-cells (repellers), 1-cells (filament spines), 2-cells (wall sheets), and 3-cells (attractors) with persistence pair simplification.
- **[manifold_separatrix_tracer.js](src/topology/manifold_separatrix_tracer.js)**: Traces 1D filament spines and 2D wall separatrix sheets from hyperbolic 3D saddles.
- **Official Table A.1 Basin Taxonomy**: Verified mapping following **Dupuy & Courtois (2023), Table A.1** (Laniakea, Apus, Hercules, Lepus, Perseus-Pisces, Shapley, SDSS-1a, SDSS-2a, SDSS-2b).

### 5. Catalog Ingestion & Kinematics (`src/data/`)
- **[remote_fits_streamer.js](src/data/remote_fits_streamer.js)**: Progressive HTTP Range-Request client streaming remote FITS primary cubes directly from IP2I Lyon and CDS Strasbourg archives.
- **IndexedDB Persistent Cache (`zrt_cosmicflows_cache_v1`)**: Automatic client-side binary blob storage for $128^3$ and $256^3$ grids, eliminating redundant network roundtrips.
- **[cf4_group_catalog.js](src/data/cf4_group_catalog.js)**: Ingestion for the 38,000 Grouped Galaxy Catalog (CF4gp), Beers biweight location/scale estimators, projected harmonic radius $R_H$, and multi-estimator Virial mass suite ($M_{\text{vir}}, M_{\text{proj}}, M_{\text{med}}, M_{\text{avg}}$).
- **[tfr_multiband_calibrator.js](src/data/tfr_multiband_calibrator.js)**: Multi-band Tully-Fisher relation calibrator (WISE W1/W2, Spitzer 3.6, SDSS $i$, 2MASS $K_s$), 21cm HI profile deprojection, dust extinction, and Malmquist bias compensators.
- **[cosmicflows_mock_generator.js](src/data/cosmicflows_mock_generator.js)**: Synthetic $\Lambda\mathrm{CDM}$ matter power spectrum $P(k)$ generator with Eisenstein & Hu transfer functions and Zel'dovich mock displacements.

### 6. Multipolar Bulk Flow & Uncertainty Quantification (`src/bulk-flow/`, `src/uncertainty/`, `src/statistics/`)
- **[multipolar_bulk_flow.js](src/bulk-flow/multipolar_bulk_flow.js)**: Radial shell decomposition ($R \in [20, 250]\,h^{-1}\mathrm{Mpc}$) with triple weighting estimators (`FIELD_VOLUME_WEIGHTED`, `CATALOG_WEIGHTED`, `INVERSE_VARIANCE_WEIGHTED`), spherical harmonic multipole expansion ($l=0$ to $4$), and window function deconvolution.
- **[bayesian_hmc_field.js](src/uncertainty/bayesian_hmc_field.js)**: CF4++ 10,000-step Hamiltonian Monte Carlo posterior mean $\bar{\mathbf{v}}(\mathbf{x})$ and RMS uncertainty field loader with conditional Gaussian realization sampling.
- **[covariance_regularizer.js](src/statistics/covariance_regularizer.js)**: Ledoit-Wolf and Oracle Approximating Shrinkage (OAS) estimators, condition number monitors, and compact Wendland support kernels.

### 7. Surface Geometry & Hydrodynamic Fluxes (`src/surfaces/`)
- **[exact_marching_tetrahedra.js](src/surfaces/exact_marching_tetrahedra.js)**: 6-tetrahedron and 5-tetrahedron simplicial cube decompositions resolving topological ambiguities, exact linear edge interpolation, and surface area $\iint dA$ / enclosed volume $\iiint dV$ calculus.
- **[watershed_manifold_mesher.js](src/surfaces/watershed_manifold_mesher.js)**: Triangulated 3D boundary interface extractor calculating hydrodynamic inter-basin mass/momentum flux integrals $\Phi_{AB} = \iint_{\partial\mathcal{B}_{AB}} (\mathbf{v} \cdot \hat{\mathbf{n}}) dA$.

### 8. Provenance, Units & Scientific Runtime (`src/export/`, `src/provenance/`, `src/runtime/`)
- **[provenance_jsonld_exporter.js](src/export/provenance_jsonld_exporter.js)**: W3C PROV-JSONLD execution bundle serializer with SHA-256 entity digests, automated LaTeX figure captions, and BibTeX citations.
- **[vo_table_fits_serializer.js](src/export/vo_table_fits_serializer.js)**: IVOA VOTable XML (v1.4) and Multi-Extension FITS (MEF) binary table serializer with Big-Endian byte packing and IEEE 32-bit checksums.
- **[memory_budget_manager.js](src/runtime/memory_budget_manager.js)**: VRAM and WebGL memory budgeting sentinel with LRU cache eviction and OOM guards for $64^3$, $128^3$, and $256^3$ grids.

---

## 📁 Repository Structure

```
zrt-cosmicflows-workbench/
├── data/                               # Local binary grids, metadata & research tables
│   ├── cf4_density_individual_64.bin   # CF4 Ungrouped Density (64³)
│   ├── cf4_velocity_individual_64.bin  # CF4 Ungrouped Velocity (64³)
│   ├── cf4_density_grouped_64.bin      # CF4 Grouped Density (64³)
│   ├── cf4_velocity_grouped_64.bin     # CF4 Grouped Velocity (64³)
│   ├── cf4_grid_meta.json              # Physical grid extents & bulk flows
│   └── cf4_research_data.json          # Table A.1 basin metrics & kinematics
├── src/                                # Core Scientific Computing Engine (ES6 Modules)
│   ├── bulk-flow/                      # Spherical harmonic multipole bulk flow estimators
│   ├── coordinates/                    # Astrometric WCS, CMB frames & dimensional guards
│   ├── data/                           # FITS parsers, group catalogs, TFR calibrator & streamer
│   ├── export/                         # W3C PROV-JSONLD, VOTable & MEF FITS exporters
│   ├── fields/                         # Okubo-Weiss, Helmholtz-Hodge, dispersion & tidal tensors
│   ├── integration/                    # Cash-Karp RK45, DOPRI5 & Yoshida symplectic integrators
│   ├── interpolation/                  # 64-point Tricubic Hermite splines & finite difference
│   ├── provenance/                     # Lineage graphs, W3C PROV-O & SHA-256 digesters
│   ├── runtime/                        # Distributed worker pools & memory budget sentinels
│   ├── statistics/                     # Covariance shrinkage, OAS, power spectra & resampling
│   ├── streamlines/                    # Poincaré sections, Lyapunov estimators & massive streamers
│   ├── surfaces/                       # Exact Marching Tetrahedra & watershed boundary meshers
│   ├── topology/                       # Betti persistent homology, Morse-Smale & Newton-Raphson
│   ├── uncertainty/                    # Bayesian HMC, NUTS sampler & Gelman-Rubin diagnostics
│   ├── units/                          # Strict physical dimension tensors & unit converters
│   ├── validation/                     # Multi-method cross-validation oracle & acceptance gates
│   └── watershed/                      # Watershed segmentation & Table A.1 classifiers
├── tests/                              # Automated Pytest & CDP Browser Test Suite (1,032 tests)
│   ├── analytic-fields/                # Analytical solutions & benchmark fields
│   ├── bulk-flow/                      # Multipole bulk flow validation tests
│   ├── coordinates/                    # WCS & reference frame transformation tests
│   ├── data/                           # FITS streamer, group catalog & TFR tests
│   ├── export/                         # PROV-JSONLD, VOTable & LaTeX export tests
│   ├── fields/                         # Okubo-Weiss, Helmholtz FFT & dispersion tensor tests
│   ├── integration/                    # RK45, DOPRI5 & symplectic integrator tests
│   ├── interpolation/                  # C^1 continuity, gradient & Hessian accuracy tests
│   ├── runtime/                        # Worker pool & memory manager tests
│   ├── statistics/                     # Covariance regularizer & moment estimator tests
│   ├── streamlines/                    # Poincaré section & Lyapunov exponent tests
│   ├── surfaces/                       # Marching Tetrahedra & flux integration tests
│   ├── topology/                       # Betti numbers, persistent homology & Morse-Smale tests
│   ├── uncertainty/                    # Bayesian HMC field & MCMC diagnostics tests
│   ├── validation/                     # Cross-validation oracle tests
│   ├── watershed/                      # Watershed classification & Table A.1 tests
│   ├── test_cosmo_time.py              # CDP browser tests: Friedmann scale factor & growth
│   └── test_erosita_overlays.py        # CDP browser tests: eROSITA WHIM bridge gas shaders
├── workers/                            # Dedicated Web Workers (Zero-copy parallel execution)
│   ├── fits_decompressor_worker.js     # Big-Endian FITS primary array decoder
│   ├── grid_field_worker.js            # Background field interpolation worker
│   └── streamline_worker.js            # Parallel streamline integration worker
├── tools/                              # Auditing and utility scripts
│   └── count-meaningful-loc.py         # Non-HTML meaningful LOC auditor (87,500+ LOC)
├── index.html                          # Interactive Three.js WebGL Research Workbench UI
└── package.json                        # ES6 module package definition
```

---

## 🚀 Quickstart & Local Hosting

### 1. Prerequisites
- **Python 3.10+** (with `pytest`, `pytest-cov`, `numpy`, `scipy`)
- **Node.js 18+** (for ES6 module execution)
- Modern web browser with WebGL 2.0 support (Chrome, Edge, Firefox, Safari)

### 2. Launch Local Workbench Server
Start the local server using Python's built-in HTTP server:
```bash
python -m http.server 8000
```
Open **[http://localhost:8000/](http://localhost:8000/)** in your browser.

---

## 🧪 Running the Test Suite

Run the full automated test suite (1,032 tests) using `pytest`:

```bash
# Run all unit and integration test suites
python -m pytest tests/ -v

# Run specific subsystem test suites
python -m pytest tests/fields/ -v
python -m pytest tests/topology/ -v
python -m pytest tests/data/ -v
python -m pytest tests/streamlines/ -v

# Run CDP headless browser visual & shader tests
python -m pytest tests/test_cosmo_time.py tests/test_erosita_overlays.py -v
```

---

## 📊 Meaningful LOC Audit

Verify the scientific codebase volume using the automated auditor:
```bash
python tools/count-meaningful-loc.py
```
Output:
```
====================================================================
ZRT COSMICFLOWS RESEARCH WORKBENCH - MEANINGFUL NON-HTML LOC AUDIT
====================================================================
Subsystem / Module Path                    |  Meaningful LOC
--------------------------------------------------------------------
src/fields                                 |           8,126
src/data                                   |           7,488
src/topology                               |           4,081
src/coordinates                            |           3,346
src/export                                 |           2,969
src/streamlines                            |           2,915
src/units                                  |           2,777
src/surfaces                               |           2,750
src/validation                             |           2,604
src/runtime                                |           2,507
src/bulk-flow                              |           2,280
src/statistics                             |           1,949
src/integration                            |           1,831
src/provenance                             |           1,487
src/uncertainty                            |           4,427
src/watershed                              |           1,074
src/interpolation                          |           1,095
tests/                                     |          30,176
workers/                                   |             819
--------------------------------------------------------------------
TOTAL MEANINGFUL NON-HTML LOC              |          87,516
TOTAL NUMBER OF CODE FILES                 |             207
====================================================================
[PASS] LOC REQUIREMENT MET: 87,516 >= 80,000
```

---

## 📜 Scientific Invariants & Attribution

All computations using official Cosmicflows-4 data products strictly follow the invariants in [`AGENTS.md`](AGENTS.md) and [`GEMINI.md`](GEMINI.md).

### Required Citations

When using or citing results produced by this workbench:

1. **Courtois et al. (2023)**:
   > Courtois, H. M., Dupuy, A., Guinet, D., et al. (2023). *Cosmicflows-4: The catalog of 56,000 galaxy distances and peculiar velocities*. Astronomy & Astrophysics, 670, L15. [DOI: 10.1051/0004-6361/202245331](https://doi.org/10.1051/0004-6361/202245331)

2. **Dupuy & Courtois (2023)**:
   > Dupuy, A., & Courtois, H. M. (2023). *Cosmicflows-4: Cosmography and Watershed Basins of Attraction*. Astronomy & Astrophysics, 678, A176. [DOI: 10.1051/0004-6361/202346802](https://doi.org/10.1051/0004-6361/202346802)

3. **Hoffman et al. (2024)**:
   > Hoffman, Y., Courtois, H. M., Tully, R. B., et al. (2024). *The Cosmicflows-4 Wiener Filter Reconstruction of the Local Universe*. Monthly Notices of the Royal Astronomical Society, 527(4), 10327–10340. [DOI: 10.1093/mnras/stad3782](https://doi.org/10.1093/mnras/stad3782)

---

## ⚖️ License

Distributed under the **MIT License**. See `LICENSE` for details.
