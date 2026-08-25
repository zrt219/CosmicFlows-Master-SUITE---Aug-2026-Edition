# Project: Cosmicflows Master Suite Cosmological Upgrade

## Architecture
The Cosmicflows Master Suite is an astrophysical cosmography visualization and simulation platform running Three.js WebGL and vector calculus engines. The upgrade incorporates three major subsystems:
1. **Cosmological Time Evolution Engine ($t \in [-13.8\text{ Gyr}, +10\text{ Gyr}]$)**:
   - Analytical $\Lambda\text{CDM}$ scale factor inversion $a(t_{\text{cosmic}}) = \left(\sqrt{\frac{\Omega_m}{\Omega_\Lambda}}\sinh\left(\frac{3}{2}H_0\sqrt{\Omega_\Lambda}t_{\text{cosmic}}\right)\right)^{2/3}$ mapping continuous cosmic time to scale factor $a$ and redshift $z$.
   - Carroll-Press-Turner / Eisenstein-Hu linear growth factor $D_+(z) = D(a)/D(1)$.
   - Zel'dovich Lagrangian displacement de-clustering for past epochs ($t < 0$, $z \to \infty$) returning galaxy point clouds toward a smooth primordial distribution.
   - Non-linear gravitational collapse for future epochs ($t > 0$, $t \to +10\text{ Gyr}$) into Shapley Supercluster Core and Great Attractor megasinks.
   - In-place single-pass CPU vectorized buffer mutation ($<0.8\text{ ms}$ for 56k points) ensuring 0 frame stutter during continuous scrubber interaction.
   - UI time scrubber slider `[-13.8, +10.0] Gyr` with telemetry HUD, preset epoch buttons, and transport play/pause.
   - Global programmatic API `window.cosmicflows.timeEngine`.
2. **Interactive Galaxy Cluster Spectroscopy Dossiers**:
   - 3D raycast selection on major galaxy clusters and superclusters (Virgo, Coma, Centaurus, Norma, Shapley, Vela, Perseus-Pisces, Fornax).
   - Scientific inspector card modal (`#spectroscopy-modal`) displaying virial mass ($M_{200}$), recessional velocity ($cz$), velocity dispersion ($\sigma_v$), peculiar motion components $(v_x, v_y, v_z)$, and primary references.
   - Interactive Gaussian velocity dispersion distribution canvas ($N(v)$) with $1\sigma / 2\sigma$ dispersion bands.
   - Dual-theme modal styling (Dark glassmorphism vs White publication vector mode) with pin/close/hover interactions.
   - Programmatic API `window.cosmicflows.spectroscopy`.
3. **eROSITA Multi-Wavelength Hot Gas Overlays**:
   - Volumetric translucent warm-hot intergalactic medium (WHIM) X-ray gas halos and filamentary bridges connecting major cluster pairs (Coma-Virgo, Shapley-Centaurus, Perseus-Pisces).
   - Dual-theme rendering: White Mode (charcoal stippled gas) and Dark Mode (luminescent emerald/violet plasma glow).
   - Smooth toggle controls (`chk-erosita-gas`, `chk-whim-bridges`, `sel-erosita-band`) and zero-leak memory management.
4. **Automated Headless CDP Regression Test Suite**:
   - Multi-tier automated test suite covering time evolution math, spectroscopy dossiers, eROSITA overlays, WebGL memory stability, and dual-theme publication exports.

---

## Feature Inventory

| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| F1 | $\Lambda\text{CDM}$ Time-to-Redshift Inversion | Exact analytical $a(t)$ and $z(t)$ calculation across $t \in [-13.8, +10]\text{ Gyr}$ | M_UPGRADE | ORIGINAL_REQUEST §R1, Survey |
| F2 | Linear Growth Factor $D_+(z)$ | Carroll-Press-Turner analytical linear growth factor $D_+(z)$ | M_UPGRADE | ORIGINAL_REQUEST §R1, Survey |
| F3 | Past Expansion De-Clustering | Zel'dovich de-clustering dissolving structure into homogeneous primordial field as $z \to \infty$ | M_UPGRADE | ORIGINAL_REQUEST §R1, Survey |
| F4 | Future Gravitational Collapse | Non-linear accretion & virial condensation into Shapley Core and Great Attractor sinks | M_UPGRADE | ORIGINAL_REQUEST §R1, Survey |
| F5 | Zero-Stutter Scrubber & UI HUD | Continuous time slider with real-time HUD telemetry, presets, transport controls, and zero frame stutter | M_UPGRADE | ORIGINAL_REQUEST §R1, Survey |
| F6 | 3D Raycasting Cluster Selection | Raycasting hit-test selection on Virgo, Coma, Centaurus, Norma, Shapley, Vela, Perseus-Pisces, Fornax | M_UPGRADE | ORIGINAL_REQUEST §R2, Survey |
| F7 | Spectroscopy Dossier Modal | Scientific inspector card modal displaying $M_{200}$, $cz$, $\sigma_v$, $(v_x, v_y, v_z)$, and citations | M_UPGRADE | ORIGINAL_REQUEST §R2, Survey |
| F8 | Gaussian Velocity Dispersion Canvas | Interactive $N(v)$ profile canvas with $1\sigma / 2\sigma$ dispersion bands | M_UPGRADE | ORIGINAL_REQUEST §R2, Survey |
| F9 | eROSITA WHIM Gas Bridges & Halos | Volumetric translucent X-ray gas halos & filamentary bridges (Coma-Virgo, Shapley-Centaurus, Perseus-Pisces) | M_UPGRADE | ORIGINAL_REQUEST §R3, Survey |
| F10 | Dual-Theme eROSITA Overlays | Charcoal stippled gas in White Mode & luminescent emerald/violet plasma in Dark Mode | M_UPGRADE | ORIGINAL_REQUEST §R3, Survey |
| F11 | Zero-Leak Disposal & Memory Mgmt | Comprehensive `disposeHierarchy` for all new meshes, geometries, materials, and textures | M_UPGRADE | ORIGINAL_REQUEST §R3, Survey |
| F12 | Programmatic API Bridges | `window.cosmicflows.timeEngine` and `window.cosmicflows.spectroscopy` and `window.cosmicflows.erosita` | M_UPGRADE | ORIGINAL_REQUEST Acceptance Criteria |
| F13 | E2E Headless CDP Test Suite | Comprehensive Pytest suite testing R1, R2, R3, memory leaks, and publication exports | E2E_Track | ORIGINAL_REQUEST Acceptance Criteria |
| F14 | Gate Verification & Audit | Multi-reviewer, challenger, and forensic auditor verification | M_VERIFY | Acceptance Criteria |

---

## Milestones

| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| E2E | E2E Test Suite Track | Build automated Pytest CDP test suite for R1, R2, R3, publish `TEST_READY.md` | none | IN_PROGRESS |
| M_UPGRADE | Cosmicflows Master Upgrade | Implement R1 Time Engine, R2 Spectroscopy Dossiers, R3 eROSITA Gas Overlays & Shaders in `index.html` | none | IN_PROGRESS |
| M_VERIFY | Multi-Agent Gate & Audit | Reviewers APPROVE, Challengers verify, Forensic Auditor CLEAN audit | E2E, M_UPGRADE | PLANNED |

---

## Interface Contracts

### Global Programmatic Bridge: `window.cosmicflows`
```javascript
window.cosmicflows = {
  version: '2026.2',
  simState: simState,
  LANDMARKS: LANDMARKS,
  ENGINES: ENGINES,
  ASTROMETRIC_FEATURES: ASTROMETRIC_FEATURES,
  timeEngine: {
    getTime: function() { return simState.cosmicTime; },
    setTime: function(tGyr) { ... },
    computeCosmology: function(tGyr) { return { a, z, D_plus, f_growth, H_z, Omega_m_z, Omega_L_z }; },
    play: function() { ... },
    pause: function() { ... },
    setSpeed: function(spd) { ... }
  },
  spectroscopy: {
    getDossier: function(clusterId) { return { name, M200, cz, sigma_v, v_pec, kBT, LX, references }; },
    openDossier: function(clusterId) { ... },
    closeDossier: function() { ... },
    getAllHubs: function() { return ['virgo_cl', 'coma_cl', 'centaurus_cl', 'ga_norma', 'shapley_core', 'vela_scl', 'perseus_cl', 'fornax_cl']; }
  },
  erosita: {
    setOverlayEnabled: function(enabled) { ... },
    setBridgesEnabled: function(enabled) { ... },
    setEnergyBand: function(band) { ... }, // 'soft_0.2_0.6' | 'mid_0.6_2.3' | 'hard_2.3_5.0'
    getOverlayState: function() { return { enabled, bridgesEnabled, band }; }
  }
};
```

---

## Code Layout
- `index.html`: Complete Master Suite application with Time Evolution Engine, Spectroscopy Modal, eROSITA Overlays, and WebGL renderer.
- `tests/test_harness.py`: Automated headless CDP test runner.
- `tests/test_cosmo_time.py`: Tier 1-4 tests for Cosmological Time Evolution Engine.
- `tests/test_spectroscopy_dossiers.py`: Tier 1-4 tests for Spectroscopy Dossiers & Raycasting.
- `tests/test_erosita_overlays.py`: Tier 1-4 tests for eROSITA Hot Gas Overlays & Dual-Theme Shaders.
- `tests/`: Full Pytest regression suite.
- `.agents/`: Agent metadata and execution journals.
