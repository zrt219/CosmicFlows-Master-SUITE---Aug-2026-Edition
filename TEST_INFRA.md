# E2E Test Infra: Cosmicflows Master Suite Upgrade

## Test Philosophy
- Opaque-box & transparent programmatic API verification via Chrome DevTools Protocol (CDP) and Pytest.
- Full requirement coverage across R1 (Time Evolution), R2 (Spectroscopy Dossiers), R3 (eROSITA Overlays), and Acceptance Criteria.

## Feature Inventory & Test Mapping
| # | Feature | Requirement | Tier 1 (Unit) | Tier 2 (Boundary) | Tier 3 (Cross) | Tier 4 (Scenario) |
|---|---------|-------------|:-------------:|:-----------------:|:--------------:|:-----------------:|
| 1 | $\Lambda\text{CDM}$ Time Mapping | R1 | ✓ ($t=0, -13.7, +10$) | ✓ ($z \to \infty, a \to 0$) | ✓ (with RK45) | ✓ (Cosmic History) |
| 2 | Linear Growth Factor $D_+(z)$ | R1 | ✓ ($D_+(0)=1, D_+(-13.7)<0.01$) | ✓ (high-z limits) | ✓ (particle scaling) | ✓ (de-clustering) |
| 3 | Past Expansion De-Clustering | R1 | ✓ (displacement) | ✓ ($t=-13.78$) | ✓ (multi-catalog) | ✓ (early universe) |
| 4 | Future Gravitational Collapse | R1 | ✓ (Shapley infall) | ✓ ($t=+10\text{ Gyr}$) | ✓ (Great Attractor) | ✓ (virialization) |
| 5 | Zero-Stutter Scrubber | R1 | ✓ (slider event) | ✓ (rapid drag) | ✓ (telemetry HUD) | ✓ (smooth 60fps) |
| 6 | 3D Raycasting Hit-Testing | R2 | ✓ (all 8 hubs) | ✓ (edge clicks) | ✓ (with orbit controls)| ✓ (cluster inspection)|
| 7 | Spectroscopy Dossier Modal | R2 | ✓ (DOM elements) | ✓ (unpinned/pinned)| ✓ (theme switch) | ✓ (scientific report)|
| 8 | Velocity Dispersion Canvas | R2 | ✓ (canvas draw) | ✓ (extreme $\sigma_v$)| ✓ (resize/theme) | ✓ (Gaussian profile)|
| 9 | eROSITA WHIM Gas Bridges | R3 | ✓ (3 bridges) | ✓ (toggle on/off) | ✓ (gas halos) | ✓ (multi-wavelength)|
| 10 | Dual-Theme eROSITA Overlays | R3 | ✓ (Dark & White) | ✓ (rapid toggle) | ✓ (UnrealBloomPass) | ✓ (publication capture)|
| 11 | Zero-Leak Disposal | Acc | ✓ (disposeHierarchy)| ✓ (rapid reload) | ✓ (theme flip) | ✓ (memory stress) |
| 12 | Headless CDP Integration | Acc | ✓ (CDP client) | ✓ (async commands) | ✓ (screenshots) | ✓ (100% pytest) |

## Test Architecture
- Test Runner: `pytest tests/ -v`
- CDP Test Harness: `tests/test_harness.py`
- Test Files:
  - `tests/test_cosmo_time.py`
  - `tests/test_spectroscopy_dossiers.py`
  - `tests/test_erosita_overlays.py`
  - `tests/test_m1.py`, `tests/test_tier1_features.py`, etc.
