# Project: CosmicFlows-4 Research Workbench (WebGL 3D Optimization & Performance Hardening)

## Architecture
The CosmicFlows-4 Research Workbench is a client-side high-performance scientific visualization and computing environment built in vanilla ES6 JavaScript, WebGL/Three.js, HTML5/CSS3, with extensive automated Python (pytest + Chrome DevTools Protocol) testing and verification suites.

### Core Subsystems
1. **60 FPS Streamline Particle Advection & 4D Cosmic Time Engine (`rebuildStreamlines`, `applyCosmicTime`)**:
   - High-throughput particle advection along 4,200 cosmic streamlines using `THREE.DynamicDrawUsage` buffers.
   - Structure-of-Arrays (SoA) TypedArray particle state (`particleProgress`, `particleSpeed`, `particleCurveIdx`) and flattened streamline vertex curves (`flatVertices`, `offsets`, `lengths`).
   - 4D Cosmic time evolution (-13.8 Gyr to +10.0 Gyr) with spatial bounding-box early-exit pruning for attractor sink gravitational decays, discrete user haptics, and throttled DOM telemetry updates.
2. **WebGL Shader Lifecycle & Zero-Leak Memory Manager (`disposeHierarchy`, `disposeMaterial`)**:
   - Comprehensive shader lifecycle management for `UnrealBloomPass`, `EffectComposer`, `LuminosityHighPassShader`, and resolution updating on resize.
   - Recursive disposal across all 15 visual subsystems (halos, ruler box, points, watersheds, ZOA, range rings, callouts, streamlines, particles, bulk dipole, quiver field, eROSITA overlays, TFR probe, beacons, compute overlays) and sweeping all 13 material texture map slots.
   - Zero GPU memory growth across 10x/30x theme switches and dataset switches with `webglcontextlost`/`webglcontextrestored` resilience.
3. **Zero-Lag Canvas Gestures & Pointer Interaction Liveness (`initKeyboardCameraNav`, pointer listeners)**:
   - OrbitControls camera manipulation with strict drag disambiguation ($\Delta r \ge 5\text{ px}$, $\Delta t \ge 350\text{ ms}$) distinguishing camera drag from discrete cluster taps.
   - Drag-aware `pointermove` listener disabling hover card raycasting during camera rotation, and `pointercancel` cleanup.
   - Comprehensive modal exclusion filter in `pointerup` covering all modal IDs (`#primer-modal`, `#video-player-modal`, `#camera-settings-modal`, `#shortcuts-modal`, `#error-dossier-modal`, `#error-hud-pill`, `#hover-card`, `#dark-telemetry-hud`).
   - Statically pre-allocated WASD navigation vectors (`_nudgeFwd`, `_nudgeRight`, `_nudgeUp`) eliminating 180 heap allocations/sec.
4. **Search & Multi-Catalog Query Optimization (`initCosmicSearchEngine`, `query56kGalaxies`)**:
   - Pre-computed normalized entity search keys and token indexes at catalog indexing time, eliminating runtime string conversions and dynamic RegExp compilation.
   - Fast $O(1)$ deterministic multi-scale coordinate and morphology hashing for all 126,000 PGC galaxies, NGC, IC, UGC, Messier, and Abell catalogs executing in < 0.25 ms with zero main-thread blocking.
5. **Quality Gates & Test Harness**:
   - 6-gate pre-commit pipeline (`scripts/pre_commit_gate.py`), KaTeX invariant validator, JS TDZ linter, synthetic canvas drag validator, headless CDP boot runner, and Pytest test suites (167 automated tests passing 100%).

## Feature Inventory
| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| F1 | Dynamic Draw Usage Buffer Optimization | Set `usage = THREE.DynamicDrawUsage` on `flowParticlesMesh.geometry.attributes.position` and `galaxyPointsMesh.geometry.attributes.position`. | M1 | Survey 1 / R1 |
| F2 | TypedArray SoA Particle State & Flattened Curves | Replace particle object arrays with flat `Float32Array` / `Int32Array` buffers and flatten streamline curves with integer offset table. | M1 | Survey 1 / R1 |
| F3 | 4D Cosmic Time Engine Optimization | Bounding-box early exit for sink attractor math, discrete haptic trigger guarding, and throttled DOM telemetry updates. | M1 | Survey 1 / R1 |
| F4 | WebGL Shader & Post-Processing Lifecycle | Resilient `UnrealBloomPass` and `EffectComposer` lifecycle with Dark/White mode conditional execution and dynamic resize handling. | M2 | Survey 2 / R2 |
| F5 | Comprehensive Resource Disposal & Theme Invariance | Recursive `disposeHierarchy()` sweeping all 15 visual subsystems and 13 material texture slots with zero VRAM leaks across dataset/theme switches. | M2 | Survey 2 / R2 |
| F6 | WebGL Context Loss & Restoration Handlers | Robust `webglcontextlost` and `webglcontextrestored` event handling reconstructing all primary scene layers. | M2 | Survey 2 / R2 |
| F7 | Strict Drag Disambiguation & Pointer Liveness | Enforce $\Delta r \ge 5\text{ px}$ / $\Delta t \ge 350\text{ ms}$ click-vs-drag threshold and add `pointercancel` cleanup. | M3 | Survey 3 / R3 |
| F8 | Drag-Aware PointerMove & Modal Filter Expansion | Disable hover raycasting during camera drags and expand `closest(...)` modal exclusion filter across all overlay modals. | M3 | Survey 3 / R3 |
| F9 | WASD Navigation GC Elimination | Pre-allocate static vectors `_nudgeFwd`, `_nudgeRight`, `_nudgeUp` in `nudgeCamera()` eliminating 180 allocations/sec. | M3 | Survey 3 / R3 |
| F10 | Pre-Computed Normalized Entity Search Index | Pre-compute normalized search strings on entity objects at index creation to eliminate runtime regex compilation. | M4 | Survey 3 / R4 |
| F11 | 126k Galaxy Catalog Sub-5ms Query Execution | Deterministic $O(1)$ coordinate and morphology hashing for all 126k astronomical catalog objects with micro-debounced input. | M4 | Survey 3 / R4 |
| F12 | Master Pre-Commit Quality Gate & Test Suite Pass | 100% pass across all 6 automated pre-commit gates (`python scripts/pre_commit_gate.py`), 167 automated Pytest cases, 0 console errors, 0 KaTeX parse errors. | M5 | Survey 3 / AC |

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| M1 | 60 FPS Streamline Particle Advection & 4D Cosmic Time Engine | Features F1, F2, F3 | none | DONE |
| M2 | WebGL Shader Lifecycle & Zero-Leak Memory Management | Features F4, F5, F6 | none | DONE |
| M3 | Zero-Lag Canvas Gestures & Pointer Interaction Liveness | Features F7, F8, F9 | none | DONE |
| M4 | Search & Multi-Catalog Performance Hardening | Features F10, F11 | none | DONE |
| M5 | Dual-Track E2E & Pre-Commit Quality Gate Pass | Feature F12 | M1, M2, M3, M4 | DONE |

## Interface Contracts

### Particle Advection & 4D Time ↔ Three.js Render Loop
- `posAttr.setUsage(THREE.DynamicDrawUsage)`
- `applyCosmicTime(tGyr: number): void`
- `rebuildStreamlines(): void`

### Memory Manager ↔ Visual Subsystems
- `disposeHierarchy(obj: THREE.Object3D): void`
- `disposeMaterial(mat: THREE.Material): void`
- `toggleTheme(targetTheme?: string): void`

### Canvas Pointer Interaction ↔ Camera / Modals
- `pointerup`: $\Delta r < 5\text{ px} \land \Delta t < 350\text{ ms} \land !\text{closest}(\text{MODAL\_SELECTOR})$ triggers cluster raycast picking.
- `pointermove`: skips raycasting when camera is being dragged.

### Search Engine ↔ Astrometric Data
- `window.cosmicflows.search.query(text: string): EntityMatch[]` (execution time $< 0.3\text{ ms}$)
- `window.cosmicflows.search.select(entityId: string): void`

## Code Layout
- `index.html`: Primary UI markup, WebGL render pipeline, astrometric datasets, search engine, tips engine, dossier renderer, video modal controller.
- `src/runtime/memory_budget_manager.js`: Enterprise multi-tier LRU cache and memory budgeting engine.
- `scripts/pre_commit_gate.py`: Master pre-commit quality gate (AST TDZ, headless CDP boot, canvas drag, dual-theme PIL variance, splash lifecycle, KaTeX parser).
- `scripts/check_katex_invariants.py`: In-browser and AST KaTeX mathematical formula validator.
- `scripts/check_js_tdz.py`: AST-based JavaScript variable scope and TDZ invariant checker.
- `scripts/check_canvas_drag.py`: Synthetic pointer event canvas drag and click disambiguation validator.
- `tests/`: Automated pytest test suites (cosmological time, coordinates, fields, streamlines, memory budget manager, search, etc.).
- `tools/count-meaningful-loc.py`: Non-HTML meaningful line-of-code counter.
