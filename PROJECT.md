# Project: CosmicFlows-4 Research Workbench (Follow-up Enhancements)

## Architecture
The CosmicFlows-4 Research Workbench is a client-side high-performance scientific visualization and computing environment built in vanilla ES6 JavaScript, WebGL/Three.js, HTML5/CSS3, with extensive automated Python (pytest + Chrome DevTools Protocol) testing and verification suites.

### Core Subsystems
1. **Search & Astrometric Catalog Engine (`initCosmicSearchEngine`)**:
   - Universal search across 38 cosmological structures, 8 CF4 supercluster basins, and comprehensive multi-catalog cross-references (Abell, ACO, Messier, NGC, PGC, UGC).
   - Space/punctuation-insensitive fuzzy matching, keyboard arrow navigation, instant category badge rendering, smooth 3D camera fly-to, expanding targeting beacon, and automatic inspection dossier triggering.
2. **Cosmography Explainer & Contextual Insights (`initCosmicTips`)**:
   - Floating glassmorphic HUD pill cycling through 10 university-grade cosmological explanations with mathematical rigor.
   - Clean Mode (`C` key) fade out, hover pause, theme resilience, and explore drawer restoration hook.
   - Universal relative distance comparators dynamically calculated against the Virgo Cluster baseline ($d_{\text{Virgo}} = 16.5\,\text{Mpc} = 53.8\,\text{Mly}$).
3. **Spectroscopy Dossier Modal & Fallback Synthesizer (`openDossier`)**:
   - Complete physical and observational profiles for all structures, including synthetic Gaussian velocity dispersion profiles, coordinates, distances, and reference media.
4. **In-App Video Modal Player (`openCosmicVideoModal`, `closeCosmicVideoModal`)**:
   - Responsive 16:9 aspect-ratio video player with `backdrop-filter: blur(20px)`, YouTube no-cookie embedding, ESC key trap, backdrop click-outside dismissal, and audio-leak prevention.
   - Verified scientific literature DOIs, Wikipedia links, and curated video documentaries.
5. **Quality Gates & Test Harness**:
   - 6-gate pre-commit pipeline (`scripts/pre_commit_gate.py`), KaTeX invariant validator, JS TDZ linter, synthetic canvas drag validator, headless CDP boot runner, and Pytest test suites (>1,000 automated tests).

## Feature Inventory
| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| F1 | Multi-Catalog Alias Expansion | Map Abell (A1656, A3558, A426, etc.), ACO (ACO 3627, ACO 3526), Messier (M87, M49, M86), NGC (NGC 4889, NGC 4874, NGC 4696), PGC, and UGC IDs across all 38 astrometric features and 8 basins. | M1 | Survey 1 / R1 |
| F2 | Normalizing Fuzzy Search Matching | Space-, hyphen-, and underscore-insensitive matching with instant search dropdown updates. | M1 | Survey 1 / R1 |
| F3 | Keyboard Arrow Navigation & Enter Fly-to | Up/Down arrow selection, auto-scrolling list, Enter key camera fly-to, 3D beacon activation, and dossier auto-open. | M1 | Survey 1 / R1 |
| F4 | Cosmography Tips Database Expansion | 10 curated university-grade explanations (Peculiar Velocity, Hubble Flow, Cosmic Watersheds, Virial Mass, Dipole Repeller, Zone of Avoidance, Tully-Fisher, Cosmic V-Web, Wiener Filter, Okubo-Weiss). | M2 | Survey 2 / R2 |
| F5 | Floating Explainer Pill Lifecycle | 16s auto-cycle, hover pause, Clean Mode fade, theme toggle resilience, and drawer restore hook. | M2 | Survey 2 / R2 |
| F6 | Universal Relative Distance Comparators | Format and display `$d\text{ Mpc} (d \times 3.26\text{ Mly}) [\sim R\times \text{farther than Virgo}]$` across all cluster inspection cards and dossiers. | M2 | Survey 2 / R2 |
| F7 | Dynamic Spectroscopy Dossier Fallback | Support full dossier inspection cards for all 38 astrometric structures and dynamically searched objects. | M2 | Survey 2 / R2 |
| F8 | In-App 16:9 Video Modal Player Polish | Responsive 16:9 container, backdrop-filter blur, YouTube no-cookie embed, ESC key listener, and backdrop click-outside dismissal. | M3 | Survey 3 / R3 |
| F9 | Reference Media & DOI Verification | Verified active study DOIs (Courtois2023, Dupuy2023, Tully2014, Hoffman2017), Wikipedia URLs, and YouTube video IDs with fallback search triggers. | M3 | Survey 3 / R3 |
| F10 | Master Pre-Commit Quality Gate & Test Suite Pass | 100% pass across all 6 automated pre-commit gates (`python scripts/pre_commit_gate.py`), >= 1,000 automated Pytest cases, 0 console errors, 0 KaTeX parse errors. | M4 | Survey 3 / AC |

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| M1 | Multi-Catalog Search & Alternate Designations | Features F1, F2, F3 | none | PLANNED |
| M2 | Interactive Cosmography Explainer & Distance Comparators | Features F4, F5, F6, F7 | M1 | PLANNED |
| M3 | In-App Video Modal Player & Reference Media Polish | Features F8, F9 | M2 | PLANNED |
| M4 | Dual-Track E2E & Pre-Commit Quality Gate Pass | Feature F10 | M1, M2, M3 | PLANNED |

## Interface Contracts

### Search Engine ↔ 3D Scene / Dossier
- `window.cosmicflows.search.query(text: string): EntityMatch[]`
  - Returns filtered array of matching structures and basins with highlighted alias matches.
- `window.cosmicflows.search.select(entityId: string): void`
  - Flies camera to `sgToThree(pos.x, pos.y, pos.z)`, creates targeting beacon, and calls `openDossier(entityId)`.

### Cosmography Explainer ↔ UI / State
- `window.cosmicflows.explainer.nextTip(): TipObject`
- `window.cosmicflows.explainer.setTip(index: number): void`
- `window.cosmicflows.explainer.getTips(): TipObject[]`
- `window.cosmicflows.explainer.toggle(visible?: boolean): void`

### Video Modal ↔ Spectroscopy Dossier
- `openCosmicVideoModal(title: string, embedUrl: string, externalUrl: string): void`
  - Sets `iframe.src = embedUrl`, displays `#video-player-modal`, sets title and external link href.
- `closeCosmicVideoModal(): void`
  - Clears `iframe.src = ''` and hides `#video-player-modal`.

## Code Layout
- `index.html`: Primary UI markup, WebGL render pipeline, astrometric datasets, search engine, tips engine, dossier renderer, video modal controller.
- `scripts/pre_commit_gate.py`: Master pre-commit quality gate (AST TDZ, headless CDP boot, canvas drag, dual-theme PIL variance, splash lifecycle, KaTeX parser).
- `scripts/check_katex_invariants.py`: In-browser and AST KaTeX mathematical formula validator.
- `scripts/check_js_tdz.py`: AST-based JavaScript variable scope and TDZ invariant checker.
- `scripts/check_canvas_drag.py`: Synthetic pointer event canvas drag and click disambiguation validator.
- `tests/`: Automated pytest test suites (coordinates, fields, streamlines, spectroscopy dossiers, etc.).
- `tools/count-meaningful-loc.py`: Non-HTML meaningful line-of-code counter.
