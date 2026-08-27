# Cosmological Field Dimensions, Scientific Semantics, LaTeX Invariants & Engineering Protocols

## 1. Physical Dimensions Are Strict
* Supergalactic Cartesian spatial positions are represented in **Mpc/h** when operating on the CosmicFlows reconstruction grids.
* Peculiar velocities are represented in **km/s**.
* Density contrast (\delta) is dimensionless.
* Never combine quantities with incompatible physical dimensions without an explicit transformation.
* Direct addition `position_Mpc_h + velocity_km_s` is strictly prohibited.
* Dimensionally invalid operations must fail loudly rather than be silently coerced.

## 2. CosmicFlows-4 Grid Ordering
For the official public CF4 density and velocity grids distributed by the CosmicFlows/IP2I project:
* **Raw file ordering:** `(SGZ, SGY, SGX)`
* **Canonical ZRT representation:** `(SGX, SGY, SGZ)`
All indexing, interpolation, derivative, visualization, streamline, watershed, and topology operations must pass through one tested canonical axis-mapping layer.

## 3. CF4 Velocity Scale Factor (Exact 52.0 Multiplier)
The official CosmicFlows project specifies that the values in the public CF4 velocity grids and associated error products must be multiplied by:
$$52.0$$
before being interpreted as physical peculiar velocities. Apply this scaling:
* only to validated CF4 velocity values and errors;
* exactly once.
Attempting to apply it twice must fail loudly.

## 4. Linear-Theory Continuity Diagnostic
Where linear perturbation theory is appropriate, evaluate consistency using:
$$\nabla \cdot \mathbf{v} \approx -H_0 f \delta$$
when $a=1$ and $H=H_0$. This is a linear-regime scientific validation diagnostic, not an exact identity in strongly non-linear collapsed regions.

## 5. Watershed IDs Are Reconstruction-Specific
Never use one universal hard-coded basin-ID table for every CosmicFlows watershed file. The official watershed files follow the ordering of **Dupuy & Courtois (2023), Table A.1**.

## 6. GitHub MathJax & KaTeX Inviolable Syntax Invariants
1. **Display Math Isolation**: Every `$$ ... $$` display equation MUST be isolated with dedicated empty blank lines before and after.
2. **Subscript / Superscript Macro Bracing**: Every command or macro in subscripts/superscripts MUST be enclosed in curly braces (`\nabla_{\mathbf{x}}`, `\sigma_{\theta}`, `\sigma_{\phi}`, `\int_0^{\infty}`, `\sigma_{\mu}`, `C^{\infty}`, `M_{\nu}`, `\Psi_{\epsilon}`, `T_{\sigma}`).
3. **Explicit Subscript Underscores**: Always include explicit underscores for variable labels: `\mathbf{v}_{\text{pot}}`, `\mathbf{v}_{\text{sol}}`, `\mathbf{v}_0`, `\mathbf{v}_1`, `E_{\text{pot}}`, `E_{\text{sol}}`, `\hat{v}_i^{\text{pot}}`, `\mathcal{P}_{ij}^{\parallel}`.
4. **Universal Matrix Syntax**: Use `\left( \begin{matrix} ... \end{matrix} \right)` or standard `\begin{matrix}` instead of `\begin{pmatrix}` with external subscripts.
5. **No Naked `aligned` Blocks in Boxes**: Avoid embedding `\begin{aligned}` directly inside `\boxed{...}` or unconfigured KaTeX display blocks without testing. Prefer clean separated display equations.
6. **Universal Operators**: Use `\mathrm{Tr}`, `\det`, `\mathrm{diag}`, `\mathrm{rank}`.

## 7. Universal Android APK Signing Protocol
1. Sideloadable Android APKs must support `minSdk = 21` and `targetSdk = 34`.
2. Every release APK must be signed with V1 (JAR) + V2 + V3 signature schemes:
   - Step 1: `jarsigner -verbose -sigalg SHA256withRSA -digestalg SHA-256 -keystore release.jks app.apk keyAlias`
   - Step 2: `zipalign -v 4 app.apk app-aligned.apk`
   - Step 3: `apksigner sign --ks release.jks --v1-signing-enabled true --v2-signing-enabled true --v3-signing-enabled true app-aligned.apk`
3. Verify with `apksigner verify --verbose` asserting `Verified using v1: true, v2: true, v3: true`.

## 8. Strict GPU Render Barrier & Splash Screen Lifecycle
1. Splash loading screens must maintain DOM presence in `<body>` and enforce a minimum display duration threshold ($\ge 1.2\,\text{s}$).
2. Always explicitly dismiss splash loading screens (`splashController.dismiss(true)`) prior to capturing automated screenshots.
3. Force explicit render passes (`composer.render()` or `renderer.render(scene, camera)`).
4. Validate non-empty WebGL point cloud buffers (`geometry.attributes.position.count > 1000`).
5. Validate image output using PIL (`size_kb > 35` and pixel std variance $\sigma > 12.0$) to guarantee zero blank or unloaded frames.

## 9. Dynamic Viewport Safe-Area & Floating HUD Clearance
1. Floating sidebars or scientific drawer panels must never overlap static header branding or promo pills.
2. Hook `resize` and `orientationchange` events to compute dynamic top offsets:
   `top = branding.offsetHeight + promoPills.offsetHeight + 24px`
3. Provide a minimum vertical clearance margin of $+35\text{px}$ across all responsive viewport widths.

## 10. Strict Deployment Target Filtering ("No Vercel" Directive)
1. When the user specifies "no vercel", "github only", or targets a single artifact, immediately terminate background deployment tasks to unwanted platforms.
2. Suppress automatic cloud deployment triggers for the affected workflow.

## 11. Windows CRLF Line Ending & Safe String Replacements
1. When modifying multi-line text files on Windows, always normalize `\r\n` to `\n` prior to string searching or multi-line regex matching.
2. Always pass replacement functions `lambda m: replacement` into `re.sub()` to prevent python template escape errors.

## 12. Scientific Provenance & Required Attribution
* Every dataset, intermediate field, and publication export must carry SHA-256 digests and W3C PROV-O lineage.
* Citing CF4 density/velocity products requires: **Courtois et al. (2023), A&A 670, L15** (DOI: `10.1051/0004-6361/202245331`).
* Citing CF4 watershed products requires: **Dupuy & Courtois (2023), A&A 678, A176** (DOI: `10.1051/0004-6361/202346802`).

## 13. JavaScript Variable Scope, TDZ Invariants & Boot Error Trapping
1. In all mathematical inversion, solver, and utility functions in JavaScript (`solveTullyFisher`, `jacobiDiagonalize`, `integrateRK4`), all intermediate quantities (`sigmaMu`, `mAbs`, `mu`, `wCorr`) MUST be declared and initialized strictly prior to any calculation that consumes them.
2. Never rely on variable hoisting or compute derived statistical error adjustments before their root parameter definitions.
3. Core WebGL bootstrap routines (`init()`, `animate()`, `rebuildStreamlines()`, `splashController.dismiss()`) must reside inside defensive try-catch error barriers so that any single non-critical formula or dataset parsing failure never blocks the 3D canvas animation loop or mouse/touch event listeners.

## 14. Mandatory Pre-Commit Zero-Exception Gate
1. Prior to committing code or publishing release builds, the codebase MUST pass the automated pre-commit quality gate (`python scripts/pre_commit_gate.py`).
2. The runtime inspection must assert:
   - **Zero Console Errors**: `console.error` count $= 0$.
   - **Zero Unhandled Exceptions**: `window.onerror` and `unhandledrejection` count $= 0$.
   - **Zero WebGL Context / Shader Failures**: No failed shader compilation or program link errors.
3. Any unhandled runtime exception, missing asset 404, or unhandled rejection encountered during page load, engine switching, theme toggle, or animation loop constitutes a blocking failure (exit code 1).

## 15. Automated Canvas Drag & Interaction Validation
1. All 3D viewport canvas interactions (Three.js OrbitControls, pointer tracking, raycast picking) must be verified by automated synthetic pointer event tests (`scripts/check_canvas_drag.py`).
2. **Drag Disambiguation Invariant**: Synthetic drag motions ($\Delta r \ge 5\,\text{px}$) must update camera view matrices without triggering raycast click events or opening inspection dossiers.
3. **Click / Tap Invariant**: Discrete taps ($\Delta r < 5\,\text{px}, \Delta t < 350\,\text{ms}$) on cluster nodes must reliably trigger raycasting and modal inspection dossiers.
4. **Interaction Liveness**: Canvas drag sequences must maintain 60 FPS performance budgets ($< 16.6\,\text{ms/frame}$) and must never leak pointer capture or lock the UI thread upon `pointerup` or `pointercancel`.

## 16. In-App Visual Error Reporting & Resilient Canvas Fallback
1. All uncaught JavaScript runtime exceptions, unhandled promise rejections, and module initialization failures must be captured by a global error registry (`window.errorRegistry`).
2. If any error occurs during initialization or interaction, the system MUST:
   - Display a floating diagnostic alert pill (`#error-hud-pill`) with one-click access to the diagnostic dossier modal.
   - Force-dismiss any active splash overlay (`pointer-events: none; display: none`).
   - Re-enable OrbitControls camera dragging and maintain active WebGL animation frames.
3. The diagnostics modal must provide a copyable JSON payload containing error stack traces, active cosmological engine, GPU context, and theme state.
## 17. Color-Blind Accessibility, High-Legibility Baseline & UI Scaling
1. **Comfortable Baseline Sizing**: All interactive buttons must maintain a minimum text size of $\ge 11\,\text{px}$ bold with explicit $1.5\,\text{px}$ solid border outlines and distinct icon badges. Never rely on subtle background hue shifts alone for state distinction.
2. **Dual-Coding Invariant**: Every actionable control must feature dual visual coding (distinct emoji/icon + high-contrast text label).
3. **Dynamic UI Size Scaling**: The application must support variable UI size scaling (`--ui-scale`: 100%, 125%, 150%). All font sizes, paddings, and header heights must scale via `calc(... * var(--ui-scale))` while automatically triggering `resize` events to preserve dynamic viewport clearance margins (Rule 9).
4. **Color-Blind Safe Palettes**: Provide accessible theme accent presets (Cyan/Blue, Solar Gold [Protan/Deutan Safe], Nebula Violet [Tritan Safe], Emerald, Monochrome).

## 18. Top Navigation Header Cluster Isolation & Scalable Drawer Grids
1. **Isolated Flex Containers**: Top navigation headers must segregate action button clusters and right-hand telemetry/model banners into dedicated flex items. Right-hand model banners must enforce `max-width` and `text-overflow: ellipsis` with `white-space: nowrap` to prevent text cut-off on narrower displays.
2. **Auto-Fit Drawer Grids**: Multi-tab scientific drawer headers must use responsive auto-fit CSS grids (`grid-template-columns: repeat(4, 1fr)` or `repeat(auto-fit, minmax(105px, 1fr))`) rather than rigid $N \times 2$ matrices, guaranteeing that adding new settings or diagnostics tabs never causes visual overflow or clipping.

## 19. Canonical Global State Synchronization Across Multi-Location Settings
1. When user preferences (e.g. Tooltip Mode, UI Scale, GPU HUD Visibility, Theme) are exposed in multiple UI locations (topbar buttons, drawer panes, modal dialogues):
   - State mutations MUST pass through one canonical setter function (e.g. `setTooltipMode()`, `setUiScale()`, `applyCustomAccent()`).
   - All canonical setters MUST be explicitly exported on `window` and sync all associated DOM `<select>` values, `<input>` ranges, text badges, and `localStorage` keys simultaneously.

## 20. Cosmological Epistemic Uncertainty Semantics & Action-Oriented Naming
1. **Action-Oriented Control Labels**: UI controls in scientific workbenches must prioritize intuitive, action-oriented plain-English names (`Choose Model`, `Show Labels`, `Filter Galaxies`, `Superclusters`, `Distances`, `Gravity Sim`, `Flow Rivers`, `Cosmic Web`, `App Settings`) with clear subtitle summaries under group headings.
2. **Epistemic Uncertainty Semantics**: In Bayesian cosmological field reconstruction, "Epistemic Code" strictly represents observational measurement uncertainty (Tully-Fisher scatter, Zone of Avoidance dust obscuration, Wiener Filter posterior variances $\sigma_{\delta}$, selection biases), distinct from physical cosmic variance. Every reconstruction product must display its epistemic status code (`● Observed`, `≈ Reconstructed`, `◌ Posterior Sample`, `◇ Simulated`, `★ Sandbox`).
