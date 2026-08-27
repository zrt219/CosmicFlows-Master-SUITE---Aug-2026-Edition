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
$$
52.0
$$
before being interpreted as physical peculiar velocities. Apply this scaling:
* only to validated CF4 velocity values and errors;
* exactly once.
Attempting to apply it twice must fail loudly.

## 4. Linear-Theory Continuity Diagnostic
Where linear perturbation theory is appropriate, evaluate consistency using:
$$
\nabla \cdot \mathbf{v} \approx -H_0 f \delta
$$
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
## 21. GitHub CommonMark Underscore, Asterisk & HTML Tag Collision Invariants
1. **Subscript Underscore Isolation**: Any mathematical equation containing multiple subscripts or integrals with underscores (e.g., $I_{jl}$, $T_{kl}$, $\mathbf{v}_{\text{pot}}$, $\mathbf{v}_{\text{sol}}$, $\bar{q}_j$, $\bar{q}_l$) MUST be rendered as an isolated display equation block `$$ ... $$` with empty lines before and after, or written with isolated single-symbol inline references so that GitHub CommonMark never converts underscores into `<em>` tags.
2. **Asterisk Superscript Disambiguation (`^{\ast}`)**: Never use raw `^*` in LaTeX math formulas; always write `^{\ast}` (e.g. `\Sigma^{\ast}`, `\lambda^{\ast}`, `v^{\ast}`, `p^{\ast}`). Literal asterisks inside inline or display math are parsed by CommonMark as markdown italics and bullet triggers, which corrupts formulas on GitHub.
3. **Subscript Superscript Double Bracing (`_{{...}}`)**: When a subscript contains a superscript (e.g. $L^2$), always wrap the outer subscript in double braces `_{{L^2}}` or write `_{L^{2}}` to prevent MathJax strict double-subscript parsing exceptions.
4. **HTML Angle Bracket Escaping**: Never write raw angle brackets directly adjacent to letters in prose or inline math (e.g., `<r`, `<d`); always insert explicit whitespace (`< r`) or use LaTeX relational symbols (`\le`, `\ge`, `\prec`).
5. **Strict Display Block Isolation**: The closing display math delimiter (`$$ ... $$`) must always reside on its own dedicated line without any trailing prose on the same line, with empty blank lines preceding and following the block.
6. **Exhaustive Dual-Mode KaTeX Invariant Gate**: Every mathematical formula across all documentation (`README.md`, `AGENTS.md`, `GEMINI.md`) and in-app UI modals/tooltips (`#primer-modal`) MUST pass the automated in-browser KaTeX linter (`python scripts/check_katex_invariants.py`) asserting:
   - Zero display math errors (`displayMode: true`)
   - Zero inline math errors (`displayMode: false`)
   - Zero unparsed LaTeX macros visible as plain text.

## 22. Comoving vs. Physical Velocity Transformation Invariants
1. Peculiar velocity $\mathbf{v}_{\text{pec}}$ and proper recession velocity $\mathbf{v}_{\text{rec}}$ must strictly adhere to cosmological expansion kinematic transformations:

$$
\mathbf{v}_{\text{rec}}(\mathbf{r}, t) = H(t)\mathbf{r} + \mathbf{v}_{\text{pec}}(\mathbf{x}, t), \quad \mathbf{v}_{\text{pec}} = a(t)\dot{\mathbf{x}}
$$

2. Direct arithmetic combinations of comoving spatial coordinates $\mathbf{x}$ and proper physical coordinates $\mathbf{r} = a(t)\mathbf{x}$ without scale factor $a(t)$ conversion are strictly prohibited.

## 23. Cosmological Growth Factor & Logarithmic Rate Invariants
1. The linear growth factor $D(a)$ and growth rate $f(a)$ must be computed via exact numerical integration of the linear perturbation ODE:

$$
\frac{\mathrm{d}^2 D}{\mathrm{d}a^2} + \left( \frac{3}{a} + \frac{1}{H(a)}\frac{\mathrm{d}H}{\mathrm{d}a} \right) \frac{\mathrm{d}D}{\mathrm{d}a} = \frac{3 \Omega_{m,0}}{2 a^5 \left[H(a)/H_0\right]^2} D(a)
$$

2. In $\Lambda\mathrm{CDM}$ cosmology, evaluate logarithmic growth rates using $f \equiv \frac{\mathrm{d}\ln D}{\mathrm{d}\ln a} \approx \Omega_m(a)^{0.55}$. Approximations must fail loudly outside $z \in [0, 10]$.

## 24. Lagrangian Zel'dovich Displacement Invariants
1. In Lagrangian perturbation theory, displacement vector fields $\mathbf{\Psi}(\mathbf{q})$ must satisfy:

$$
\mathbf{x}(\mathbf{q}, t) = \mathbf{q} - D(t)\nabla_{\mathbf{q}} \Phi_0(\mathbf{q}), \quad \mathbf{\Psi}(\mathbf{q}) \equiv -\nabla_{\mathbf{q}}\Phi_0(\mathbf{q})
$$

where $\nabla_{\mathbf{q}} \cdot \mathbf{\Psi}(\mathbf{q}) = -\delta_0(\mathbf{q})$.
2. Displacement vectors must maintain strict units of comoving $\text{Mpc}/h$.

## 25. Rate-of-Strain & Vorticity Tensor Orthogonality Invariant
1. The symmetric rate-of-strain $S_{ij} = \frac{1}{2}(\partial_j v_i + \partial_i v_j)$ and anti-symmetric vorticity tensor $\Omega_{ij} = \frac{1}{2}(\partial_j v_i - \partial_i v_j)$ must satisfy strict Frobenius inner-product orthogonality:

$$
\mathrm{Tr}(S\Omega) = S_{ij}\Omega_{ij} \equiv 0
$$

2. Any numerical algorithm computing the Okubo-Weiss invariant $Q \equiv \mathrm{Tr}(S^2) - \mathrm{Tr}(\Omega\Omega^{\top})$ must enforce this contraction identity to machine precision ($< 10^{-14}$).

## 26. Cosmic Web Deformation Eigenvalue Ordering Invariant
1. All V-web and gravitational deformation tensor analyses must enforce descending eigenvalue sorting:

$$
\lambda_1 \ge \lambda_2 \ge \lambda_3
$$

2. Threshold classification functions must never accept unsorted eigenvalues and must validate threshold parameter $\gamma_{\text{th}} > 0$.

## 27. Halo Virial Mass & Concentration Invariants
1. Dark matter halo virial quantities must be evaluated at spherical overdensity $\Delta = 200$ relative to critical density $\rho_{\text{crit}}(z) = \frac{3 H^2(z)}{8\pi G}$:

$$
M_{200} = \frac{4\pi}{3} 200 \rho_{\text{crit}}(z) R_{200}^3, \quad c_{200} \equiv \frac{R_{200}}{r_s}
$$

2. Direct addition of mass $M_{200}$ and radius $R_{200}$ is dimensionally prohibited; mass must be represented in solar masses $M_{\odot}$ or $10^{14}\,h^{-1}M_{\odot}$.

## 28. Wiener Filter Reconstruction Matrix Invariants
1. The minimum-variance linear Wiener Filter field reconstruction must satisfy:

$$
\mathbf{s}_{\text{WF}} = \mathbf{S} \mathbf{R}^{\dagger} \left( \mathbf{R} \mathbf{S} \mathbf{R}^{\dagger} + \mathbf{N} \right)^{-1} \mathbf{d}
$$

where $\mathbf{S} \equiv \langle \mathbf{s}\mathbf{s}^{\dagger} \rangle$ is the signal cosmological covariance matrix, $\mathbf{N} \equiv \langle \mathbf{n}\mathbf{n}^{\dagger} \rangle$ is the noise covariance matrix, and $\mathbf{R}$ is the survey selection operator.
2. Numerical matrix inversions must employ Cholesky or LDL decomposition with condition number checking.

## 29. Symplectic Integrator Phase-Space Volume Preservation
1. Hamiltonian Monte Carlo trajectory integrators must use symplectic Störmer-Verlet leapfrog time-stepping.
2. Trajectory integrators must preserve phase-space volume:

$$
\det J = \left| \det \frac{\partial (v(t+\epsilon), p(t+\epsilon))}{\partial (v(t), p(t))} \right| \equiv 1
$$

3. Energy drift $\Delta \mathcal{H} = |\mathcal{H}(t) - \mathcal{H}(0)|$ exceeding $10^{-2}$ per step must trigger adaptive step-size reduction.

## 30. MCMC Convergence Diagnostic Threshold Invariants
1. Posterior chain convergence must be validated against strict statistical thresholds:
   - Rank-split potential scale reduction factor: $\hat{R} \le 1.01$
   - Multivariate potential scale reduction factor: $\text{MPSRF} \le 1.05$
   - Minimum effective sample size: $\text{ESS}_{\text{bulk}} \ge 400$, $\text{ESS}_{\text{tail}} \ge 200$
   - Geweke stationarity $Z$-score: $|Z| \le 1.96$ ($p > 0.05$)

## 31. Optimal Covariance Shrinkage Intensity Invariants
1. Ledoit-Wolf and Oracle Approximating Shrinkage (OAS) estimators must bound the shrinkage intensity $\hat{\lambda}^{\ast} \in [0, 1]$ and preserve positive-definiteness:

$$
\Sigma^{\ast} = (1 - \lambda^{\ast}) S + \lambda^{\ast} \mu I_p, \quad \lambda_{\min}(\Sigma^{\ast}) > 0
$$

## 32. Malmquist Bias Sign & Correction Formula Invariants
1. Homogeneous Malmquist distance modulus corrections must strictly apply a negative shift to observed distance moduli:

$$
\Delta\mu_{\text{hom}} = -\frac{3\ln 10}{5}\sigma_{\mu}^2 \approx -1.38155\,\sigma_{\mu}^2, \quad d_{\text{corr}} = d_{\text{obs}} \cdot 10^{\frac{\Delta\mu_{\text{hom}}}{5}}
$$

2. Positive homogeneous Malmquist corrections are unphysical and must fail loudly.

## 33. Minimum-Variance Bulk Flow Weighting Invariants
1. Galaxy peculiar velocity weights $w_n$ must incorporate both observational error $\sigma_n$ and 1D thermal cosmic velocity dispersion $\sigma_v \approx 187\,\text{km}/\text{s}$:

$$
w_n = \frac{1}{\sigma_n^2 + \sigma_v^2}
$$

2. Weights must never be zero or negative ($\sigma_n^2 + \sigma_v^2 > 0$).

## 34. Homology Boundary Operator Nilpotency Invariant
1. All cubical and simplicial persistent homology engines must satisfy the fundamental chain complex identity:

$$
\partial_k \circ \partial_{k+1} \equiv 0 \pmod 2 \quad (\partial^2 = 0)
$$

2. Every cell filtration $K(\delta_{\text{th}})$ must verify nilpotency before computing Betti numbers $\beta_k = \dim(\ker \partial_k) - \dim(\mathrm{im}\,\partial_{k+1})$.

## 35. Morse Critical Point Index Invariant
1. 3D stationary point classification on smooth velocity fields must map the Hessian/Jacobian eigenvalue signature directly to Morse index $\mu \in \{0, 1, 2, 3\}$:
   - $\mu = 0$: Local minimum / Repeller Void ($\lambda_1, \lambda_2, \lambda_3 > 0$)
   - $\mu = 1$: 1-Saddle / Cosmic Wall Hub ($\lambda_1 < 0 < \lambda_2 \le \lambda_3$)
   - $\mu = 2$: 2-Saddle / Cosmic Filament Hub ($\lambda_1 \le \lambda_2 < 0 < \lambda_3$)
   - $\mu = 3$: Local maximum / Attractor Node ($\lambda_1, \lambda_2, \lambda_3 < 0$)

## 36. Simplicial Surface Gauss-Bonnet Invariant
1. All 3D watershed basin boundaries extracted via Marching Tetrahedra must form closed, orientable 2-manifolds satisfying the discrete Gauss-Bonnet angular defect theorem:

$$
\sum_{v \in V} K_v = \sum_{v \in V} \left( 2\pi - \sum_{f \in F(v)} \theta_f(v) \right) = 2\pi \chi(\mathcal{M}) = 4\pi
$$

2. Non-zero boundary gaps or self-intersecting non-manifold facets constitute critical geometry generation failures.

## 37. Tomita-Gott Gaussian Random Field Morphometry Invariant
1. The theoretical Euler characteristic density $V_3(\nu) = N (\nu^2 - 1) e^{-\nu^2/2}$ must evaluate to zero at $\nu = \pm 1$ and exhibit the universal Gaussian peak-to-trough asymmetry constant:

$$
A_{\text{GRF}} \equiv \frac{|V_3(0)| - V_3(\sqrt{3})}{|V_3(0)| + V_3(\sqrt{3})} = \frac{1 - 2 e^{-3/2}}{1 + 2 e^{-3/2}} \approx 0.38288
$$

## 38. Inter-Basin Hydrodynamic Flux Facet Orientation Invariant
1. Inter-basin mass flux integrals across triangulated separatrix surfaces must compute outward unit normals via standard counter-clockwise vertex ordering:

$$
\hat{\mathbf{n}}_k = \frac{(\mathbf{x}_{k,1} - \mathbf{x}_{k,0}) \times (\mathbf{x}_{k,2} - \mathbf{x}_{k,0})}{\|(\mathbf{x}_{k,1} - \mathbf{x}_{k,0}) \times (\mathbf{x}_{k,2} - \mathbf{x}_{k,0})\|}, \quad \Phi_{AB} = \sum_{k=1}^{N_{\text{tri}}} \rho_k (\mathbf{v}_k \cdot \hat{\mathbf{n}}_k) A_k
$$

2. Inconsistent normal flips across adjacent facets are strictly forbidden.

## 39. Explicit WebGL Memory Disposal & Zero VRAM Leak Invariant
1. When destroying scenes, swapping reconstruction datasets, or rebuilding streamlines:
   - Geometries: `geometry.dispose()`
   - Materials: `material.dispose()`
   - Textures: `texture.dispose()`
   - Render Targets: `renderTarget.dispose()`
2. Continuous dataset switching must demonstrate zero GPU memory growth.

## 40. Shader Floating-Point Defensive Clamping Invariant
1. Custom GLSL shaders must specify `precision highp float;` and guard all division, square roots, and logarithms against `NaN` or `Infinity`:
   - `float safeDiv(float n, float d) { return n / max(abs(d), 1e-7); }`
   - `float safeLog(float x) { return log(max(x, 1e-7)); }`

## 41. UnrealBloomPass & Tone-Mapping HDR Clamping
1. Luminance thresholding in post-processing bloom passes must use half-float render targets (`THREE.HalfFloatType`) and enforce Reinhard or ACESFilmic tone-mapping to prevent white clipping artifacts.

## 42. Instanced Double-Buffering for 60 FPS Particle Advection
1. Particle advection across thousands of streamlines must use instanced buffer attributes with double-buffered vertex buffers to eliminate garbage collection pauses and avoid WebGL buffer reallocation during continuous time scrubbing.

## 43. Raycast Click vs. Drag Disambiguation Invariant
1. Canvas pointer events must distinguish between discrete selection taps ($\Delta r < 5\,\text{px}, \Delta t < 350\,\text{ms}$) and continuous OrbitControls camera dragging ($\Delta r \ge 5\,\text{px}$).
2. Drag gestures must never trigger inspection modals or cluster picking.

## 44. WebGL Context Loss Recovery & Fallback Invariant
1. All 3D WebGL renderers must listen for `webglcontextlost` and `webglcontextrestored` events on canvas elements, cleanly pausing the animation loop on loss and re-initializing shaders/buffers upon context restoration.

## 45. Pure CSS Variable Theming & Dynamic Rem/Calc Scaling
1. All layout spacing, font sizes, border radii, and accent colors must use CSS variables (`--ui-scale`, `--accent-primary`, `--bg-card`) with dynamic `calc()` expressions.
2. Hardcoded pixel font sizes in UI components are forbidden.

## 46. Focus-Visible Keyboard Navigation Invariant
1. Every interactive button, dropdown, slider, and checkbox must display a high-contrast focus ring (`outline: 2px solid var(--accent-primary)`) when navigated via Tab / arrow keys.
2. All modal dialogues must dismiss cleanly on Escape key press.

## 47. Live Telemetry & Epistemic Uncertainty Badge Synchronization
1. Physical and epistemological status badges must synchronize in real time across the primary topbar, floating HUD, and drawer inspection panels when astronomical datasets or active engines are toggled.

## 48. Subsystem Global Namespace Encapsulation
1. All JavaScript classes, mathematical solvers, data loaders, and controllers must reside inside the single top-level namespace `window.cosmicflows` or export clean ES6 modules.
2. Unscoped global variables on `window` are strictly prohibited.

## 49. Mandatory Dual-Mode KaTeX In-Browser CDP Verification Gate
1. Every pull request or commit must pass `scripts/check_katex_invariants.py`, asserting 0 parse errors across both `displayMode: true` and `displayMode: false` across all documentation and in-app templates.

## 50. Astronomical Data File Format & Byte-Order Integrity
1. FITS and CSV data ingestion routines must parse big-endian binary float arrays correctly and assert data shape and non-NaN values upon loading.

## 51. Zero Console Error & Zero Unhandled Exception Production Gate
1. The master pre-commit gate (`python scripts/pre_commit_gate.py`) must pass 100% across all automated suites (AST TDZ, headless CDP boot exceptions, canvas drag, image variance, splash lifecycle, KaTeX linter) with zero warnings or errors before publishing builds.
