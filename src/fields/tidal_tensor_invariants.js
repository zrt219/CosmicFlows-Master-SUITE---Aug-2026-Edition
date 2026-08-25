/**
 * @file tidal_tensor_invariants.js
 * @description Comprehensive, production-grade implementation of 3D Gravitational Tidal Tensor Invariants,
 * Zel'dovich Web Morphological Classification, Cosmological Anisotropy Metrics, Singularity Collapse Time
 * Estimators, Tidal Torque Theory (TTT) Angular Momentum Acquisition, and High-Order Spectral/Grid Solvers
 * for the CosmicFlows-4 Research Workbench.
 *
 * Mathematical and Astrophysical Foundations:
 * -------------------------------------------
 * 1. Gravitational Potential and Poisson Equation:
 *    \nabla^2 \Phi(\mathbf{x}) = 4 \pi G \bar{\rho} a^2 \delta(\mathbf{x}) = S(\mathbf{x})
 *    where \delta(\mathbf{x}) = (\rho(\mathbf{x}) - \bar{\rho}) / \bar{\rho} is the density contrast,
 *    a is the cosmological scale factor, and \bar{\rho} is the mean matter density.
 *
 * 2. Gravitational Deformation Tensor (Tidal Hessian):
 *    \Psi_{ij}(\mathbf{x}) = \frac{\partial^2 \Phi}{\partial x_i \partial x_j}(\mathbf{x})  (i, j \in \{x, y, z\})
 *    Normalized dimensionless deformation tensor:
 *    \tilde{\Psi}_{ij} = \frac{\Psi_{ij}}{4 \pi G \bar{\rho} a^2}
 *    \text{Tr}(\tilde{\Psi}) = \tilde{\Psi}_{xx} + \tilde{\Psi}_{yy} + \tilde{\Psi}_{zz} = \delta(\mathbf{x})
 *
 * 3. Traceless Gravitational Tidal Tensor:
 *    T_ij(\mathbf{x}) = \Psi_{ij}(\mathbf{x}) - \frac{1}{3} \nabla^2 \Phi(\mathbf{x}) \delta_{ij}
 *                     = \Psi_{ij}(\mathbf{x}) - \frac{1}{3} \text{Tr}(\Psi) \delta_{ij}
 *    By construction: \text{Tr}(T) = T_{xx} + T_{yy} + T_{zz} = 0.
 *    Normalized traceless tidal shear tensor:
 *    t_ij(\mathbf{x}) = \tilde{\Psi}_{ij}(\mathbf{x}) - \frac{1}{3} \delta(\mathbf{x}) \delta_{ij}
 *
 * 4. Eigenvalue Decomposition and Principal Invariants:
 *    Sorted real eigenvalues: \lambda_1 \ge \lambda_2 \ge \lambda_3.
 *    Orthonormal eigenvectors: \mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3 forming SO(3) principal frame.
 *    - First Principal Invariant:
 *      I_1 = \text{Tr}(\Psi) = \lambda_1 + \lambda_2 + \lambda_3 = \delta
 *      (For traceless tidal tensor T: I_1(T) = \text{Tr}(T) = 0)
 *    - Second Principal Invariant:
 *      I_2 = \frac{1}{2} [(\text{Tr} \Psi)^2 - \text{Tr}(\Psi^2)] = \lambda_1 \lambda_2 + \lambda_2 \lambda_3 + \lambda_3 \lambda_1
 *      For traceless T: J_2 = -\frac{1}{2} \text{Tr}(T^2) = \lambda_1 \lambda_2 + \lambda_2 \lambda_3 + \lambda_3 \lambda_1 \le 0
 *      \text{Tr}(T^2) = \frac{1}{3} [(\lambda_1 - \lambda_2)^2 + (\lambda_2 - \lambda_3)^2 + (\lambda_3 - \lambda_1)^2]
 *    - Third Principal Invariant:
 *      I_3 = \det(\Psi) = \lambda_1 \lambda_2 \lambda_3
 *      For traceless T: J_3 = \det(T) = \lambda_1 \lambda_2 \lambda_3
 *    - Characteristic Polynomial & Cayley-Hamilton Identity:
 *      \det(\lambda I - \Psi) = \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3 = 0
 *      \Psi^3 - I_1 \Psi^2 + I_2 \Psi - I_3 I_{3\times 3} = 0_{3\times 3}
 *
 * 5. Cosmic Web Classification (Hahn et al. 2007; Forero-Romero et al. 2009):
 *    Classifies morphology based on the number of deformation eigenvalues exceeding threshold \gamma_{\text{th}}:
 *    - Type 0: Void      (\lambda_1 \le \gamma_{\text{th}})               -> 0 axes collapsing, 3 expanding
 *    - Type 1: Sheet     (\lambda_1 > \gamma_{\text{th}} \ge \lambda_2)   -> 1 axis collapsing (\mathbf{e}_1 normal), 2 expanding
 *    - Type 2: Filament  (\lambda_2 > \gamma_{\text{th}} \ge \lambda_3)   -> 2 axes collapsing, 1 expanding (\mathbf{e}_3 spine)
 *    - Type 3: Knot/Node (\lambda_3 > \gamma_{\text{th}})               -> 3 axes collapsing (gravitational cluster peak)
 *
 * 6. Cosmological Anisotropy & Morphology Metrics:
 *    - Anisotropy Parameter: \alpha = \frac{\lambda_1 - \lambda_3}{2 \sum_{i} \lambda_i} = \frac{\lambda_1 - \lambda_3}{2 \text{Tr}(\Psi)}
 *    - Ellipticity (Doroshkevich 1970; BBKS 1986): e = \frac{\lambda_1 - \lambda_3}{2 (\lambda_1 + \lambda_2 + \lambda_3)}
 *    - Prolateness: p = \frac{\lambda_1 + \lambda_3 - 2 \lambda_2}{2 (\lambda_1 + \lambda_2 + \lambda_3)} (-e \le p \le e)
 *      * p > 0: Prolate / filamentary deformation (\lambda_1 - \lambda_2 > \lambda_2 - \lambda_3)
 *      * p < 0: Oblate / pancake deformation (\lambda_1 - \lambda_2 < \lambda_2 - \lambda_3)
 *      * p = 0: Triaxial symmetric (e \ne 0) or spherical (e = 0)
 *    - Tidal Shear Magnitude: q = \sqrt{\frac{1}{2} [(\lambda_1-\lambda_2)^2 + (\lambda_2-\lambda_3)^2 + (\lambda_3-\lambda_1)^2]} = \sqrt{\frac{3}{2} \text{Tr}(T^2)}
 *    - Triaxiality Parameter: T_{\text{triax}} = \frac{\lambda_1^2 - \lambda_2^2}{\lambda_1^2 - \lambda_3^2} \in [0, 1]
 *    - Sphericity Index: S = 1 - 3e
 *
 * 7. Zel'dovich Approximation & Singularity Collapse Estimators:
 *    - Eulerian mapping: \mathbf{x}(\mathbf{q}, t) = \mathbf{q} - D(t) \nabla_\mathbf{q} \Phi_0(\mathbf{q})
 *    - Deformation Jacobian: J(\mathbf{q}, t) = \det(\delta_{ij} - D(t) \Psi_{ij}) = (1 - D(t)\lambda_1)(1 - D(t)\lambda_2)(1 - D(t)\lambda_3)
 *    - Density: \rho(\mathbf{x}, t) = \bar{\rho} / J(\mathbf{q}, t)
 *    - Collapse condition along i-th principal axis: D(t_{\text{col}, i}) = 1 / \lambda_i (for \lambda_i > 0)
 *      * 1st collapse (Pancake caustic / sheet formation): D(t_{\text{col}, 1}) = 1 / \lambda_1
 *      * 2nd collapse (Filament formation): D(t_{\text{col}, 2}) = 1 / \lambda_2
 *      * 3rd collapse (Cluster / Halo formation): D(t_{\text{col}, 3}) = 1 / \lambda_3
 *    - EdS model: D(z) = 1/(1+z), t_{\text{col}, i} = t_0 / (1 + z_{\text{col}, i})^{3/2} = t_0 / \lambda_i^{3/2}
 *    - Flat \Lambda\text{CDM}: Exact numerical growth factor D(z) and cosmic time t(z).
 *
 * 8. Tidal Torque Theory (TTT) Angular Momentum Acquisition:
 *    - Angular momentum acquired by proto-halo of Lagrangian volume V_L:
 *      L_i(t) = a^2(t) \dot{D}(t) \sum_{j,k,l} \epsilon_{ijk} T_{jl}(\mathbf{q}_{\text{cm}}) I_{lk}
 *    - Torque Vector: \tau_i = \sum_{j,k,l} \epsilon_{ijk} T_{jl} I_{lk} = \frac{1}{2} \sum_{j,k} \epsilon_{ijk} [T, I]_{jk}
 *    - Proto-halo Inertia Tensor: I_{lk} = \int_{V_L} \rho(\mathbf{q}) (q_l - \bar{q}_l)(q_k - \bar{q}_k) d^3q
 *    - Non-zero spin requires misalignment / non-commutation: [T, I] = T I - I T \ne 0.
 *    - Dimensionless Peebles Spin Parameter: \lambda_{\text{spin}} = \frac{L |E|^{1/2}}{G M^{5/2}}
 *    - Bullock Modified Spin Parameter: \lambda' = \frac{L}{\sqrt{2} M V_{\text{vir}} R_{\text{vir}}}
 *
 * References:
 * - Zel'dovich, Ya. B. (1970). "Gravitational instability: An approximate theory for large density perturbations." A&A, 5, 84.
 * - Doroshkevich, A. G. (1970). "Spatial structure of perturbations and the origin of rotation of galaxies." Astrofizika, 6, 581.
 * - Bardeen, J. M., Bond, J. R., Kaiser, N., & Szalay, A. S. (1986). "The statistics of peaks of Gaussian random fields." ApJ, 304, 15 (BBKS).
 * - White, S. D. M. (1984). "Angular momentum growth in protogalaxies." ApJ, 286, 38.
 * - Catelan, P., & Theuns, T. (1996). "Evolution of the angular momentum of protogalaxies from tidal torques." MNRAS, 282, 436.
 * - Porciani, C., Dekel, A., & Hoffman, Y. (2002). "Testing tidal torque theory - I. Spin amplitude and direction." MNRAS, 332, 325.
 * - Hahn, O., Porciani, C., Carollo, C. M., & Dekel, A. (2007). "Properties of cosmic web structures." MNRAS, 375, 489.
 * - Forero-Romero, J. E., et al. (2009). "A dynamical classification of the cosmic web." MNRAS, 396, 1815.
 * - Hoffman, Y., Metuki, O., Yepes, G., et al. (2012). "A kinematic classification of the cosmic web." MNRAS, 425, 2049.
 * - Libeskind, N. I., et al. (2018). "The velocity shear tensor: tracer of filamentary cosmic web." MNRAS, 473, 1195.
 * - Courtois, H. M., et al. (2023). "Cosmicflows-4: The Velocity Field and Cosmography." ApJ, 944, 94.
 *
 * @module fields/tidal_tensor_invariants
 */

import { GridIndexer, BoundaryMode } from './grid_indexer.js';
import { ScalarField3D } from './scalar_field_3d.js';
import { VectorField3D } from './vector_field_3d.js';

// ============================================================================
// ENUMERATIONS AND CONSTANTS
// ============================================================================

/**
 * Standard Cosmic Web morphological classification types.
 * Defined by Hahn et al. (2007) and Forero-Romero et al. (2009).
 * @readonly
 * @enum {number}
 */
export const CosmicWebType = Object.freeze({
  VOID: 0,      // 0 eigenvalues > gamma_th: 3D volumetric expansion
  SHEET: 1,     // 1 eigenvalue  > gamma_th: 1D collapse (wall / pancake)
  FILAMENT: 2,  // 2 eigenvalues > gamma_th: 2D collapse (filament spine)
  KNOT: 3       // 3 eigenvalues > gamma_th: 3D collapse (cluster / halo peak)
});

/**
 * Descriptive human-readable labels for Cosmic Web types.
 * @readonly
 */
export const CosmicWebLabels = Object.freeze({
  [CosmicWebType.VOID]: 'Void (3D Expansion)',
  [CosmicWebType.SHEET]: 'Sheet / Wall (1D Collapse)',
  [CosmicWebType.FILAMENT]: 'Filament (2D Collapse)',
  [CosmicWebType.KNOT]: 'Knot / Cluster (3D Collapse)'
});

/**
 * Standard Cosmological World Models.
 * @readonly
 * @enum {string}
 */
export const CosmologicalModel = Object.freeze({
  EDS: 'eds',                 // Einstein-de Sitter: Omega_m = 1.0, Omega_Lambda = 0.0
  FLAT_LCDM: 'flat_lcdm',     // Flat LambdaCDM: Omega_m + Omega_Lambda = 1.0
  OPEN_CDM: 'open_cdm'        // Open CDM: Omega_m < 1.0, Omega_Lambda = 0.0
});

/**
 * Spatial smoothing kernels for multiscale tidal analysis.
 * @readonly
 * @enum {string}
 */
export const TidalSmoothingKernel = Object.freeze({
  GAUSSIAN: 'gaussian',   // W(k, R) = exp(-k^2 R^2 / 2)
  TOP_HAT: 'top_hat',     // W(k, R) = 3 * (sin(kR) - kR cos(kR)) / (kR)^3
  SHARP_K: 'sharp_k',     // W(k, R) = Theta(1 - kR)
  NONE: 'none'
});

/**
 * Finite-difference differentiation stencil order.
 * @readonly
 * @enum {number}
 */
export const DifferentiationOrder = Object.freeze({
  SECOND: 2, // 3-point stencil (O(h^2))
  FOURTH: 4, // 5-point stencil (O(h^4))
  SIXTH: 6   // 7-point stencil (O(h^6))
});

/**
 * Standard physical and cosmological constants in astrophysical units.
 * Units: Distance in Mpc/h, Velocity in km/s, Mass in M_sun/h, Time in Gyr/h.
 * @readonly
 */
export const CosmologicalConstants = Object.freeze({
  /** Gravitational constant G in (km/s)^2 * (Mpc/h) / (10^10 M_sun/h) */
  G_ASTRO: 4.30091e-4,
  /** Default Hubble constant H_0 = 100 h km/s/Mpc (so h = 1.0 in h-units) */
  H0_DEFAULT: 100.0,
  /** Default matter density parameter Omega_m (Planck 2018 / CF4 baseline) */
  OMEGA_M_DEFAULT: 0.315,
  /** Default dark energy density parameter Omega_Lambda */
  OMEGA_LAMBDA_DEFAULT: 0.685,
  /** Default radiation density parameter Omega_r */
  OMEGA_R_DEFAULT: 0.0,
  /** Critical density rho_crit in h^2 M_sun / Mpc^3 */
  RHO_CRIT_0: 2.77536627e11,
  /** Speed of light c in km/s */
  C_KMS: 299792.458,
  /** Default 4 * pi * G * rho_bar normalizer for Poisson equation */
  FOUR_PI_G_RHO_BAR: 1.0
});

// ============================================================================
// DATA STRUCTURES & VALUE OBJECTS
// ============================================================================

/**
 * Encapsulates the complete set of algebraic, morphological, and dynamical invariants
 * of a 3x3 Gravitational Deformation / Tidal Tensor at a point.
 */
export class TidalTensorInvariants {
  /**
   * Constructs a TidalTensorInvariants instance.
   * 
   * @param {Object} params Configuration parameters.
   * @param {Float64Array} params.tensor 3x3 flat array of the tensor (row-major).
   * @param {Float64Array} params.tracelessTensor 3x3 flat array of traceless tidal shear T_ij.
   * @param {number} params.trace First invariant I_1 = Tr(Psi) = lambda_1 + lambda_2 + lambda_3.
   * @param {number} params.secondInvariant Second invariant I_2 = lambda_1*lambda_2 + lambda_2*lambda_3 + lambda_3*lambda_1.
   * @param {number} params.thirdInvariant Third invariant I_3 = det(Psi) = lambda_1*lambda_2*lambda_3.
   * @param {number} params.shearSecondInvariant Traceless shear invariant J_2 = -0.5 * Tr(T^2).
   * @param {number} params.shearThirdInvariant Traceless shear invariant J_3 = det(T).
   * @param {number} params.discriminant Polynomial discriminant Delta = 18 I_1 I_2 I_3 - 4 I_1^3 I_3 + I_1^2 I_2^2 - 4 I_2^3 - 27 I_3^2.
   * @param {Float64Array} params.eigenvalues Sorted eigenvalues [lambda_1, lambda_2, lambda_3] (lambda_1 >= lambda_2 >= lambda_3).
   * @param {Array<Float64Array>} params.eigenvectors Orthonormal eigenvectors [e_1, e_2, e_3] corresponding to eigenvalues.
   * @param {number} params.anisotropy Anisotropy parameter alpha = (lambda_1 - lambda_3) / (2 * sum(lambda_i)).
   * @param {number} params.ellipticity Ellipticity e = (lambda_1 - lambda_3) / (2 * sum(lambda_i)).
   * @param {number} params.prolateness Prolateness p = (lambda_1 + lambda_3 - 2*lambda_2) / (2 * sum(lambda_i)).
   * @param {number} params.shearMagnitude Tidal shear magnitude q = sqrt(3 * J_2_pos) = sqrt(0.5 * sum((li - lj)^2)).
   * @param {number} params.triaxiality Triaxiality parameter T_triax = (lambda_1^2 - lambda_2^2) / (lambda_1^2 - lambda_3^2).
   * @param {number} params.sphericity Sphericity index S = 1 - 3*e.
   * @param {number} params.webType CosmicWebType enum value.
   * @param {number} params.gammaThreshold Classification threshold used.
   * @param {Float64Array} [params.collapseRedshifts] Collapse redshifts [z_col_1, z_col_2, z_col_3].
   * @param {Float64Array} [params.collapseTimesGyr] Collapse cosmic times in Gyr [t_col_1, t_col_2, t_col_3].
   */
  constructor(params) {
    this.tensor = new Float64Array(params.tensor);
    this.tracelessTensor = new Float64Array(params.tracelessTensor);
    this.trace = Number(params.trace);
    this.secondInvariant = Number(params.secondInvariant);
    this.thirdInvariant = Number(params.thirdInvariant);
    this.shearSecondInvariant = Number(params.shearSecondInvariant);
    this.shearThirdInvariant = Number(params.shearThirdInvariant);
    this.discriminant = Number(params.discriminant);
    this.eigenvalues = new Float64Array(params.eigenvalues);
    this.eigenvectors = [
      new Float64Array(params.eigenvectors[0]),
      new Float64Array(params.eigenvectors[1]),
      new Float64Array(params.eigenvectors[2])
    ];
    this.anisotropy = Number(params.anisotropy);
    this.ellipticity = Number(params.ellipticity);
    this.prolateness = Number(params.prolateness);
    this.shearMagnitude = Number(params.shearMagnitude);
    this.triaxiality = Number(params.triaxiality);
    this.sphericity = Number(params.sphericity);
    this.webType = Number(params.webType);
    this.gammaThreshold = Number(params.gammaThreshold);

    this.collapseRedshifts = params.collapseRedshifts
      ? new Float64Array(params.collapseRedshifts)
      : new Float64Array([NaN, NaN, NaN]);
    this.collapseTimesGyr = params.collapseTimesGyr
      ? new Float64Array(params.collapseTimesGyr)
      : new Float64Array([NaN, NaN, NaN]);
  }

  /**
   * First sorted eigenvalue (fastest collapsing direction).
   * @returns {number}
   */
  get lambda1() {
    return this.eigenvalues[0];
  }

  /**
   * Second sorted eigenvalue (intermediate collapsing direction).
   * @returns {number}
   */
  get lambda2() {
    return this.eigenvalues[1];
  }

  /**
   * Third sorted eigenvalue (slowest collapsing / fastest expanding direction).
   * @returns {number}
   */
  get lambda3() {
    return this.eigenvalues[2];
  }

  /**
   * Normal vector to sheet / pancake (principal axis e_1).
   * @returns {Float64Array}
   */
  get sheetNormal() {
    return this.eigenvectors[0];
  }

  /**
   * Tangent spine vector to filament (principal axis e_3).
   * @returns {Float64Array}
   */
  get filamentSpine() {
    return this.eigenvectors[2];
  }

  /**
   * Human-readable classification string.
   * @returns {string}
   */
  get webTypeName() {
    return CosmicWebLabels[this.webType] || 'Unknown';
  }

  /**
   * Whether this point has experienced pancake caustic formation (1st collapse).
   * @param {number} [currentGrowth=1.0] Current linear growth factor D(t).
   * @returns {boolean}
   */
  hasPancakeCollapsed(currentGrowth = 1.0) {
    return this.lambda1 > 0.0 && (currentGrowth * this.lambda1 >= 1.0);
  }

  /**
   * Whether this point has experienced filamentary caustic collapse (2nd collapse).
   * @param {number} [currentGrowth=1.0] Current linear growth factor D(t).
   * @returns {boolean}
   */
  hasFilamentCollapsed(currentGrowth = 1.0) {
    return this.lambda2 > 0.0 && (currentGrowth * this.lambda2 >= 1.0);
  }

  /**
   * Whether this point has undergone complete 3D cluster virialization (3rd collapse).
   * @param {number} [currentGrowth=1.0] Current linear growth factor D(t).
   * @returns {boolean}
   */
  hasClusterCollapsed(currentGrowth = 1.0) {
    return this.lambda3 > 0.0 && (currentGrowth * this.lambda3 >= 1.0);
  }

  /**
   * Serializes the invariants to a plain JavaScript object.
   * @returns {Object}
   */
  toJSON() {
    return {
      tensor: Array.from(this.tensor),
      tracelessTensor: Array.from(this.tracelessTensor),
      trace: this.trace,
      secondInvariant: this.secondInvariant,
      thirdInvariant: this.thirdInvariant,
      shearSecondInvariant: this.shearSecondInvariant,
      shearThirdInvariant: this.shearThirdInvariant,
      discriminant: this.discriminant,
      eigenvalues: Array.from(this.eigenvalues),
      eigenvectors: this.eigenvectors.map(e => Array.from(e)),
      anisotropy: this.anisotropy,
      ellipticity: this.ellipticity,
      prolateness: this.prolateness,
      shearMagnitude: this.shearMagnitude,
      triaxiality: this.triaxiality,
      sphericity: this.sphericity,
      webType: this.webType,
      webTypeName: this.webTypeName,
      gammaThreshold: this.gammaThreshold,
      collapseRedshifts: Array.from(this.collapseRedshifts),
      collapseTimesGyr: Array.from(this.collapseTimesGyr)
    };
  }
}

/**
 * Encapsulates the kinematic and dynamical state in the Zel'dovich Approximation.
 */
export class ZeldovichState {
  /**
   * @param {Object} params
   * @param {Float64Array} params.lagrangianPos Initial Lagrangian coordinate q [x, y, z].
   * @param {Float64Array} params.eulerianPos Mapped Eulerian coordinate x(q, t) [x, y, z].
   * @param {Float64Array} params.displacement Displacement vector Psi(q) = -grad Phi_0(q).
   * @param {number} params.growthFactor Linear growth factor D(t).
   * @param {number} params.jacobianDet Jacobian determinant J(q, t) = det(dx / dq).
   * @param {number} params.density Contrast (rho - rho_bar) / rho_bar or physical rho.
   * @param {number} params.streamCount Number of intersecting velocity streams (1, 3, or more).
   * @param {boolean} params.isCaustic True if J <= 0 (post-shell-crossing caustic).
   */
  constructor(params) {
    this.lagrangianPos = new Float64Array(params.lagrangianPos);
    this.eulerianPos = new Float64Array(params.eulerianPos);
    this.displacement = new Float64Array(params.displacement);
    this.growthFactor = Number(params.growthFactor);
    this.jacobianDet = Number(params.jacobianDet);
    this.density = Number(params.density);
    this.streamCount = Number(params.streamCount || (this.jacobianDet <= 0 ? 3 : 1));
    this.isCaustic = Boolean(params.isCaustic || this.jacobianDet <= 0);
  }

  toJSON() {
    return {
      lagrangianPos: Array.from(this.lagrangianPos),
      eulerianPos: Array.from(this.eulerianPos),
      displacement: Array.from(this.displacement),
      growthFactor: this.growthFactor,
      jacobianDet: this.jacobianDet,
      density: this.density,
      streamCount: this.streamCount,
      isCaustic: this.isCaustic
    };
  }
}

/**
 * Encapsulates Tidal Torque Theory (TTT) Angular Momentum Acquisition and Misalignment Diagnostics.
 */
export class TidalTorqueState {
  /**
   * @param {Object} params
   * @param {Float64Array} params.angularMomentum Acquired angular momentum vector L = [Lx, Ly, Lz].
   * @param {Float64Array} params.torque Vector torque tau = [tau_x, tau_y, tau_z] = eps_ijk T_jl I_lk.
   * @param {Float64Array} params.tidalTensor 3x3 Tidal tensor T_jl.
   * @param {Float64Array} params.inertiaTensor 3x3 Proto-halo inertia tensor I_lk.
   * @param {Float64Array} params.commutator 3x3 Commutator matrix [T, I] = T*I - I*T.
   * @param {number} params.misalignmentAngleRad Angle between major principal axis of T and I (radians).
   * @param {number} params.peeblesSpin Dimensionless Peebles spin parameter lambda = L * |E|^0.5 / (G * M^2.5).
   * @param {number} params.bullockSpin Dimensionless Bullock spin parameter lambda' = L / (sqrt(2) * M * V_vir * R_vir).
   * @param {number} params.haloMass Total halo mass M (M_sun/h).
   */
  constructor(params) {
    this.angularMomentum = new Float64Array(params.angularMomentum);
    this.torque = new Float64Array(params.torque);
    this.tidalTensor = new Float64Array(params.tidalTensor);
    this.inertiaTensor = new Float64Array(params.inertiaTensor);
    this.commutator = new Float64Array(params.commutator);
    this.misalignmentAngleRad = Number(params.misalignmentAngleRad);
    this.peeblesSpin = Number(params.peeblesSpin);
    this.bullockSpin = Number(params.bullockSpin);
    this.haloMass = Number(params.haloMass);
  }

  /**
   * Angular momentum vector magnitude |L|.
   * @returns {number}
   */
  get angularMomentumMagnitude() {
    const lx = this.angularMomentum[0];
    const ly = this.angularMomentum[1];
    const lz = this.angularMomentum[2];
    return Math.sqrt(lx * lx + ly * ly + lz * lz);
  }

  /**
   * Torque vector magnitude |tau|.
   * @returns {number}
   */
  get torqueMagnitude() {
    const tx = this.torque[0];
    const ty = this.torque[1];
    const tz = this.torque[2];
    return Math.sqrt(tx * tx + ty * ty + tz * tz);
  }

  /**
   * Misalignment angle in degrees.
   * @returns {number}
   */
  get misalignmentAngleDeg() {
    return (this.misalignmentAngleRad * 180.0) / Math.PI;
  }

  toJSON() {
    return {
      angularMomentum: Array.from(this.angularMomentum),
      torque: Array.from(this.torque),
      angularMomentumMagnitude: this.angularMomentumMagnitude,
      torqueMagnitude: this.torqueMagnitude,
      misalignmentAngleRad: this.misalignmentAngleRad,
      misalignmentAngleDeg: this.misalignmentAngleDeg,
      peeblesSpin: this.peeblesSpin,
      bullockSpin: this.bullockSpin,
      haloMass: this.haloMass
    };
  }
}

/**
 * Encapsulates global cosmological web morphological statistics across a 3D grid.
 */
export class CosmicWebSummary {
  /**
   * @param {Object} params
   * @param {number} params.totalVoxels Total number of grid voxels analyzed.
   * @param {number} params.gammaThreshold Threshold gamma_th used for classification.
   * @param {Float64Array} params.volumeFractions Volume filling fractions [f_V_void, f_V_sheet, f_V_filament, f_V_knot].
   * @param {Float64Array} params.massFractions Mass filling fractions [f_M_void, f_M_sheet, f_M_filament, f_M_knot].
   * @param {Float64Array} params.voxelCounts Raw voxel counts per category.
   * @param {Object} params.meanInvariantsPerType Average invariants (trace, shear, ellipticity) per web type.
   * @param {Array<Object>} [params.thresholdScan] Multi-threshold scan results across gamma_th.
   */
  constructor(params) {
    this.totalVoxels = Number(params.totalVoxels);
    this.gammaThreshold = Number(params.gammaThreshold);
    this.volumeFractions = new Float64Array(params.volumeFractions);
    this.massFractions = new Float64Array(params.massFractions);
    this.voxelCounts = new Float64Array(params.voxelCounts);
    this.meanInvariantsPerType = params.meanInvariantsPerType || {};
    this.thresholdScan = params.thresholdScan || [];
  }

  toJSON() {
    return {
      totalVoxels: this.totalVoxels,
      gammaThreshold: this.gammaThreshold,
      volumeFractions: {
        void: this.volumeFractions[CosmicWebType.VOID],
        sheet: this.volumeFractions[CosmicWebType.SHEET],
        filament: this.volumeFractions[CosmicWebType.FILAMENT],
        knot: this.volumeFractions[CosmicWebType.KNOT]
      },
      massFractions: {
        void: this.massFractions[CosmicWebType.VOID],
        sheet: this.massFractions[CosmicWebType.SHEET],
        filament: this.massFractions[CosmicWebType.FILAMENT],
        knot: this.massFractions[CosmicWebType.KNOT]
      },
      voxelCounts: {
        void: this.voxelCounts[CosmicWebType.VOID],
        sheet: this.voxelCounts[CosmicWebType.SHEET],
        filament: this.voxelCounts[CosmicWebType.FILAMENT],
        knot: this.voxelCounts[CosmicWebType.KNOT]
      },
      meanInvariantsPerType: this.meanInvariantsPerType,
      thresholdScan: this.thresholdScan
    };
  }
}

// ============================================================================
// MATHEMATICAL & NUMERICAL CORE ENGINES
// ============================================================================

/**
 * Pure mathematical algorithms for 3x3 symmetric matrix eigensystems, Cardano/Vieta roots,
 * Principal Invariants, Zel'dovich dynamics, and Tidal Torque algebra.
 */
export class TidalTensorMath {
  /**
   * Decomposes a 3x3 deformation Hessian matrix Psi_ij into its trace (density contrast delta)
   * and traceless tidal shear tensor T_ij = Psi_ij - (1/3) * Tr(Psi) * delta_ij.
   * 
   * @param {ArrayLike<number>} Psi 3x3 flat array [xx, xy, xz, yx, yy, yz, zx, zy, zz].
   * @returns {{ traceless: Float64Array, trace: number, normalizedTrace: number }}
   */
  static decomposeTidalTensor(Psi) {
    if (!Psi || Psi.length < 9) {
      throw new Error('TidalTensorMath.decomposeTidalTensor: Input matrix Psi must have length >= 9.');
    }

    const xx = Psi[0], xy = Psi[1], xz = Psi[2];
    const yx = Psi[3], yy = Psi[4], yz = Psi[5];
    const zx = Psi[6], zy = Psi[7], zz = Psi[8];

    const trace = xx + yy + zz;
    const thirdTrace = trace / 3.0;

    const traceless = new Float64Array(9);
    traceless[0] = xx - thirdTrace;
    traceless[1] = 0.5 * (xy + yx);
    traceless[2] = 0.5 * (xz + zx);

    traceless[3] = traceless[1];
    traceless[4] = yy - thirdTrace;
    traceless[5] = 0.5 * (yz + zy);

    traceless[6] = traceless[2];
    traceless[7] = traceless[5];
    traceless[8] = zz - thirdTrace;

    return {
      traceless,
      trace,
      normalizedTrace: trace
    };
  }

  /**
   * Computes the 3 Principal Invariants (I_1, I_2, I_3) and characteristic polynomial discriminant
   * of a 3x3 symmetric matrix A.
   * 
   * Characteristic Equation: det(lambda * I - A) = lambda^3 - I_1 * lambda^2 + I_2 * lambda - I_3 = 0.
   * 
   * @param {ArrayLike<number>} A 3x3 flat matrix (row-major).
   * @returns {{ I1: number, I2: number, I3: number, J2: number, J3: number, discriminant: number }}
   */
  static computeInvariants(A) {
    if (!A || A.length < 9) {
      throw new Error('TidalTensorMath.computeInvariants: Matrix A must have length >= 9.');
    }

    const xx = A[0], xy = A[1], xz = A[2];
    const yx = A[3], yy = A[4], yz = A[5];
    const zx = A[6], zy = A[7], zz = A[8];

    // First Invariant: Trace
    const I1 = xx + yy + zz;

    // Second Invariant: Sum of principal minors
    // I_2 = (A_xx*A_yy - A_xy*A_yx) + (A_yy*A_zz - A_yz*A_zy) + (A_xx*A_zz - A_xz*A_zx)
    const I2 = (xx * yy - xy * yx) + (yy * zz - yz * zy) + (xx * zz - xz * zx);

    // Third Invariant: Determinant
    const I3 = xx * (yy * zz - yz * zy) - xy * (yx * zz - yz * zx) + xz * (yx * zy - yy * zx);

    // Deviatoric / Traceless invariants
    // J_2 = -0.5 * Tr(dev(A)^2) = I_2 - (1/3)*I_1^2
    const J2 = I2 - (I1 * I1) / 3.0;

    // J_3 = det(dev(A)) = I_3 - (1/3)*I_1*I_2 + (2/27)*I_1^3
    const J3 = I3 - (I1 * I2) / 3.0 + (2.0 * I1 * I1 * I1) / 27.0;

    // Polynomial Discriminant Delta:
    // Delta = 18*I_1*I_2*I_3 - 4*I_1^3*I_3 + I_1^2*I_2^2 - 4*I_2^3 - 27*I_3^2
    // For depressed cubic: Delta_depressed = -4*p^3 - 27*q^2 where p = J_2, q = -J_3
    const discriminant = 18.0 * I1 * I2 * I3 - 4.0 * (I1 ** 3) * I3 + (I1 ** 2) * (I2 ** 2) - 4.0 * (I2 ** 3) - 27.0 * (I3 ** 2);

    return {
      I1,
      I2,
      I3,
      J2,
      J3,
      discriminant
    };
  }

  /**
   * Computes the 3 real eigenvalues of a 3x3 symmetric matrix A using analytical trigonometric Vieta/Cardano formulas.
   * Guaranteed to return sorted eigenvalues: lambda1 >= lambda2 >= lambda3.
   * 
   * @param {ArrayLike<number>} A 3x3 flat array.
   * @returns {Float64Array} Sorted eigenvalues [lambda1, lambda2, lambda3].
   */
  static cardanoEigenvalues(A) {
    if (!A || A.length < 9) {
      throw new Error('TidalTensorMath.cardanoEigenvalues: Matrix A must have length >= 9.');
    }

    const xx = A[0], xy = 0.5 * (A[1] + A[3]), xz = 0.5 * (A[2] + A[6]);
    const yy = A[4], yz = 0.5 * (A[5] + A[7]), zz = A[8];

    // Compute trace and shifted matrix B = A - (Tr(A)/3)*I
    const q = (xx + yy + zz) / 3.0;

    const b00 = xx - q;
    const b11 = yy - q;
    const b22 = zz - q;

    // Frobenius norm squared of deviatoric part: Tr(B^2) = sum B_ij^2
    const pSq = (b00 * b00 + b11 * b11 + b22 * b22 + 2.0 * (xy * xy + xz * xz + yz * yz)) / 6.0;

    if (pSq <= 1e-30) {
      // Degenerate scalar isotropic matrix
      return new Float64Array([q, q, q]);
    }

    const p = Math.sqrt(pSq);

    // Determinant of normalized deviatoric matrix B / p
    // det(B) = b00*(b11*b22 - yz^2) - xy*(xy*b22 - yz*xz) + xz*(xy*yz - b11*xz)
    const detB = b00 * (b11 * b22 - yz * yz) - xy * (xy * b22 - yz * xz) + xz * (xy * yz - b11 * xz);
    const r = detB / (2.0 * p * p * p);

    // Clamp r to [-1, 1] to prevent NaN from rounding errors
    const rClamped = Math.max(-1.0, Math.min(1.0, r));

    const phi = Math.acos(rClamped) / 3.0;

    // Three roots in depressed form
    const eig1 = q + 2.0 * p * Math.cos(phi);
    const eig3 = q + 2.0 * p * Math.cos(phi + (2.0 * Math.PI) / 3.0);
    const eig2 = 3.0 * q - eig1 - eig3; // Using Tr(A) = eig1 + eig2 + eig3

    // Sort descending: lambda1 >= lambda2 >= lambda3
    const roots = [eig1, eig2, eig3];
    roots.sort((a, b) => b - a);

    return new Float64Array(roots);
  }

  /**
   * Computes the complete eigensystem (eigenvalues and orthonormal eigenvectors) of a 3x3 symmetric
   * matrix using the cyclic Jacobi rotation method.
   * 
   * @param {ArrayLike<number>} A 3x3 symmetric matrix (flat array).
   * @param {number} [maxIter=100] Maximum Jacobi sweep iterations.
   * @param {number} [epsilon=1e-15] Convergence tolerance on off-diagonal Frobenius norm.
   * @returns {{ eigenvalues: Float64Array, eigenvectors: Array<Float64Array>, iterations: number }}
   */
  static jacobiDiagonalization(A, maxIter = 100, epsilon = 1e-15) {
    if (!A || A.length < 9) {
      throw new Error('TidalTensorMath.jacobiDiagonalization: Matrix A must have length >= 9.');
    }

    // Copy symmetric matrix to working buffer
    const D = new Float64Array(9);
    D[0] = A[0];
    D[1] = 0.5 * (A[1] + A[3]);
    D[2] = 0.5 * (A[2] + A[6]);

    D[3] = D[1];
    D[4] = A[4];
    D[5] = 0.5 * (A[5] + A[7]);

    D[6] = D[2];
    D[7] = D[5];
    D[8] = A[8];

    // Initialize eigenvector matrix V as 3x3 identity
    const V = new Float64Array([
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0
    ]);

    let iter = 0;
    for (; iter < maxIter; iter++) {
      // Off-diagonal norm
      const offDiag = Math.sqrt(D[1] * D[1] + D[2] * D[2] + D[5] * D[5]);
      if (offDiag < epsilon) {
        break;
      }

      // Sweep through pairs (0,1), (0,2), (1,2)
      const pairs = [[0, 1], [0, 2], [1, 2]];
      for (const [p, q] of pairs) {
        const idxPQ = p * 3 + q;
        const app = D[p * 3 + p];
        const aqq = D[q * 3 + q];
        const apq = D[idxPQ];

        if (Math.abs(apq) < 1e-18) continue;

        const theta = (aqq - app) / (2.0 * apq);
        let t;
        if (Math.abs(theta) > 1e10) {
          t = 1.0 / (2.0 * theta);
        } else {
          t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta * theta));
          if (theta < 0.0) t = -t;
        }

        const c = 1.0 / Math.sqrt(1.0 + t * t);
        const s = t * c;
        const tau = s / (1.0 + c);

        // Update diagonal elements
        D[p * 3 + p] = app - t * apq;
        D[q * 3 + q] = aqq + t * apq;
        D[idxPQ] = 0.0;
        D[q * 3 + p] = 0.0;

        // Update remaining off-diagonal elements
        for (let r = 0; r < 3; r++) {
          if (r !== p && r !== q) {
            const arp = D[r * 3 + p];
            const arq = D[r * 3 + q];
            const newArp = arp - s * (arq + tau * arp);
            const newArq = arq + s * (arp - tau * arq);

            D[r * 3 + p] = newArp;
            D[p * 3 + r] = newArp;
            D[r * 3 + q] = newArq;
            D[q * 3 + r] = newArq;
          }
        }

        // Accumulate rotation into V
        for (let r = 0; r < 3; r++) {
          const vrp = V[r * 3 + p];
          const vrq = V[r * 3 + q];
          V[r * 3 + p] = vrp - s * (vrq + tau * vrp);
          V[r * 3 + q] = vrq + s * (vrp - tau * vrq);
        }
      }
    }

    // Extract eigenvalues and eigenvector column vectors
    const eigValues = [D[0], D[4], D[8]];
    const eigVectors = [
      new Float64Array([V[0], V[3], V[6]]), // Col 0
      new Float64Array([V[1], V[4], V[7]]), // Col 1
      new Float64Array([V[2], V[5], V[8]])  // Col 2
    ];

    // Sort eigenvalues descending (lambda1 >= lambda2 >= lambda3) and permute eigenvectors
    const indices = [0, 1, 2];
    indices.sort((a, b) => eigValues[b] - eigValues[a]);

    const sortedValues = new Float64Array(3);
    const sortedVectors = [];

    for (let k = 0; k < 3; k++) {
      const idx = indices[k];
      sortedValues[k] = eigValues[idx];

      // Normalize eigenvector
      const vec = eigVectors[idx];
      let norm = Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
      if (norm < 1e-12) norm = 1.0;

      sortedVectors.push(new Float64Array([
        vec[0] / norm,
        vec[1] / norm,
        vec[2] / norm
      ]));
    }

    // Ensure right-handed orthonormal SO(3) coordinate system: det(R) = +1
    const e1 = sortedVectors[0];
    const e2 = sortedVectors[1];
    const e3 = sortedVectors[2];

    const crossX = e1[1] * e2[2] - e1[2] * e2[1];
    const crossY = e1[2] * e2[0] - e1[0] * e2[2];
    const crossZ = e1[0] * e2[1] - e1[1] * e2[0];

    const dot = crossX * e3[0] + crossY * e3[1] + crossZ * e3[2];
    if (dot < 0.0) {
      // Flip e3 to make right-handed
      e3[0] = -e3[0];
      e3[1] = -e3[1];
      e3[2] = -e3[2];
    }

    return {
      eigenvalues: sortedValues,
      eigenvectors: sortedVectors,
      iterations: iter
    };
  }

  /**
   * Computes the dimensionless cosmological anisotropy, ellipticity, prolateness,
   * tidal shear magnitude, and triaxiality parameters from sorted eigenvalues.
   * 
   * @param {ArrayLike<number>} eigenvalues Sorted eigenvalues [lambda1, lambda2, lambda3].
   * @returns {{ anisotropy: number, ellipticity: number, prolateness: number, shearMagnitude: number, triaxiality: number, sphericity: number }}
   */
  static computeAnisotropyParameters(eigenvalues) {
    const l1 = eigenvalues[0];
    const l2 = eigenvalues[1];
    const l3 = eigenvalues[2];

    const trace = l1 + l2 + l3;
    const absTrace = Math.abs(trace);
    const denom = absTrace > 1e-12 ? 2.0 * trace : 2.0 * 1e-12;

    // Anisotropy parameter alpha = (lambda1 - lambda3) / (2 * sum(lambda_i))
    const anisotropy = (l1 - l3) / denom;

    // Ellipticity e = (lambda1 - lambda3) / (2 * sum(lambda_i))
    const ellipticity = (l1 - l3) / denom;

    // Prolateness p = (lambda1 + lambda3 - 2*lambda2) / (2 * sum(lambda_i))
    const prolateness = (l1 + l3 - 2.0 * l2) / denom;

    // Tidal shear magnitude q = sqrt(0.5 * [(l1-l2)^2 + (l2-l3)^2 + (l3-l1)^2])
    const diff12 = l1 - l2;
    const diff23 = l2 - l3;
    const diff31 = l3 - l1;
    const shearMagnitude = Math.sqrt(0.5 * (diff12 * diff12 + diff23 * diff23 + diff31 * diff31));

    // Triaxiality parameter T_triax = (lambda1^2 - lambda2^2) / (lambda1^2 - lambda3^2)
    const l1Sq = l1 * l1;
    const l2Sq = l2 * l2;
    const l3Sq = l3 * l3;
    const triaxDenom = l1Sq - l3Sq;
    const triaxiality = Math.abs(triaxDenom) > 1e-12
      ? Math.max(0.0, Math.min(1.0, (l1Sq - l2Sq) / triaxDenom))
      : 0.5;

    // Sphericity index S = 1 - 3*e
    const sphericity = 1.0 - 3.0 * Math.abs(ellipticity);

    return {
      anisotropy,
      ellipticity,
      prolateness,
      shearMagnitude,
      triaxiality,
      sphericity
    };
  }

  /**
   * Classifies the Cosmic Web morphology according to the Hahn et al. (2007) / Forero-Romero et al. (2009)
   * threshold criterion based on sorted eigenvalues lambda1 >= lambda2 >= lambda3.
   * 
   * @param {ArrayLike<number>} eigenvalues Sorted eigenvalues [lambda1, lambda2, lambda3].
   * @param {number} [gammaThreshold=0.0] Eigenvalue threshold gamma_th.
   * @returns {number} CosmicWebType (VOID=0, SHEET=1, FILAMENT=2, KNOT=3).
   */
  static classifyCosmicWeb(eigenvalues, gammaThreshold = 0.0) {
    const l1 = eigenvalues[0];
    const l2 = eigenvalues[1];
    const l3 = eigenvalues[2];

    if (l3 > gammaThreshold) {
      return CosmicWebType.KNOT;      // 3 eigenvalues > gamma_th (all 3 axes collapsing)
    } else if (l2 > gammaThreshold) {
      return CosmicWebType.FILAMENT;  // 2 eigenvalues > gamma_th (2 collapsing, 1 expanding)
    } else if (l1 > gammaThreshold) {
      return CosmicWebType.SHEET;     // 1 eigenvalue > gamma_th (1 collapsing, 2 expanding)
    } else {
      return CosmicWebType.VOID;      // 0 eigenvalues > gamma_th (all 3 axes expanding)
    }
  }

  /**
   * Verifies the Cayley-Hamilton Theorem for a 3x3 matrix:
   * A^3 - I_1 * A^2 + I_2 * A - I_3 * I = 0.
   * Returns the Frobenius residual norm ||A^3 - I_1 A^2 + I_2 A - I_3 I||_F.
   * 
   * @param {ArrayLike<number>} A 3x3 flat array.
   * @param {Object} invariants Result from computeInvariants(A).
   * @returns {number} Frobenius norm of the residual matrix.
   */
  static verifyCayleyHamilton(A, invariants) {
    const { I1, I2, I3 } = invariants;

    // A^2 = A * A
    const A2 = new Float64Array(9);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let sum = 0.0;
        for (let k = 0; k < 3; k++) {
          sum += A[i * 3 + k] * A[k * 3 + j];
        }
        A2[i * 3 + j] = sum;
      }
    }

    // A^3 = A^2 * A
    const A3 = new Float64Array(9);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let sum = 0.0;
        for (let k = 0; k < 3; k++) {
          sum += A2[i * 3 + k] * A[k * 3 + j];
        }
        A3[i * 3 + j] = sum;
      }
    }

    // Residual = A^3 - I_1 * A^2 + I_2 * A - I_3 * I
    let frobeniusSq = 0.0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const delta_ij = (i === j) ? 1.0 : 0.0;
        const res = A3[i * 3 + j] - I1 * A2[i * 3 + j] + I2 * A[i * 3 + j] - I3 * delta_ij;
        frobeniusSq += res * res;
      }
    }

    return Math.sqrt(frobeniusSq);
  }

  /**
   * Computes the Tidal Torque Theory (TTT) torque and angular momentum acquisition:
   * tau_i = sum_{j,k,l} eps_ijk T_jl I_lk
   * L_i(t) = a^2(t) * \dot{D}(t) * tau_i
   * 
   * @param {ArrayLike<number>} T_3x3 3x3 Tidal tensor at halo center of mass (flat array).
   * @param {ArrayLike<number>} I_3x3 3x3 Proto-halo inertia tensor (flat array).
   * @param {Object} [options]
   * @param {number} [options.scaleFactor=1.0] Cosmological scale factor a(t).
   * @param {number} [options.growthRateDot=1.0] Time derivative of growth factor \dot{D}(t).
   * @param {number} [options.haloMass=1e12] Halo total mass M in M_sun/h.
   * @param {number} [options.haloRadius=1.0] Virial radius R_vir in Mpc/h.
   * @param {number} [options.circularVelocity=200.0] Virial velocity V_vir in km/s.
   * @param {number} [options.bindingEnergy=1e58] Absolute gravitational binding energy |E|.
   * @returns {TidalTorqueState}
   */
  static computeTidalTorque(T_3x3, I_3x3, options = {}) {
    const a = options.scaleFactor ?? 1.0;
    const Ddot = options.growthRateDot ?? 1.0;
    const mass = options.haloMass ?? 1e12;
    const rVir = options.haloRadius ?? 1.0;
    const vVir = options.circularVelocity ?? 200.0;
    const energy = options.bindingEnergy ?? 1e58;

    // Matrix product M = T * I (3x3)
    const TI = new Float64Array(9);
    const IT = new Float64Array(9);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let sumTI = 0.0;
        let sumIT = 0.0;
        for (let k = 0; k < 3; k++) {
          sumTI += T_3x3[i * 3 + k] * I_3x3[k * 3 + j];
          sumIT += I_3x3[i * 3 + k] * T_3x3[k * 3 + j];
        }
        TI[i * 3 + j] = sumTI;
        IT[i * 3 + j] = sumIT;
      }
    }

    // Commutator [T, I] = T*I - I*T
    const commutator = new Float64Array(9);
    for (let i = 0; i < 9; i++) {
      commutator[i] = TI[i] - IT[i];
    }

    // Torque vector tau_i = eps_ijk T_jl I_lk = eps_ijk (TI)_jk
    // tau_0 = (TI)_12 - (TI)_21
    // tau_1 = (TI)_20 - (TI)_02
    // tau_2 = (TI)_01 - (TI)_10
    const tau_x = TI[1 * 3 + 2] - TI[2 * 3 + 1];
    const tau_y = TI[2 * 3 + 0] - TI[0 * 3 + 2];
    const tau_z = TI[0 * 3 + 1] - TI[1 * 3 + 0];

    const torque = new Float64Array([tau_x, tau_y, tau_z]);

    // Angular momentum L = a^2 * \dot{D} * tau (or a^3 in some conventions)
    const factor = (a * a) * Ddot;
    const L_x = factor * tau_x;
    const L_y = factor * tau_y;
    const L_z = factor * tau_z;
    const angularMomentum = new Float64Array([L_x, L_y, L_z]);
    const L_mag = Math.sqrt(L_x * L_x + L_y * L_y + L_z * L_z);

    // Compute principal axes of T and I to find misalignment angle
    const diagT = TidalTensorMath.jacobiDiagonalization(T_3x3);
    const diagI = TidalTensorMath.jacobiDiagonalization(I_3x3);

    // Major principal axis of T: e_T1
    const eT1 = diagT.eigenvectors[0];
    // Major principal axis of I: e_I1
    const eI1 = diagI.eigenvectors[0];

    const dotTI = Math.abs(eT1[0] * eI1[0] + eT1[1] * eI1[1] + eT1[2] * eI1[2]);
    const misalignmentAngleRad = Math.acos(Math.min(1.0, Math.max(0.0, dotTI)));

    // Peebles spin parameter: lambda = L * |E|^0.5 / (G * M^2.5)
    const G = CosmologicalConstants.G_ASTRO;
    const peeblesSpin = (L_mag * Math.sqrt(Math.abs(energy))) / (G * (mass ** 2.5) + 1e-30);

    // Bullock spin parameter: lambda' = L / (sqrt(2) * M * V_vir * R_vir)
    const bullockSpin = L_mag / (Math.SQRT2 * mass * vVir * rVir + 1e-30);

    return new TidalTorqueState({
      angularMomentum,
      torque,
      tidalTensor: new Float64Array(T_3x3),
      inertiaTensor: new Float64Array(I_3x3),
      commutator,
      misalignmentAngleRad,
      peeblesSpin,
      bullockSpin,
      haloMass: mass
    });
  }

  /**
   * Computes the linear cosmological growth factor D(z) normalized such that D(z=0) = 1.0.
   * 
   * In Einstein-de Sitter: D(z) = 1 / (1 + z).
   * In flat LambdaCDM: uses Carroll, Press & Turner (1992) / Eisenstein & Hu (1998) exact formula:
   * D(z) = (g(z) / (1 + z)) / g(0)
   * where g(z) = 2.5 * Omega_m(z) / [Omega_m(z)^(4/7) - Omega_Lambda(z) + (1 + Omega_m(z)/2)*(1 + Omega_Lambda(z)/70)].
   * 
   * @param {number} z Redshift z >= 0.
   * @param {Object} [cosmoParams]
   * @param {string} [cosmoParams.model='flat_lcdm'] Cosmological model.
   * @param {number} [cosmoParams.omegaM=0.315] Matter density parameter.
   * @param {number} [cosmoParams.omegaLambda=0.685] Dark energy density parameter.
   * @returns {number} Normalized growth factor D(z).
   */
  static linearGrowthFactor(z, cosmoParams = {}) {
    const model = cosmoParams.model || CosmologicalModel.FLAT_LCDM;
    const omegaM0 = cosmoParams.omegaM ?? CosmologicalConstants.OMEGA_M_DEFAULT;
    const omegaLambda0 = cosmoParams.omegaLambda ?? CosmologicalConstants.OMEGA_LAMBDA_DEFAULT;

    if (model === CosmologicalModel.EDS || Math.abs(omegaM0 - 1.0) < 1e-5) {
      return 1.0 / (1.0 + Math.max(-0.999, z));
    }

    // Helper for g(Omega_m, Omega_Lambda)
    const computeG = (omM, omL) => {
      const denom = (omM ** (4.0 / 7.0)) - omL + (1.0 + omM / 2.0) * (1.0 + omL / 70.0);
      return (2.5 * omM) / (denom + 1e-30);
    };

    const g0 = computeG(omegaM0, omegaLambda0);

    const onePlusZ = 1.0 + Math.max(-0.999, z);
    const onePlusZ3 = onePlusZ ** 3;
    const E_zSq = omegaM0 * onePlusZ3 + omegaLambda0;
    const omM_z = (omegaM0 * onePlusZ3) / E_zSq;
    const omL_z = omegaLambda0 / E_zSq;

    const gz = computeG(omM_z, omL_z);
    return (gz / (onePlusZ * g0));
  }

  /**
   * Computes the dimensionless growth rate f(z) = d ln D / d ln a \approx \Omega_m(z)^{0.55}.
   * 
   * @param {number} z Redshift.
   * @param {Object} [cosmoParams]
   * @returns {number} f(z) = d ln D / d ln a.
   */
  static linearGrowthRate(z, cosmoParams = {}) {
    const omegaM0 = cosmoParams.omegaM ?? CosmologicalConstants.OMEGA_M_DEFAULT;
    const omegaLambda0 = cosmoParams.omegaLambda ?? CosmologicalConstants.OMEGA_LAMBDA_DEFAULT;

    const onePlusZ = 1.0 + Math.max(-0.999, z);
    const onePlusZ3 = onePlusZ ** 3;
    const E_zSq = omegaM0 * onePlusZ3 + omegaLambda0;
    const omM_z = (omegaM0 * onePlusZ3) / E_zSq;

    // Linder (2005) gamma index = 0.55 for standard GR LambdaCDM
    return omM_z ** 0.55;
  }

  /**
   * Computes the Zel'dovich singularity collapse redshifts and cosmic times along the 3 principal axes:
   * D(z_{\text{col}, i}) = 1 / \lambda_i (for \lambda_i > 0).
   * 
   * @param {ArrayLike<number>} eigenvalues Sorted eigenvalues [lambda1, lambda2, lambda3].
   * @param {Object} [cosmoParams]
   * @param {number} [cosmoParams.H0=100.0] Hubble constant.
   * @param {number} [cosmoParams.omegaM=0.315] Matter density.
   * @param {number} [cosmoParams.omegaLambda=0.685] Dark energy density.
   * @returns {{ collapseRedshifts: Float64Array, collapseTimesGyr: Float64Array }}
   */
  static computeZeldovichCollapse(eigenvalues, cosmoParams = {}) {
    const H0 = cosmoParams.H0 ?? CosmologicalConstants.H0_DEFAULT;
    const omM = cosmoParams.omegaM ?? CosmologicalConstants.OMEGA_M_DEFAULT;
    const omL = cosmoParams.omegaLambda ?? CosmologicalConstants.OMEGA_LAMBDA_DEFAULT;

    // Cosmic age at z=0 in Gyr: t_0 \approx (2 / (3 * H_0 * sqrt(Omega_Lambda))) * asinh(sqrt(Omega_Lambda / Omega_m))
    const H0_in_per_Gyr = (H0 * 1.02271e-3); // km/s/Mpc -> 1/Gyr
    let t0_Gyr;
    if (omL > 1e-4) {
      t0_Gyr = (2.0 / (3.0 * H0_in_per_Gyr * Math.sqrt(omL))) * Math.asinh(Math.sqrt(omL / omM));
    } else {
      t0_Gyr = 2.0 / (3.0 * H0_in_per_Gyr);
    }

    const zCols = new Float64Array([NaN, NaN, NaN]);
    const tCols = new Float64Array([NaN, NaN, NaN]);

    for (let i = 0; i < 3; i++) {
      const lam = eigenvalues[i];
      if (lam <= 1e-6) {
        // Expanding or neutral axis; never collapses in linear regime
        continue;
      }

      const targetD = 1.0 / lam;

      // Invert D(z) = targetD using robust 1D bisection / Newton-Raphson
      let zLow = -0.99;
      let zHigh = 2000.0;

      // Check boundaries
      const dHigh = TidalTensorMath.linearGrowthFactor(zHigh, cosmoParams);
      const dLow = TidalTensorMath.linearGrowthFactor(zLow, cosmoParams);

      if (targetD < dHigh) {
        zCols[i] = zHigh;
      } else if (targetD > dLow) {
        zCols[i] = zLow;
      } else {
        // Bisection to 1e-7 accuracy
        for (let iter = 0; iter < 60; iter++) {
          const zMid = 0.5 * (zLow + zHigh);
          const dMid = TidalTensorMath.linearGrowthFactor(zMid, cosmoParams);
          if (Math.abs(dMid - targetD) < 1e-9) {
            zLow = zMid;
            break;
          }
          if (dMid > targetD) {
            // targetD is smaller, meaning higher redshift
            zLow = zMid;
          } else {
            zHigh = zMid;
          }
        }
        zCols[i] = 0.5 * (zLow + zHigh);
      }

      // Convert redshift z_col to cosmic time t(z_col)
      const zCol = zCols[i];
      if (zCol >= -0.99) {
        const onePz = 1.0 + zCol;
        if (omL > 1e-4) {
          const arg = Math.sqrt((omL / omM) * (onePz ** -3));
          tCols[i] = (2.0 / (3.0 * H0_in_per_Gyr * Math.sqrt(omL))) * Math.asinh(arg);
        } else {
          tCols[i] = t0_Gyr / (onePz ** 1.5);
        }
      }
    }

    return {
      collapseRedshifts: zCols,
      collapseTimesGyr: tCols
    };
  }

  /**
   * Computes the Zel'dovich map state at a specific Lagrangian point q given displacement vector
   * and growth factor D(t).
   * 
   * @param {ArrayLike<number>} q Lagrangian coordinate [qx, qy, qz].
   * @param {ArrayLike<number>} gradPhi Displacement field Psi(q) = grad Phi(q).
   * @param {ArrayLike<number>} deformationTensor 3x3 deformation Hessian Psi_ij(q).
   * @param {number} growthFactor Linear growth factor D(t).
   * @returns {ZeldovichState}
   */
  static evaluateZeldovichState(q, gradPhi, deformationTensor, growthFactor) {
    const x = new Float64Array([
      q[0] - growthFactor * gradPhi[0],
      q[1] - growthFactor * gradPhi[1],
      q[2] - growthFactor * gradPhi[2]
    ]);

    const eigenvalues = TidalTensorMath.cardanoEigenvalues(deformationTensor);
    const l1 = eigenvalues[0];
    const l2 = eigenvalues[1];
    const l3 = eigenvalues[2];

    const j1 = 1.0 - growthFactor * l1;
    const j2 = 1.0 - growthFactor * l2;
    const j3 = 1.0 - growthFactor * l3;
    const jacobianDet = j1 * j2 * j3;

    // Density rho = rho_bar / |J| (with caustic cap)
    const density = Math.abs(jacobianDet) > 1e-5 ? 1.0 / Math.abs(jacobianDet) : 1e5;

    let streamCount = 1;
    if (j1 <= 0) streamCount += 2;
    if (j2 <= 0) streamCount += 2;
    if (j3 <= 0) streamCount += 2;

    return new ZeldovichState({
      lagrangianPos: new Float64Array(q),
      eulerianPos: x,
      displacement: new Float64Array(gradPhi),
      growthFactor,
      jacobianDet,
      density,
      streamCount,
      isCaustic: jacobianDet <= 0.0
    });
  }
}

// ============================================================================
// FIELD-LEVEL ANALYZER & SOLVER
// ============================================================================

/**
 * High-performance Cosmological Gravitational Tidal Tensor Field Analyzer.
 * Supports:
 * 1. Finite-difference differentiation of potential grids (2nd, 4th, 6th order stencils).
 * 2. Spectral 3D FFT Poisson solver for computing potential Phi and tidal shear T_ij directly from density delta(x).
 * 3. Multiscale spatial Gaussian smoothing filters.
 * 4. Full-grid eigenvalue mapping, cosmic web classification, ellipticity/prolateness fields.
 * 5. Volumetric and mass-weighted Cosmic Web morphology statistics.
 */
export class TidalTensorAnalyzer {
  /**
   * Constructs a TidalTensorAnalyzer instance.
   * 
   * @param {ScalarField3D|Object} field Source potential field or density field.
   * @param {Object} [options] Configuration options.
   * @param {boolean} [options.isDensity=false] True if input field is density contrast delta(x); false if potential Phi(x).
   * @param {number} [options.order=DifferentiationOrder.FOURTH] Finite difference stencil order (2, 4, or 6).
   * @param {number} [options.gammaThreshold=0.0] Web classification threshold gamma_th.
   * @param {number} [options.smoothingScale=0.0] Gaussian smoothing scale R_smooth in physical units (0 = unsmoothed).
   * @param {string} [options.boundaryMode=BoundaryMode.PERIODIC] Boundary handling strategy.
   * @param {Object} [options.cosmology] Cosmological parameters (omegaM, omegaLambda, H0).
   */
  constructor(field, options = {}) {
    if (!field || !field.grid) {
      throw new TypeError('TidalTensorAnalyzer: Input field must have a valid GridIndexer.');
    }

    this.grid = field.grid;
    this.field = field;
    this.isDensity = Boolean(options.isDensity);
    this.order = options.order ?? DifferentiationOrder.FOURTH;
    this.gammaThreshold = options.gammaThreshold ?? 0.0;
    this.smoothingScale = options.smoothingScale ?? 0.0;
    this.boundaryMode = options.boundaryMode || this.grid.boundaryMode || BoundaryMode.PERIODIC;
    this.cosmology = Object.freeze({
      omegaM: options.cosmology?.omegaM ?? CosmologicalConstants.OMEGA_M_DEFAULT,
      omegaLambda: options.cosmology?.omegaLambda ?? CosmologicalConstants.OMEGA_LAMBDA_DEFAULT,
      H0: options.cosmology?.H0 ?? CosmologicalConstants.H0_DEFAULT,
      fourPiGRhoBar: options.cosmology?.fourPiGRhoBar ?? CosmologicalConstants.FOUR_PI_G_RHO_BAR
    });

    // Cached potential field and precomputed tensor components
    this._potentialField = null;
    this._tensorFields = null; // [Txx, Txy, Txz, Tyy, Tyz, Tzz]
  }

  /**
   * Lazily computes or retrieves the gravitational potential field Phi(x).
   * If input is density delta(x), solves the Poisson equation via spectral 3D FFT.
   * 
   * @returns {ScalarField3D}
   */
  get potentialField() {
    if (!this._potentialField) {
      if (this.isDensity) {
        this._potentialField = this._solvePoissonSpectral(this.field);
      } else {
        this._potentialField = this.field;
      }
    }
    return this._potentialField;
  }

  /**
   * 1D In-Place Radix-2 Cooley-Tukey FFT.
   * @private
   */
  static _fft1D(real, imag, inverse = false) {
    const n = real.length;
    if ((n & (n - 1)) !== 0) {
      throw new Error(`TidalTensorAnalyzer._fft1D: Length must be power of 2, received ${n}`);
    }

    // Bit-reversal permutation
    let j = 0;
    for (let i = 0; i < n - 1; i++) {
      if (i < j) {
        const tempR = real[i];
        real[i] = real[j];
        real[j] = tempR;

        const tempI = imag[i];
        imag[i] = imag[j];
        imag[j] = tempI;
      }
      let k = n >> 1;
      while (k <= j) {
        j -= k;
        k >>= 1;
      }
      j += k;
    }

    // Butterfly passes
    const sign = inverse ? 1.0 : -1.0;
    for (let len = 2; len <= n; len <<= 1) {
      const halfLen = len >> 1;
      const angle = (sign * 2.0 * Math.PI) / len;
      const wStepR = Math.cos(angle);
      const wStepI = Math.sin(angle);

      for (let i = 0; i < n; i += len) {
        let wr = 1.0;
        let wi = 0.0;
        for (let k = 0; k < halfLen; k++) {
          const idxU = i + k;
          const idxV = i + k + halfLen;

          const vr = real[idxV];
          const vi = imag[idxV];

          const twiddleR = vr * wr - vi * wi;
          const twiddleI = vr * wi + vi * wr;

          const ur = real[idxU];
          const ui = imag[idxU];

          real[idxU] = ur + twiddleR;
          imag[idxU] = ui + twiddleI;

          real[idxV] = ur - twiddleR;
          imag[idxV] = ui - twiddleI;

          const nextWr = wr * wStepR - wi * wStepI;
          const nextWi = wr * wStepI + wi * wStepR;
          wr = nextWr;
          wi = nextWi;
        }
      }
    }

    if (inverse) {
      const invN = 1.0 / n;
      for (let i = 0; i < n; i++) {
        real[i] *= invN;
        imag[i] *= invN;
      }
    }
  }

  /**
   * 3D In-Place FFT on flat grid buffers.
   * @private
   */
  static _fft3D(real, imag, nx, ny, nz, inverse = false) {
    const lineX_R = new Float64Array(nx);
    const lineX_I = new Float64Array(nx);
    const lineY_R = new Float64Array(ny);
    const lineY_I = new Float64Array(ny);
    const lineZ_R = new Float64Array(nz);
    const lineZ_I = new Float64Array(nz);

    // 1. Transform along X (fastest axis)
    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        const offset = (iz * ny + iy) * nx;
        for (let ix = 0; ix < nx; ix++) {
          lineX_R[ix] = real[offset + ix];
          lineX_I[ix] = imag[offset + ix];
        }
        TidalTensorAnalyzer._fft1D(lineX_R, lineX_I, inverse);
        for (let ix = 0; ix < nx; ix++) {
          real[offset + ix] = lineX_R[ix];
          imag[offset + ix] = lineX_I[ix];
        }
      }
    }

    // 2. Transform along Y (middle axis)
    for (let iz = 0; iz < nz; iz++) {
      for (let ix = 0; ix < nx; ix++) {
        for (let iy = 0; iy < ny; iy++) {
          const idx = (iz * ny + iy) * nx + ix;
          lineY_R[iy] = real[idx];
          lineY_I[iy] = imag[idx];
        }
        TidalTensorAnalyzer._fft1D(lineY_R, lineY_I, inverse);
        for (let iy = 0; iy < ny; iy++) {
          const idx = (iz * ny + iy) * nx + ix;
          real[idx] = lineY_R[iy];
          imag[idx] = lineY_I[iy];
        }
      }
    }

    // 3. Transform along Z (slowest axis)
    for (let iy = 0; iy < ny; iy++) {
      for (let ix = 0; ix < nx; ix++) {
        for (let iz = 0; iz < nz; iz++) {
          const idx = (iz * ny + iy) * nx + ix;
          lineZ_R[iz] = real[idx];
          lineZ_I[iz] = imag[idx];
        }
        TidalTensorAnalyzer._fft1D(lineZ_R, lineZ_I, inverse);
        for (let iz = 0; iz < nz; iz++) {
          const idx = (iz * ny + iy) * nx + ix;
          real[idx] = lineZ_R[iz];
          imag[idx] = lineZ_I[iz];
        }
      }
    }
  }

  /**
   * Solves the cosmological Poisson equation \nabla^2 \Phi = 4 \pi G \bar{\rho} a^2 \delta
   * using a periodic 3D spectral FFT Green's function solver.
   * 
   * @param {ScalarField3D} densityField
   * @returns {ScalarField3D}
   * @private
   */
  _solvePoissonSpectral(densityField) {
    const grid = this.grid;
    const nx = grid.nx;
    const ny = grid.ny;
    const nz = grid.nz;
    const total = grid.totalCells;

    const real = new Float64Array(total);
    const imag = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      real[i] = densityField.data[i];
      imag[i] = 0.0;
    }

    // Forward 3D FFT: delta(x) -> delta_hat(k)
    TidalTensorAnalyzer._fft3D(real, imag, nx, ny, nz, false);

    const Lx = grid.boxSize[0];
    const Ly = grid.boxSize[1];
    const Lz = grid.boxSize[2];

    const dkx = (2.0 * Math.PI) / Lx;
    const dky = (2.0 * Math.PI) / Ly;
    const dkz = (2.0 * Math.PI) / Lz;

    const sourceFactor = this.cosmology.fourPiGRhoBar;

    for (let iz = 0; iz < nz; iz++) {
      const kzIdx = iz <= nz / 2 ? iz : iz - nz;
      const kz = kzIdx * dkz;

      for (let iy = 0; iy < ny; iy++) {
        const kyIdx = iy <= ny / 2 ? iy : iy - ny;
        const ky = kyIdx * dky;

        for (let ix = 0; ix < nx; ix++) {
          const kxIdx = ix <= nx / 2 ? ix : ix - nx;
          const kx = kxIdx * dkx;

          const idx = (iz * ny + iy) * nx + ix;
          const kSq = kx * kx + ky * ky + kz * kz;

          if (kSq < 1e-14) {
            // Monopole / zero-mode: gauge fix Phi(k=0) = 0
            real[idx] = 0.0;
            imag[idx] = 0.0;
          } else {
            // Phi_hat(k) = - (4 * pi * G * rho_bar) / k^2 * delta_hat(k)
            const green = -sourceFactor / kSq;
            real[idx] *= green;
            imag[idx] *= green;
          }
        }
      }
    }

    // Inverse 3D FFT: Phi_hat(k) -> Phi(x)
    TidalTensorAnalyzer._fft3D(real, imag, nx, ny, nz, true);

    return new ScalarField3D(grid, real, 'gravitational_potential', 'km^2/s^2');
  }

  /**
   * Computes the 6 unique symmetric components of the 3x3 deformation Hessian
   * Psi_ij = d^2 Phi / (dx_i dx_j) across the entire grid.
   * 
   * Uses high-order finite differences (2nd, 4th, or 6th order central stencils).
   * 
   * @returns {{ xx: Float64Array, xy: Float64Array, xz: Float64Array, yy: Float64Array, yz: Float64Array, zz: Float64Array }}
   */
  computeGridDeformationComponents() {
    if (this._tensorFields) {
      return this._tensorFields;
    }

    const phi = this.potentialField.data;
    const grid = this.grid;
    const nx = grid.nx;
    const ny = grid.ny;
    const nz = grid.nz;
    const total = grid.totalCells;

    const dx = grid.dx;
    const dy = grid.dy;
    const dz = grid.dz;

    const xx = new Float64Array(total);
    const xy = new Float64Array(total);
    const xz = new Float64Array(total);
    const yy = new Float64Array(total);
    const yz = new Float64Array(total);
    const zz = new Float64Array(total);

    const getPhi = (ix, iy, iz) => {
      let gx = ix, gy = iy, gz = iz;
      if (this.boundaryMode === BoundaryMode.PERIODIC) {
        gx = (ix % nx + nx) % nx;
        gy = (iy % ny + ny) % ny;
        gz = (iz % nz + nz) % nz;
      } else {
        gx = Math.max(0, Math.min(nx - 1, ix));
        gy = Math.max(0, Math.min(ny - 1, iy));
        gz = Math.max(0, Math.min(nz - 1, iz));
      }
      return phi[(gz * ny + gy) * nx + gx];
    };

    const order = this.order;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const idx = (iz * ny + iy) * nx + ix;
          const pC = getPhi(ix, iy, iz);

          let d2xx, d2yy, d2zz;
          let d2xy, d2xz, d2yz;

          if (order === DifferentiationOrder.SECOND) {
            // Standard 2nd order 3-point stencils
            const pxP = getPhi(ix + 1, iy, iz);
            const pxM = getPhi(ix - 1, iy, iz);
            d2xx = (pxP - 2.0 * pC + pxM) / (dx * dx);

            const pyP = getPhi(ix, iy + 1, iz);
            const pyM = getPhi(ix, iy - 1, iz);
            d2yy = (pyP - 2.0 * pC + pyM) / (dy * dy);

            const pzP = getPhi(ix, iy, iz + 1);
            const pzM = getPhi(ix, iy, iz - 1);
            d2zz = (pzP - 2.0 * pC + pzM) / (dz * dz);

            // Cross derivatives (2nd order)
            const pXP_YP = getPhi(ix + 1, iy + 1, iz);
            const pXP_YM = getPhi(ix + 1, iy - 1, iz);
            const pXM_YP = getPhi(ix - 1, iy + 1, iz);
            const pXM_YM = getPhi(ix - 1, iy - 1, iz);
            d2xy = (pXP_YP - pXP_YM - pXM_YP + pXM_YM) / (4.0 * dx * dy);

            const pXP_ZP = getPhi(ix + 1, iy, iz + 1);
            const pXP_ZM = getPhi(ix + 1, iy, iz - 1);
            const pXM_ZP = getPhi(ix - 1, iy, iz + 1);
            const pXM_ZM = getPhi(ix - 1, iy, iz - 1);
            d2xz = (pXP_ZP - pXP_ZM - pXM_ZP + pXM_ZM) / (4.0 * dx * dz);

            const pYP_ZP = getPhi(ix, iy + 1, iz + 1);
            const pYP_ZM = getPhi(ix, iy + 1, iz - 1);
            const pYM_ZP = getPhi(ix, iy - 1, iz + 1);
            const pYM_ZM = getPhi(ix, iy - 1, iz - 1);
            d2yz = (pYP_ZP - pYP_ZM - pYM_ZP + pYM_ZM) / (4.0 * dy * dz);
          } else if (order === DifferentiationOrder.SIXTH) {
            // 6th-order 7-point central difference stencil for 2nd derivative:
            // d2f/dx2 = (2*f(x-3) - 27*f(x-2) + 270*f(x-1) - 490*f(x) + 270*f(x+1) - 27*f(x+2) + 2*f(x+3)) / (180*h^2)
            const pxP1 = getPhi(ix + 1, iy, iz), pxM1 = getPhi(ix - 1, iy, iz);
            const pxP2 = getPhi(ix + 2, iy, iz), pxM2 = getPhi(ix - 2, iy, iz);
            const pxP3 = getPhi(ix + 3, iy, iz), pxM3 = getPhi(ix - 3, iy, iz);
            d2xx = (2.0 * (pxP3 + pxM3) - 27.0 * (pxP2 + pxM2) + 270.0 * (pxP1 + pxM1) - 490.0 * pC) / (180.0 * dx * dx);

            const pyP1 = getPhi(ix, iy + 1, iz), pyM1 = getPhi(ix, iy - 1, iz);
            const pyP2 = getPhi(ix, iy + 2, iz), pyM2 = getPhi(ix, iy - 2, iz);
            const pyP3 = getPhi(ix, iy + 3, iz), pyM3 = getPhi(ix, iy - 3, iz);
            d2yy = (2.0 * (pyP3 + pyM3) - 27.0 * (pyP2 + pyM2) + 270.0 * (pyP1 + pyM1) - 490.0 * pC) / (180.0 * dy * dy);

            const pzP1 = getPhi(ix, iy, iz + 1), pzM1 = getPhi(ix, iy, iz - 1);
            const pzP2 = getPhi(ix, iy, iz + 2), pzM2 = getPhi(ix, iy, iz - 2);
            const pzP3 = getPhi(ix, iy, iz + 3), pzM3 = getPhi(ix, iy, iz - 3);
            d2zz = (2.0 * (pzP3 + pzM3) - 27.0 * (pzP2 + pzM2) + 270.0 * (pzP1 + pzM1) - 490.0 * pC) / (180.0 * dz * dz);

            // 4th order cross derivatives
            const df_dy_at = (xOffset) => {
              const pY_P2 = getPhi(ix + xOffset, iy + 2, iz);
              const pY_P1 = getPhi(ix + xOffset, iy + 1, iz);
              const pY_M1 = getPhi(ix + xOffset, iy - 1, iz);
              const pY_M2 = getPhi(ix + xOffset, iy - 2, iz);
              return (-pY_P2 + 8.0 * pY_P1 - 8.0 * pY_M1 + pY_M2) / (12.0 * dy);
            };
            d2xy = (-df_dy_at(2) + 8.0 * df_dy_at(1) - 8.0 * df_dy_at(-1) + df_dy_at(-2)) / (12.0 * dx);

            const df_dz_at_x = (xOffset) => {
              const pZ_P2 = getPhi(ix + xOffset, iy, iz + 2);
              const pZ_P1 = getPhi(ix + xOffset, iy, iz + 1);
              const pZ_M1 = getPhi(ix + xOffset, iy, iz - 1);
              const pZ_M2 = getPhi(ix + xOffset, iy, iz - 2);
              return (-pZ_P2 + 8.0 * pZ_P1 - 8.0 * pZ_M1 + pZ_M2) / (12.0 * dz);
            };
            d2xz = (-df_dz_at_x(2) + 8.0 * df_dz_at_x(1) - 8.0 * df_dz_at_x(-1) + df_dz_at_x(-2)) / (12.0 * dx);

            const df_dz_at_y = (yOffset) => {
              const pZ_P2 = getPhi(ix, iy + yOffset, iz + 2);
              const pZ_P1 = getPhi(ix, iy + yOffset, iz + 1);
              const pZ_M1 = getPhi(ix, iy + yOffset, iz - 1);
              const pZ_M2 = getPhi(ix, iy + yOffset, iz - 2);
              return (-pZ_P2 + 8.0 * pZ_P1 - 8.0 * pZ_M1 + pZ_M2) / (12.0 * dz);
            };
            d2yz = (-df_dz_at_y(2) + 8.0 * df_dz_at_y(1) - 8.0 * df_dz_at_y(-1) + df_dz_at_y(-2)) / (12.0 * dy);
          } else {
            // Default: 4th-order 5-point central difference stencil:
            // d2f/dx2 = (-f(x-2) + 16*f(x-1) - 30*f(x) + 16*f(x+1) - f(x+2)) / (12*h^2)
            const pxP1 = getPhi(ix + 1, iy, iz), pxM1 = getPhi(ix - 1, iy, iz);
            const pxP2 = getPhi(ix + 2, iy, iz), pxM2 = getPhi(ix - 2, iy, iz);
            d2xx = (-pxP2 + 16.0 * pxP1 - 30.0 * pC + 16.0 * pxM1 - pxM2) / (12.0 * dx * dx);

            const pyP1 = getPhi(ix, iy + 1, iz), pyM1 = getPhi(ix, iy - 1, iz);
            const pyP2 = getPhi(ix, iy + 2, iz), pyM2 = getPhi(ix, iy - 2, iz);
            d2yy = (-pyP2 + 16.0 * pyP1 - 30.0 * pC + 16.0 * pyM1 - pyM2) / (12.0 * dy * dy);

            const pzP1 = getPhi(ix, iy, iz + 1), pzM1 = getPhi(ix, iy, iz - 1);
            const pzP2 = getPhi(ix, iy, iz + 2), pzM2 = getPhi(ix, iy, iz - 2);
            d2zz = (-pzP2 + 16.0 * pzP1 - 30.0 * pC + 16.0 * pzM1 - pzM2) / (12.0 * dz * dz);

            // 4th order cross derivatives: d2f/dxdy
            const df_dy_at = (xOffset) => {
              const pY_P2 = getPhi(ix + xOffset, iy + 2, iz);
              const pY_P1 = getPhi(ix + xOffset, iy + 1, iz);
              const pY_M1 = getPhi(ix + xOffset, iy - 1, iz);
              const pY_M2 = getPhi(ix + xOffset, iy - 2, iz);
              return (-pY_P2 + 8.0 * pY_P1 - 8.0 * pY_M1 + pY_M2) / (12.0 * dy);
            };
            d2xy = (-df_dy_at(2) + 8.0 * df_dy_at(1) - 8.0 * df_dy_at(-1) + df_dy_at(-2)) / (12.0 * dx);

            const df_dz_at_x = (xOffset) => {
              const pZ_P2 = getPhi(ix + xOffset, iy, iz + 2);
              const pZ_P1 = getPhi(ix + xOffset, iy, iz + 1);
              const pZ_M1 = getPhi(ix + xOffset, iy, iz - 1);
              const pZ_M2 = getPhi(ix + xOffset, iy, iz - 2);
              return (-pZ_P2 + 8.0 * pZ_P1 - 8.0 * pZ_M1 + pZ_M2) / (12.0 * dz);
            };
            d2xz = (-df_dz_at_x(2) + 8.0 * df_dz_at_x(1) - 8.0 * df_dz_at_x(-1) + df_dz_at_x(-2)) / (12.0 * dx);

            const df_dz_at_y = (yOffset) => {
              const pZ_P2 = getPhi(ix, iy + yOffset, iz + 2);
              const pZ_P1 = getPhi(ix, iy + yOffset, iz + 1);
              const pZ_M1 = getPhi(ix, iy + yOffset, iz - 1);
              const pZ_M2 = getPhi(ix, iy + yOffset, iz - 2);
              return (-pZ_P2 + 8.0 * pZ_P1 - 8.0 * pZ_M1 + pZ_M2) / (12.0 * dz);
            };
            d2yz = (-df_dz_at_y(2) + 8.0 * df_dz_at_y(1) - 8.0 * df_dz_at_y(-1) + df_dz_at_y(-2)) / (12.0 * dy);
          }

          xx[idx] = d2xx;
          xy[idx] = d2xy;
          xz[idx] = d2xz;
          yy[idx] = d2yy;
          yz[idx] = d2yz;
          zz[idx] = d2zz;
        }
      }
    }

    this._tensorFields = { xx, xy, xz, yy, yz, zz };
    return this._tensorFields;
  }

  /**
   * Evaluates the complete 3x3 Tidal Tensor Invariants at a specific grid index (ix, iy, iz).
   * 
   * @param {number} ix Grid X index.
   * @param {number} iy Grid Y index.
   * @param {number} iz Grid Z index.
   * @param {Object} [options]
   * @param {number} [options.gammaThreshold] Override classification threshold.
   * @returns {TidalTensorInvariants}
   */
  evaluateAtVoxel(ix, iy, iz, options = {}) {
    const comp = this.computeGridDeformationComponents();
    const grid = this.grid;
    const nx = grid.nx;
    const ny = grid.ny;
    const nz = grid.nz;

    let gx = ix, gy = iy, gz = iz;
    if (this.boundaryMode === BoundaryMode.PERIODIC) {
      gx = (ix % nx + nx) % nx;
      gy = (iy % ny + ny) % ny;
      gz = (iz % nz + nz) % nz;
    } else {
      gx = Math.max(0, Math.min(nx - 1, ix));
      gy = Math.max(0, Math.min(ny - 1, iy));
      gz = Math.max(0, Math.min(nz - 1, iz));
    }

    const idx = (gz * ny + gy) * nx + gx;

    const xx = comp.xx[idx];
    const xy = comp.xy[idx];
    const xz = comp.xz[idx];
    const yy = comp.yy[idx];
    const yz = comp.yz[idx];
    const zz = comp.zz[idx];

    const Psi = new Float64Array([
      xx, xy, xz,
      xy, yy, yz,
      xz, yz, zz
    ]);

    const { traceless, trace } = TidalTensorMath.decomposeTidalTensor(Psi);
    const invariants = TidalTensorMath.computeInvariants(Psi);
    const shearInvariants = TidalTensorMath.computeInvariants(traceless);

    // Compute exact eigenvalues and eigenvectors
    const diag = TidalTensorMath.jacobiDiagonalization(Psi);
    const eigenvalues = diag.eigenvalues;
    const eigenvectors = diag.eigenvectors;

    // Anisotropy and morphology
    const aniso = TidalTensorMath.computeAnisotropyParameters(eigenvalues);

    const gammaTh = options.gammaThreshold ?? this.gammaThreshold;
    const webType = TidalTensorMath.classifyCosmicWeb(eigenvalues, gammaTh);

    // Zel'dovich collapse metrics
    const collapse = TidalTensorMath.computeZeldovichCollapse(eigenvalues, this.cosmology);

    return new TidalTensorInvariants({
      tensor: Psi,
      tracelessTensor: traceless,
      trace,
      secondInvariant: invariants.I2,
      thirdInvariant: invariants.I3,
      shearSecondInvariant: shearInvariants.I2,
      shearThirdInvariant: shearInvariants.I3,
      discriminant: invariants.discriminant,
      eigenvalues,
      eigenvectors,
      anisotropy: aniso.anisotropy,
      ellipticity: aniso.ellipticity,
      prolateness: aniso.prolateness,
      shearMagnitude: aniso.shearMagnitude,
      triaxiality: aniso.triaxiality,
      sphericity: aniso.sphericity,
      webType,
      gammaThreshold: gammaTh,
      collapseRedshifts: collapse.collapseRedshifts,
      collapseTimesGyr: collapse.collapseTimesGyr
    });
  }

  /**
   * Evaluates the continuous Tidal Tensor Invariants at an arbitrary physical coordinate (x, y, z)
   * in Mpc/h using trilinear interpolation on the 6 independent Hessian components.
   * 
   * @param {number} x Physical X position.
   * @param {number} y Physical Y position.
   * @param {number} z Physical Z position.
   * @param {Object} [options]
   * @returns {TidalTensorInvariants}
   */
  evaluateAtCoordinate(x, y, z, options = {}) {
    const comp = this.computeGridDeformationComponents();
    const grid = this.grid;

    const [gx, gy, gz] = grid.coordToGridIndex(x, y, z);
    const nx = grid.nx;
    const ny = grid.ny;
    const nz = grid.nz;

    const i0 = Math.floor(gx);
    const j0 = Math.floor(gy);
    const k0 = Math.floor(gz);

    const tx = gx - i0;
    const ty = gy - j0;
    const tz = gz - k0;

    const interpBuffer = (buffer) => {
      let val = 0.0;
      for (let dk = 0; dk <= 1; dk++) {
        const kz = k0 + dk;
        const wz = dk === 0 ? (1.0 - tz) : tz;

        let gzIdx = kz;
        if (this.boundaryMode === BoundaryMode.PERIODIC) {
          gzIdx = (kz % nz + nz) % nz;
        } else {
          gzIdx = Math.max(0, Math.min(nz - 1, kz));
        }

        for (let dj = 0; dj <= 1; dj++) {
          const jy = j0 + dj;
          const wy = dj === 0 ? (1.0 - ty) : ty;

          let gyIdx = jy;
          if (this.boundaryMode === BoundaryMode.PERIODIC) {
            gyIdx = (jy % ny + ny) % ny;
          } else {
            gyIdx = Math.max(0, Math.min(ny - 1, jy));
          }

          for (let di = 0; di <= 1; di++) {
            const ix = i0 + di;
            const wx = di === 0 ? (1.0 - tx) : tx;

            let gxIdx = ix;
            if (this.boundaryMode === BoundaryMode.PERIODIC) {
              gxIdx = (ix % nx + nx) % nx;
            } else {
              gxIdx = Math.max(0, Math.min(nx - 1, ix));
            }

            const idx = (gzIdx * ny + gyIdx) * nx + gxIdx;
            val += wx * wy * wz * buffer[idx];
          }
        }
      }
      return val;
    };

    const xx = interpBuffer(comp.xx);
    const xy = interpBuffer(comp.xy);
    const xz = interpBuffer(comp.xz);
    const yy = interpBuffer(comp.yy);
    const yz = interpBuffer(comp.yz);
    const zz = interpBuffer(comp.zz);

    const Psi = new Float64Array([
      xx, xy, xz,
      xy, yy, yz,
      xz, yz, zz
    ]);

    const { traceless, trace } = TidalTensorMath.decomposeTidalTensor(Psi);
    const invariants = TidalTensorMath.computeInvariants(Psi);
    const shearInvariants = TidalTensorMath.computeInvariants(traceless);

    const diag = TidalTensorMath.jacobiDiagonalization(Psi);
    const eigenvalues = diag.eigenvalues;
    const eigenvectors = diag.eigenvectors;

    const aniso = TidalTensorMath.computeAnisotropyParameters(eigenvalues);
    const gammaTh = options.gammaThreshold ?? this.gammaThreshold;
    const webType = TidalTensorMath.classifyCosmicWeb(eigenvalues, gammaTh);
    const collapse = TidalTensorMath.computeZeldovichCollapse(eigenvalues, this.cosmology);

    return new TidalTensorInvariants({
      tensor: Psi,
      tracelessTensor: traceless,
      trace,
      secondInvariant: invariants.I2,
      thirdInvariant: invariants.I3,
      shearSecondInvariant: shearInvariants.I2,
      shearThirdInvariant: shearInvariants.I3,
      discriminant: invariants.discriminant,
      eigenvalues,
      eigenvectors,
      anisotropy: aniso.anisotropy,
      ellipticity: aniso.ellipticity,
      prolateness: aniso.prolateness,
      shearMagnitude: aniso.shearMagnitude,
      triaxiality: aniso.triaxiality,
      sphericity: aniso.sphericity,
      webType,
      gammaThreshold: gammaTh,
      collapseRedshifts: collapse.collapseRedshifts,
      collapseTimesGyr: collapse.collapseTimesGyr
    });
  }

  /**
   * Generates continuous 3D ScalarField3D maps for eigenvalues (lambda1, lambda2, lambda3),
   * cosmic web classification, ellipticity, prolateness, and tidal shear magnitude.
   * 
   * @param {Object} [options]
   * @param {number} [options.gammaThreshold] Override classification threshold.
   * @returns {{
   *   lambda1: ScalarField3D,
   *   lambda2: ScalarField3D,
   *   lambda3: ScalarField3D,
   *   cosmicWeb: ScalarField3D,
   *   ellipticity: ScalarField3D,
   *   prolateness: ScalarField3D,
   *   shearMagnitude: ScalarField3D,
   *   collapseRedshift1: ScalarField3D
   * }}
   */
  generateFieldMaps(options = {}) {
    const comp = this.computeGridDeformationComponents();
    const grid = this.grid;
    const total = grid.totalCells;

    const l1Buf = new Float64Array(total);
    const l2Buf = new Float64Array(total);
    const l3Buf = new Float64Array(total);
    const webBuf = new Float64Array(total);
    const eBuf = new Float64Array(total);
    const pBuf = new Float64Array(total);
    const qBuf = new Float64Array(total);
    const z1Buf = new Float64Array(total);

    const gammaTh = options.gammaThreshold ?? this.gammaThreshold;
    const PsiLocal = new Float64Array(9);

    for (let i = 0; i < total; i++) {
      const xx = comp.xx[i];
      const xy = comp.xy[i];
      const xz = comp.xz[i];
      const yy = comp.yy[i];
      const yz = comp.yz[i];
      const zz = comp.zz[i];

      PsiLocal[0] = xx; PsiLocal[1] = xy; PsiLocal[2] = xz;
      PsiLocal[3] = xy; PsiLocal[4] = yy; PsiLocal[5] = yz;
      PsiLocal[6] = xz; PsiLocal[7] = yz; PsiLocal[8] = zz;

      const eigenvalues = TidalTensorMath.cardanoEigenvalues(PsiLocal);
      const l1 = eigenvalues[0];
      const l2 = eigenvalues[1];
      const l3 = eigenvalues[2];

      l1Buf[i] = l1;
      l2Buf[i] = l2;
      l3Buf[i] = l3;

      const webType = TidalTensorMath.classifyCosmicWeb(eigenvalues, gammaTh);
      webBuf[i] = webType;

      const aniso = TidalTensorMath.computeAnisotropyParameters(eigenvalues);
      eBuf[i] = aniso.ellipticity;
      pBuf[i] = aniso.prolateness;
      qBuf[i] = aniso.shearMagnitude;

      // 1st collapse redshift (pancake)
      if (l1 > 1e-4) {
        z1Buf[i] = l1 - 1.0; // EdS benchmark or 1st order
      } else {
        z1Buf[i] = -1.0;
      }
    }

    return {
      lambda1: new ScalarField3D(grid, l1Buf, 'eigenvalue_lambda1', 'km^2/s^2/(Mpc/h)^2'),
      lambda2: new ScalarField3D(grid, l2Buf, 'eigenvalue_lambda2', 'km^2/s^2/(Mpc/h)^2'),
      lambda3: new ScalarField3D(grid, l3Buf, 'eigenvalue_lambda3', 'km^2/s^2/(Mpc/h)^2'),
      cosmicWeb: new ScalarField3D(grid, webBuf, 'cosmic_web_classification', 'type_index'),
      ellipticity: new ScalarField3D(grid, eBuf, 'tidal_ellipticity', 'dimensionless'),
      prolateness: new ScalarField3D(grid, pBuf, 'tidal_prolateness', 'dimensionless'),
      shearMagnitude: new ScalarField3D(grid, qBuf, 'tidal_shear_magnitude', 'km^2/s^2/(Mpc/h)^2'),
      collapseRedshift1: new ScalarField3D(grid, z1Buf, 'pancake_collapse_redshift', 'redshift')
    };
  }

  /**
   * Computes the global cosmological volume and mass filling fractions and morphological statistics
   * across the grid for a given threshold gamma_th, with optional threshold scanning.
   * 
   * @param {Object} [options]
   * @param {number} [options.gammaThreshold] Override classification threshold.
   * @param {Array<number>} [options.thresholdScanRange] Optional list of thresholds to evaluate.
   * @param {ScalarField3D} [options.densityContrastField] Optional density field for mass-weighted fractions.
   * @returns {CosmicWebSummary}
   */
  computeCosmicWebSummary(options = {}) {
    const gammaTh = options.gammaThreshold ?? this.gammaThreshold;
    const comp = this.computeGridDeformationComponents();
    const grid = this.grid;
    const total = grid.totalCells;

    const densityData = options.densityContrastField?.data || (this.isDensity ? this.field.data : null);

    const counts = new Float64Array(4); // Void, Sheet, Filament, Knot
    const massCounts = new Float64Array(4);
    let totalMass = 0.0;

    const sumTrace = new Float64Array(4);
    const sumEllipticity = new Float64Array(4);
    const sumProlateness = new Float64Array(4);
    const sumShear = new Float64Array(4);

    const PsiLocal = new Float64Array(9);

    for (let i = 0; i < total; i++) {
      const xx = comp.xx[i];
      const xy = comp.xy[i];
      const xz = comp.xz[i];
      const yy = comp.yy[i];
      const yz = comp.yz[i];
      const zz = comp.zz[i];

      PsiLocal[0] = xx; PsiLocal[1] = xy; PsiLocal[2] = xz;
      PsiLocal[3] = xy; PsiLocal[4] = yy; PsiLocal[5] = yz;
      PsiLocal[6] = xz; PsiLocal[7] = yz; PsiLocal[8] = zz;

      const eigenvalues = TidalTensorMath.cardanoEigenvalues(PsiLocal);
      const webType = TidalTensorMath.classifyCosmicWeb(eigenvalues, gammaTh);
      const aniso = TidalTensorMath.computeAnisotropyParameters(eigenvalues);

      counts[webType] += 1.0;

      const delta = densityData ? densityData[i] : (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
      const massWeight = Math.max(0.0, 1.0 + delta);
      massCounts[webType] += massWeight;
      totalMass += massWeight;

      sumTrace[webType] += (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
      sumEllipticity[webType] += aniso.ellipticity;
      sumProlateness[webType] += aniso.prolateness;
      sumShear[webType] += aniso.shearMagnitude;
    }

    const volumeFractions = new Float64Array(4);
    const massFractions = new Float64Array(4);

    for (let k = 0; k < 4; k++) {
      volumeFractions[k] = counts[k] / total;
      massFractions[k] = totalMass > 0 ? (massCounts[k] / totalMass) : volumeFractions[k];
    }

    const meanInvariantsPerType = {};
    for (let k = 0; k < 4; k++) {
      const n = counts[k] > 0 ? counts[k] : 1.0;
      meanInvariantsPerType[k] = {
        meanTrace: sumTrace[k] / n,
        meanEllipticity: sumEllipticity[k] / n,
        meanProlateness: sumProlateness[k] / n,
        meanShear: sumShear[k] / n
      };
    }

    // Optional scan over a series of thresholds
    const scanRange = options.thresholdScanRange || null;
    const thresholdScan = [];

    if (scanRange && Array.isArray(scanRange)) {
      for (const th of scanRange) {
        const scCounts = new Float64Array(4);
        for (let i = 0; i < total; i++) {
          const xx = comp.xx[i];
          const xy = comp.xy[i];
          const xz = comp.xz[i];
          const yy = comp.yy[i];
          const yz = comp.yz[i];
          const zz = comp.zz[i];

          PsiLocal[0] = xx; PsiLocal[1] = xy; PsiLocal[2] = xz;
          PsiLocal[3] = xy; PsiLocal[4] = yy; PsiLocal[5] = yz;
          PsiLocal[6] = xz; PsiLocal[7] = yz; PsiLocal[8] = zz;

          const eigenvalues = TidalTensorMath.cardanoEigenvalues(PsiLocal);
          const type = TidalTensorMath.classifyCosmicWeb(eigenvalues, th);
          scCounts[type] += 1.0;
        }

        thresholdScan.push({
          gammaThreshold: th,
          voidFraction: scCounts[0] / total,
          sheetFraction: scCounts[1] / total,
          filamentFraction: scCounts[2] / total,
          knotFraction: scCounts[3] / total
        });
      }
    }

    return new CosmicWebSummary({
      totalVoxels: total,
      gammaThreshold: gammaTh,
      volumeFractions,
      massFractions,
      voxelCounts: counts,
      meanInvariantsPerType,
      thresholdScan
    });
  }
}

// ============================================================================
// ANALYTICAL BENCHMARK & TEST FIELD GENERATORS
// ============================================================================

/**
 * Exact analytical solutions for gravitational potentials, deformation Hessians,
 * traceless tidal tensors, and eigenvalue distributions for rigorous unit testing.
 */
export class AnalyticTidalFields {
  /**
   * Keplerian Point Mass Potential:
   * Phi(r) = - G * M / r
   * d_i d_j Phi(r) = (G * M / r^3) * [delta_ij - 3 * (x_i * x_j / r^2)]
   * Tr(Psi) = 0 for r > 0.
   * Eigenvalues: lambda_parallel = - 2 * G * M / r^3, lambda_perp = + G * M / r^3 (multiplicity 2).
   * 
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {number} [GM=1.0] G * M product.
   * @returns {{ potential: number, tensor: Float64Array, traceless: Float64Array, eigenvalues: Float64Array }}
   */
  static keplerianMonopole(x, y, z, GM = 1.0) {
    const rSq = x * x + y * y + z * z;
    const r = Math.sqrt(rSq);
    if (r < 1e-12) {
      throw new Error('AnalyticTidalFields.keplerianMonopole: Singularity at origin r = 0.');
    }

    const r3 = r * rSq;
    const r5 = r3 * rSq;
    const pot = -GM / r;

    const tensor = new Float64Array(9);
    // Psi_ij = GM * (delta_ij / r^3 - 3 * x_i * x_j / r^5)
    tensor[0] = GM * (1.0 / r3 - (3.0 * x * x) / r5);
    tensor[1] = GM * (-3.0 * x * y / r5);
    tensor[2] = GM * (-3.0 * x * z / r5);

    tensor[3] = tensor[1];
    tensor[4] = GM * (1.0 / r3 - (3.0 * y * y) / r5);
    tensor[5] = GM * (-3.0 * y * z / r5);

    tensor[6] = tensor[2];
    tensor[7] = tensor[5];
    tensor[8] = GM * (1.0 / r3 - (3.0 * z * z) / r5);

    // Exact sorted eigenvalues: lambda1 = GM/r^3, lambda2 = GM/r^3, lambda3 = -2*GM/r^3
    const eigPerp = GM / r3;
    const eigPar = -2.0 * GM / r3;
    const eigenvalues = new Float64Array([eigPerp, eigPerp, eigPar]);

    return {
      potential: pot,
      tensor,
      traceless: new Float64Array(tensor), // Already traceless
      eigenvalues
    };
  }

  /**
   * Triaxial Quadrupole Potential:
   * Phi(x, y, z) = 0.5 * (A * x^2 + B * y^2 + C * z^2)
   * Psi_ij = diag(A, B, C) (constant uniform tidal tensor).
   * Tr(Psi) = A + B + C.
   * 
   * @param {number} A
   * @param {number} B
   * @param {number} C
   * @returns {{ tensor: Float64Array, traceless: Float64Array, eigenvalues: Float64Array }}
   */
  static triaxialQuadrupole(A, B, C) {
    const tensor = new Float64Array([
      A, 0.0, 0.0,
      0.0, B, 0.0,
      0.0, 0.0, C
    ]);

    const trace = A + B + C;
    const thirdTrace = trace / 3.0;

    const traceless = new Float64Array([
      A - thirdTrace, 0.0, 0.0,
      0.0, B - thirdTrace, 0.0,
      0.0, 0.0, C - thirdTrace
    ]);

    const sorted = [A, B, C].sort((a, b) => b - a);
    const eigenvalues = new Float64Array(sorted);

    return {
      tensor,
      traceless,
      eigenvalues
    };
  }

  /**
   * Plummer Sphere Potential:
   * Phi(r) = - G * M / sqrt(r^2 + b^2)
   * d_i d_j Phi(r) = (G * M / (r^2 + b^2)^(3/2)) * delta_ij - 3 * G * M * x_i * x_j / (r^2 + b^2)^(5/2)
   * 
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {number} [GM=1.0]
   * @param {number} [b=1.0] Core softening radius.
   * @returns {{ potential: number, tensor: Float64Array, traceless: Float64Array }}
   */
  static plummerSphere(x, y, z, GM = 1.0, b = 1.0) {
    const rSq = x * x + y * y + z * z;
    const denom1 = Math.sqrt(rSq + b * b);
    const denom3 = denom1 * (rSq + b * b);
    const denom5 = denom3 * (rSq + b * b);

    const pot = -GM / denom1;

    const tensor = new Float64Array(9);
    const term1 = GM / denom3;
    const factor5 = (3.0 * GM) / denom5;

    tensor[0] = term1 - factor5 * x * x;
    tensor[1] = -factor5 * x * y;
    tensor[2] = -factor5 * x * z;

    tensor[3] = tensor[1];
    tensor[4] = term1 - factor5 * y * y;
    tensor[5] = -factor5 * y * z;

    tensor[6] = tensor[2];
    tensor[7] = tensor[5];
    tensor[8] = term1 - factor5 * z * z;

    const { traceless } = TidalTensorMath.decomposeTidalTensor(tensor);

    return {
      potential: pot,
      tensor,
      traceless
    };
  }

  /**
   * 3D Fourier Plane Wave Standing Potential:
   * Phi(x, y, z) = Phi_0 * cos(kx * x + ky * y + kz * z)
   * Psi_ij = - ki * kj * Phi_0 * cos(k \cdot x)
   * 
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {ArrayLike<number>} k Wavevector [kx, ky, kz].
   * @param {number} [Phi0=1.0] Amplitude.
   * @returns {{ potential: number, tensor: Float64Array }}
   */
  static fourierPlaneWave(x, y, z, k, Phi0 = 1.0) {
    const kx = k[0], ky = k[1], kz = k[2];
    const phase = kx * x + ky * y + kz * z;
    const cosVal = Math.cos(phase);

    const pot = Phi0 * cosVal;
    const factor = -Phi0 * cosVal;

    const tensor = new Float64Array([
      factor * kx * kx, factor * kx * ky, factor * kx * kz,
      factor * ky * kx, factor * ky * ky, factor * ky * kz,
      factor * kz * kx, factor * kz * ky, factor * kz * kz
    ]);

    return {
      potential: pot,
      tensor
    };
  }
}
