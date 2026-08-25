/**
 * @file velocity_dispersion_tensor.js
 * @description Comprehensive, production-grade mathematical implementation of the 3D Velocity Dispersion Tensor,
 * Cosmic Thermal and Kinetic Pressure Fields, Spherical Velocity Anisotropy Parameter beta(r),
 * Jeans Equation Cosmological Mass Estimator, and Supercluster Multiscale Virial Diagnostics
 * for the CosmicFlows-4 Computational Research Workbench.
 *
 * Mathematical and Scientific Foundations:
 * ---------------------------------------
 * 1. 3D Velocity Dispersion Tensor:
 *    \sigma_{ij}^2(\mathbf{x}) = \langle v_i v_j \rangle(\mathbf{x}) - \langle v_i \rangle(\mathbf{x}) \langle v_j \rangle(\mathbf{x})
 *    where spatial averaging is performed with cosmological smoothing kernels (Gaussian, Top-Hat, SPH Cubic Spline, Wendland C^4).
 *    Tensor properties: Real, symmetric (3x3), positive semi-definite.
 *    Principal frame eigenvalues: \lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0.
 *    Scalar 1D velocity dispersion: \sigma_{1D} = \sqrt{\frac{1}{3} \text{Tr}(\sigma^2)}.
 *    Scalar 3D velocity dispersion: \sigma_{3D} = \sqrt{\text{Tr}(\sigma^2)} = \sqrt{3} \sigma_{1D}.
 *    Triaxiality parameter: T = \frac{\lambda_1 - \lambda_2}{\lambda_1 - \lambda_3} \in [0, 1].
 *    Ellipticity: e = \frac{\lambda_1 - \lambda_3}{2(\lambda_1 + \lambda_2 + \lambda_3)}.
 *    Prolateness: p = \frac{\lambda_1 - 2\lambda_2 + \lambda_3}{2(\lambda_1 + \lambda_2 + \lambda_3)}.
 *
 * 2. Spherical Halo Velocity Anisotropy Parameter \beta(\mathbf{x}):
 *    In spherical coordinates centered on a cluster/halo (\mathbf{x}_c):
 *    \mathbf{r} = \mathbf{x} - \mathbf{x}_c, \quad r = \|\mathbf{r}\|.
 *    Basis vectors: \hat{\mathbf{e}}_r, \hat{\mathbf{e}}_\theta, \hat{\mathbf{e}}_\phi.
 *    Spherical dispersion tensor: \mathbf{\Sigma}_{\text{sph}} = \mathbf{R} \mathbf{\Sigma}_{\text{cart}} \mathbf{R}^T
 *    Radial dispersion: \sigma_r^2 = \hat{\mathbf{e}}_r^T \mathbf{\Sigma} \hat{\mathbf{e}}_r
 *    Tangential dispersions: \sigma_\theta^2 = \hat{\mathbf{e}}_\theta^T \mathbf{\Sigma} \hat{\mathbf{e}}_\theta, \quad \sigma_\phi^2 = \hat{\mathbf{e}}_\phi^T \mathbf{\Sigma} \hat{\mathbf{e}}_\phi
 *    Combined tangential dispersion: \sigma_t^2 = \frac{\sigma_\theta^2 + \sigma_\phi^2}{2}
 *    Binney (1980) Anisotropy Parameter:
 *      \beta(\mathbf{x}) = 1 - \frac{\sigma_\theta^2 + \sigma_\phi^2}{2 \sigma_r^2} = 1 - \frac{\sigma_t^2}{\sigma_r^2}
 *    Regimes:
 *      - \beta = 0: Isotropic orbit distribution (\sigma_r = \sigma_\theta = \sigma_\phi)
 *      - 0 < \beta \le 1: Radially biased velocity dispersion (pure radial plunge \beta = 1)
 *      - \beta < 0: Tangentially biased velocity dispersion (circular orbits \beta \to -\infty)
 *      - Symmetrized anisotropy: \gamma_\beta = \frac{2\beta}{2-\beta} = \frac{\sigma_r^2 - \sigma_t^2}{\sigma_r^2 + \sigma_t^2} \in [-1, 1]
 *
 * 3. Kinetic Pressure Tensor, Cosmic Thermal Pressure Field & Sound Speed:
 *    Local matter density: \rho(\mathbf{x}) = \bar{\rho}_m (1 + \delta(\mathbf{x}))
 *    Kinetic Pressure Tensor: P_{ij}(\mathbf{x}) = \rho(\mathbf{x}) \sigma_{ij}^2(\mathbf{x})
 *    Scalar Isotropic Pressure: P(\mathbf{x}) = \frac{1}{3} \text{Tr}(P_{ij}) = \rho(\mathbf{x}) \sigma_{1D}^2(\mathbf{x})
 *    Anisotropic Stress Tensor: \Pi_{ij}(\mathbf{x}) = P_{ij}(\mathbf{x}) - P(\mathbf{x}) \delta_{ij}
 *    Effective Cosmic Sound Speed:
 *      c_s(\mathbf{x}) = \sqrt{\gamma \frac{P(\mathbf{x})}{\rho(\mathbf{x})}} = \sqrt{\gamma \sigma_{1D}^2(\mathbf{x})} = \sqrt{\frac{\gamma}{3} \text{Tr}(\sigma^2)}
 *      where \gamma = 5/3 (monoatomic ideal gas) or \gamma = 1 (isothermal).
 *    Cosmic Mach Number: \mathcal{M}(\mathbf{x}) = \frac{\|\mathbf{v}(\mathbf{x})\|}{c_s(\mathbf{x})}
 *    Cosmic Temperature Proxy: T_K(\mathbf{x}) = \frac{\mu m_p}{k_B} \sigma_{1D}^2(\mathbf{x})
 *
 * 4. Spherical Jeans Equation Cosmological Mass Estimator:
 *    From the steady-state, collisionless Boltzmann equation in spherical symmetry:
 *      M_{\text{Jeans}}(<r) = -\frac{r \sigma_r^2(r)}{G} \left( \frac{d \ln \rho}{d \ln r} + \frac{d \ln \sigma_r^2}{d \ln r} + 2\beta(r) \right)
 *    Let \gamma_\rho(r) = -\frac{d \ln \rho}{d \ln r} and \gamma_\sigma(r) = -\frac{d \ln \sigma_r^2}{d \ln r}:
 *      M_{\text{Jeans}}(<r) = \frac{r \sigma_r^2(r)}{G} \left( \gamma_\rho(r) + \gamma_\sigma(r) - 2\beta(r) \right)
 *    Cosmological units:
 *      G = 4.30091727 \times 10^{-9} (\text{km/s})^2 \text{Mpc} M_\odot^{-1}
 *      In h-scaled units: G_h = 4.30091727 \times 10^{-9} h (\text{km/s})^2 (\text{Mpc}/h) M_\odot^{-1}
 *
 * 5. Supercluster Multiscale Virial Ratio 2K / |W|:
 *    Kinetic Energy: K = K_{\text{bulk}} + K_{\text{disp}}
 *      K_{\text{bulk}} = \frac{1}{2} \int_V \rho(\mathbf{x}) \|\mathbf{v}(\mathbf{x}) - \mathbf{v}_{\text{CM}}\|^2 d^3x
 *      K_{\text{disp}} = \frac{1}{2} \int_V \rho(\mathbf{x}) \text{Tr}(\sigma_{ij}^2(\mathbf{x})) d^3x
 *    Gravitational Potential Energy:
 *      W = \frac{1}{2} \int_V \rho(\mathbf{x}) \Phi(\mathbf{x}) d^3x \quad \text{or} \quad W = -4\pi G \int_0^R \rho(r) M(<r) r dr
 *    Virial Ratio: \mathcal{V} = \frac{2K}{|W|}
 *      - \mathcal{V} \approx 1: Dynamical Virial Equilibrium
 *      - \mathcal{V} < 1: Super-virial / Gravitationally collapsing structure
 *      - \mathcal{V} > 1: Sub-virial / Unbound / Dispersing supercluster
 *
 * References:
 * - Binney, J., & Tremaine, S. (2008). "Galactic Dynamics" (2nd ed.). Princeton University Press.
 * - Courtois, H. M., et al. (2023). "Cosmicflows-4: The Velocity Field and Cosmography." ApJ, 944, 94.
 * - Tully, R. B., et al. (2014). "The Laniakea supercluster of galaxies." Nature, 513(7516), 71-73.
 * - Peebles, P. J. E. (1980). "The Large-Scale Structure of the Universe." Princeton University Press.
 * - Monaghan, J. J. (1992). "Smoothed Particle Hydrodynamics." Annu. Rev. Astron. Astrophys., 30, 543.
 * - Navarro, J. F., Frenk, C. S., & White, S. D. (1996). "The Structure of Cold Dark Matter Halos." ApJ, 462, 563.
 *
 * @module fields/velocity_dispersion_tensor
 */

import { GridIndexer, BoundaryMode } from './grid_indexer.js';
import { ScalarField3D } from './scalar_field_3d.js';
import { VectorField3D } from './vector_field_3d.js';
import { VelocityField } from './velocity_field.js';
import { DensityField } from './density_field.js';
import {
  SupergalacticPosition,
  VelocityVector,
  EigenSystem3D,
  ScientificUnits,
  sha256Hex
} from '../coordinates/scientific_types.js';

// ============================================================================
// CONSTANTS & PHYSICAL PARAMETERS
// ============================================================================

/**
 * Universal Gravitational Constant G in standard astronomical/cosmological units:
 * (km/s)^2 * Mpc / M_sun
 * Exact conversion:
 * G = 6.67430e-11 m^3 kg^-1 s^-2
 * 1 Mpc = 3.085677581491367e22 m
 * 1 M_sun = 1.988409870698051e30 kg
 * 1 km/s = 1000 m/s
 * G_astro = G * M_sun / (Mpc * (1000 m/s)^2) = 4.3009172706e-9 (km/s)^2 Mpc / M_sun
 */
export const G_COSMO_MPC_MSUN = 4.3009172706e-9;

/**
 * Critical Density of the Universe at z=0 for H0 = 100 h km/s/Mpc:
 * rho_crit,0 = 3 H0^2 / (8 pi G) = 2.77536627e11 h^2 M_sun / Mpc^3
 * In (M_sun/h) / (Mpc/h)^3: rho_crit = 2.77536627e11
 */
export const RHO_CRIT_0_MSUN_MPC3 = 2.77536627e11;

/**
 * Default Planck/CF4 fiducial matter density parameter Omega_m.
 */
export const DEFAULT_OMEGA_M = 0.3111;

/**
 * Default Hubble constant in km/s/Mpc.
 */
export const DEFAULT_H0 = 74.6;

/**
 * Default dimensionless Hubble parameter h = H0 / 100.
 */
export const DEFAULT_LITTLE_H = 0.746;

/**
 * Default monoatomic adiabatic index gamma.
 */
export const DEFAULT_GAMMA_ADIABATIC = 5.0 / 3.0;

/**
 * Proton mass in kg.
 */
export const PROTON_MASS_KG = 1.67262192369e-27;

/**
 * Boltzmann constant in J/K.
 */
export const BOLTZMANN_CONSTANT_J_PER_K = 1.380649e-23;

/**
 * Mean molecular weight for fully ionized cosmological primordial plasma (X=0.76, Y=0.24).
 * mu = 1 / (2*X + 3/4*Y) = 1 / (1.52 + 0.18) = 0.588235
 */
export const MEAN_MOLECULAR_WEIGHT_PRIMORDIAL = 0.588235;

/**
 * Small numerical epsilon to prevent division by zero in singular radial regimes.
 */
export const NUMERICAL_EPSILON = 1e-15;

/**
 * Conversion factor from (km/s)^2 to Kelvin for gas temperature proxy:
 * T = (mu * m_p / k_B) * sigma_1D^2
 * with sigma_1D in km/s (so 10^6 m^2/s^2).
 * T_factor = (0.588235 * 1.67262192e-27 / 1.380649e-23) * 1.0e6 = 71.264 K / (km/s)^2
 */
export const VELOCITY_DISPERSION_TO_KELVIN_FACTOR = (MEAN_MOLECULAR_WEIGHT_PRIMORDIAL * PROTON_MASS_KG / BOLTZMANN_CONSTANT_J_PER_K) * 1.0e6;

/**
 * Supported spatial smoothing kernel types for dispersion estimation.
 * @readonly
 * @enum {string}
 */
export const SmoothingKernelType = Object.freeze({
  GAUSSIAN: 'gaussian',
  TOPHAT: 'tophat',
  CUBIC_SPLINE: 'cubic_spline',
  WENDLAND_C4: 'wendland_c4'
});

/**
 * Velocity orbital anisotropy regime classification.
 * @readonly
 * @enum {string}
 */
export const AnisotropyRegime = Object.freeze({
  PURE_RADIAL: 'pure_radial',             // beta = 1.0 (or beta > 0.8)
  RADIAL_BIASED: 'radial_biased',         // 0.2 < beta <= 0.8
  ISOTROPIC: 'isotropic',                 // -0.2 <= beta <= 0.2
  TANGENTIAL_BIASED: 'tangential_biased', // -2.0 <= beta < -0.2
  CIRCULAR_TANGENTIAL: 'circular_tangential' // beta < -2.0 (pure circular orbits)
});

/**
 * Cosmological dynamical virial equilibrium state classification.
 * @readonly
 * @enum {string}
 */
export const VirialState = Object.freeze({
  VIRIALIZED_EQUILIBRIUM: 'virialized_equilibrium',     // 0.85 <= 2K/|W| <= 1.15
  COLLAPSING_SUPER_VIRIAL: 'collapsing_super_virial',   // 2K/|W| < 0.85 (gravitational binding dominates)
  UNBOUND_SUB_VIRIAL: 'unbound_sub_virial',             // 2K/|W| > 1.15 (kinetic energy dominates / dispersing)
  HIGHLY_UNBOUND_EXPANDING: 'highly_unbound_expanding'   // 2K/|W| > 2.0
});

// ============================================================================
// 1. EXACT MATHEMATICAL UTILITIES & EIGENSOLVER
// ============================================================================

/**
 * Exact 3D Jacobi eigenvalue diagonalization for symmetric 3x3 matrices.
 * Iteratively applies Givens orthogonal similarity transformations: A' = P^T A P.
 * Guarantees quadratic convergence and machine-precision orthogonality.
 *
 * @param {Float64Array|Array<number>} mat 9-element row-major 3x3 symmetric matrix.
 * @param {number} [maxSweeps=50] Maximum number of Jacobi sweeps.
 * @param {number} [tolerance=1e-15] Convergence threshold on off-diagonal Frobenius norm.
 * @returns {{ eigenvalues: Float64Array, eigenvectors: Float64Array, sweepCount: number }}
 */
export function diagonalizeSymmetric3x3(mat, maxSweeps = 50, tolerance = 1e-15) {
  if (!mat || mat.length < 9) {
    throw new TypeError('diagonalizeSymmetric3x3: Input matrix must have at least 9 elements.');
  }

  // Symmetrize input matrix to avoid roundoff asymmetry
  const a = new Float64Array([
    mat[0], 0.5 * (mat[1] + mat[3]), 0.5 * (mat[2] + mat[6]),
    0.5 * (mat[1] + mat[3]), mat[4], 0.5 * (mat[5] + mat[7]),
    0.5 * (mat[2] + mat[6]), 0.5 * (mat[5] + mat[7]), mat[8]
  ]);

  // Initialize eigenvector matrix V as 3x3 Identity
  const v = new Float64Array([
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0
  ]);

  let sweep = 0;
  for (; sweep < maxSweeps; sweep++) {
    // Measure off-diagonal Frobenius norm
    const offDiagNorm = Math.sqrt(2.0 * (a[1] * a[1] + a[2] * a[2] + a[5] * a[5]));
    if (offDiagNorm < tolerance) {
      break;
    }

    // Iterate over 3 upper-triangular pairs: (0,1), (0,2), (1,2)
    const pairs = [[0, 1], [0, 2], [1, 2]];
    for (let p = 0; p < 3; p++) {
      const [i, j] = pairs[p];
      const a_ij = a[i * 3 + j];
      if (Math.abs(a_ij) < tolerance * 0.1) continue;

      const a_ii = a[i * 3 + i];
      const a_jj = a[j * 3 + j];
      const diff = a_jj - a_ii;

      let t;
      if (Math.abs(diff) < 1e-30) {
        t = a_ij >= 0 ? 1.0 : -1.0;
      } else {
        const phi = diff / (2.0 * a_ij);
        t = 1.0 / (Math.abs(phi) + Math.sqrt(phi * phi + 1.0));
        if (phi < 0.0) t = -t;
      }

      const c = 1.0 / Math.sqrt(t * t + 1.0);
      const s = t * c;
      const tau = s / (1.0 + c);

      // Update diagonal elements
      a[i * 3 + i] -= t * a_ij;
      a[j * 3 + j] += t * a_ij;
      a[i * 3 + j] = 0.0;
      a[j * 3 + i] = 0.0;

      // Update remaining matrix elements
      for (let k = 0; k < 3; k++) {
        if (k !== i && k !== j) {
          const a_ik = a[k * 3 + i];
          const a_jk = a[k * 3 + j];
          const new_ik = a_ik - s * (a_jk + tau * a_ik);
          const new_jk = a_jk + s * (a_ik - tau * a_jk);
          a[k * 3 + i] = new_ik;
          a[i * 3 + k] = new_ik;
          a[k * 3 + j] = new_jk;
          a[j * 3 + k] = new_jk;
        }
      }

      // Update eigenvectors
      for (let k = 0; k < 3; k++) {
        const v_ki = v[k * 3 + i];
        const v_kj = v[k * 3 + j];
        v[k * 3 + i] = v_ki - s * (v_kj + tau * v_ki);
        v[k * 3 + j] = v_kj + s * (v_ki - tau * v_kj);
      }
    }
  }

  // Extract eigenvalues
  const rawEvals = [
    { val: a[0], col: 0 },
    { val: a[4], col: 1 },
    { val: a[8], col: 2 }
  ];

  // Sort eigenvalues in descending order lambda_1 >= lambda_2 >= lambda_3
  rawEvals.sort((x, y) => y.val - x.val);

  const eigenvalues = new Float64Array([rawEvals[0].val, rawEvals[1].val, rawEvals[2].val]);
  const eigenvectors = new Float64Array(9);

  // Column k of eigenvectors corresponds to k-th sorted eigenvalue
  for (let k = 0; k < 3; k++) {
    const srcCol = rawEvals[k].col;
    let normSq = 0.0;
    for (let row = 0; row < 3; row++) {
      const comp = v[row * 3 + srcCol];
      eigenvectors[row * 3 + k] = comp;
      normSq += comp * comp;
    }
    const invNorm = normSq > 0 ? 1.0 / Math.sqrt(normSq) : 1.0;
    for (let row = 0; row < 3; row++) {
      eigenvectors[row * 3 + k] *= invNorm;
    }
  }

  // Ensure right-handed orthonormal coordinate frame: det(R) = +1
  const detV = eigenvectors[0] * (eigenvectors[4] * eigenvectors[8] - eigenvectors[5] * eigenvectors[7])
             - eigenvectors[1] * (eigenvectors[3] * eigenvectors[8] - eigenvectors[5] * eigenvectors[6])
             + eigenvectors[2] * (eigenvectors[3] * eigenvectors[7] - eigenvectors[4] * eigenvectors[6]);

  if (detV < 0.0) {
    // Flip third column to guarantee right-handed orientation
    eigenvectors[2] = -eigenvectors[2];
    eigenvectors[5] = -eigenvectors[5];
    eigenvectors[8] = -eigenvectors[8];
  }

  return { eigenvalues, eigenvectors, sweepCount: sweep };
}

/**
 * Analytical Cardano cubic eigensolver for a 3x3 real symmetric matrix.
 * Used for ultra-fast validation and direct eigenvalue checks.
 *
 * @param {Float64Array|Array<number>} mat 9-element 3x3 symmetric matrix.
 * @returns {Float64Array} Eigenvalues sorted in descending order [\lambda_1, \lambda_2, \lambda_3].
 */
export function cardanoEigenvaluesSymmetric3x3(mat) {
  if (!mat || mat.length < 9) {
    throw new TypeError('cardanoEigenvaluesSymmetric3x3: Matrix must have at least 9 elements.');
  }

  const a11 = mat[0], a12 = 0.5 * (mat[1] + mat[3]), a13 = 0.5 * (mat[2] + mat[6]);
  const a22 = mat[4], a23 = 0.5 * (mat[5] + mat[7]);
  const a33 = mat[8];

  // Invariants of 3x3 matrix:
  // p1 = a12^2 + a13^2 + a23^2
  const p1 = a12 * a12 + a13 * a13 + a23 * a23;
  const tr = a11 + a22 + a33;
  const q = tr / 3.0;

  if (p1 < 1e-30) {
    // Diagonal matrix
    const vals = [a11, a22, a33].sort((a, b) => b - a);
    return new Float64Array(vals);
  }

  const b11 = a11 - q, b22 = a22 - q, b33 = a33 - q;
  const p2 = b11 * b11 + b22 * b22 + b33 * b33 + 2.0 * p1;
  const p = Math.sqrt(p2 / 6.0);

  // det(B / p)
  const invP = 1.0 / p;
  const c11 = b11 * invP, c22 = b22 * invP, c33 = b33 * invP;
  const c12 = a12 * invP, c13 = a13 * invP, c23 = a23 * invP;

  const detB_over_p3 = c11 * (c22 * c33 - c23 * c23)
                     - c12 * (c12 * c33 - c23 * c13)
                     + c13 * (c12 * c23 - c22 * c13);

  const r = 0.5 * detB_over_p3;
  const phi = Math.acos(Math.max(-1.0, Math.min(1.0, r))) / 3.0;

  const eig1 = q + 2.0 * p * Math.cos(phi);
  const eig3 = q + 2.0 * p * Math.cos(phi + (2.0 * Math.PI / 3.0));
  const eig2 = 3.0 * q - eig1 - eig3; // tr = eig1 + eig2 + eig3

  const res = [eig1, eig2, eig3].sort((a, b) => b - a);
  return new Float64Array(res);
}

/**
 * Evaluates smoothing kernel weight W(r, R).
 *
 * @param {number} r Distance from kernel center.
 * @param {number} R Smoothing scale radius.
 * @param {string} [kernelType=SmoothingKernelType.GAUSSIAN] Kernel type.
 * @returns {number} Normalized kernel weight.
 */
export function evaluateSmoothingKernel(r, R, kernelType = SmoothingKernelType.GAUSSIAN) {
  if (R <= 0) {
    throw new RangeError('evaluateSmoothingKernel: Smoothing scale R must be strictly positive.');
  }

  const q = r / R;

  switch (kernelType) {
    case SmoothingKernelType.GAUSSIAN: {
      // 3D Normalized Gaussian: W(r, R) = 1 / ((2 pi)^(3/2) R^3) * exp(-r^2 / (2 R^2))
      const norm = 1.0 / (Math.pow(2.0 * Math.PI, 1.5) * R * R * R);
      return norm * Math.exp(-0.5 * q * q);
    }
    case SmoothingKernelType.TOPHAT: {
      // 3D Normalized Top-Hat: W(r, R) = 3 / (4 pi R^3) for r <= R, 0 otherwise
      if (r <= R) {
        return 3.0 / (4.0 * Math.PI * R * R * R);
      }
      return 0.0;
    }
    case SmoothingKernelType.CUBIC_SPLINE: {
      // Standard Monaghan (1992) SPH Cubic Spline kernel normalized in 3D
      // Support extends to 2*R
      const norm = 1.0 / (Math.PI * R * R * R);
      if (q <= 1.0) {
        return norm * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
      } else if (q <= 2.0) {
        const term = 2.0 - q;
        return norm * 0.25 * term * term * term;
      }
      return 0.0;
    }
    case SmoothingKernelType.WENDLAND_C4: {
      // Wendland C^4 kernel normalized in 3D, compact support [0, 2*R]
      const norm = 21.0 / (2.0 * Math.PI * Math.pow(2.0 * R, 3));
      if (q <= 2.0) {
        const u = 1.0 - 0.5 * q;
        const u4 = u * u * u * u;
        return norm * u4 * (1.0 + 2.0 * q);
      }
      return 0.0;
    }
    default:
      throw new TypeError(`evaluateSmoothingKernel: Unknown kernel type '${kernelType}'.`);
  }
}

// ============================================================================
// 2. DISPERSION TENSOR POINT & EIGENSYSTEM
// ============================================================================

/**
 * Representation of a single local 3D velocity dispersion tensor point \sigma_{ij}^2(\mathbf{x}).
 */
export class DispersionTensorPoint {
  /**
   * Constructs a DispersionTensorPoint.
   *
   * @param {number} sxx \sigma_{xx}^2 component (km/s)^2.
   * @param {number} sxy \sigma_{xy}^2 component (km/s)^2.
   * @param {number} sxz \sigma_{xz}^2 component (km/s)^2.
   * @param {number} syy \sigma_{yy}^2 component (km/s)^2.
   * @param {number} syz \sigma_{yz}^2 component (km/s)^2.
   * @param {number} szz \sigma_{zz}^2 component (km/s)^2.
   * @param {Array<number>} [position=[0,0,0]] Spatial position coordinates [x, y, z] in Mpc/h.
   */
  constructor(sxx, sxy, sxz, syy, syz, szz, position = [0.0, 0.0, 0.0]) {
    this.sxx = Number(sxx);
    this.sxy = Number(sxy);
    this.sxz = Number(sxz);
    this.syy = Number(syy);
    this.syz = Number(syz);
    this.szz = Number(szz);
    this.position = new Float64Array(position);

    this._eigenCached = null;
  }

  /**
   * Returns the 3x3 matrix in row-major Float64Array.
   * @returns {Float64Array}
   */
  toMatrix3x3() {
    return new Float64Array([
      this.sxx, this.sxy, this.sxz,
      this.sxy, this.syy, this.syz,
      this.sxz, this.syz, this.szz
    ]);
  }

  /**
   * Trace of the velocity dispersion tensor: Tr(\sigma^2) = \sigma_{xx}^2 + \sigma_{yy}^2 + \sigma_{zz}^2.
   * Units: (km/s)^2.
   * @returns {number}
   */
  get trace() {
    return this.sxx + this.syy + this.szz;
  }

  /**
   * 1D scalar velocity dispersion: \sigma_{1D} = \sqrt{\frac{1}{3} \text{Tr}(\sigma^2)}.
   * Units: km/s.
   * @returns {number}
   */
  get sigma1D() {
    const tr = Math.max(0.0, this.trace);
    return Math.sqrt(tr / 3.0);
  }

  /**
   * 3D scalar velocity dispersion: \sigma_{3D} = \sqrt{\text{Tr}(\sigma^2)}.
   * Units: km/s.
   * @returns {number}
   */
  get sigma3D() {
    const tr = Math.max(0.0, this.trace);
    return Math.sqrt(tr);
  }

  /**
   * Determinant of the 3x3 velocity dispersion tensor.
   * det(\sigma^2) = sxx(syy*szz - syz^2) - sxy(sxy*szz - syz*sxz) + sxz(sxy*syz - syy*sxz).
   * @returns {number}
   */
  get determinant() {
    return this.sxx * (this.syy * this.szz - this.syz * this.syz)
         - this.sxy * (this.sxy * this.szz - this.syz * this.sxz)
         + this.sxz * (this.sxy * this.syz - this.syy * this.sxz);
  }

  /**
   * Second Principal Invariant I_2 = 0.5 * ((Tr \sigma^2)^2 - Tr((\sigma^2)^2)).
   * @returns {number}
   */
  get secondInvariant() {
    const tr1 = this.trace;
    const tr2 = (this.sxx * this.sxx + this.syy * this.syy + this.szz * this.szz)
              + 2.0 * (this.sxy * this.sxy + this.sxz * this.sxz + this.syz * this.syz);
    return 0.5 * (tr1 * tr1 - tr2);
  }

  /**
   * Eigensystem decomposition with eigenvalues sorted in descending order: \lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0.
   * @returns {{ eigenvalues: Float64Array, eigenvectors: Float64Array }}
   */
  get eigensystem() {
    if (!this._eigenCached) {
      const mat = this.toMatrix3x3();
      const res = diagonalizeSymmetric3x3(mat);
      // Guard against tiny negative eigenvalues due to machine precision roundoff
      const l1 = Math.max(0.0, res.eigenvalues[0]);
      const l2 = Math.max(0.0, res.eigenvalues[1]);
      const l3 = Math.max(0.0, res.eigenvalues[2]);
      this._eigenCached = {
        eigenvalues: new Float64Array([l1, l2, l3]),
        eigenvectors: res.eigenvectors
      };
    }
    return this._eigenCached;
  }

  /**
   * Triaxiality parameter T = (\lambda_1 - \lambda_2) / (\lambda_1 - \lambda_3).
   * T \in [0, 1]: 0 = oblate (\lambda_1 = \lambda_2 > \lambda_3), 1 = prolate (\lambda_1 > \lambda_2 = \lambda_3).
   * @returns {number}
   */
  get triaxiality() {
    const ev = this.eigensystem.eigenvalues;
    const denom = ev[0] - ev[2];
    if (denom < NUMERICAL_EPSILON) {
      return 0.5; // Perfectly spherical dispersion
    }
    return Math.max(0.0, Math.min(1.0, (ev[0] - ev[1]) / denom));
  }

  /**
   * Ellipticity e = (\lambda_1 - \lambda_3) / (2 * (\lambda_1 + \lambda_2 + \lambda_3)).
   * @returns {number}
   */
  get ellipticity() {
    const ev = this.eigensystem.eigenvalues;
    const sum = ev[0] + ev[1] + ev[2];
    if (sum < NUMERICAL_EPSILON) return 0.0;
    return (ev[0] - ev[2]) / (2.0 * sum);
  }

  /**
   * Prolateness p = (\lambda_1 - 2*\lambda_2 + \lambda_3) / (2 * (\lambda_1 + \lambda_2 + \lambda_3)).
   * @returns {number}
   */
  get prolateness() {
    const ev = this.eigensystem.eigenvalues;
    const sum = ev[0] + ev[1] + ev[2];
    if (sum < NUMERICAL_EPSILON) return 0.0;
    return (ev[0] - 2.0 * ev[1] + ev[2]) / (2.0 * sum);
  }

  /**
   * Condition number \kappa = \lambda_1 / \max(\lambda_3, \epsilon).
   * @returns {number}
   */
  get conditionNumber() {
    const ev = this.eigensystem.eigenvalues;
    return ev[0] / Math.max(NUMERICAL_EPSILON, ev[2]);
  }

  /**
   * Checks if the dispersion tensor is positive semi-definite.
   * @returns {boolean}
   */
  isPositiveSemiDefinite() {
    const ev = this.eigensystem.eigenvalues;
    return ev[0] >= -1e-12 && ev[1] >= -1e-12 && ev[2] >= -1e-12;
  }
}

// ============================================================================
// 3. 3D VELOCITY DISPERSION TENSOR FIELD CONTAINER
// ============================================================================

/**
 * VelocityDispersionTensor3D encapsulates a continuous or discrete 3D cosmological grid
 * containing the 6 independent components of the symmetric velocity dispersion tensor:
 * \sigma_{xx}^2, \sigma_{xy}^2, \sigma_{xz}^2, \sigma_{yy}^2, \sigma_{yz}^2, \sigma_{zz}^2.
 */
export class VelocityDispersionTensor3D {
  /**
   * Constructs a VelocityDispersionTensor3D.
   *
   * @param {GridIndexer} gridIndexer Associated grid geometry and coordinate mappings.
   * @param {Float64Array|Float32Array} sxxBuffer Array of \sigma_{xx}^2 values.
   * @param {Float64Array|Float32Array} sxyBuffer Array of \sigma_{xy}^2 values.
   * @param {Float64Array|Float32Array} sxzBuffer Array of \sigma_{xz}^2 values.
   * @param {Float64Array|Float32Array} syyBuffer Array of \sigma_{yy}^2 values.
   * @param {Float64Array|Float32Array} syzBuffer Array of \sigma_{yz}^2 values.
   * @param {Float64Array|Float32Array} szzBuffer Array of \sigma_{zz}^2 values.
   * @param {Object} [options] Metadata and configuration options.
   * @param {number} [options.smoothingRadius=0.0] Smoothing filter scale in Mpc/h.
   * @param {string} [options.kernelType='gaussian'] Kernel type used for estimation.
   */
  constructor(gridIndexer, sxxBuffer, sxyBuffer, sxzBuffer, syyBuffer, syzBuffer, szzBuffer, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('VelocityDispersionTensor3D: gridIndexer must be an instance of GridIndexer.');
    }

    const total = gridIndexer.totalCells;
    const buffers = [sxxBuffer, sxyBuffer, sxzBuffer, syyBuffer, syzBuffer, szzBuffer];
    const names = ['sxx', 'sxy', 'sxz', 'syy', 'syz', 'szz'];

    for (let i = 0; i < 6; i++) {
      if (!buffers[i] || buffers[i].length !== total) {
        throw new Error(`VelocityDispersionTensor3D: Buffer '${names[i]}' size (${buffers[i]?.length}) must match grid total cells (${total}).`);
      }
    }

    this.grid = gridIndexer;
    this.smoothingRadius = Number(options.smoothingRadius || 0.0);
    this.kernelType = options.kernelType || SmoothingKernelType.GAUSSIAN;

    // Allocate continuous Float64Array internal buffers
    this.sxx = new Float64Array(total);
    this.sxy = new Float64Array(total);
    this.sxz = new Float64Array(total);
    this.syy = new Float64Array(total);
    this.syz = new Float64Array(total);
    this.szz = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      this.sxx[i] = Math.max(0.0, sxxBuffer[i]);
      this.sxy[i] = sxyBuffer[i];
      this.sxz[i] = sxzBuffer[i];
      this.syy[i] = Math.max(0.0, syyBuffer[i]);
      this.syz[i] = syzBuffer[i];
      this.szz[i] = Math.max(0.0, szzBuffer[i]);
    }
  }

  /**
   * Factory method: Computes the 3D Velocity Dispersion Tensor field from a VelocityField
   * via separable spatial Gaussian convolution or local neighborhood kernel smoothing:
   * \langle v_i \rangle = W * v_i, \quad \langle v_i v_j \rangle = W * (v_i v_j)
   * \sigma_{ij}^2 = \langle v_i v_j \rangle - \langle v_i \rangle \langle v_j \rangle
   *
   * @param {VelocityField} velocityField Source peculiar velocity field (vx, vy, vz in km/s).
   * @param {number} smoothingRadius Smoothing kernel radius R in Mpc/h.
   * @param {Object} [options] Options for kernel weighting and boundary conditions.
   * @param {string} [options.kernelType=SmoothingKernelType.GAUSSIAN] Smoothing kernel.
   * @returns {VelocityDispersionTensor3D} Computed velocity dispersion tensor field.
   */
  static computeFromVelocityField(velocityField, smoothingRadius, options = {}) {
    if (!(velocityField instanceof VelocityField)) {
      throw new TypeError('VelocityDispersionTensor3D.computeFromVelocityField: velocityField must be an instance of VelocityField.');
    }
    if (typeof smoothingRadius !== 'number' || smoothingRadius <= 0.0) {
      throw new RangeError('VelocityDispersionTensor3D.computeFromVelocityField: smoothingRadius must be a positive number in Mpc/h.');
    }

    const grid = velocityField.grid;
    const { nx, ny, nz } = grid;
    const total = grid.totalCells;
    const dx = grid.dx, dy = grid.dy, dz = grid.dz;
    const kernelType = options.kernelType || SmoothingKernelType.GAUSSIAN;

    // Kernel truncation radius in grid voxels (3-sigma for Gaussian, 2-support for splines)
    const spanScale = kernelType === SmoothingKernelType.GAUSSIAN ? 3.0 : 2.0;
    const rx = Math.max(1, Math.ceil(spanScale * smoothingRadius / dx));
    const ry = Math.max(1, Math.ceil(spanScale * smoothingRadius / dy));
    const rz = Math.max(1, Math.ceil(spanScale * smoothingRadius / dz));

    // Allocate output buffers
    const sxxOut = new Float64Array(total);
    const sxyOut = new Float64Array(total);
    const sxzOut = new Float64Array(total);
    const syyOut = new Float64Array(total);
    const syzOut = new Float64Array(total);
    const szzOut = new Float64Array(total);

    // Perform exact local spatial kernel convolution
    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const centerIdx = grid.getLinearIndex(ix, iy, iz);

          let sumW = 0.0;
          let sumVx = 0.0, sumVy = 0.0, sumVz = 0.0;
          let sumVxx = 0.0, sumVxy = 0.0, sumVxz = 0.0;
          let sumVyy = 0.0, sumVyz = 0.0, sumVzz = 0.0;

          const minKz = Math.max(0, iz - rz);
          const maxKz = Math.min(nz - 1, iz + rz);
          const minKy = Math.max(0, iy - ry);
          const maxKy = Math.min(ny - 1, iy + ry);
          const minKx = Math.max(0, ix - rx);
          const maxKx = Math.min(nx - 1, ix + rx);

          for (let kz = minKz; kz <= maxKz; kz++) {
            const dzPhys = (kz - iz) * dz;
            for (let ky = minKy; ky <= maxKy; ky++) {
              const dyPhys = (ky - iy) * dy;
              for (let kx = minKx; kx <= maxKx; kx++) {
                const dxPhys = (kx - ix) * dx;
                const dist = Math.sqrt(dxPhys * dxPhys + dyPhys * dyPhys + dzPhys * dzPhys);

                const w = evaluateSmoothingKernel(dist, smoothingRadius, kernelType);
                if (w <= 0.0) continue;

                const nIdx = grid.getLinearIndex(kx, ky, kz);
                const vx_k = velocityField.vx[nIdx];
                const vy_k = velocityField.vy[nIdx];
                const vz_k = velocityField.vz[nIdx];

                sumW += w;
                sumVx += w * vx_k;
                sumVy += w * vy_k;
                sumVz += w * vz_k;

                sumVxx += w * vx_k * vx_k;
                sumVxy += w * vx_k * vy_k;
                sumVxz += w * vx_k * vz_k;
                sumVyy += w * vy_k * vy_k;
                sumVyz += w * vy_k * vz_k;
                sumVzz += w * vz_k * vz_k;
              }
            }
          }

          if (sumW > 0.0) {
            const invW = 1.0 / sumW;
            const meanVx = sumVx * invW;
            const meanVy = sumVy * invW;
            const meanVz = sumVz * invW;

            const varXx = Math.max(0.0, (sumVxx * invW) - (meanVx * meanVx));
            const covXy = (sumVxy * invW) - (meanVx * meanVy);
            const covXz = (sumVxz * invW) - (meanVx * meanVz);
            const varYy = Math.max(0.0, (sumVyy * invW) - (meanVy * meanVy));
            const covYz = (sumVyz * invW) - (meanVy * meanVz);
            const varZz = Math.max(0.0, (sumVzz * invW) - (meanVz * meanVz));

            sxxOut[centerIdx] = varXx;
            sxyOut[centerIdx] = covXy;
            sxzOut[centerIdx] = covXz;
            syyOut[centerIdx] = varYy;
            syzOut[centerIdx] = covYz;
            szzOut[centerIdx] = varZz;
          }
        }
      }
    }

    return new VelocityDispersionTensor3D(grid, sxxOut, sxyOut, sxzOut, syyOut, syzOut, szzOut, {
      smoothingRadius,
      kernelType
    });
  }

  /**
   * Factory method: Computes the 3D Velocity Dispersion Tensor field from a discrete particle/galaxy catalog.
   *
   * @param {Array<Array<number>>|Float64Array} positions Array of particle positions [x, y, z] in Mpc/h.
   * @param {Array<Array<number>>|Float64Array} velocities Array of particle velocities [vx, vy, vz] in km/s.
   * @param {GridIndexer} gridIndexer Target grid geometry.
   * @param {number} smoothingRadius Kernel smoothing scale R in Mpc/h.
   * @param {Object} [options] Particle weights and kernel options.
   * @param {Array<number>|Float64Array} [options.weights=null] Optional mass or observational weights.
   * @param {string} [options.kernelType=SmoothingKernelType.GAUSSIAN] Kernel type.
   * @returns {VelocityDispersionTensor3D}
   */
  static computeFromParticleCatalog(positions, velocities, gridIndexer, smoothingRadius, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('VelocityDispersionTensor3D.computeFromParticleCatalog: gridIndexer must be an instance of GridIndexer.');
    }
    if (!positions || !velocities || positions.length !== velocities.length) {
      throw new Error('VelocityDispersionTensor3D.computeFromParticleCatalog: Positions and velocities count must match.');
    }
    if (typeof smoothingRadius !== 'number' || smoothingRadius <= 0.0) {
      throw new RangeError('VelocityDispersionTensor3D.computeFromParticleCatalog: smoothingRadius must be strictly positive.');
    }

    const nPart = positions.length;
    const total = gridIndexer.totalCells;
    const { nx, ny, nz } = gridIndexer;
    const kernelType = options.kernelType || SmoothingKernelType.GAUSSIAN;
    const weights = options.weights || null;

    const sxxOut = new Float64Array(total);
    const sxyOut = new Float64Array(total);
    const sxzOut = new Float64Array(total);
    const syyOut = new Float64Array(total);
    const syzOut = new Float64Array(total);
    const szzOut = new Float64Array(total);

    // Compute at each grid voxel
    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const [gx, gy, gz] = gridIndexer.getNodeCoord(ix, iy, iz);
          const gridIdx = gridIndexer.getLinearIndex(ix, iy, iz);

          let sumW = 0.0;
          let sumVx = 0.0, sumVy = 0.0, sumVz = 0.0;
          let sumVxx = 0.0, sumVxy = 0.0, sumVxz = 0.0;
          let sumVyy = 0.0, sumVyz = 0.0, sumVzz = 0.0;

          for (let p = 0; p < nPart; p++) {
            const px = positions[p][0];
            const py = positions[p][1];
            const pz = positions[p][2];

            const dx = px - gx;
            const dy = py - gy;
            const dz = pz - gz;
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

            const kw = evaluateSmoothingKernel(dist, smoothingRadius, kernelType);
            if (kw <= 0.0) continue;

            const pw = weights ? weights[p] : 1.0;
            const w = kw * pw;

            const vx = velocities[p][0];
            const vy = velocities[p][1];
            const vz = velocities[p][2];

            sumW += w;
            sumVx += w * vx;
            sumVy += w * vy;
            sumVz += w * vz;

            sumVxx += w * vx * vx;
            sumVxy += w * vx * vy;
            sumVxz += w * vx * vz;
            sumVyy += w * vy * vy;
            sumVyz += w * vy * vz;
            sumVzz += w * vz * vz;
          }

          if (sumW > 0.0) {
            const invW = 1.0 / sumW;
            const meanVx = sumVx * invW;
            const meanVy = sumVy * invW;
            const meanVz = sumVz * invW;

            sxxOut[gridIdx] = Math.max(0.0, (sumVxx * invW) - (meanVx * meanVx));
            sxyOut[gridIdx] = (sumVxy * invW) - (meanVx * meanVy);
            sxzOut[gridIdx] = (sumVxz * invW) - (meanVx * meanVz);
            syyOut[gridIdx] = Math.max(0.0, (sumVyy * invW) - (meanVy * meanVy));
            syzOut[gridIdx] = (sumVyz * invW) - (meanVy * meanVz);
            szzOut[gridIdx] = Math.max(0.0, (sumVzz * invW) - (meanVz * meanVz));
          }
        }
      }
    }

    return new VelocityDispersionTensor3D(gridIndexer, sxxOut, sxyOut, sxzOut, syyOut, syzOut, szzOut, {
      smoothingRadius,
      kernelType
    });
  }

  /**
   * Factory method: Constructs an analytic velocity dispersion tensor field from a function.
   *
   * @param {GridIndexer} gridIndexer Grid indexer.
   * @param {Function} tensorFn Function (x, y, z) => [sxx, sxy, sxz, syy, syz, szz].
   * @returns {VelocityDispersionTensor3D}
   */
  static fromAnalyticFunction(gridIndexer, tensorFn) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('VelocityDispersionTensor3D.fromAnalyticFunction: gridIndexer must be an instance of GridIndexer.');
    }
    if (typeof tensorFn !== 'function') {
      throw new TypeError('VelocityDispersionTensor3D.fromAnalyticFunction: tensorFn must be a function.');
    }

    const total = gridIndexer.totalCells;
    const sxx = new Float64Array(total);
    const sxy = new Float64Array(total);
    const sxz = new Float64Array(total);
    const syy = new Float64Array(total);
    const syz = new Float64Array(total);
    const szz = new Float64Array(total);

    for (let iz = 0; iz < gridIndexer.nz; iz++) {
      for (let iy = 0; iy < gridIndexer.ny; iy++) {
        for (let ix = 0; ix < gridIndexer.nx; ix++) {
          const idx = gridIndexer.getLinearIndex(ix, iy, iz);
          const [x, y, z] = gridIndexer.getNodeCoord(ix, iy, iz);
          const [xx, xy, xz, yy, yz, zz] = tensorFn(x, y, z);
          sxx[idx] = Math.max(0.0, xx);
          sxy[idx] = xy;
          sxz[idx] = xz;
          syy[idx] = Math.max(0.0, yy);
          syz[idx] = yz;
          szz[idx] = Math.max(0.0, zz);
        }
      }
    }

    return new VelocityDispersionTensor3D(gridIndexer, sxx, sxy, sxz, syy, syz, szz);
  }

  /**
   * Gets the local dispersion tensor point at discrete integer grid indices (ix, iy, iz).
   *
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {DispersionTensorPoint}
   */
  getTensorAt(ix, iy, iz) {
    const idx = this.grid.getLinearIndex(ix, iy, iz);
    const pos = this.grid.getNodeCoord(ix, iy, iz);
    return new DispersionTensorPoint(
      this.sxx[idx],
      this.sxy[idx],
      this.sxz[idx],
      this.syy[idx],
      this.syz[idx],
      this.szz[idx],
      pos
    );
  }

  /**
   * Continuous trilinear interpolation of the 6 dispersion tensor components at any physical (x, y, z) in Mpc/h.
   *
   * @param {number} x Supergalactic X (Mpc/h).
   * @param {number} y Supergalactic Y (Mpc/h).
   * @param {number} z Supergalactic Z (Mpc/h).
   * @returns {DispersionTensorPoint}
   */
  evaluateAt(x, y, z) {
    const [u, v, w] = this.grid.coordToGridIndex(x, y, z);
    const { nx, ny, nz } = this.grid;

    // Floor indices for the 8 surrounding voxels
    const i0 = Math.max(0, Math.min(nx - 2, Math.floor(u)));
    const j0 = Math.max(0, Math.min(ny - 2, Math.floor(v)));
    const k0 = Math.max(0, Math.min(nz - 2, Math.floor(w)));

    const i1 = i0 + 1;
    const j1 = j0 + 1;
    const k1 = k0 + 1;

    // Trilinear fractional weights
    const fu = Math.max(0.0, Math.min(1.0, u - i0));
    const fv = Math.max(0.0, Math.min(1.0, v - j0));
    const fw = Math.max(0.0, Math.min(1.0, w - k0));

    const w000 = (1.0 - fu) * (1.0 - fv) * (1.0 - fw);
    const w100 = fu * (1.0 - fv) * (1.0 - fw);
    const w010 = (1.0 - fu) * fv * (1.0 - fw);
    const w110 = fu * fv * (1.0 - fw);
    const w001 = (1.0 - fu) * (1.0 - fv) * fw;
    const w101 = fu * (1.0 - fv) * fw;
    const w011 = (1.0 - fu) * fv * fw;
    const w111 = fu * fv * fw;

    const idx000 = this.grid.getLinearIndex(i0, j0, k0);
    const idx100 = this.grid.getLinearIndex(i1, j0, k0);
    const idx010 = this.grid.getLinearIndex(i0, j1, k0);
    const idx110 = this.grid.getLinearIndex(i1, j1, k0);
    const idx001 = this.grid.getLinearIndex(i0, j0, k1);
    const idx101 = this.grid.getLinearIndex(i1, j0, k1);
    const idx011 = this.grid.getLinearIndex(i0, j1, k1);
    const idx111 = this.grid.getLinearIndex(i1, j1, k1);

    const interp = (buf) => {
      return w000 * buf[idx000] + w100 * buf[idx100] +
             w010 * buf[idx010] + w110 * buf[idx110] +
             w001 * buf[idx001] + w101 * buf[idx101] +
             w011 * buf[idx011] + w111 * buf[idx111];
    };

    return new DispersionTensorPoint(
      Math.max(0.0, interp(this.sxx)),
      interp(this.sxy),
      interp(this.sxz),
      Math.max(0.0, interp(this.syy)),
      interp(this.syz),
      Math.max(0.0, interp(this.szz)),
      [x, y, z]
    );
  }

  /**
   * Generates a ScalarField3D containing the 1D scalar velocity dispersion \sigma_{1D}(\mathbf{x}) = \sqrt{\frac{1}{3} \text{Tr}(\sigma^2)}.
   * Units: km/s.
   * @returns {ScalarField3D}
   */
  getScalarSigma1DField() {
    const total = this.grid.totalCells;
    const data = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      const tr = Math.max(0.0, this.sxx[i] + this.syy[i] + this.szz[i]);
      data[i] = Math.sqrt(tr / 3.0);
    }
    return new ScalarField3D(this.grid, data, 'sigma_1d', ScientificUnits.KM_PER_S);
  }

  /**
   * Generates a ScalarField3D containing the 3D scalar velocity dispersion \sigma_{3D}(\mathbf{x}) = \sqrt{\text{Tr}(\sigma^2)}.
   * Units: km/s.
   * @returns {ScalarField3D}
   */
  getScalarSigma3DField() {
    const total = this.grid.totalCells;
    const data = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      const tr = Math.max(0.0, this.sxx[i] + this.syy[i] + this.szz[i]);
      data[i] = Math.sqrt(tr);
    }
    return new ScalarField3D(this.grid, data, 'sigma_3d', ScientificUnits.KM_PER_S);
  }

  /**
   * Generates a ScalarField3D containing the Triaxiality parameter T(\mathbf{x}) \in [0, 1].
   * @returns {ScalarField3D}
   */
  getTriaxialityField() {
    const total = this.grid.totalCells;
    const data = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      const pt = new DispersionTensorPoint(
        this.sxx[i], this.sxy[i], this.sxz[i],
        this.syy[i], this.syz[i], this.szz[i]
      );
      data[i] = pt.triaxiality;
    }
    return new ScalarField3D(this.grid, data, 'triaxiality', ScientificUnits.DIMENSIONLESS);
  }

  /**
   * Generates a ScalarField3D containing the Ellipticity parameter e(\mathbf{x}).
   * @returns {ScalarField3D}
   */
  getEllipticityField() {
    const total = this.grid.totalCells;
    const data = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      const pt = new DispersionTensorPoint(
        this.sxx[i], this.sxy[i], this.sxz[i],
        this.syy[i], this.syz[i], this.szz[i]
      );
      data[i] = pt.ellipticity;
    }
    return new ScalarField3D(this.grid, data, 'ellipticity', ScientificUnits.DIMENSIONLESS);
  }

  /**
   * Returns a VectorField3D of the primary principal dispersion axis (eigenvector \mathbf{e}_1).
   * @returns {VectorField3D}
   */
  getPrincipalAxis1Field() {
    const total = this.grid.totalCells;
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const pt = new DispersionTensorPoint(
        this.sxx[i], this.sxy[i], this.sxz[i],
        this.syy[i], this.syz[i], this.szz[i]
      );
      const evec = pt.eigensystem.eigenvectors;
      // Column 0 is first principal axis
      vx[i] = evec[0];
      vy[i] = evec[3];
      vz[i] = evec[6];
    }

    return new VectorField3D(this.grid, vx, vy, vz, 'principal_axis_1', ScientificUnits.DIMENSIONLESS);
  }
}

// ============================================================================
// 4. SPHERICAL HALO ANISOTROPY ANALYZER (BETA PARAMETER)
// ============================================================================

/**
 * Result container for spherical velocity dispersion transformation and anisotropy evaluation.
 */
export class SphericalAnisotropyPoint {
  /**
   * Constructs a SphericalAnisotropyPoint.
   *
   * @param {number} r Distance from center (Mpc/h).
   * @param {number} sigmaR2 Radial velocity dispersion squared \sigma_r^2 (km/s)^2.
   * @param {number} sigmaTheta2 Polar/theta velocity dispersion squared \sigma_\theta^2 (km/s)^2.
   * @param {number} sigmaPhi2 Azimuthal/phi velocity dispersion squared \sigma_\phi^2 (km/s)^2.
   * @param {Array<number>} [position=[0,0,0]] Physical coordinates.
   */
  constructor(r, sigmaR2, sigmaTheta2, sigmaPhi2, position = [0.0, 0.0, 0.0]) {
    this.r = Number(r);
    this.sigmaR2 = Math.max(0.0, Number(sigmaR2));
    this.sigmaTheta2 = Math.max(0.0, Number(sigmaTheta2));
    this.sigmaPhi2 = Math.max(0.0, Number(sigmaPhi2));
    this.position = new Float64Array(position);
  }

  /**
   * Combined tangential velocity dispersion squared: \sigma_t^2 = \frac{\sigma_\theta^2 + \sigma_\phi^2}{2}.
   * @returns {number}
   */
  get sigmaTan2() {
    return 0.5 * (this.sigmaTheta2 + this.sigmaPhi2);
  }

  /**
   * Total tangential dispersion squared: \sigma_{\text{tan,tot}}^2 = \sigma_\theta^2 + \sigma_\phi^2.
   * @returns {number}
   */
  get sigmaTanTot2() {
    return this.sigmaTheta2 + this.sigmaPhi2;
  }

  /**
   * 1D Radial velocity dispersion: \sigma_r (km/s).
   * @returns {number}
   */
  get sigmaR() {
    return Math.sqrt(this.sigmaR2);
  }

  /**
   * 1D Polar velocity dispersion: \sigma_\theta (km/s).
   * @returns {number}
   */
  get sigmaTheta() {
    return Math.sqrt(this.sigmaTheta2);
  }

  /**
   * 1D Azimuthal velocity dispersion: \sigma_\phi (km/s).
   * @returns {number}
   */
  get sigmaPhi() {
    return Math.sqrt(this.sigmaPhi2);
  }

  /**
   * 1D Tangential velocity dispersion: \sigma_t (km/s).
   * @returns {number}
   */
  get sigmaTan() {
    return Math.sqrt(this.sigmaTan2);
  }

  /**
   * Binney (1980) Velocity Anisotropy Parameter:
   * \beta = 1 - \frac{\sigma_\theta^2 + \sigma_\phi^2}{2 \sigma_r^2} = 1 - \frac{\sigma_t^2}{\sigma_r^2}
   * @returns {number}
   */
  get beta() {
    if (this.sigmaR2 < NUMERICAL_EPSILON) {
      if (this.sigmaTan2 < NUMERICAL_EPSILON) return 0.0; // Both zero -> isotropic zero dispersion
      return -Infinity; // Pure circular tangential motion
    }
    return 1.0 - (this.sigmaTan2 / this.sigmaR2);
  }

  /**
   * Symmetrized anisotropy parameter:
   * \gamma_\beta = \frac{2\beta}{2 - \beta} = \frac{\sigma_r^2 - \sigma_t^2}{\sigma_r^2 + \sigma_t^2} \in [-1, 1]
   * @returns {number}
   */
  get symmetrizedBeta() {
    const sum = this.sigmaR2 + this.sigmaTan2;
    if (sum < NUMERICAL_EPSILON) return 0.0;
    return (this.sigmaR2 - this.sigmaTan2) / sum;
  }

  /**
   * Classifies orbital anisotropy regime into standard dynamical categories.
   * @returns {string}
   */
  classifyRegime() {
    const b = this.beta;
    if (b > 0.8) return AnisotropyRegime.PURE_RADIAL;
    if (b > 0.2) return AnisotropyRegime.RADIAL_BIASED;
    if (b >= -0.2) return AnisotropyRegime.ISOTROPIC;
    if (b >= -2.0) return AnisotropyRegime.TANGENTIAL_BIASED;
    return AnisotropyRegime.CIRCULAR_TANGENTIAL;
  }
}

/**
 * SphericalAnisotropyAnalyzer performs coordinate transformations from Cartesian
 * velocity dispersion tensors to spherical halo-centric frames, computing \beta(r) profiles.
 */
export class SphericalAnisotropyAnalyzer {
  /**
   * Transforms a Cartesian velocity dispersion tensor \mathbf{\Sigma}_{\text{cart}} into spherical
   * components (\sigma_r^2, \sigma_\theta^2, \sigma_\phi^2) at position \mathbf{x} relative to halo center \mathbf{x}_c.
   *
   * @param {DispersionTensorPoint} tensorPoint Local Cartesian dispersion tensor point.
   * @param {Array<number>} haloCenter Center coordinates [xc, yc, zc] in Mpc/h.
   * @returns {SphericalAnisotropyPoint}
   */
  static transformToSpherical(tensorPoint, haloCenter) {
    if (!(tensorPoint instanceof DispersionTensorPoint)) {
      throw new TypeError('SphericalAnisotropyAnalyzer.transformToSpherical: tensorPoint must be a DispersionTensorPoint.');
    }
    if (!haloCenter || haloCenter.length < 3) {
      throw new TypeError('SphericalAnisotropyAnalyzer.transformToSpherical: haloCenter must have 3 coordinates.');
    }

    const rx = tensorPoint.position[0] - haloCenter[0];
    const ry = tensorPoint.position[1] - haloCenter[1];
    const rz = tensorPoint.position[2] - haloCenter[2];

    const r2 = rx * rx + ry * ry + rz * rz;
    const r = Math.sqrt(r2);
    const rxy2 = rx * rx + ry * ry;
    const rxy = Math.sqrt(rxy2);

    // Center point singularity handling
    if (r < NUMERICAL_EPSILON) {
      const avgDiag = tensorPoint.trace / 3.0;
      return new SphericalAnisotropyPoint(0.0, avgDiag, avgDiag, avgDiag, tensorPoint.position);
    }

    // Spherical basis unit vectors
    // 1. Radial unit vector e_r = (rx, ry, rz) / r
    const er_x = rx / r;
    const er_y = ry / r;
    const er_z = rz / r;

    // 2. Polar & azimuthal unit vectors with pole singularity protection
    let etheta_x, etheta_y, etheta_z;
    let ephi_x, ephi_y, ephi_z;

    if (rxy < NUMERICAL_EPSILON) {
      // Exactly on the z-axis (north or south pole)
      etheta_x = 1.0; etheta_y = 0.0; etheta_z = 0.0;
      ephi_x = 0.0; ephi_y = 1.0; ephi_z = 0.0;
    } else {
      // e_theta = (cos(theta) cos(phi), cos(theta) sin(phi), -sin(theta))
      // cos(theta) = rz/r, sin(theta) = rxy/r, cos(phi) = rx/rxy, sin(phi) = ry/rxy
      etheta_x = (rz * rx) / (r * rxy);
      etheta_y = (rz * ry) / (r * rxy);
      etheta_z = -rxy / r;

      // e_phi = (-sin(phi), cos(phi), 0)
      ephi_x = -ry / rxy;
      ephi_y = rx / rxy;
      ephi_z = 0.0;
    }

    const { sxx, sxy, sxz, syy, syz, szz } = tensorPoint;

    // Quadratic form for radial dispersion: \sigma_r^2 = \hat{\mathbf{e}}_r^T \mathbf{\Sigma} \hat{\mathbf{e}}_r
    const sigmaR2 = er_x * (sxx * er_x + sxy * er_y + sxz * er_z) +
                    er_y * (sxy * er_x + syy * er_y + syz * er_z) +
                    er_z * (sxz * er_x + syz * er_y + szz * er_z);

    // Quadratic form for theta dispersion: \sigma_\theta^2 = \hat{\mathbf{e}}_\theta^T \mathbf{\Sigma} \hat{\mathbf{e}}_\theta
    const sigmaTheta2 = etheta_x * (sxx * etheta_x + sxy * etheta_y + sxz * etheta_z) +
                        etheta_y * (sxy * etheta_x + syy * etheta_y + syz * etheta_z) +
                        etheta_z * (sxz * etheta_x + syz * etheta_y + szz * etheta_z);

    // Quadratic form for phi dispersion: \sigma_\phi^2 = \hat{\mathbf{e}}_\phi^T \mathbf{\Sigma} \hat{\mathbf{e}}_\phi
    const sigmaPhi2 = ephi_x * (sxx * ephi_x + sxy * ephi_y + sxz * ephi_z) +
                      ephi_y * (sxy * ephi_x + syy * ephi_y + syz * ephi_z) +
                      ephi_z * (sxz * ephi_x + syz * ephi_y + szz * ephi_z);

    return new SphericalAnisotropyPoint(
      r,
      Math.max(0.0, sigmaR2),
      Math.max(0.0, sigmaTheta2),
      Math.max(0.0, sigmaPhi2),
      tensorPoint.position
    );
  }

  /**
   * Evaluates theoretical Osipkov-Merritt anisotropy profile:
   * \beta(r) = \frac{r^2}{r^2 + r_a^2}
   *
   * @param {number} r Radial distance in Mpc/h.
   * @param {number} anisotropyRadius Anisotropy radius r_a in Mpc/h.
   * @returns {number} Anisotropy parameter \beta(r) \in [0, 1].
   */
  static osipkovMerrittBeta(r, anisotropyRadius) {
    if (r < 0.0 || anisotropyRadius <= 0.0) {
      throw new RangeError('osipkovMerrittBeta: r >= 0 and anisotropyRadius > 0 required.');
    }
    const r2 = r * r;
    const ra2 = anisotropyRadius * anisotropyRadius;
    return r2 / (r2 + ra2);
  }

  /**
   * Computes radial profile of spherical velocity dispersion components and anisotropy parameter \beta(r)
   * in concentric spherical shells around a halo/cluster center.
   *
   * @param {VelocityDispersionTensor3D} tensorField 3D dispersion tensor field.
   * @param {Array<number>} haloCenter Center coordinates [xc, yc, zc] in Mpc/h.
   * @param {Object} [options] Binning options.
   * @param {number} [options.rMin=0.5] Minimum radial distance in Mpc/h.
   * @param {number} [options.rMax=20.0] Maximum radial distance in Mpc/h.
   * @param {number} [options.numBins=30] Number of radial shell bins.
   * @param {boolean} [options.logSpacing=false] Whether bins are logarithmically spaced.
   * @returns {Array<{ r: number, rInner: number, rOuter: number, count: number, sigmaR: number, sigmaTheta: number, sigmaPhi: number, sigmaTan: number, sigmaR2: number, sigmaTan2: number, beta: number, symmetrizedBeta: number, regime: string }>}
   */
  static computeRadialAnisotropyProfile(tensorField, haloCenter, options = {}) {
    if (!(tensorField instanceof VelocityDispersionTensor3D)) {
      throw new TypeError('SphericalAnisotropyAnalyzer.computeRadialAnisotropyProfile: tensorField must be a VelocityDispersionTensor3D.');
    }

    const rMin = Number(options.rMin ?? 0.5);
    const rMax = Number(options.rMax ?? 20.0);
    const numBins = Math.max(2, parseInt(options.numBins ?? 30, 10));
    const logSpacing = Boolean(options.logSpacing);

    if (rMin <= 0.0 || rMax <= rMin) {
      throw new RangeError('computeRadialAnisotropyProfile: Invalid radial range [rMin, rMax].');
    }

    // Setup bin edges
    const binEdges = new Float64Array(numBins + 1);
    if (logSpacing) {
      const logMin = Math.log(rMin);
      const logMax = Math.log(rMax);
      const dLog = (logMax - logMin) / numBins;
      for (let b = 0; b <= numBins; b++) {
        binEdges[b] = Math.exp(logMin + b * dLog);
      }
    } else {
      const dr = (rMax - rMin) / numBins;
      for (let b = 0; b <= numBins; b++) {
        binEdges[b] = rMin + b * dr;
      }
    }

    // Accumulators for each bin
    const binCounts = new Int32Array(numBins);
    const sumR = new Float64Array(numBins);
    const sumSigmaR2 = new Float64Array(numBins);
    const sumSigmaTheta2 = new Float64Array(numBins);
    const sumSigmaPhi2 = new Float64Array(numBins);

    const grid = tensorField.grid;
    const { nx, ny, nz } = grid;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const pt = tensorField.getTensorAt(ix, iy, iz);
          const sph = SphericalAnisotropyAnalyzer.transformToSpherical(pt, haloCenter);

          const r = sph.r;
          if (r < binEdges[0] || r >= binEdges[numBins]) continue;

          // Find bin index via binary search
          let low = 0, high = numBins - 1, binIdx = -1;
          while (low <= high) {
            const mid = (low + high) >> 1;
            if (r >= binEdges[mid] && r < binEdges[mid + 1]) {
              binIdx = mid;
              break;
            } else if (r < binEdges[mid]) {
              high = mid - 1;
            } else {
              low = mid + 1;
            }
          }

          if (binIdx >= 0) {
            binCounts[binIdx]++;
            sumR[binIdx] += r;
            sumSigmaR2[binIdx] += sph.sigmaR2;
            sumSigmaTheta2[binIdx] += sph.sigmaTheta2;
            sumSigmaPhi2[binIdx] += sph.sigmaPhi2;
          }
        }
      }
    }

    // Compile profile results
    const profile = [];
    for (let b = 0; b < numBins; b++) {
      const count = binCounts[b];
      const rInner = binEdges[b];
      const rOuter = binEdges[b + 1];
      const rMid = count > 0 ? sumR[b] / count : 0.5 * (rInner + rOuter);

      if (count === 0) {
        profile.push({
          r: rMid,
          rInner,
          rOuter,
          count: 0,
          sigmaR: 0.0,
          sigmaTheta: 0.0,
          sigmaPhi: 0.0,
          sigmaTan: 0.0,
          sigmaR2: 0.0,
          sigmaTan2: 0.0,
          beta: 0.0,
          symmetrizedBeta: 0.0,
          regime: AnisotropyRegime.ISOTROPIC
        });
        continue;
      }

      const meanSr2 = sumSigmaR2[b] / count;
      const meanStheta2 = sumSigmaTheta2[b] / count;
      const meanSphi2 = sumSigmaPhi2[b] / count;
      const meanStan2 = 0.5 * (meanStheta2 + meanSphi2);

      const sphPt = new SphericalAnisotropyPoint(rMid, meanSr2, meanStheta2, meanSphi2);

      profile.push({
        r: rMid,
        rInner,
        rOuter,
        count,
        sigmaR: sphPt.sigmaR,
        sigmaTheta: sphPt.sigmaTheta,
        sigmaPhi: sphPt.sigmaPhi,
        sigmaTan: sphPt.sigmaTan,
        sigmaR2: meanSr2,
        sigmaTan2: meanStan2,
        beta: sphPt.beta,
        symmetrizedBeta: sphPt.symmetrizedBeta,
        regime: sphPt.classifyRegime()
      });
    }

    return profile;
  }
}

// ============================================================================
// 5. COSMIC THERMAL & KINETIC PRESSURE FIELDS & SOUND SPEED
// ============================================================================

/**
 * CosmicPressureField encapsulates the 3D Kinetic and Thermal Pressure Tensor:
 * P_{ij}(\mathbf{x}) = \rho(\mathbf{x}) \sigma_{ij}^2(\mathbf{x})
 * Scalar Isotropic Pressure: P(\mathbf{x}) = \frac{1}{3} \text{Tr}(P_{ij}) = \rho(\mathbf{x}) \sigma_{1D}^2(\mathbf{x})
 * Effective Cosmic Sound Speed: c_s(\mathbf{x}) = \sqrt{\gamma \frac{P}{\rho}} = \sqrt{\gamma \sigma_{1D}^2}
 */
export class CosmicPressureField {
  /**
   * Constructs a CosmicPressureField.
   *
   * @param {VelocityDispersionTensor3D} dispersionField Velocity dispersion tensor field.
   * @param {DensityField|ScalarField3D} [densityField=null] Matter density field (\rho or 1+\delta).
   * @param {Object} [options] Configuration options.
   * @param {number} [options.meanDensity=RHO_CRIT_0_MSUN_MPC3 * DEFAULT_OMEGA_M] Background density \bar{\rho}_m in M_sun / (Mpc/h)^3.
   * @param {number} [options.gammaAdiabatic=DEFAULT_GAMMA_ADIABATIC] Adiabatic index \gamma (default 5/3).
   * @param {boolean} [options.isDensityContrast=true] True if density field represents \delta = (\rho - \bar{\rho})/\bar{\rho}.
   */
  constructor(dispersionField, densityField = null, options = {}) {
    if (!(dispersionField instanceof VelocityDispersionTensor3D)) {
      throw new TypeError('CosmicPressureField: dispersionField must be an instance of VelocityDispersionTensor3D.');
    }

    const grid = dispersionField.grid;
    this.dispersionField = dispersionField;
    this.densityField = densityField;
    this.grid = grid;

    this.meanDensity = Number(options.meanDensity ?? (RHO_CRIT_0_MSUN_MPC3 * DEFAULT_OMEGA_M));
    this.gammaAdiabatic = Number(options.gammaAdiabatic ?? DEFAULT_GAMMA_ADIABATIC);
    this.isDensityContrast = options.isDensityContrast !== false;

    const total = grid.totalCells;

    // Allocate 6 pressure tensor components: P_ij = rho * sigma_ij^2
    this.pxx = new Float64Array(total);
    this.pxy = new Float64Array(total);
    this.pxz = new Float64Array(total);
    this.pyy = new Float64Array(total);
    this.pyz = new Float64Array(total);
    this.pzz = new Float64Array(total);
    this.pIsotropic = new Float64Array(total);
    this.soundSpeed = new Float64Array(total);
    this.temperatureProxy = new Float64Array(total);

    // Compute pressure fields
    for (let i = 0; i < total; i++) {
      let rhoLocal = this.meanDensity;
      if (densityField) {
        const rawVal = densityField.data ? densityField.data[i] : (densityField.density ? densityField.density[i] : 0.0);
        if (this.isDensityContrast) {
          rhoLocal = Math.max(1e-6 * this.meanDensity, this.meanDensity * (1.0 + rawVal));
        } else {
          rhoLocal = Math.max(1e-6 * this.meanDensity, rawVal);
        }
      }

      const sxx = dispersionField.sxx[i];
      const sxy = dispersionField.sxy[i];
      const sxz = dispersionField.sxz[i];
      const syy = dispersionField.syy[i];
      const syz = dispersionField.syz[i];
      const szz = dispersionField.szz[i];

      this.pxx[i] = rhoLocal * sxx;
      this.pxy[i] = rhoLocal * sxy;
      this.pxz[i] = rhoLocal * sxz;
      this.pyy[i] = rhoLocal * syy;
      this.pyz[i] = rhoLocal * syz;
      this.pzz[i] = rhoLocal * szz;

      const pScalar = (this.pxx[i] + this.pyy[i] + this.pzz[i]) / 3.0;
      this.pIsotropic[i] = pScalar;

      // c_s = sqrt(gamma * P / rho) = sqrt(gamma * sigma_1D^2)
      const sig1d2 = (sxx + syy + szz) / 3.0;
      this.soundSpeed[i] = Math.sqrt(Math.max(0.0, this.gammaAdiabatic * sig1d2));
      this.temperatureProxy[i] = Math.max(0.0, sig1d2 * VELOCITY_DISPERSION_TO_KELVIN_FACTOR);
    }
  }

  /**
   * Returns a ScalarField3D of the isotropic thermal/kinetic pressure field P(\mathbf{x}).
   * Units: (M_sun/h) / (Mpc/h)^3 * (km/s)^2.
   * @returns {ScalarField3D}
   */
  getIsotropicPressureField() {
    return new ScalarField3D(this.grid, this.pIsotropic, 'isotropic_pressure', 'M_sun*(km/s)^2/(Mpc/h)^3');
  }

  /**
   * Returns a ScalarField3D of the effective sound speed field c_s(\mathbf{x}).
   * Units: km/s.
   * @returns {ScalarField3D}
   */
  getSoundSpeedField() {
    return new ScalarField3D(this.grid, this.soundSpeed, 'sound_speed', ScientificUnits.KM_PER_S);
  }

  /**
   * Returns a ScalarField3D of the cosmic gas temperature proxy T_K(\mathbf{x}).
   * Units: Kelvin (K).
   * @returns {ScalarField3D}
   */
  getTemperatureField() {
    return new ScalarField3D(this.grid, this.temperatureProxy, 'gas_temperature', 'K');
  }

  /**
   * Computes the cosmic Mach number field \mathcal{M}(\mathbf{x}) = \|\mathbf{v}(\mathbf{x})\| / c_s(\mathbf{x}).
   *
   * @param {VelocityField} velocityField Cosmic peculiar velocity field.
   * @returns {ScalarField3D}
   */
  computeMachNumberField(velocityField) {
    if (!(velocityField instanceof VelocityField)) {
      throw new TypeError('CosmicPressureField.computeMachNumberField: velocityField must be an instance of VelocityField.');
    }

    const total = this.grid.totalCells;
    const machData = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const vx = velocityField.vx[i];
      const vy = velocityField.vy[i];
      const vz = velocityField.vz[i];
      const vMag = Math.sqrt(vx * vx + vy * vy + vz * vz);
      const cs = this.soundSpeed[i];

      machData[i] = cs > NUMERICAL_EPSILON ? vMag / cs : 0.0;
    }

    return new ScalarField3D(this.grid, machData, 'mach_number', ScientificUnits.DIMENSIONLESS);
  }
}

// ============================================================================
// 6. JEANS EQUATION COSMOLOGICAL MASS ESTIMATOR
// ============================================================================

/**
 * JeansMassEstimator solves the Spherical Jeans Equation for dynamical mass profiles:
 * M_{\text{Jeans}}(<r) = -\frac{r \sigma_r^2(r)}{G} \left( \frac{d \ln \rho}{d \ln r} + \frac{d \ln \sigma_r^2}{d \ln r} + 2\beta(r) \right)
 *
 * Let \gamma_\rho = -\frac{d \ln \rho}{d \ln r} (density logarithmic slope)
 * and \gamma_\sigma = -\frac{d \ln \sigma_r^2}{d \ln r} (dispersion logarithmic slope):
 * M_{\text{Jeans}}(<r) = \frac{r \sigma_r^2(r)}{G} \left( \gamma_\rho(r) + \gamma_\sigma(r) - 2\beta(r) \right)
 */
export class JeansMassEstimator {
  /**
   * Computes the regularized logarithmic derivative d ln(f) / d ln(r) using central finite differences on log scale.
   *
   * @param {Float64Array|Array<number>} rArr Radial distances (must be strictly positive and strictly increasing).
   * @param {Float64Array|Array<number>} fArr Function values (must be strictly positive).
   * @returns {Float64Array} Logarithmic derivatives at each radial point.
   */
  static computeLogarithmicSlope(rArr, fArr) {
    const n = rArr.length;
    if (n < 2) {
      throw new Error('JeansMassEstimator.computeLogarithmicSlope: Array length must be at least 2.');
    }

    const slopes = new Float64Array(n);
    const logR = new Float64Array(n);
    const logF = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      if (rArr[i] <= 0.0) throw new RangeError(`rArr[${i}] must be strictly positive.`);
      logR[i] = Math.log(rArr[i]);
      logF[i] = Math.log(Math.max(NUMERICAL_EPSILON, fArr[i]));
    }

    // Boundary at i=0 (forward difference)
    slopes[0] = (logF[1] - logF[0]) / (logR[1] - logR[0]);

    // Interior points (second-order central difference)
    for (let i = 1; i < n - 1; i++) {
      const h1 = logR[i] - logR[i - 1];
      const h2 = logR[i + 1] - logR[i];
      // Non-uniform grid 3-point central derivative formula
      const df1 = logF[i] - logF[i - 1];
      const df2 = logF[i + 1] - logF[i];
      slopes[i] = (df1 / h1 * h2 + df2 / h2 * h1) / (h1 + h2);
    }

    // Boundary at i=n-1 (backward difference)
    slopes[n - 1] = (logF[n - 1] - logF[n - 2]) / (logR[n - 1] - logR[n - 2]);

    return slopes;
  }

  /**
   * Estimates the enclosed dynamical mass profile M_{\text{Jeans}}(<r) from radial profiles of
   * density \rho(r), radial velocity dispersion \sigma_r^2(r), and anisotropy parameter \beta(r).
   *
   * @param {Array<number>|Float64Array} radii Radial distances r in Mpc/h.
   * @param {Array<number>|Float64Array} densityProfile Matter density profile \rho(r) in M_sun / (Mpc/h)^3.
   * @param {Array<number>|Float64Array} sigmaR2Profile Radial dispersion squared \sigma_r^2(r) in (km/s)^2.
   * @param {Array<number>|Float64Array} [betaProfile=null] Anisotropy profile \beta(r) (defaults to isotropic \beta=0).
   * @param {Object} [options] Configuration options.
   * @param {number} [options.G=G_COSMO_MPC_MSUN] Gravitational constant.
   * @param {boolean} [options.enforceMonotonicity=true] Whether to enforce non-decreasing mass dM/dr >= 0.
   * @returns {Array<{ r: number, mass: number, dlnRho_dlnR: number, dlnSigma_dlnR: number, beta: number, vCirc: number }>}
   */
  static estimateEnclosedMassProfile(radii, densityProfile, sigmaR2Profile, betaProfile = null, options = {}) {
    const n = radii.length;
    if (densityProfile.length !== n || sigmaR2Profile.length !== n) {
      throw new Error('JeansMassEstimator.estimateEnclosedMassProfile: Profile arrays must have identical length.');
    }

    const G = Number(options.G ?? G_COSMO_MPC_MSUN);
    const enforceMonotonicity = options.enforceMonotonicity !== false;

    // Compute logarithmic derivatives
    const dlnRho = JeansMassEstimator.computeLogarithmicSlope(radii, densityProfile);
    const dlnSigma = JeansMassEstimator.computeLogarithmicSlope(radii, sigmaR2Profile);

    const result = [];
    let prevMass = 0.0;

    for (let i = 0; i < n; i++) {
      const r = radii[i];
      const sr2 = Math.max(0.0, sigmaR2Profile[i]);
      const beta = betaProfile ? betaProfile[i] : 0.0;

      // Jeans bracket: -(dlnRho/dlnR + dlnSigma2/dlnR + 2*beta)
      // or: gamma_rho + gamma_sigma - 2*beta
      const gammaRho = -dlnRho[i];
      const gammaSigma = -dlnSigma[i];
      const bracket = gammaRho + gammaSigma - 2.0 * beta;

      let mass = (r * sr2 / G) * bracket;
      mass = Math.max(0.0, mass);

      if (enforceMonotonicity && i > 0 && mass < prevMass) {
        mass = prevMass;
      }
      prevMass = mass;

      // Circular velocity v_c(r) = sqrt(G M / r)
      const vCirc = Math.sqrt(Math.max(0.0, G * mass / r));

      result.push({
        r,
        mass,
        dlnRho_dlnR: dlnRho[i],
        dlnSigma_dlnR: dlnSigma[i],
        beta,
        vCirc
      });
    }

    return result;
  }

  /**
   * Exact analytical Jeans mass verification for a Singular Isothermal Sphere (SIS):
   * \rho(r) = \frac{\sigma^2}{2 \pi G r^2} \implies \frac{d \ln \rho}{d \ln r} = -2
   * \sigma_r^2(r) = \sigma^2 = \text{const} \implies \frac{d \ln \sigma_r^2}{d \ln r} = 0
   * Isotropic orbits \implies \beta = 0
   * Exact Solution: M(<r) = \frac{2 \sigma^2 r}{G}
   *
   * @param {number} r Radial distance in Mpc/h.
   * @param {number} sigma1D 1D isothermal velocity dispersion in km/s.
   * @param {number} [G=G_COSMO_MPC_MSUN] Gravitational constant.
   * @returns {number} Exact SIS enclosed mass in M_sun.
   */
  static singularIsothermalSphereExactMass(r, sigma1D, G = G_COSMO_MPC_MSUN) {
    if (r <= 0.0 || sigma1D <= 0.0) {
      throw new RangeError('singularIsothermalSphereExactMass: r and sigma1D must be strictly positive.');
    }
    return (2.0 * sigma1D * sigma1D * r) / G;
  }

  /**
   * Exact analytical Jeans mass for a Hernquist profile:
   * \rho(r) = \frac{M_{\text{tot}} a}{2 \pi r (r + a)^3}
   * M_{\text{exact}}(<r) = M_{\text{tot}} \frac{r^2}{(r + a)^2}
   *
   * @param {number} r Radius in Mpc/h.
   * @param {number} totalMass Total mass M_tot in M_sun.
   * @param {number} scaleRadius Hernquist scale radius a in Mpc/h.
   * @returns {number} Exact Hernquist enclosed mass in M_sun.
   */
  static hernquistExactMass(r, totalMass, scaleRadius) {
    if (r <= 0.0 || totalMass <= 0.0 || scaleRadius <= 0.0) {
      throw new RangeError('hernquistExactMass: Parameters must be strictly positive.');
    }
    const term = r / (r + scaleRadius);
    return totalMass * term * term;
  }

  /**
   * Exact analytical NFW profile enclosed mass:
   * M(<r) = 4 \pi \rho_0 r_s^3 \left[ \ln\left(1 + \frac{r}{r_s}\right) - \frac{r / r_s}{1 + r / r_s} \right]
   *
   * @param {number} r Radius in Mpc/h.
   * @param {number} rho0 Characteristic scale density \rho_0 in M_sun / (Mpc/h)^3.
   * @param {number} rs Scale radius r_s in Mpc/h.
   * @returns {number} Exact NFW enclosed mass in M_sun.
   */
  static nfwExactMass(r, rho0, rs) {
    if (r <= 0.0 || rho0 <= 0.0 || rs <= 0.0) {
      throw new RangeError('nfwExactMass: Parameters must be strictly positive.');
    }
    const x = r / rs;
    const factor = 4.0 * Math.PI * rho0 * rs * rs * rs;
    return factor * (Math.log(1.0 + x) - (x / (1.0 + x)));
  }
}

// ============================================================================
// 7. SUPERCLUSTER MULTISCALE VIRIAL RATIO (2K / |W|)
// ============================================================================

/**
 * SuperclusterVirialAnalyzer computes the global and multiscale Virial Ratio 2K / |W|
 * for cosmological supercluster structures (e.g. Laniakea, Shapley, Perseus-Pisces):
 *
 * Kinetic Energy:
 *   K_{\text{bulk}} = \frac{1}{2} \int_V \rho(\mathbf{x}) \|\mathbf{v}(\mathbf{x}) - \mathbf{v}_{\text{CM}}\|^2 d^3x
 *   K_{\text{disp}} = \frac{1}{2} \int_V \rho(\mathbf{x}) \text{Tr}(\sigma^2(\mathbf{x})) d^3x
 *   K_{\text{tot}} = K_{\text{bulk}} + K_{\text{disp}}
 *
 * Potential Energy:
 *   W = \frac{1}{2} \int_V \rho(\mathbf{x}) \Phi(\mathbf{x}) d^3x \quad \text{or} \quad W = -4\pi G \int_0^R \rho(r) M(<r) r dr
 *
 * Virial Ratio:
 *   \mathcal{V} = \frac{2 K_{\text{tot}}}{|W|}
 */
export class SuperclusterVirialAnalyzer {
  /**
   * Computes the global kinetic energy (bulk flow + internal dispersion) across a grid volume.
   *
   * @param {GridIndexer} grid Grid indexer.
   * @param {VelocityField} velocityField Peculiar velocity field (vx, vy, vz).
   * @param {VelocityDispersionTensor3D} dispersionField Velocity dispersion tensor field.
   * @param {ScalarField3D|DensityField} [densityField=null] Matter density field.
   * @param {Object} [options] Configuration options.
   * @param {number} [options.meanDensity=RHO_CRIT_0_MSUN_MPC3 * DEFAULT_OMEGA_M] Background density.
   * @returns {{ kBulk: number, kDisp: number, kTotal: number, centerOfMassVelocity: Float64Array, totalMass: number }}
   */
  static computeKineticEnergy(grid, velocityField, dispersionField, densityField = null, options = {}) {
    if (!(grid instanceof GridIndexer)) {
      throw new TypeError('computeKineticEnergy: grid must be an instance of GridIndexer.');
    }
    if (!(velocityField instanceof VelocityField) || !(dispersionField instanceof VelocityDispersionTensor3D)) {
      throw new TypeError('computeKineticEnergy: velocityField and dispersionField must be valid field instances.');
    }

    const total = grid.totalCells;
    const dV = grid.dx * grid.dy * grid.dz; // (Mpc/h)^3
    const meanRho = Number(options.meanDensity ?? (RHO_CRIT_0_MSUN_MPC3 * DEFAULT_OMEGA_M));

    // First pass: Compute total mass and center-of-mass bulk velocity
    let totalMass = 0.0;
    let sumPx = 0.0, sumPy = 0.0, sumPz = 0.0;

    for (let i = 0; i < total; i++) {
      let rho = meanRho;
      if (densityField) {
        const dVal = densityField.data ? densityField.data[i] : (densityField.density ? densityField.density[i] : 0.0);
        rho = Math.max(0.0, meanRho * (1.0 + dVal));
      }

      const cellMass = rho * dV;
      totalMass += cellMass;

      sumPx += cellMass * velocityField.vx[i];
      sumPy += cellMass * velocityField.vy[i];
      sumPz += cellMass * velocityField.vz[i];
    }

    const vcm_x = totalMass > 0.0 ? sumPx / totalMass : 0.0;
    const vcm_y = totalMass > 0.0 ? sumPy / totalMass : 0.0;
    const vcm_z = totalMass > 0.0 ? sumPz / totalMass : 0.0;

    // Second pass: Integrate bulk kinetic energy and internal dispersion kinetic energy
    let kBulk = 0.0;
    let kDisp = 0.0;

    for (let i = 0; i < total; i++) {
      let rho = meanRho;
      if (densityField) {
        const dVal = densityField.data ? densityField.data[i] : (densityField.density ? densityField.density[i] : 0.0);
        rho = Math.max(0.0, meanRho * (1.0 + dVal));
      }

      const cellMass = rho * dV;

      // Relative peculiar velocity relative to COM
      const dvx = velocityField.vx[i] - vcm_x;
      const dvy = velocityField.vy[i] - vcm_y;
      const dvz = velocityField.vz[i] - vcm_z;
      const vRelSq = dvx * dvx + dvy * dvy + dvz * dvz;

      kBulk += 0.5 * cellMass * vRelSq;

      // Internal velocity dispersion kinetic energy: 0.5 * m * Tr(\sigma^2)
      const trSigma2 = dispersionField.sxx[i] + dispersionField.syy[i] + dispersionField.szz[i];
      kDisp += 0.5 * cellMass * trSigma2;
    }

    return {
      kBulk,
      kDisp,
      kTotal: kBulk + kDisp,
      centerOfMassVelocity: new Float64Array([vcm_x, vcm_y, vcm_z]),
      totalMass
    };
  }

  /**
   * Computes gravitational potential energy W from density and gravitational potential fields:
   * W = \frac{1}{2} \int_V \rho(\mathbf{x}) \Phi(\mathbf{x}) d^3x
   *
   * @param {GridIndexer} grid Grid indexer.
   * @param {ScalarField3D|DensityField} densityField Density field \rho.
   * @param {ScalarField3D} potentialField Gravitational potential field \Phi in (km/s)^2.
   * @param {Object} [options] Configuration options.
   * @returns {{ potentialEnergy: number, absPotentialEnergy: number }}
   */
  static computePotentialEnergyFromGrid(grid, densityField, potentialField, options = {}) {
    if (!(grid instanceof GridIndexer)) {
      throw new TypeError('computePotentialEnergyFromGrid: grid must be an instance of GridIndexer.');
    }
    if (!(potentialField instanceof ScalarField3D)) {
      throw new TypeError('computePotentialEnergyFromGrid: potentialField must be an instance of ScalarField3D.');
    }

    const total = grid.totalCells;
    const dV = grid.dx * grid.dy * grid.dz;
    const meanRho = Number(options.meanDensity ?? (RHO_CRIT_0_MSUN_MPC3 * DEFAULT_OMEGA_M));

    let W = 0.0;
    for (let i = 0; i < total; i++) {
      let rho = meanRho;
      if (densityField) {
        const dVal = densityField.data ? densityField.data[i] : (densityField.density ? densityField.density[i] : 0.0);
        rho = Math.max(0.0, meanRho * (1.0 + dVal));
      }

      const cellMass = rho * dV;
      const phi = potentialField.data[i];

      W += 0.5 * cellMass * phi;
    }

    return {
      potentialEnergy: W,
      absPotentialEnergy: Math.abs(W)
    };
  }

  /**
   * Computes spherical gravitational potential energy by radial shell integration:
   * W = -4\pi G \int_0^R \rho(r) M(<r) r dr
   *
   * @param {Array<number>|Float64Array} radii Radial bin centers in Mpc/h.
   * @param {Array<number>|Float64Array} densityProfile Density profile \rho(r) in M_sun / (Mpc/h)^3.
   * @param {Array<number>|Float64Array} enclosedMassProfile Enclosed mass profile M(<r) in M_sun.
   * @param {number} [G=G_COSMO_MPC_MSUN] Gravitational constant.
   * @returns {number} Gravitational potential energy W in M_sun * (km/s)^2.
   */
  static computeSphericalPotentialEnergy(radii, densityProfile, enclosedMassProfile, G = G_COSMO_MPC_MSUN) {
    const n = radii.length;
    if (n < 2 || densityProfile.length !== n || enclosedMassProfile.length !== n) {
      throw new Error('computeSphericalPotentialEnergy: Invalid profile lengths.');
    }

    let W = 0.0;
    for (let i = 0; i < n - 1; i++) {
      const r1 = radii[i];
      const r2 = radii[i + 1];
      const dr = r2 - r1;
      const rMid = 0.5 * (r1 + r2);
      const rhoMid = 0.5 * (densityProfile[i] + densityProfile[i + 1]);
      const mMid = 0.5 * (enclosedMassProfile[i] + enclosedMassProfile[i + 1]);

      // Shell potential contribution: dW = - (G * M(<r) / r) * dM = - (G * M(<r) / r) * (4 pi r^2 rho dr)
      const dW = -4.0 * Math.PI * G * rhoMid * mMid * rMid * dr;
      W += dW;
    }

    return W;
  }

  /**
   * Evaluates the global Virial Ratio \mathcal{V} = 2K / |W| and classifies the dynamical state.
   *
   * @param {number} kineticEnergy Total kinetic energy K (M_sun * (km/s)^2).
   * @param {number} potentialEnergy Gravitational potential energy W (M_sun * (km/s)^2).
   * @param {Object} [options] Thresholds.
   * @param {number} [options.tolerance=0.15] Relative tolerance for equilibrium (|2K/|W| - 1| <= tolerance).
   * @returns {{ virialRatio: number, isVirialized: boolean, virialState: string, kineticEnergy: number, potentialEnergy: number, virialDeficit: number }}
   */
  static evaluateVirialRatio(kineticEnergy, potentialEnergy, options = {}) {
    const K = Number(kineticEnergy);
    const W = Number(potentialEnergy);
    const absW = Math.abs(W);
    const tol = Number(options.tolerance ?? 0.15);

    if (absW < NUMERICAL_EPSILON) {
      return {
        virialRatio: Infinity,
        isVirialized: false,
        virialState: VirialState.HIGHLY_UNBOUND_EXPANDING,
        kineticEnergy: K,
        potentialEnergy: W,
        virialDeficit: 2.0 * K
      };
    }

    const virialRatio = (2.0 * K) / absW;
    const virialDeficit = 2.0 * K + W; // 2K + W (should be 0 in exact equilibrium)
    const isVirialized = Math.abs(virialRatio - 1.0) <= tol;

    let virialState;
    if (virialRatio > 2.0) {
      virialState = VirialState.HIGHLY_UNBOUND_EXPANDING;
    } else if (virialRatio > 1.0 + tol) {
      virialState = VirialState.UNBOUND_SUB_VIRIAL;
    } else if (virialRatio < 1.0 - tol) {
      virialState = VirialState.COLLAPSING_SUPER_VIRIAL;
    } else {
      virialState = VirialState.VIRIALIZED_EQUILIBRIUM;
    }

    return {
      virialRatio,
      isVirialized,
      virialState,
      kineticEnergy: K,
      potentialEnergy: W,
      virialDeficit
    };
  }
}

// ============================================================================
// 8. DIAGNOSTIC AUDIT & W3C PROV-O LINEAGE MANIFEST
// ============================================================================

/**
 * Generates an exhaustive cryptographic diagnostic audit manifest for the velocity dispersion
 * and cosmic pressure field calculations, adhering to W3C PROV-O provenance standards.
 *
 * @param {VelocityDispersionTensor3D} dispersionField Velocity dispersion tensor field.
 * @param {CosmicPressureField} [pressureField=null] Optional cosmic pressure field.
 * @param {Object} [customMetadata={}] Additional run parameters.
 * @returns {Object} Full cryptographic audit report.
 */
export function generateDispersionAuditReport(dispersionField, pressureField = null, customMetadata = {}) {
  if (!(dispersionField instanceof VelocityDispersionTensor3D)) {
    throw new TypeError('generateDispersionAuditReport: dispersionField must be a VelocityDispersionTensor3D.');
  }

  const grid = dispersionField.grid;
  const total = grid.totalCells;

  // Measure distribution statistics of 1D dispersion
  let minSig1D = Infinity, maxSig1D = -Infinity, sumSig1D = 0.0, sumSig1DSq = 0.0;
  let minTriax = Infinity, maxTriax = -Infinity, sumTriax = 0.0;

  for (let i = 0; i < total; i++) {
    const tr = Math.max(0.0, dispersionField.sxx[i] + dispersionField.syy[i] + dispersionField.szz[i]);
    const sig1d = Math.sqrt(tr / 3.0);

    if (sig1d < minSig1D) minSig1D = sig1d;
    if (sig1d > maxSig1D) maxSig1D = sig1d;
    sumSig1D += sig1d;
    sumSig1DSq += sig1d * sig1d;

    const pt = new DispersionTensorPoint(
      dispersionField.sxx[i], dispersionField.sxy[i], dispersionField.sxz[i],
      dispersionField.syy[i], dispersionField.syz[i], dispersionField.szz[i]
    );
    const triax = pt.triaxiality;
    if (triax < minTriax) minTriax = triax;
    if (triax > maxTriax) maxTriax = triax;
    sumTriax += triax;
  }

  const meanSig1D = sumSig1D / total;
  const stdSig1D = Math.sqrt(Math.max(0.0, (sumSig1DSq / total) - (meanSig1D * meanSig1D)));
  const meanTriax = sumTriax / total;

  const payload = {
    activity: 'CosmicFlows-4:VelocityDispersionAndCosmicPressureAnalysis',
    timestampIso: new Date().toISOString(),
    grid: {
      nx: grid.nx,
      ny: grid.ny,
      nz: grid.nz,
      origin: grid.origin,
      boxSize: grid.boxSize,
      spacing: [grid.dx, grid.dy, grid.dz],
      totalCells: total
    },
    dispersionField: {
      smoothingRadius: dispersionField.smoothingRadius,
      kernelType: dispersionField.kernelType,
      sigma1D: {
        min: minSig1D,
        max: maxSig1D,
        mean: meanSig1D,
        std: stdSig1D,
        unit: ScientificUnits.KM_PER_S
      },
      triaxiality: {
        min: minTriax,
        max: maxTriax,
        mean: meanTriax
      }
    },
    pressureField: pressureField ? {
      meanDensity: pressureField.meanDensity,
      gammaAdiabatic: pressureField.gammaAdiabatic,
      isDensityContrast: pressureField.isDensityContrast
    } : null,
    constants: {
      G: G_COSMO_MPC_MSUN,
      rhoCrit0: RHO_CRIT_0_MSUN_MPC3,
      omegaM: DEFAULT_OMEGA_M,
      h0: DEFAULT_H0
    },
    customMetadata
  };

  const payloadJson = JSON.stringify(payload);
  const digestSha256 = sha256Hex(payloadJson);

  return {
    ...payload,
    provenance: {
      specification: 'W3C PROV-O / RFC 6234 SHA-256 Lineage',
      integrityDigest: digestSha256
    }
  };
}
