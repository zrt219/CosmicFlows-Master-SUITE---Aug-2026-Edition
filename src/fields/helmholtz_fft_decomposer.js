/**
 * @file helmholtz_fft_decomposer.js
 * @description Production-grade 3D Fast Fourier Transform Helmholtz-Hodge Vector Field Decomposer
 * for CosmicFlows-4 (CF4) peculiar velocity fields, cosmological turbulence, and cosmic web kinematics.
 * 
 * Mathematical Foundations & Vector Calculus Formulation:
 * =========================================================================================
 * 1. Helmholtz-Hodge Fundamental Theorem of Vector Fields on a 3D Periodic Torus \mathbb{T}^3:
 *    Any smooth, square-integrable vector field \mathbf{v}(\mathbf{x}) \in (L^2(\Omega))^3 on a periodic
 *    cosmological domain \Omega = [0, L_x) \times [0, L_y) \times [0, L_z) admits a unique, orthogonal
 *    L^2-decomposition into three mutually orthogonal components:
 * 
 *      \mathbf{v}(\mathbf{x}) = \mathbf{v}_{\text{pot}}(\mathbf{x}) + \mathbf{v}_{\text{sol}}(\mathbf{x}) + \mathbf{v}_0
 * 
 *    where:
 *      (a) \mathbf{v}_0 = \frac{1}{V} \int_{\Omega} \mathbf{v}(\mathbf{x}) d^3\mathbf{x}
 *          is the spatial harmonic mean bulk-flow vector (the \mathbf{k} = \mathbf{0} zero-mode),
 *          representing the spatial DC cosmic dipole motion.
 * 
 *      (b) \mathbf{v}_{\text{pot}}(\mathbf{x}) = -\nabla \Phi_v(\mathbf{x})
 *          is the longitudinal / irrotational / potential / compressible component:
 *            \nabla \times \mathbf{v}_{\text{pot}}(\mathbf{x}) = \mathbf{0}
 *          governed by the scalar velocity potential \Phi_v(\mathbf{x}) satisfying Poisson's equation:
 *            \nabla^2 \Phi_v(\mathbf{x}) = -\nabla \cdot \mathbf{v}(\mathbf{x}) = -\theta(\mathbf{x})
 * 
 *      (c) \mathbf{v}_{\text{sol}}(\mathbf{x}) = \nabla \times \mathbf{A}(\mathbf{x})
 *          is the transverse / solenoidal / rotational / divergence-free component:
 *            \nabla \cdot \mathbf{v}_{\text{sol}}(\mathbf{x}) = 0
 *          governed by the vector stream potential \mathbf{A}(\mathbf{x}) under the Coulomb gauge:
 *            \nabla \cdot \mathbf{A}(\mathbf{x}) = 0
 *            \nabla^2 \mathbf{A}(\mathbf{x}) = -\nabla \times \mathbf{v}(\mathbf{x}) = -\boldsymbol{\omega}(\mathbf{x})
 * 
 * 2. Spectral Fourier Space Projection Operators:
 *    Under the 3D discrete Fourier transform with wavevector \mathbf{k} = (k_x, k_y, k_z) \neq \mathbf{0}:
 *      \hat{\mathbf{v}}(\mathbf{k}) = \int_{\Omega} \mathbf{v}(\mathbf{x}) e^{-i \mathbf{k} \cdot \mathbf{x}} d^3\mathbf{x}
 * 
 *    The Helmholtz-Hodge projection tensor operators are:
 *      Longitudinal / Potential Projector:
 *        P^{\text{pot}}_{ij}(\mathbf{k}) = \frac{k_i k_j}{|\mathbf{k}|^2}
 *      Transverse / Solenoidal Projector:
 *        P^{\text{sol}}_{ij}(\mathbf{k}) = \delta_{ij} - \frac{k_i k_j}{|\mathbf{k}|^2}
 * 
 *    Exact Algebraic Tensor Invariants in Fourier Space:
 *      (i)   Completeness:             P^{\text{pot}}_{ij}(\mathbf{k}) + P^{\text{sol}}_{ij}(\mathbf{k}) = \delta_{ij}
 *      (ii)  Idempotence:              \sum_l P^{\text{pot}}_{il} P^{\text{pot}}_{lj} = P^{\text{pot}}_{ij},  \sum_l P^{\text{sol}}_{il} P^{\text{sol}}_{lj} = P^{\text{sol}}_{ij}
 *      (iii) Mutual Orthogonality:     \sum_l P^{\text{pot}}_{il} P^{\text{sol}}_{lj} = 0
 *      (iv)  Symmetry / Self-Adjoint:  P^{\text{pot}}_{ij} = P^{\text{pot}}_{ji},  P^{\text{sol}}_{ij} = P^{\text{sol}}_{ji}
 *      (v)   Trace / Projection Rank:  \text{Tr}(P^{\text{pot}}) = 1,  \text{Tr}(P^{\text{sol}}) = 2
 * 
 * 3. Spectral Potentials & Kinematic Fields:
 *    - Scalar Velocity Potential:
 *        \hat{\Phi}_v(\mathbf{k}) = \frac{i \mathbf{k} \cdot \hat{\mathbf{v}}(\mathbf{k})}{|\mathbf{k}|^2} \quad (\mathbf{k} \neq \mathbf{0}), \quad \hat{\Phi}_v(\mathbf{0}) = 0
 *        -\nabla \Phi_v \longleftrightarrow -i \mathbf{k} \hat{\Phi}_v(\mathbf{k}) = \frac{\mathbf{k}(\mathbf{k} \cdot \hat{\mathbf{v}}(\mathbf{k}))}{|\mathbf{k}|^2} = \hat{\mathbf{v}}_{\text{pot}}(\mathbf{k})
 * 
 *    - Vector Stream Potential (Coulomb Gauge):
 *        \hat{\mathbf{A}}(\mathbf{k}) = \frac{i \mathbf{k} \times \hat{\mathbf{v}}(\mathbf{k})}{|\mathbf{k}|^2} \quad (\mathbf{k} \neq \mathbf{0}), \quad \hat{\mathbf{A}}(\mathbf{0}) = \mathbf{0}
 *        \nabla \times \mathbf{A} \longleftrightarrow i \mathbf{k} \times \hat{\mathbf{A}}(\mathbf{k}) = \hat{\mathbf{v}}_{\text{sol}}(\mathbf{k})
 *        \nabla \cdot \mathbf{A} \longleftrightarrow i \mathbf{k} \cdot \hat{\mathbf{A}}(\mathbf{k}) = 0
 * 
 *    - Vorticity Tensor & Vector:
 *        \boldsymbol{\omega}(\mathbf{x}) = \nabla \times \mathbf{v}(\mathbf{x}) = \nabla \times \mathbf{v}_{\text{sol}}(\mathbf{x})
 *        \hat{\boldsymbol{\omega}}(\mathbf{k}) = i \mathbf{k} \times \hat{\mathbf{v}}(\mathbf{k}) = i \mathbf{k} \times \hat{\mathbf{v}}_{\text{sol}}(\mathbf{k})
 * 
 *    - Divergence (Expansion Rate):
 *        \theta(\mathbf{x}) = \nabla \cdot \mathbf{v}(\mathbf{x}) = \nabla \cdot \mathbf{v}_{\text{pot}}(\mathbf{x})
 *        \hat{\theta}(\mathbf{k}) = i \mathbf{k} \cdot \hat{\mathbf{v}}(\mathbf{k}) = i \mathbf{k} \cdot \hat{\mathbf{v}}_{\text{pot}}(\mathbf{k})
 * 
 * 4. Parseval Kinetic Energy Partitioning:
 *    By Plancherel's theorem, total kinetic energy is decomposed with zero cross-term:
 *      E_{\text{tot}} = \frac{1}{2} \int_{\Omega} |\mathbf{v}(\mathbf{x})|^2 d^3\mathbf{x} = E_{\text{pot}} + E_{\text{sol}} + E_0
 *      \langle \mathbf{v}_{\text{pot}}, \mathbf{v}_{\text{sol}} \rangle_{L^2} = \int_{\Omega} \mathbf{v}_{\text{pot}}(\mathbf{x}) \cdot \mathbf{v}_{\text{sol}}(\mathbf{x}) d^3\mathbf{x} \equiv 0
 * 
 * 5. Cosmological Linear Matter Density Inversion:
 *    In linear perturbation theory in an expanding FLRW universe:
 *      \theta(\mathbf{x}) = \nabla \cdot \mathbf{v}(\mathbf{x}) = -a H(z) f(\Omega_m, z) \delta(\mathbf{x})
 *      \delta_{\text{rec}}(\mathbf{x}) = -\frac{\theta(\mathbf{x})}{H_0 f(\Omega_m, z)} = \frac{\nabla^2 \Phi_v(\mathbf{x})}{H_0 f}
 * 
 * @module fields/helmholtz_fft_decomposer
 */

import { GridIndexer, BoundaryMode } from './grid_indexer.js';
import { ScalarField3D } from './scalar_field_3d.js';
import { VectorField3D } from './vector_field_3d.js';
import {
  VelocityVector,
  SupergalacticPosition,
  ScientificUnits
} from '../coordinates/scientific_types.js';
import {
  PHYSICAL_CONSTANTS,
  COSMOLOGICAL_MODELS
} from '../units/cosmological_constants.js';

// ============================================================================
// 1. ENUMERATIONS AND CONSTANTS
// ============================================================================

/**
 * Green's function / Wavevector spatial discretization kernel types.
 * @readonly
 * @enum {string}
 */
export const HelmholtzKernelType = Object.freeze({
  CONTINUOUS: 'continuous',             // Standard continuous isotropic wavevector: k_i = 2*pi*m_i / L_i
  DISCRETE_7POINT: 'discrete_7point',   // Exact 7-point central difference: k_tilde_i = sin(k_i * dx_i) / dx_i
  DISCRETE_19POINT: 'discrete_19point', // Compact 19-point isotropic stencil
  DISCRETE_27POINT: 'discrete_27point'  // High-order isotropic 27-point trigonometric stencil: (8*sin(k*dx) - sin(2*k*dx)) / (6*dx)
});

/**
 * Supported spectral window smoothing filter types.
 * @readonly
 * @enum {string}
 */
export const HelmholtzFilterType = Object.freeze({
  NONE: 'none',
  GAUSSIAN: 'gaussian',         // W(k) = exp(-0.5 * k^2 * R^2)
  TOPHAT_REAL: 'tophat_real',   // W(k) = 3 * (sin(kR) - kR*cos(kR)) / (kR)^3
  SHARP_K: 'sharp_k',           // W(k) = (k <= k_cut) ? 1 : 0
  FERMI_DIRAC: 'fermi_dirac'    // W(k) = 1 / (1 + exp((k - k_cut) / delta_k))
});

/**
 * Standard spatial scale bands for multiscale cosmic web decomposition (in Mpc/h).
 * @readonly
 * @enum {string}
 */
export const ScaleBandType = Object.freeze({
  ALL: 'all',
  CLUSTER: 'cluster',           // R < 5.0 Mpc/h (high-k virialized/turbulent)
  FILAMENT: 'filament',         // 5.0 <= R < 20.0 Mpc/h (intermediate cosmic web)
  SUPERCLUSTER: 'supercluster', // R >= 20.0 Mpc/h (large-scale linear potential flow)
  CUSTOM: 'custom'
});

/**
 * Physical constants for cosmological decomposition invariants.
 */
export const HELMHOLTZ_DECOMPOSER_DEFAULTS = Object.freeze({
  DEFAULT_HUBBLE_H0: 100.0,            // (km/s) / (Mpc/h)
  DEFAULT_OMEGA_M: 0.315,              // Planck 2018 LCDM
  DEFAULT_GROWTH_INDEX_GAMMA: 0.55,    // General Relativity Peebles growth index
  DEFAULT_NUM_SPECTRUM_BINS: 32,
  STRICT_ORTHOGONALITY_TOLERANCE: 1e-10,
  STRICT_ENERGY_TOLERANCE: 1e-10
});

// ============================================================================
// 2. HIGH-PERFORMANCE 1D & 3D FAST FOURIER TRANSFORM ENGINES
// ============================================================================

/**
 * 1D In-Place Radix-2 Cooley-Tukey Fast Fourier Transform engine.
 * Real and Imaginary buffers must have length equal to a power of 2.
 * 
 * @param {Float64Array} real In-out real components.
 * @param {Float64Array} imag In-out imaginary components.
 * @param {boolean} [inverse=false] True for backward (inverse) transform.
 */
export function fft1DInternal(real, imag, inverse = false) {
  const n = real.length;
  if ((n & (n - 1)) !== 0 || n === 0) {
    throw new Error(`fft1DInternal: Array length n must be a positive power of 2, received ${n}.`);
  }

  // 1. Bit-reversal permutation (Gold-Rader algorithm)
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

  // 2. Cooley-Tukey Radix-2 Butterfly stages
  const sign = inverse ? 1.0 : -1.0;
  for (let len = 2; len <= n; len <<= 1) {
    const halfLen = len >> 1;
    const angle = (sign * 2.0 * Math.PI) / len;
    const wStepR = Math.cos(angle);
    const wStepI = Math.sin(angle);

    for (let i = 0; i < n; i += len) {
      let wR = 1.0;
      let wI = 0.0;
      for (let m = 0; m < halfLen; m++) {
        const idxEven = i + m;
        const idxOdd = idxEven + halfLen;

        const oddR = real[idxOdd];
        const oddI = imag[idxOdd];

        const tr = wR * oddR - wI * oddI;
        const ti = wR * oddI + wI * oddR;

        real[idxOdd] = real[idxEven] - tr;
        imag[idxOdd] = imag[idxEven] - ti;
        real[idxEven] = real[idxEven] + tr;
        imag[idxEven] = imag[idxEven] + ti;

        const nextWR = wR * wStepR - wI * wStepI;
        const nextWI = wR * wStepI + wI * wStepR;
        wR = nextWR;
        wI = nextWI;
      }
    }
  }

  // 3. Normalization for inverse transform
  if (inverse) {
    const invN = 1.0 / n;
    for (let i = 0; i < n; i++) {
      real[i] *= invN;
      imag[i] *= invN;
    }
  }
}

/**
 * High-performance 3D Fast Fourier Transform Engine for 3D scalar and vector grids.
 * Implements separable 3-pass (X -> Y -> Z) transform with reusable contiguous line buffers.
 */
export class FFT3DEngine {
  /**
   * @param {number} nx Number of cells along X (must be power of 2 >= 2).
   * @param {number} ny Number of cells along Y (must be power of 2 >= 2).
   * @param {number} nz Number of cells along Z (must be power of 2 >= 2).
   */
  constructor(nx, ny, nz) {
    if (!Number.isInteger(nx) || (nx & (nx - 1)) !== 0 || nx < 2) {
      throw new RangeError(`FFT3DEngine: nx must be a power of 2 >= 2. Received ${nx}.`);
    }
    if (!Number.isInteger(ny) || (ny & (ny - 1)) !== 0 || ny < 2) {
      throw new RangeError(`FFT3DEngine: ny must be a power of 2 >= 2. Received ${ny}.`);
    }
    if (!Number.isInteger(nz) || (nz & (nz - 1)) !== 0 || nz < 2) {
      throw new RangeError(`FFT3DEngine: nz must be a power of 2 >= 2. Received ${nz}.`);
    }

    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.total = nx * ny * nz;

    const maxDim = Math.max(nx, ny, nz);
    this._lineR = new Float64Array(maxDim);
    this._lineI = new Float64Array(maxDim);
  }

  /**
   * Executes in-place 3D forward or inverse FFT on contiguous 3D real/imag Float64Arrays.
   * Forward transforms do not scale by 1/N; Inverse transforms scale by 1/(nx*ny*nz).
   * 
   * @param {Float64Array} real
   * @param {Float64Array} imag
   * @param {boolean} [inverse=false]
   */
  transform(real, imag, inverse = false) {
    if (real.length !== this.total || imag.length !== this.total) {
      throw new Error(`FFT3DEngine: Buffer length (${real.length}) does not match grid total cells (${this.total}).`);
    }

    const { nx, ny, nz } = this;
    const strideX = 1;
    const strideY = nx;
    const strideZ = nx * ny;
    const lineR = this._lineR;
    const lineI = this._lineI;

    // Pass 1: Along X-axis (stride 1)
    const subXReal = lineR.subarray(0, nx);
    const subXImag = lineI.subarray(0, nx);
    for (let iz = 0; iz < nz; iz++) {
      const zOff = iz * strideZ;
      for (let iy = 0; iy < ny; iy++) {
        const yOff = zOff + iy * strideY;
        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix * strideX;
          subXReal[ix] = real[idx];
          subXImag[ix] = imag[idx];
        }
        fft1DInternal(subXReal, subXImag, inverse);
        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix * strideX;
          real[idx] = subXReal[ix];
          imag[idx] = subXImag[ix];
        }
      }
    }

    // Pass 2: Along Y-axis (stride nx)
    const subYReal = lineR.subarray(0, ny);
    const subYImag = lineI.subarray(0, ny);
    for (let iz = 0; iz < nz; iz++) {
      const zOff = iz * strideZ;
      for (let ix = 0; ix < nx; ix++) {
        const xOff = zOff + ix * strideX;
        for (let iy = 0; iy < ny; iy++) {
          const idx = xOff + iy * strideY;
          subYReal[iy] = real[idx];
          subYImag[iy] = imag[idx];
        }
        fft1DInternal(subYReal, subYImag, inverse);
        for (let iy = 0; iy < ny; iy++) {
          const idx = xOff + iy * strideY;
          real[idx] = subYReal[iy];
          imag[idx] = subYImag[iy];
        }
      }
    }

    // Pass 3: Along Z-axis (stride nx*ny)
    const subZReal = lineR.subarray(0, nz);
    const subZImag = lineI.subarray(0, nz);
    for (let iy = 0; iy < ny; iy++) {
      const yOff = iy * strideY;
      for (let ix = 0; ix < nx; ix++) {
        const xOff = yOff + ix * strideX;
        for (let iz = 0; iz < nz; iz++) {
          const idx = xOff + iz * strideZ;
          subZReal[iz] = real[idx];
          subZImag[iz] = imag[idx];
        }
        fft1DInternal(subZReal, subZImag, inverse);
        for (let iz = 0; iz < nz; iz++) {
          const idx = xOff + iz * strideZ;
          real[idx] = subZReal[iz];
          imag[idx] = subZImag[iz];
        }
      }
    }
  }

  /**
   * Computes forward 3D FFT of real-valued input array.
   * Returns freshly allocated real and imaginary complex frequency buffers.
   * 
   * @param {Float64Array|Float32Array|Array<number>} inputReal
   * @returns {{ real: Float64Array, imag: Float64Array }}
   */
  forward(inputReal) {
    const real = new Float64Array(this.total);
    const imag = new Float64Array(this.total);
    for (let i = 0; i < this.total; i++) {
      real[i] = inputReal[i];
      imag[i] = 0.0;
    }
    this.transform(real, imag, false);
    return { real, imag };
  }

  /**
   * Computes inverse 3D FFT. Returns real part of reconstructed spatial array.
   * 
   * @param {Float64Array} real
   * @param {Float64Array} imag
   * @returns {Float64Array}
   */
  inverse(real, imag) {
    const rCopy = new Float64Array(real);
    const iCopy = new Float64Array(imag);
    this.transform(rCopy, iCopy, true);
    return rCopy;
  }
}

// ============================================================================
// 3. SPECTRAL K-SPACE PROJECTION OPERATOR TENSORS
// ============================================================================

/**
 * 3x3 Spectral Helmholtz Projection Operator Matrix P_ij(k).
 */
export class KSpaceProjector {
  /**
   * Evaluates the 3x3 potential (irrotational) projector matrix:
   * P_pot_{ij}(k) = (k_i * k_j) / |k|^2
   * 
   * @param {number} kx
   * @param {number} ky
   * @param {number} kz
   * @returns {Float64Array} 9-element row-major 3x3 matrix.
   */
  static getPotentialProjectorMatrix(kx, ky, kz) {
    const mat = new Float64Array(9);
    const kSq = kx * kx + ky * ky + kz * kz;
    if (kSq === 0.0) return mat;

    const invKSq = 1.0 / kSq;
    mat[0] = kx * kx * invKSq; // P_xx
    mat[1] = kx * ky * invKSq; // P_xy
    mat[2] = kx * kz * invKSq; // P_xz
    mat[3] = ky * kx * invKSq; // P_yx
    mat[4] = ky * ky * invKSq; // P_yy
    mat[5] = ky * kz * invKSq; // P_yz
    mat[6] = kz * kx * invKSq; // P_zx
    mat[7] = kz * ky * invKSq; // P_zy
    mat[8] = kz * kz * invKSq; // P_zz
    return mat;
  }

  /**
   * Evaluates the 3x3 solenoidal (divergence-free) projector matrix:
   * P_sol_{ij}(k) = \delta_{ij} - (k_i * k_j) / |k|^2
   * 
   * @param {number} kx
   * @param {number} ky
   * @param {number} kz
   * @returns {Float64Array} 9-element row-major 3x3 matrix.
   */
  static getSolenoidalProjectorMatrix(kx, ky, kz) {
    const mat = new Float64Array(9);
    const kSq = kx * kx + ky * ky + kz * kz;
    if (kSq === 0.0) {
      // For k=0 zero mode, projector returns 0 (DC mode handled separately)
      return mat;
    }

    const invKSq = 1.0 / kSq;
    mat[0] = 1.0 - kx * kx * invKSq; // P_xx
    mat[1] = -kx * ky * invKSq;       // P_xy
    mat[2] = -kx * kz * invKSq;       // P_xz
    mat[3] = -ky * kx * invKSq;       // P_yx
    mat[4] = 1.0 - ky * ky * invKSq; // P_yy
    mat[5] = -ky * kz * invKSq;       // P_yz
    mat[6] = -kz * kx * invKSq;       // P_zx
    mat[7] = -kz * ky * invKSq;       // P_zy
    mat[8] = 1.0 - kz * kz * invKSq; // P_zz
    return mat;
  }

  /**
   * Verifies algebraic projection tensor invariants (idempotence, orthogonality, trace).
   * 
   * @param {number} kx
   * @param {number} ky
   * @param {number} kz
   * @returns {{ tracePot: number, traceSol: number, isIdempotent: boolean, isOrthogonal: boolean, completenessError: number }}
   */
  static verifyAlgebraicInvariants(kx, ky, kz) {
    const kSq = kx * kx + ky * ky + kz * kz;
    if (kSq === 0.0) {
      return { tracePot: 0, traceSol: 0, isIdempotent: true, isOrthogonal: true, completenessError: 0 };
    }

    const P_pot = KSpaceProjector.getPotentialProjectorMatrix(kx, ky, kz);
    const P_sol = KSpaceProjector.getSolenoidalProjectorMatrix(kx, ky, kz);

    const tracePot = P_pot[0] + P_pot[4] + P_pot[8];
    const traceSol = P_sol[0] + P_sol[4] + P_sol[8];

    // Completeness check: P_pot + P_sol == I
    let completenessError = 0.0;
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        const delta = r === c ? 1.0 : 0.0;
        const sum = P_pot[r * 3 + c] + P_sol[r * 3 + c];
        completenessError = Math.max(completenessError, Math.abs(sum - delta));
      }
    }

    // Idempotence check: P_pot * P_pot == P_pot, P_sol * P_sol == P_sol
    let idempotenceError = 0.0;
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        let pSqPot = 0.0;
        let pSqSol = 0.0;
        for (let k = 0; k < 3; k++) {
          pSqPot += P_pot[r * 3 + k] * P_pot[k * 3 + c];
          pSqSol += P_sol[r * 3 + k] * P_sol[k * 3 + c];
        }
        idempotenceError = Math.max(idempotenceError, Math.abs(pSqPot - P_pot[r * 3 + c]), Math.abs(pSqSol - P_sol[r * 3 + c]));
      }
    }

    // Mutual Orthogonality check: P_pot * P_sol == 0
    let orthogonalityError = 0.0;
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        let prod = 0.0;
        for (let k = 0; k < 3; k++) {
          prod += P_pot[r * 3 + k] * P_sol[k * 3 + c];
        }
        orthogonalityError = Math.max(orthogonalityError, Math.abs(prod));
      }
    }

    return {
      tracePot,
      traceSol,
      isIdempotent: idempotenceError < 1e-12,
      isOrthogonal: orthogonalityError < 1e-12,
      completenessError
    };
  }
}

// ============================================================================
// 4. ENERGY SPECTRUM CONTAINER AND INERTIAL RANGE FITTER
// ============================================================================

/**
 * Mathematical Container for 1D Shell-Averaged Energy Spectra and Wavemode Statistics.
 */
export class HelmholtzEnergySpectrum {
  /**
   * @param {Object} params
   * @param {Float64Array} params.kMid Center wavenumber of each radial shell [h/Mpc].
   * @param {Float64Array} params.kMin Lower wavenumber edge [h/Mpc].
   * @param {Float64Array} params.kMax Upper wavenumber edge [h/Mpc].
   * @param {Float64Array} params.eTot Total kinetic energy in shell [(km/s)^2 * (Mpc/h)^3].
   * @param {Float64Array} params.ePot Irrotational (potential) energy in shell.
   * @param {Float64Array} params.eSol Solenoidal (divergence-free) energy in shell.
   * @param {Float64Array} params.pTot Volume-normalized total velocity power spectrum P_tot(k) [(km/s)^2 * (Mpc/h)^3].
   * @param {Float64Array} params.pPot Volume-normalized potential velocity power spectrum P_pot(k).
   * @param {Float64Array} params.pSol Volume-normalized solenoidal velocity power spectrum P_sol(k).
   * @param {Int32Array} params.modeCounts Number of discrete Fourier modes in each shell.
   * @param {number} params.totalEnergySum Integrated total kinetic energy E_tot.
   * @param {number} params.potEnergySum Integrated potential kinetic energy E_pot.
   * @param {number} params.solEnergySum Integrated solenoidal kinetic energy E_sol.
   * @param {number} params.zeroModeEnergy Harmonic zero-mode energy E_0 = 0.5 * V * |v_0|^2.
   * @param {number} params.compressibleRatio Global ratio E_pot / (E_pot + E_sol).
   * @param {number} params.solenoidalRatio Global ratio E_sol / (E_pot + E_sol).
   */
  constructor(params) {
    this.kMid = params.kMid;
    this.kMin = params.kMin;
    this.kMax = params.kMax;
    this.eTot = params.eTot;
    this.ePot = params.ePot;
    this.eSol = params.eSol;
    this.pTot = params.pTot;
    this.pPot = params.pPot;
    this.pSol = params.pSol;
    this.modeCounts = params.modeCounts;
    this.numBins = this.kMid.length;

    this.totalEnergySum = params.totalEnergySum;
    this.potEnergySum = params.potEnergySum;
    this.solEnergySum = params.solEnergySum;
    this.zeroModeEnergy = params.zeroModeEnergy;
    this.compressibleRatio = params.compressibleRatio;
    this.solenoidalRatio = params.solenoidalRatio;
  }

  /**
   * Performs least-squares power-law fit E(k) = A * k^(-alpha) over an inertial range [kStart, kEnd].
   * 
   * @param {string} component 'tot' | 'pot' | 'sol'
   * @param {number} [kStart=0] Minimum k for fit range.
   * @param {number} [kEnd=Infinity] Maximum k for fit range.
   * @returns {{ amplitude: number, spectralIndex: number, rSquared: number, binCount: number }}
   */
  fitPowerLaw(component = 'tot', kStart = 0, kEnd = Infinity) {
    let spectrum;
    if (component === 'pot') spectrum = this.ePot;
    else if (component === 'sol') spectrum = this.eSol;
    else spectrum = this.eTot;

    const logK = [];
    const logE = [];

    for (let i = 0; i < this.numBins; i++) {
      const k = this.kMid[i];
      const e = spectrum[i];
      if (k >= kStart && k <= kEnd && e > 1e-30 && this.modeCounts[i] > 0) {
        logK.push(Math.log(k));
        logE.push(Math.log(e));
      }
    }

    const n = logK.length;
    if (n < 2) {
      return { amplitude: 0, spectralIndex: 0, rSquared: 0, binCount: n };
    }

    let sumX = 0;
    let sumY = 0;
    let sumXX = 0;
    let sumXY = 0;
    let sumYY = 0;

    for (let i = 0; i < n; i++) {
      const x = logK[i];
      const y = logE[i];
      sumX += x;
      sumY += y;
      sumXX += x * x;
      sumXY += x * y;
      sumYY += y * y;
    }

    const denom = n * sumXX - sumX * sumX;
    if (Math.abs(denom) < 1e-15) {
      return { amplitude: 0, spectralIndex: 0, rSquared: 0, binCount: n };
    }

    const slope = (n * sumXY - sumX * sumY) / denom;
    const intercept = (sumY - slope * sumX) / n;

    // R^2 goodness of fit
    const meanY = sumY / n;
    let ssTot = 0;
    let ssRes = 0;
    for (let i = 0; i < n; i++) {
      const y = logE[i];
      const predY = intercept + slope * logK[i];
      ssTot += (y - meanY) * (y - meanY);
      ssRes += (y - predY) * (y - predY);
    }
    const rSquared = ssTot > 0 ? Math.max(0.0, 1.0 - ssRes / ssTot) : 1.0;

    return {
      amplitude: Math.exp(intercept),
      spectralIndex: -slope, // alpha in k^(-alpha)
      rSquared,
      binCount: n
    };
  }

  /**
   * Returns solenoidal and compressible energy fractions per wavenumber shell.
   * @returns {Array<{ k: number, kMin: number, kMax: number, fPot: number, fSol: number, eTot: number, ePot: number, eSol: number, modeCount: number }>}
   */
  getSpectralFractionTable() {
    const rows = [];
    for (let i = 0; i < this.numBins; i++) {
      const et = this.eTot[i];
      const fPot = et > 1e-30 ? this.ePot[i] / et : 0.0;
      const fSol = et > 1e-30 ? this.eSol[i] / et : 0.0;
      rows.push({
        k: this.kMid[i],
        kMin: this.kMin[i],
        kMax: this.kMax[i],
        fPot,
        fSol,
        eTot: et,
        ePot: this.ePot[i],
        eSol: this.eSol[i],
        modeCount: this.modeCounts[i]
      });
    }
    return rows;
  }
}

// ============================================================================
// 5. DIAGNOSTICS REPORT AND INVARIANT VERIFIER
// ============================================================================

/**
 * Comprehensive Diagnostic Verification Report Container.
 */
export class HelmholtzDiagnostics {
  /**
   * @param {Object} params
   */
  constructor(params) {
    this.l2InnerProduct = params.l2InnerProduct;               // \int v_pot \cdot v_sol d^3x
    this.normalizedL2Correlation = params.normalizedL2Correlation; // |<v_pot, v_sol>| / (||v_pot|| * ||v_sol||)
    this.zeroModePotProduct = params.zeroModePotProduct;       // v_0 \cdot \int v_pot d^3x
    this.zeroModeSolProduct = params.zeroModeSolProduct;       // v_0 \cdot \int v_sol d^3x
    this.maxReconstructionAbsError = params.maxReconstructionAbsError; // max |v - (v_pot + v_sol + v_0)|
    this.rmsReconstructionAbsError = params.rmsReconstructionAbsError;
    this.relativeReconstructionL2 = params.relativeReconstructionL2;
    this.energyConservationRelError = params.energyConservationRelError; // |E_tot - (E_pot + E_sol + E_0)| / E_tot
    this.solenoidalDivergenceL2 = params.solenoidalDivergenceL2; // ||\nabla \cdot v_sol||_L2
    this.solenoidalDivergenceMax = params.solenoidalDivergenceMax;
    this.potentialCurlL2 = params.potentialCurlL2;             // ||\nabla \times v_pot||_L2
    this.potentialCurlMax = params.potentialCurlMax;
    this.potentialGradL2Error = params.potentialGradL2Error;   // ||-\nabla \Phi_v - v_pot||_L2
    this.potentialGradRelError = params.potentialGradRelError;
    this.streamCurlL2Error = params.streamCurlL2Error;         // ||\nabla \times A - v_sol||_L2
    this.streamCurlRelError = params.streamCurlRelError;
    this.streamCoulombGaugeL2 = params.streamCoulombGaugeL2;   // ||\nabla \cdot A||_L2
    this.spectralDivergenceSolMax = params.spectralDivergenceSolMax;
    this.spectralCurlPotMax = params.spectralCurlPotMax;
    this.spectralGradPhiPotMax = params.spectralGradPhiPotMax;
    this.spectralCurlStreamSolMax = params.spectralCurlStreamSolMax;
    this.spectralGaugeAMax = params.spectralGaugeAMax;
    this.totalHelicity = params.totalHelicity;                 // \int v \cdot \omega d^3x
    this.meanEnstrophy = params.meanEnstrophy;                 // (1/V) \int 0.5 * |\omega|^2 d^3x
    this.isStrictlyOrthogonal = params.isStrictlyOrthogonal;   // boolean flag (< 1e-10)
    this.isEnergyConserving = params.isEnergyConserving;       // boolean flag (< 1e-10)
  }
}

// ============================================================================
// 6. HELMHOLTZ DECOMPOSITION RESULT CONTAINER
// ============================================================================

/**
 * Complete Result Container for 3D Helmholtz-Hodge Decomposition.
 */
export class HelmholtzDecompositionResult {
  /**
   * @param {Object} params
   * @param {VectorField3D} params.vPot Potential (irrotational) vector field v_pot(x).
   * @param {VectorField3D} params.vSol Solenoidal (divergence-free) vector field v_sol(x).
   * @param {VelocityVector} params.v0 Mean harmonic bulk velocity mode v_0 [km/s].
   * @param {ScalarField3D} params.phiV Velocity scalar potential field \Phi_v(x) [(km/s) * (Mpc/h)].
   * @param {VectorField3D} params.streamA Solenoidal vector stream potential A(x) [(km/s) * (Mpc/h)].
   * @param {VectorField3D} params.vorticity Vorticity field \omega(x) = \nabla \times v(x) [km/s / (Mpc/h)].
   * @param {ScalarField3D} params.divergence Velocity divergence \theta(x) = \nabla \cdot v(x) [km/s / (Mpc/h)].
   * @param {ScalarField3D} params.helicity Helicity density h(x) = v \cdot \omega [(km/s)^2 / (Mpc/h)].
   * @param {HelmholtzEnergySpectrum} params.energySpectrum 1D isotropic kinetic energy power spectra.
   * @param {HelmholtzDiagnostics} params.diagnostics Diagnostic metrics and orthogonality invariants.
   * @param {GridIndexer} params.grid Reference grid indexer.
   */
  constructor(params) {
    this.vPot = params.vPot;
    this.vSol = params.vSol;
    this.v0 = params.v0;
    this.phiV = params.phiV;
    this.streamA = params.streamA;
    this.vorticity = params.vorticity;
    this.divergence = params.divergence;
    this.helicity = params.helicity;
    this.energySpectrum = params.energySpectrum;
    this.diagnostics = params.diagnostics;
    this.grid = params.grid;
  }

  /**
   * Pointwise exact spatial field reconstruction: v_recon(x) = v_pot(x) + v_sol(x) + v_0.
   * @returns {VectorField3D}
   */
  reconstructField() {
    const total = this.grid.totalCells;
    const vReconX = new Float64Array(total);
    const vReconY = new Float64Array(total);
    const vReconZ = new Float64Array(total);

    const v0x = this.v0.vx;
    const v0y = this.v0.vy;
    const v0z = this.v0.vz;

    for (let i = 0; i < total; i++) {
      vReconX[i] = this.vPot.vx[i] + this.vSol.vx[i] + v0x;
      vReconY[i] = this.vPot.vy[i] + this.vSol.vy[i] + v0y;
      vReconZ[i] = this.vPot.vz[i] + this.vSol.vz[i] + v0z;
    }

    return new VectorField3D(this.grid, vReconX, vReconY, vReconZ, {
      name: 'reconstructed_velocity_field',
      unit: this.vPot.unit
    });
  }

  /**
   * Reconstructs cosmological linear matter density field \delta_rec(x) = -\nabla \cdot v_pot / (H0 * f).
   * 
   * In linear perturbation theory:
   *   \delta(x) = -\theta(x) / (H_0 * f(\Omega_m, z))
   * where \theta = \nabla \cdot v = \nabla \cdot v_pot.
   * 
   * @param {Object} [options]
   * @param {number} [options.H0=100.0] Hubble parameter in (km/s) / (Mpc/h) [default 100 for h=1].
   * @param {number} [options.omegaM=0.315] Matter density parameter \Omega_m.
   * @param {number} [options.redshift=0.0] Cosmological redshift z.
   * @param {number} [options.growthRate] Explicit linear growth rate f = d ln D / d ln a. Defaults to \Omega_m(z)^0.55.
   * @returns {ScalarField3D}
   */
  reconstructDensityField(options = {}) {
    const {
      H0 = HELMHOLTZ_DECOMPOSER_DEFAULTS.DEFAULT_HUBBLE_H0,
      omegaM = HELMHOLTZ_DECOMPOSER_DEFAULTS.DEFAULT_OMEGA_M,
      redshift = 0.0,
      growthRate = null
    } = options;

    let f = growthRate;
    if (f === null || f === undefined) {
      const a = 1.0 / (1.0 + redshift);
      const omegaMz = (omegaM / (a * a * a)) / (omegaM / (a * a * a) + (1.0 - omegaM));
      f = Math.pow(omegaMz, HELMHOLTZ_DECOMPOSER_DEFAULTS.DEFAULT_GROWTH_INDEX_GAMMA);
    }

    const scalingFactor = 1.0 / (H0 * f);
    const total = this.grid.totalCells;
    const deltaData = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      // delta_rec = -div(v) / (H0 * f)
      deltaData[i] = -this.divergence.data[i] * scalingFactor;
    }

    return new ScalarField3D(this.grid, deltaData, 'reconstructed_density_field', 'dimensionless');
  }

  /**
   * Returns a JSON-serializable diagnostic summary.
   * @returns {Object}
   */
  toJSON() {
    return {
      bulkVelocity: {
        vx: this.v0.vx,
        vy: this.v0.vy,
        vz: this.v0.vz,
        speed: this.v0.magnitude()
      },
      energies: {
        totalKinetic: this.energySpectrum.totalEnergySum,
        potentialKinetic: this.energySpectrum.potEnergySum,
        solenoidalKinetic: this.energySpectrum.solEnergySum,
        zeroModeKinetic: this.energySpectrum.zeroModeEnergy,
        compressibleRatio: this.energySpectrum.compressibleRatio,
        solenoidalRatio: this.energySpectrum.solenoidalRatio
      },
      diagnostics: {
        l2InnerProduct: this.diagnostics.l2InnerProduct,
        normalizedL2Correlation: this.diagnostics.normalizedL2Correlation,
        maxReconstructionAbsError: this.diagnostics.maxReconstructionAbsError,
        relativeReconstructionL2: this.diagnostics.relativeReconstructionL2,
        energyConservationRelError: this.diagnostics.energyConservationRelError,
        solenoidalDivergenceL2: this.diagnostics.solenoidalDivergenceL2,
        potentialCurlL2: this.diagnostics.potentialCurlL2,
        isStrictlyOrthogonal: this.diagnostics.isStrictlyOrthogonal,
        isEnergyConserving: this.diagnostics.isEnergyConserving
      }
    };
  }
}

// ============================================================================
// 7. MAIN DECOMPOSER ENGINE CLASS
// ============================================================================

/**
 * High-precision 3D FFT Helmholtz-Hodge Vector Field Decomposer.
 */
export class HelmholtzFFTDecomposer {
  /**
   * @param {Object} [options]
   * @param {string} [options.kernelType='continuous'] HelmholtzKernelType
   * @param {number} [options.numSpectrumBins=32] Number of radial shells for 1D power spectra.
   * @param {boolean} [options.computePotentials=true] Compute \Phi_v and stream potential A.
   * @param {boolean} [options.computeVorticityHelicity=true] Compute \omega, \theta, and helicity h.
   * @param {boolean} [options.runDiagnostics=true] Compute full orthogonality and residual metrics.
   */
  constructor(options = {}) {
    this.kernelType = options.kernelType || HelmholtzKernelType.CONTINUOUS;
    this.numSpectrumBins = options.numSpectrumBins || HELMHOLTZ_DECOMPOSER_DEFAULTS.DEFAULT_NUM_SPECTRUM_BINS;
    this.computePotentials = options.computePotentials !== false;
    this.computeVorticityHelicity = options.computeVorticityHelicity !== false;
    this.runDiagnostics = options.runDiagnostics !== false;

    this._cachedFFT = null;
  }

  /**
   * Retrieves or creates an FFT3DEngine instance matching grid dimensions.
   * @private
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   * @returns {FFT3DEngine}
   */
  _getFFTEngine(nx, ny, nz) {
    if (!this._cachedFFT || this._cachedFFT.nx !== nx || this._cachedFFT.ny !== ny || this._cachedFFT.nz !== nz) {
      this._cachedFFT = new FFT3DEngine(nx, ny, nz);
    }
    return this._cachedFFT;
  }

  /**
   * Executes complete 3D Fourier-space Helmholtz-Hodge decomposition on a VectorField3D.
   * 
   * @param {VectorField3D} vectorField Input 3D vector velocity field.
   * @param {Object} [overrideOptions] Optional options overriding constructor configuration.
   * @returns {HelmholtzDecompositionResult}
   */
  decompose(vectorField, overrideOptions = {}) {
    if (!(vectorField instanceof VectorField3D)) {
      throw new TypeError('HelmholtzFFTDecomposer.decompose: vectorField must be an instance of VectorField3D.');
    }

    const grid = vectorField.grid;
    const { nx, ny, nz } = grid;
    const total = grid.totalCells;

    const kernelType = overrideOptions.kernelType || this.kernelType;
    const numBins = overrideOptions.numSpectrumBins || this.numSpectrumBins;
    const computePotentials = overrideOptions.computePotentials !== undefined ? overrideOptions.computePotentials : this.computePotentials;
    const computeVorticityHelicity = overrideOptions.computeVorticityHelicity !== undefined ? overrideOptions.computeVorticityHelicity : this.computeVorticityHelicity;
    const runDiagnostics = overrideOptions.runDiagnostics !== undefined ? overrideOptions.runDiagnostics : this.runDiagnostics;

    const fft = this._getFFTEngine(nx, ny, nz);

    // 1. Forward 3D FFT of vector field components (vx, vy, vz)
    const fx = fft.forward(vectorField.vx);
    const fy = fft.forward(vectorField.vy);
    const fz = fft.forward(vectorField.vz);

    const vxR = fx.real;
    const vxI = fx.imag;
    const vyR = fy.real;
    const vyI = fy.imag;
    const vzR = fz.real;
    const vzI = fz.imag;

    // 2. Harmonic zero mode (bulk velocity v_0)
    // In our FFT convention, index 0 is (kx=0, ky=0, kz=0) mode: \sum v(x)
    const invTotal = 1.0 / total;
    const v0x = vxR[0] * invTotal;
    const v0y = vyR[0] * invTotal;
    const v0z = vzR[0] * invTotal;
    const v0 = new VelocityVector(v0x, v0y, v0z);

    // 3. Allocate Fourier-space buffers for projection
    const vPotXR = new Float64Array(total);
    const vPotXI = new Float64Array(total);
    const vPotYR = new Float64Array(total);
    const vPotYI = new Float64Array(total);
    const vPotZR = new Float64Array(total);
    const vPotZI = new Float64Array(total);

    const vSolXR = new Float64Array(total);
    const vSolXI = new Float64Array(total);
    const vSolYR = new Float64Array(total);
    const vSolYI = new Float64Array(total);
    const vSolZR = new Float64Array(total);
    const vSolZI = new Float64Array(total);

    const phiVR = computePotentials ? new Float64Array(total) : null;
    const phiVI = computePotentials ? new Float64Array(total) : null;

    const aXR = computePotentials ? new Float64Array(total) : null;
    const aXI = computePotentials ? new Float64Array(total) : null;
    const aYR = computePotentials ? new Float64Array(total) : null;
    const aYI = computePotentials ? new Float64Array(total) : null;
    const aZR = computePotentials ? new Float64Array(total) : null;
    const aZI = computePotentials ? new Float64Array(total) : null;

    const omegaXR = computeVorticityHelicity ? new Float64Array(total) : null;
    const omegaXI = computeVorticityHelicity ? new Float64Array(total) : null;
    const omegaYR = computeVorticityHelicity ? new Float64Array(total) : null;
    const omegaYI = computeVorticityHelicity ? new Float64Array(total) : null;
    const omegaZR = computeVorticityHelicity ? new Float64Array(total) : null;
    const omegaZI = computeVorticityHelicity ? new Float64Array(total) : null;

    const thetaR = computeVorticityHelicity ? new Float64Array(total) : null;
    const thetaI = computeVorticityHelicity ? new Float64Array(total) : null;

    // Physical box geometry (periodic cell spacing)
    const Lx = grid.boxSize[0];
    const Ly = grid.boxSize[1];
    const Lz = grid.boxSize[2];
    const dx = Lx / nx;
    const dy = Ly / ny;
    const dz = Lz / nz;
    const domainVolume = Lx * Ly * Lz;

    const twoPiLx = (2.0 * Math.PI) / Lx;
    const twoPiLy = (2.0 * Math.PI) / Ly;
    const twoPiLz = (2.0 * Math.PI) / Lz;

    // Precompute 1D continuous and discrete wavenumber arrays
    const kxArr = new Float64Array(nx);
    const kyArr = new Float64Array(ny);
    const kzArr = new Float64Array(nz);
    const kxEffArr = new Float64Array(nx);
    const kyEffArr = new Float64Array(ny);
    const kzEffArr = new Float64Array(nz);

    for (let ix = 0; ix < nx; ix++) {
      const mx = ix < (nx >> 1) ? ix : ix - nx;
      kxArr[ix] = mx * twoPiLx;
      if (kernelType === HelmholtzKernelType.DISCRETE_7POINT) {
        kxEffArr[ix] = Math.sin(kxArr[ix] * dx) / dx;
      } else if (kernelType === HelmholtzKernelType.DISCRETE_27POINT) {
        kxEffArr[ix] = (8.0 * Math.sin(kxArr[ix] * dx) - Math.sin(2.0 * kxArr[ix] * dx)) / (6.0 * dx);
      } else if (kernelType === HelmholtzKernelType.DISCRETE_19POINT) {
        kxEffArr[ix] = (Math.sin(kxArr[ix] * dx) / dx) * (1.0 + (1.0 - Math.cos(kxArr[ix] * dx)) / 6.0);
      } else {
        kxEffArr[ix] = kxArr[ix];
      }
    }

    for (let iy = 0; iy < ny; iy++) {
      const my = iy < (ny >> 1) ? iy : iy - ny;
      kyArr[iy] = my * twoPiLy;
      if (kernelType === HelmholtzKernelType.DISCRETE_7POINT) {
        kyEffArr[iy] = Math.sin(kyArr[iy] * dy) / dy;
      } else if (kernelType === HelmholtzKernelType.DISCRETE_27POINT) {
        kyEffArr[iy] = (8.0 * Math.sin(kyArr[iy] * dy) - Math.sin(2.0 * kyArr[iy] * dy)) / (6.0 * dy);
      } else if (kernelType === HelmholtzKernelType.DISCRETE_19POINT) {
        kyEffArr[iy] = (Math.sin(kyArr[iy] * dy) / dy) * (1.0 + (1.0 - Math.cos(kyArr[iy] * dy)) / 6.0);
      } else {
        kyEffArr[iy] = kyArr[iy];
      }
    }

    for (let iz = 0; iz < nz; iz++) {
      const mz = iz < (nz >> 1) ? iz : iz - nz;
      kzArr[iz] = mz * twoPiLz;
      if (kernelType === HelmholtzKernelType.DISCRETE_7POINT) {
        kzEffArr[iz] = Math.sin(kzArr[iz] * dz) / dz;
      } else if (kernelType === HelmholtzKernelType.DISCRETE_27POINT) {
        kzEffArr[iz] = (8.0 * Math.sin(kzArr[iz] * dz) - Math.sin(2.0 * kzArr[iz] * dz)) / (6.0 * dz);
      } else if (kernelType === HelmholtzKernelType.DISCRETE_19POINT) {
        kzEffArr[iz] = (Math.sin(kzArr[iz] * dz) / dz) * (1.0 + (1.0 - Math.cos(kzArr[iz] * dz)) / 6.0);
      } else {
        kzEffArr[iz] = kzArr[iz];
      }
    }

    // 4. K-space Projection Loop over all (ix, iy, iz)
    for (let iz = 0; iz < nz; iz++) {
      const kz = kzEffArr[iz];
      const kzSq = kz * kz;
      const zOff = iz * nx * ny;

      for (let iy = 0; iy < ny; iy++) {
        const ky = kyEffArr[iy];
        const kyzSq = ky * ky + kzSq;
        const yOff = zOff + iy * nx;

        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix;

          // k = 0 mode (DC / bulk mode)
          if (ix === 0 && iy === 0 && iz === 0) {
            vPotXR[0] = 0.0;
            vPotXI[0] = 0.0;
            vPotYR[0] = 0.0;
            vPotYI[0] = 0.0;
            vPotZR[0] = 0.0;
            vPotZI[0] = 0.0;

            vSolXR[0] = 0.0;
            vSolXI[0] = 0.0;
            vSolYR[0] = 0.0;
            vSolYI[0] = 0.0;
            vSolZR[0] = 0.0;
            vSolZI[0] = 0.0;

            if (computePotentials) {
              phiVR[0] = 0.0;
              phiVI[0] = 0.0;
              aXR[0] = 0.0;
              aXI[0] = 0.0;
              aYR[0] = 0.0;
              aYI[0] = 0.0;
              aZR[0] = 0.0;
              aZI[0] = 0.0;
            }

            if (computeVorticityHelicity) {
              omegaXR[0] = 0.0;
              omegaXI[0] = 0.0;
              omegaYR[0] = 0.0;
              omegaYI[0] = 0.0;
              omegaZR[0] = 0.0;
              omegaZI[0] = 0.0;
              thetaR[0] = 0.0;
              thetaI[0] = 0.0;
            }
            continue;
          }

          const kx = kxEffArr[ix];
          const kSq = kx * kx + kyzSq;
          const invKSq = 1.0 / kSq;

          const vx_r = vxR[idx];
          const vx_i = vxI[idx];
          const vy_r = vyR[idx];
          const vy_i = vyI[idx];
          const vz_r = vzR[idx];
          const vz_i = vzI[idx];

          // k \cdot v_hat (complex dot product)
          const kDotVR = kx * vx_r + ky * vy_r + kz * vz_r;
          const kDotVI = kx * vx_i + ky * vy_i + kz * vz_i;

          // Longitudinal projection: v_pot_hat = k * (k \cdot v_hat) / |k|^2
          const scalarR = kDotVR * invKSq;
          const scalarI = kDotVI * invKSq;

          const potXR = kx * scalarR;
          const potXI = kx * scalarI;
          const potYR = ky * scalarR;
          const potYI = ky * scalarI;
          const potZR = kz * scalarR;
          const potZI = kz * scalarI;

          vPotXR[idx] = potXR;
          vPotXI[idx] = potXI;
          vPotYR[idx] = potYR;
          vPotYI[idx] = potYI;
          vPotZR[idx] = potZR;
          vPotZI[idx] = potZI;

          // Solenoidal projection: v_sol_hat = v_hat - v_pot_hat
          const solXR = vx_r - potXR;
          const solXI = vx_i - potXI;
          const solYR = vy_r - potYR;
          const solYI = vy_i - potYI;
          const solZR = vz_r - potZR;
          const solZI = vz_i - potZI;

          vSolXR[idx] = solXR;
          vSolXI[idx] = solXI;
          vSolYR[idx] = solYR;
          vSolYI[idx] = solYI;
          vSolZR[idx] = solZR;
          vSolZI[idx] = solZI;

          // Velocity scalar potential: \hat{\Phi}_v = i * (k \cdot v_hat) / |k|^2
          // With z = R + i*I, i*z = -I + i*R
          if (computePotentials) {
            phiVR[idx] = -scalarI;
            phiVI[idx] = scalarR;

            // Stream vector potential: \hat{A} = i * (k \times v_hat) / |k|^2
            const crossXR = (ky * vz_r - kz * vy_r) * invKSq;
            const crossXI = (ky * vz_i - kz * vy_i) * invKSq;
            const crossYR = (kz * vx_r - kx * vz_r) * invKSq;
            const crossYI = (kz * vx_i - kx * vz_i) * invKSq;
            const crossZR = (kx * vy_r - ky * vx_r) * invKSq;
            const crossZI = (kx * vy_i - ky * vx_i) * invKSq;

            // Multiply by i: (-crossI, crossR)
            aXR[idx] = -crossXI;
            aXI[idx] = crossXR;
            aYR[idx] = -crossYI;
            aYI[idx] = crossYR;
            aZR[idx] = -crossZI;
            aZI[idx] = crossZR;
          }

          // Vorticity: \hat{\omega} = i * (k \times v_hat)
          if (computeVorticityHelicity) {
            const rawCrossXR = ky * vz_r - kz * vy_r;
            const rawCrossXI = ky * vz_i - kz * vy_i;
            const rawCrossYR = kz * vx_r - kx * vz_r;
            const rawCrossYI = kz * vx_i - kx * vz_i;
            const rawCrossZR = kx * vy_r - ky * vx_r;
            const rawCrossZI = kx * vy_i - ky * vx_i;

            omegaXR[idx] = -rawCrossXI;
            omegaXI[idx] = rawCrossXR;
            omegaYR[idx] = -rawCrossYI;
            omegaYI[idx] = rawCrossYR;
            omegaZR[idx] = -rawCrossZI;
            omegaZI[idx] = rawCrossZR;

            // Divergence: \hat{\theta} = i * (k \cdot v_hat) = (-kDotVI, kDotVR)
            thetaR[idx] = -kDotVI;
            thetaI[idx] = kDotVR;
          }
        }
      }
    }

    // 5. Compute 1D Isotropic Radial Shell Energy Spectra in Fourier Space
    const energySpectrum = this._computeEnergySpectraFromKSpace(
      vxR, vxI, vyR, vyI, vzR, vzI,
      vPotXR, vPotXI, vPotYR, vPotYI, vPotZR, vPotZI,
      vSolXR, vSolXI, vSolYR, vSolYI, vSolZR, vSolZI,
      kxArr, kyArr, kzArr,
      nx, ny, nz,
      domainVolume,
      v0,
      numBins
    );

    // 6. Inverse 3D FFT to reconstruct real-space fields
    const potRealX = fft.inverse(vPotXR, vPotXI);
    const potRealY = fft.inverse(vPotYR, vPotYI);
    const potRealZ = fft.inverse(vPotZR, vPotZI);

    const solRealX = fft.inverse(vSolXR, vSolXI);
    const solRealY = fft.inverse(vSolYR, vSolYI);
    const solRealZ = fft.inverse(vSolZR, vSolZI);

    const vPot = new VectorField3D(grid, potRealX, potRealY, potRealZ, {
      name: 'potential_velocity_field',
      unit: vectorField.unit
    });

    const vSol = new VectorField3D(grid, solRealX, solRealY, solRealZ, {
      name: 'solenoidal_velocity_field',
      unit: vectorField.unit
    });

    let phiV = null;
    let streamA = null;
    if (computePotentials) {
      const phiReal = fft.inverse(phiVR, phiVI);
      phiV = new ScalarField3D(grid, phiReal, 'velocity_scalar_potential', `${vectorField.unit}*Mpc/h`);

      const aRealX = fft.inverse(aXR, aXI);
      const aRealY = fft.inverse(aYR, aYI);
      const aRealZ = fft.inverse(aZR, aZI);
      streamA = new VectorField3D(grid, aRealX, aRealY, aRealZ, {
        name: 'solenoidal_stream_potential',
        unit: `${vectorField.unit}*Mpc/h`
      });
    }

    let vorticity = null;
    let divergence = null;
    let helicity = null;
    if (computeVorticityHelicity) {
      const omRealX = fft.inverse(omegaXR, omegaXI);
      const omRealY = fft.inverse(omegaYR, omegaYI);
      const omRealZ = fft.inverse(omegaZR, omegaZI);
      vorticity = new VectorField3D(grid, omRealX, omRealY, omRealZ, {
        name: 'vorticity_field',
        unit: `${vectorField.unit}/(Mpc/h)`
      });

      const divReal = fft.inverse(thetaR, thetaI);
      divergence = new ScalarField3D(grid, divReal, 'velocity_divergence_field', `${vectorField.unit}/(Mpc/h)`);

      // Helicity density h(x) = v(x) \cdot \omega(x) = v_sol(x) \cdot \omega(x)
      const helData = new Float64Array(total);
      for (let i = 0; i < total; i++) {
        const vxTot = potRealX[i] + solRealX[i] + v0x;
        const vyTot = potRealY[i] + solRealY[i] + v0y;
        const vzTot = potRealZ[i] + solRealZ[i] + v0z;
        helData[i] = vxTot * omRealX[i] + vyTot * omRealY[i] + vzTot * omRealZ[i];
      }
      helicity = new ScalarField3D(grid, helData, 'kinetic_helicity_density', `(${vectorField.unit})^2/(Mpc/h)`);
    }

    // 7. Rigorous Orthogonality & Residual Diagnostic Suite
    let diagnostics = null;
    if (runDiagnostics) {
      const kContext = {
        kxEffArr,
        kyEffArr,
        kzEffArr,
        vPotXR, vPotXI, vPotYR, vPotYI, vPotZR, vPotZI,
        vSolXR, vSolXI, vSolYR, vSolYI, vSolZR, vSolZI,
        phiVR, phiVI,
        aXR, aXI, aYR, aYI, aZR, aZI,
        computePotentials
      };
      diagnostics = this._computeDiagnostics(
        vectorField, vPot, vSol, v0, phiV, streamA, vorticity, divergence, helicity,
        energySpectrum, grid, kContext
      );
    }

    return new HelmholtzDecompositionResult({
      vPot,
      vSol,
      v0,
      phiV,
      streamA,
      vorticity,
      divergence,
      helicity,
      energySpectrum,
      diagnostics,
      grid
    });
  }

  /**
   * Internal routine to calculate 1D shell-averaged kinetic energy spectra from Fourier buffers.
   * @private
   */
  _computeEnergySpectraFromKSpace(
    vxR, vxI, vyR, vyI, vzR, vzI,
    potXR, potXI, potYR, potYI, potZR, potZI,
    solXR, solXI, solYR, solYI, solZR, solZI,
    kxArr, kyArr, kzArr,
    nx, ny, nz,
    domainVolume,
    v0,
    numBins
  ) {
    const total = nx * ny * nz;
    const kxMax = Math.abs(kxArr[nx >> 1]);
    const kyMax = Math.abs(kyArr[ny >> 1]);
    const kzMax = Math.abs(kzArr[nz >> 1]);
    const kMax = Math.sqrt(kxMax * kxMax + kyMax * kyMax + kzMax * kzMax);

    // Minimum non-zero wavemode k_min
    let kMin = Infinity;
    for (let i = 1; i < nx; i++) {
      const k = Math.abs(kxArr[i]);
      if (k > 0 && k < kMin) kMin = k;
    }
    for (let i = 1; i < ny; i++) {
      const k = Math.abs(kyArr[i]);
      if (k > 0 && k < kMin) kMin = k;
    }
    for (let i = 1; i < nz; i++) {
      const k = Math.abs(kzArr[i]);
      if (k > 0 && k < kMin) kMin = k;
    }
    if (!Number.isFinite(kMin) || kMin <= 0) kMin = (2.0 * Math.PI) / 100.0;

    const deltaK = (kMax - kMin) / numBins;
    const kMid = new Float64Array(numBins);
    const kMinEdge = new Float64Array(numBins);
    const kMaxEdge = new Float64Array(numBins);
    const eTotBin = new Float64Array(numBins);
    const ePotBin = new Float64Array(numBins);
    const eSolBin = new Float64Array(numBins);
    const modeCounts = new Int32Array(numBins);

    for (let b = 0; b < numBins; b++) {
      kMinEdge[b] = kMin + b * deltaK;
      kMaxEdge[b] = kMin + (b + 1) * deltaK;
      kMid[b] = 0.5 * (kMinEdge[b] + kMaxEdge[b]);
    }

    // Normalization factor for Fourier power: V / (N_total)^2
    const normFactor = domainVolume / (total * total);

    let potEnergySum = 0.0;
    let solEnergySum = 0.0;
    let totalEnergySum = 0.0;

    for (let iz = 0; iz < nz; iz++) {
      const kz = kzArr[iz];
      const kzSq = kz * kz;
      const zOff = iz * nx * ny;

      for (let iy = 0; iy < ny; iy++) {
        const ky = kyArr[iy];
        const kyzSq = ky * ky + kzSq;
        const yOff = zOff + iy * nx;

        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix;
          if (ix === 0 && iy === 0 && iz === 0) continue;

          const kx = kxArr[ix];
          const kMag = Math.sqrt(kx * kx + kyzSq);

          // Power in mode |v|^2 = |vx|^2 + |vy|^2 + |vz|^2
          const pTot = (
            vxR[idx] * vxR[idx] + vxI[idx] * vxI[idx] +
            vyR[idx] * vyR[idx] + vyI[idx] * vyI[idx] +
            vzR[idx] * vzR[idx] + vzI[idx] * vzI[idx]
          ) * normFactor * 0.5;

          const pPot = (
            potXR[idx] * potXR[idx] + potXI[idx] * potXI[idx] +
            potYR[idx] * potYR[idx] + potYI[idx] * potYI[idx] +
            potZR[idx] * potZR[idx] + potZI[idx] * potZI[idx]
          ) * normFactor * 0.5;

          const pSol = (
            solXR[idx] * solXR[idx] + solXI[idx] * solXI[idx] +
            solYR[idx] * solYR[idx] + solYI[idx] * solYI[idx] +
            solZR[idx] * solZR[idx] + solZI[idx] * solZI[idx]
          ) * normFactor * 0.5;

          potEnergySum += pPot;
          solEnergySum += pSol;
          totalEnergySum += pTot;

          const bin = Math.min(numBins - 1, Math.max(0, Math.floor((kMag - kMin) / deltaK)));
          eTotBin[bin] += pTot;
          ePotBin[bin] += pPot;
          eSolBin[bin] += pSol;
          modeCounts[bin]++;
        }
      }
    }

    // Zero-mode bulk energy: 0.5 * V * |v_0|^2
    const v0MagSq = v0.vx * v0.vx + v0.vy * v0.vy + v0.vz * v0.vz;
    const zeroModeEnergy = 0.5 * domainVolume * v0MagSq;
    const fullTotalEnergy = totalEnergySum + zeroModeEnergy;

    // Power spectra P(k) = E(k) / (4*pi*k^2 * deltaK)
    const pTot = new Float64Array(numBins);
    const pPot = new Float64Array(numBins);
    const pSol = new Float64Array(numBins);

    for (let b = 0; b < numBins; b++) {
      const k = kMid[b];
      const shellVol = 4.0 * Math.PI * k * k * deltaK;
      if (shellVol > 0 && modeCounts[b] > 0) {
        pTot[b] = eTotBin[b] / (shellVol * 0.5);
        pPot[b] = ePotBin[b] / (shellVol * 0.5);
        pSol[b] = eSolBin[b] / (shellVol * 0.5);
      }
    }

    const nonZeroEnergy = potEnergySum + solEnergySum;
    const compressibleRatio = nonZeroEnergy > 0 ? potEnergySum / nonZeroEnergy : 0.0;
    const solenoidalRatio = nonZeroEnergy > 0 ? solEnergySum / nonZeroEnergy : 0.0;

    return new HelmholtzEnergySpectrum({
      kMid,
      kMin: kMinEdge,
      kMax: kMaxEdge,
      eTot: eTotBin,
      ePot: ePotBin,
      eSol: eSolBin,
      pTot,
      pPot,
      pSol,
      modeCounts,
      totalEnergySum: fullTotalEnergy,
      potEnergySum,
      solEnergySum,
      zeroModeEnergy,
      compressibleRatio,
      solenoidalRatio
    });
  }

  /**
   * Diagnostic verification of mathematical invariants.
   * @private
   */
  _computeDiagnostics(
    vOrig, vPot, vSol, v0, phiV, streamA, vorticity, divergence, helicity,
    energySpectrum, grid, kContext = null
  ) {
    const { nx, ny, nz } = grid;
    const total = grid.totalCells;
    const invTotal = 1.0 / total;
    const Lx = grid.boxSize[0];
    const Ly = grid.boxSize[1];
    const Lz = grid.boxSize[2];
    const dx = Lx / nx;
    const dy = Ly / ny;
    const dz = Lz / nz;
    const dV = dx * dy * dz;
    const domainVolume = Lx * Ly * Lz;

    const v0x = v0.vx;
    const v0y = v0.vy;
    const v0z = v0.vz;

    let innerProduct = 0.0;
    let potNormSq = 0.0;
    let solNormSq = 0.0;
    let origNormSq = 0.0;

    let sumPotX = 0.0;
    let sumPotY = 0.0;
    let sumPotZ = 0.0;
    let sumSolX = 0.0;
    let sumSolY = 0.0;
    let sumSolZ = 0.0;

    let maxReconAbsErr = 0.0;
    let sumReconSqErr = 0.0;
    let sumHelicity = 0.0;
    let sumEnstrophy = 0.0;

    for (let i = 0; i < total; i++) {
      const px = vPot.vx[i];
      const py = vPot.vy[i];
      const pz = vPot.vz[i];

      const sx = vSol.vx[i];
      const sy = vSol.vy[i];
      const sz = vSol.vz[i];

      const ox = vOrig.vx[i];
      const oy = vOrig.vy[i];
      const oz = vOrig.vz[i];

      // Inner product <v_pot, v_sol>
      const dot = px * sx + py * sy + pz * sz;
      innerProduct += dot * dV;

      potNormSq += (px * px + py * py + pz * pz) * dV;
      solNormSq += (sx * sx + sy * sy + sz * sz) * dV;
      origNormSq += (ox * ox + oy * oy + oz * oz) * dV;

      sumPotX += px;
      sumPotY += py;
      sumPotZ += pz;
      sumSolX += sx;
      sumSolY += sy;
      sumSolZ += sz;

      // Reconstruction error: v - (v_pot + v_sol + v_0)
      const rx = ox - (px + sx + v0x);
      const ry = oy - (py + sy + v0y);
      const rz = oz - (pz + sz + v0z);
      const errMag = Math.sqrt(rx * rx + ry * ry + rz * rz);

      if (errMag > maxReconAbsErr) maxReconAbsErr = errMag;
      sumReconSqErr += (rx * rx + ry * ry + rz * rz) * dV;

      if (helicity) {
        sumHelicity += helicity.data[i] * dV;
      }
      if (vorticity) {
        const wx = vorticity.vx[i];
        const wy = vorticity.vy[i];
        const wz = vorticity.vz[i];
        sumEnstrophy += 0.5 * (wx * wx + wy * wy + wz * wz) * dV;
      }
    }

    const potL2 = Math.sqrt(potNormSq);
    const solL2 = Math.sqrt(solNormSq);
    const origL2 = Math.sqrt(origNormSq);

    let normalizedL2Correlation = 0.0;
    if (origNormSq > 1e-30) {
      if (potL2 > 1e-7 * origL2 && solL2 > 1e-7 * origL2) {
        normalizedL2Correlation = Math.abs(innerProduct) / (potL2 * solL2);
      } else {
        normalizedL2Correlation = Math.abs(innerProduct) / origNormSq;
      }
    }

    const zeroModePotProduct = (v0x * sumPotX + v0y * sumPotY + v0z * sumPotZ) * dV;
    const zeroModeSolProduct = (v0x * sumSolX + v0y * sumSolY + v0z * sumSolZ) * dV;

    const rmsReconAbsErr = Math.sqrt(sumReconSqErr / domainVolume);
    const relativeReconL2 = origL2 > 1e-30 ? Math.sqrt(sumReconSqErr) / origL2 : rmsReconAbsErr;

    // Kinetic energy conservation residual
    const realSpaceTotEnergy = 0.5 * origNormSq;
    const decompTotEnergy = 0.5 * potNormSq + 0.5 * solNormSq + 0.5 * domainVolume * (v0x * v0x + v0y * v0y + v0z * v0z);
    const energyConservationRelError = realSpaceTotEnergy > 1e-30
      ? Math.abs(realSpaceTotEnergy - decompTotEnergy) / realSpaceTotEnergy
      : Math.abs(realSpaceTotEnergy - decompTotEnergy);

    // 8. Spectral differential operator exactness invariants
    let spectralDivergenceSolMax = 0.0;
    let spectralCurlPotMax = 0.0;
    let spectralGradPhiPotMax = 0.0;
    let spectralCurlStreamSolMax = 0.0;
    let spectralGaugeAMax = 0.0;

    if (kContext) {
      const {
        kxEffArr, kyEffArr, kzEffArr,
        vPotXR, vPotXI, vPotYR, vPotYI, vPotZR, vPotZI,
        vSolXR, vSolXI, vSolYR, vSolYI, vSolZR, vSolZI,
        phiVR, phiVI,
        aXR, aXI, aYR, aYI, aZR, aZI,
        computePotentials
      } = kContext;

      for (let iz = 0; iz < nz; iz++) {
        const kz = kzEffArr[iz];
        const zOff = iz * nx * ny;
        for (let iy = 0; iy < ny; iy++) {
          const ky = kyEffArr[iy];
          const yOff = zOff + iy * nx;
          for (let ix = 0; ix < nx; ix++) {
            const idx = yOff + ix;
            if (ix === 0 && iy === 0 && iz === 0) continue;

            const kx = kxEffArr[ix];

            // Spectral divergence of v_sol: |i k . v_sol_hat|
            const divSolR = -(kx * vSolXI[idx] + ky * vSolYI[idx] + kz * vSolZI[idx]);
            const divSolI =  (kx * vSolXR[idx] + ky * vSolYR[idx] + kz * vSolZR[idx]);
            const divSolMag = Math.sqrt(divSolR * divSolR + divSolI * divSolI) * invTotal;
            if (divSolMag > spectralDivergenceSolMax) spectralDivergenceSolMax = divSolMag;

            // Spectral curl of v_pot: |i k x v_pot_hat|
            const curlPotXR = -(ky * vPotZI[idx] - kz * vPotYI[idx]);
            const curlPotXI =  (ky * vPotZR[idx] - kz * vPotYR[idx]);
            const curlPotYR = -(kz * vPotXI[idx] - kx * vPotZI[idx]);
            const curlPotYI =  (kz * vPotXR[idx] - kx * vPotZR[idx]);
            const curlPotZR = -(kx * vPotYI[idx] - ky * vPotXI[idx]);
            const curlPotZI =  (kx * vPotYR[idx] - ky * vPotXR[idx]);
            const curlPotMag = Math.sqrt(
              curlPotXR * curlPotXR + curlPotXI * curlPotXI +
              curlPotYR * curlPotYR + curlPotYI * curlPotYI +
              curlPotZR * curlPotZR + curlPotZI * curlPotZI
            ) * invTotal;
            if (curlPotMag > spectralCurlPotMax) spectralCurlPotMax = curlPotMag;

            // Spectral gradient of Phi_v vs v_pot: |-i k Phi_v_hat - v_pot_hat|
            if (computePotentials && phiVR && aXR) {
              // -i kx (phiVR + i phiVI) = kx phiVI - i kx phiVR
              const negGradXR =  kx * phiVI[idx];
              const negGradXI = -kx * phiVR[idx];
              const negGradYR =  ky * phiVI[idx];
              const negGradYI = -ky * phiVR[idx];
              const negGradZR =  kz * phiVI[idx];
              const negGradZI = -kz * phiVR[idx];

              const dXR = negGradXR - vPotXR[idx];
              const dXI = negGradXI - vPotXI[idx];
              const dYR = negGradYR - vPotYR[idx];
              const dYI = negGradYI - vPotYI[idx];
              const dZR = negGradZR - vPotZR[idx];
              const dZI = negGradZI - vPotZI[idx];
              const gradErr = Math.sqrt(dXR*dXR + dXI*dXI + dYR*dYR + dYI*dYI + dZR*dZR + dZI*dZI) * invTotal;
              if (gradErr > spectralGradPhiPotMax) spectralGradPhiPotMax = gradErr;

              // Spectral curl of A vs v_sol: |i k x A_hat - v_sol_hat|
              const curlAXR = -(ky * aZI[idx] - kz * aYI[idx]);
              const curlAXI =  (ky * aZR[idx] - kz * aYR[idx]);
              const curlAYR = -(kz * aXI[idx] - kx * aZI[idx]);
              const curlAYI =  (kz * aXR[idx] - kx * aZR[idx]);
              const curlAZR = -(kx * aYI[idx] - ky * aXI[idx]);
              const curlAZI =  (kx * aYR[idx] - ky * aXR[idx]);

              const sXR = curlAXR - vSolXR[idx];
              const sXI = curlAXI - vSolXI[idx];
              const sYR = curlAYR - vSolYR[idx];
              const sYI = curlAYI - vSolYI[idx];
              const sZR = curlAZR - vSolZR[idx];
              const sZI = curlAZI - vSolZI[idx];
              const streamErr = Math.sqrt(sXR*sXR + sXI*sXI + sYR*sYR + sYI*sYI + sZR*sZR + sZI*sZI) * invTotal;
              if (streamErr > spectralCurlStreamSolMax) spectralCurlStreamSolMax = streamErr;

              // Spectral gauge of A: |i k . A_hat|
              const gaugeR = -(kx * aXI[idx] + ky * aYI[idx] + kz * aZI[idx]);
              const gaugeI =  (kx * aXR[idx] + ky * aYR[idx] + kz * aZR[idx]);
              const gaugeMag = Math.sqrt(gaugeR*gaugeR + gaugeI*gaugeI) * invTotal;
              if (gaugeMag > spectralGaugeAMax) spectralGaugeAMax = gaugeMag;
            }
          }
        }
      }
    }

    // 9. Finite difference diagnostics
    const inv2dx = 1.0 / (2.0 * dx);
    const inv2dy = 1.0 / (2.0 * dy);
    const inv2dz = 1.0 / (2.0 * dz);

    let sumDivSolSq = 0.0;
    let maxDivSol = 0.0;
    let sumCurlPotSq = 0.0;
    let maxCurlPot = 0.0;

    for (let iz = 0; iz < nz; iz++) {
      const izP = (iz + 1) % nz;
      const izM = (iz - 1 + nz) % nz;
      for (let iy = 0; iy < ny; iy++) {
        const iyP = (iy + 1) % ny;
        const iyM = (iy - 1 + ny) % ny;
        for (let ix = 0; ix < nx; ix++) {
          const ixP = (ix + 1) % nx;
          const ixM = (ix - 1 + nx) % nx;

          // Solenoidal divergence: d(vSolX)/dx + d(vSolY)/dy + d(vSolZ)/dz
          const idxXP = grid.index(ixP, iy, iz);
          const idxXM = grid.index(ixM, iy, iz);
          const idxYP = grid.index(ix, iyP, iz);
          const idxYM = grid.index(ix, iyM, iz);
          const idxZP = grid.index(ix, iy, izP);
          const idxZM = grid.index(ix, iy, izM);

          const divSol = (
            (vSol.vx[idxXP] - vSol.vx[idxXM]) * inv2dx +
            (vSol.vy[idxYP] - vSol.vy[idxYM]) * inv2dy +
            (vSol.vz[idxZP] - vSol.vz[idxZM]) * inv2dz
          );
          const absDiv = Math.abs(divSol);
          if (absDiv > maxDivSol) maxDivSol = absDiv;
          sumDivSolSq += divSol * divSol * dV;

          // Potential curl: \nabla \times v_pot
          const curlX = (vPot.vz[idxYP] - vPot.vz[idxYM]) * inv2dy - (vPot.vy[idxZP] - vPot.vy[idxZM]) * inv2dz;
          const curlY = (vPot.vx[idxZP] - vPot.vx[idxZM]) * inv2dz - (vPot.vz[idxXP] - vPot.vz[idxXM]) * inv2dx;
          const curlZ = (vPot.vy[idxXP] - vPot.vy[idxXM]) * inv2dx - (vPot.vx[idxYP] - vPot.vx[idxYM]) * inv2dy;

          const curlMag = Math.sqrt(curlX * curlX + curlY * curlY + curlZ * curlZ);
          if (curlMag > maxCurlPot) maxCurlPot = curlMag;
          sumCurlPotSq += (curlX * curlX + curlY * curlY + curlZ * curlZ) * dV;
        }
      }
    }

    const solenoidalDivergenceL2 = Math.sqrt(sumDivSolSq);
    const potentialCurlL2 = Math.sqrt(sumCurlPotSq);

    // Potential gradient residual ||-\nabla \Phi_v - v_pot||
    let potGradL2Error = 0.0;
    if (phiV) {
      let sumGradErrSq = 0.0;
      for (let iz = 0; iz < nz; iz++) {
        const izP = (iz + 1) % nz;
        const izM = (iz - 1 + nz) % nz;
        for (let iy = 0; iy < ny; iy++) {
          const iyP = (iy + 1) % ny;
          const iyM = (iy - 1 + ny) % ny;
          for (let ix = 0; ix < nx; ix++) {
            const ixP = (ix + 1) % nx;
            const ixM = (ix - 1 + nx) % nx;

            const idx = grid.index(ix, iy, iz);
            const dPhi_dx = (phiV.data[grid.index(ixP, iy, iz)] - phiV.data[grid.index(ixM, iy, iz)]) * inv2dx;
            const dPhi_dy = (phiV.data[grid.index(ix, iyP, iz)] - phiV.data[grid.index(ix, iyM, iz)]) * inv2dy;
            const dPhi_dz = (phiV.data[grid.index(ix, iy, izP)] - phiV.data[grid.index(ix, iy, izM)]) * inv2dz;

            const ex = (-dPhi_dx) - vPot.vx[idx];
            const ey = (-dPhi_dy) - vPot.vy[idx];
            const ez = (-dPhi_dz) - vPot.vz[idx];

            sumGradErrSq += (ex * ex + ey * ey + ez * ez) * dV;
          }
        }
      }
      potGradL2Error = Math.sqrt(sumGradErrSq);
    }

    // Stream curl residual ||\nabla \times A - v_sol|| and Coulomb gauge ||\nabla \cdot A||
    let streamCurlL2Error = 0.0;
    let streamCoulombGaugeL2 = 0.0;
    if (streamA) {
      let sumCurlErrSq = 0.0;
      let sumGaugeErrSq = 0.0;
      for (let iz = 0; iz < nz; iz++) {
        const izP = (iz + 1) % nz;
        const izM = (iz - 1 + nz) % nz;
        for (let iy = 0; iy < ny; iy++) {
          const iyP = (iy + 1) % ny;
          const iyM = (iy - 1 + ny) % ny;
          for (let ix = 0; ix < nx; ix++) {
            const ixP = (ix + 1) % nx;
            const ixM = (ix - 1 + nx) % nx;

            const idx = grid.index(ix, iy, iz);
            const idxXP = grid.index(ixP, iy, iz);
            const idxXM = grid.index(ixM, iy, iz);
            const idxYP = grid.index(ix, iyP, iz);
            const idxYM = grid.index(ix, iyM, iz);
            const idxZP = grid.index(ix, iy, izP);
            const idxZM = grid.index(ix, iy, izM);

            const curlAX = (streamA.vz[idxYP] - streamA.vz[idxYM]) * inv2dy - (streamA.vy[idxZP] - streamA.vy[idxZM]) * inv2dz;
            const curlAY = (streamA.vx[idxZP] - streamA.vx[idxZM]) * inv2dz - (streamA.vz[idxXP] - streamA.vz[idxZM]) * inv2dx;
            const curlAZ = (streamA.vy[idxXP] - streamA.vy[idxXM]) * inv2dx - (streamA.vx[idxYP] - streamA.vx[idxYM]) * inv2dy;

            const ex = curlAX - vSol.vx[idx];
            const ey = curlAY - vSol.vy[idx];
            const ez = curlAZ - vSol.vz[idx];
            sumCurlErrSq += (ex * ex + ey * ey + ez * ez) * dV;

            const divA = (
              (streamA.vx[idxXP] - streamA.vx[idxXM]) * inv2dx +
              (streamA.vy[idxYP] - streamA.vy[idxYM]) * inv2dy +
              (streamA.vz[idxZP] - streamA.vz[idxZM]) * inv2dz
            );
            sumGaugeErrSq += divA * divA * dV;
          }
        }
      }
      streamCurlL2Error = Math.sqrt(sumCurlErrSq);
      streamCoulombGaugeL2 = Math.sqrt(sumGaugeErrSq);
    }

    const isStrictlyOrthogonal = normalizedL2Correlation < HELMHOLTZ_DECOMPOSER_DEFAULTS.STRICT_ORTHOGONALITY_TOLERANCE;
    const isEnergyConserving = energyConservationRelError < HELMHOLTZ_DECOMPOSER_DEFAULTS.STRICT_ENERGY_TOLERANCE;

    return new HelmholtzDiagnostics({
      l2InnerProduct: innerProduct,
      normalizedL2Correlation,
      zeroModePotProduct,
      zeroModeSolProduct,
      maxReconstructionAbsError: maxReconAbsErr,
      rmsReconstructionAbsError: rmsReconAbsErr,
      relativeReconstructionL2: relativeReconL2,
      energyConservationRelError,
      solenoidalDivergenceL2,
      solenoidalDivergenceMax: maxDivSol,
      potentialCurlL2,
      potentialCurlMax: maxCurlPot,
      potentialGradL2Error: potGradL2Error,
      potentialGradRelError: potL2 > 1e-30 ? potGradL2Error / potL2 : 0.0,
      streamCurlL2Error,
      streamCurlRelError: solL2 > 1e-30 ? streamCurlL2Error / solL2 : 0.0,
      streamCoulombGaugeL2,
      spectralDivergenceSolMax,
      spectralCurlPotMax,
      spectralGradPhiPotMax,
      spectralCurlStreamSolMax,
      spectralGaugeAMax,
      totalHelicity: sumHelicity,
      meanEnstrophy: sumEnstrophy / domainVolume,
      isStrictlyOrthogonal,
      isEnergyConserving
    });
  }

  /**
   * Static helper: Performs one-shot decomposition.
   * 
   * @param {VectorField3D} vectorField
   * @param {Object} [options]
   * @returns {HelmholtzDecompositionResult}
   */
  static decompose(vectorField, options = {}) {
    const decomposer = new HelmholtzFFTDecomposer(options);
    return decomposer.decompose(vectorField, options);
  }

  /**
   * Static helper: Computes energy spectra of a vector field.
   * 
   * @param {VectorField3D} vectorField
   * @param {Object} [options]
   * @returns {HelmholtzEnergySpectrum}
   */
  static computeEnergySpectra(vectorField, options = {}) {
    const res = HelmholtzFFTDecomposer.decompose(vectorField, {
      ...options,
      computePotentials: false,
      computeVorticityHelicity: false,
      runDiagnostics: false
    });
    return res.energySpectrum;
  }

  /**
   * Factory: Generates an analytical pure Potential (irrotational) vector field.
   * v(x) = -\nabla \Phi(x) where \Phi(x) = \sum A_m \cos(k_m \cdot x + \phi_m)
   * 
   * @param {GridIndexer} grid
   * @param {Array<{ k: Array<number>, amplitude: number, phase?: number }>} modes
   * @returns {VectorField3D}
   */
  static createSyntheticPurePotential(grid, modes) {
    if (!(grid instanceof GridIndexer)) {
      throw new TypeError('createSyntheticPurePotential: grid must be an instance of GridIndexer.');
    }
    const total = grid.totalCells;
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    const { nx, ny, nz, xMin, yMin, zMin } = grid;
    const dx = grid.boxSize[0] / nx;
    const dy = grid.boxSize[1] / ny;
    const dz = grid.boxSize[2] / nz;

    for (let iz = 0; iz < nz; iz++) {
      const z = zMin + iz * dz;
      for (let iy = 0; iy < ny; iy++) {
        const y = yMin + iy * dy;
        for (let ix = 0; ix < nx; ix++) {
          const x = xMin + ix * dx;
          const idx = grid.index(ix, iy, iz);

          let vXTotal = 0.0;
          let vYTotal = 0.0;
          let vZTotal = 0.0;

          for (let m = 0; m < modes.length; m++) {
            const mode = modes[m];
            const kx = mode.k[0];
            const ky = mode.k[1];
            const kz = mode.k[2];
            const amp = mode.amplitude;
            const phase = mode.phase || 0.0;

            const kDotX = kx * x + ky * y + kz * z + phase;
            // Phi = amp * cos(k.x + phase) => -\nabla Phi = amp * sin(k.x + phase) * k
            const sinVal = Math.sin(kDotX);
            vXTotal += amp * kx * sinVal;
            vYTotal += amp * ky * sinVal;
            vZTotal += amp * kz * sinVal;
          }

          vx[idx] = vXTotal;
          vy[idx] = vYTotal;
          vz[idx] = vZTotal;
        }
      }
    }

    return new VectorField3D(grid, vx, vy, vz, { name: 'synthetic_pure_potential', unit: 'km/s' });
  }

  /**
   * Factory: Generates an analytical pure Solenoidal ABC (Arnold-Beltrami-Childress) helical vector field.
   * Exact divergence is zero: \nabla \cdot v = 0.
   * Vorticity is collinear: \nabla \times v = k_0 * v.
   * 
   * v_x = A * sin(k_0 * z) + C * cos(k_0 * y)
   * v_y = B * sin(k_0 * x) + A * cos(k_0 * z)
   * v_z = C * sin(k_0 * y) + B * cos(k_0 * x)
   * 
   * @param {GridIndexer} grid
   * @param {Object} params
   * @param {number} [params.A=100.0] Parameter A [km/s]
   * @param {number} [params.B=100.0] Parameter B [km/s]
   * @param {number} [params.C=100.0] Parameter C [km/s]
   * @param {number} [params.k0] Wavemode wavenumber [h/Mpc]. Defaults to 2*pi / L.
   * @returns {VectorField3D}
   */
  static createSyntheticPureSolenoidalABC(grid, params = {}) {
    if (!(grid instanceof GridIndexer)) {
      throw new TypeError('createSyntheticPureSolenoidalABC: grid must be an instance of GridIndexer.');
    }
    const {
      A = 100.0,
      B = 100.0,
      C = 100.0,
      k0 = (2.0 * Math.PI) / grid.boxSize[0]
    } = params;

    const total = grid.totalCells;
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    const { nx, ny, nz, xMin, yMin, zMin } = grid;
    const dx = grid.boxSize[0] / nx;
    const dy = grid.boxSize[1] / ny;
    const dz = grid.boxSize[2] / nz;

    for (let iz = 0; iz < nz; iz++) {
      const z = zMin + iz * dz;
      const sinK0Z = Math.sin(k0 * z);
      const cosK0Z = Math.cos(k0 * z);

      for (let iy = 0; iy < ny; iy++) {
        const y = yMin + iy * dy;
        const sinK0Y = Math.sin(k0 * y);
        const cosK0Y = Math.cos(k0 * y);

        for (let ix = 0; ix < nx; ix++) {
          const x = xMin + ix * dx;
          const sinK0X = Math.sin(k0 * x);
          const cosK0X = Math.cos(k0 * x);

          const idx = grid.index(ix, iy, iz);
          vx[idx] = A * sinK0Z + C * cosK0Y;
          vy[idx] = B * sinK0X + A * cosK0Z;
          vz[idx] = C * sinK0Y + B * cosK0X;
        }
      }
    }

    return new VectorField3D(grid, vx, vy, vz, { name: 'synthetic_pure_solenoidal_abc', unit: 'km/s' });
  }

  /**
   * Factory: Generates an analytical 3D Taylor-Green vortex solenoidal field.
   * 
   * v_x = U_0 * sin(k * x) * cos(k * y) * cos(k * z)
   * v_y = -U_0 * cos(k * x) * sin(k * y) * cos(k * z)
   * v_z = 0
   * 
   * @param {GridIndexer} grid
   * @param {Object} params
   * @param {number} [params.U0=200.0] Peak velocity [km/s]
   * @param {number} [params.k] Wavevector magnitude [h/Mpc]
   * @returns {VectorField3D}
   */
  static createSyntheticTaylorGreen(grid, params = {}) {
    if (!(grid instanceof GridIndexer)) {
      throw new TypeError('createSyntheticTaylorGreen: grid must be an instance of GridIndexer.');
    }
    const {
      U0 = 200.0,
      k = (2.0 * Math.PI) / grid.boxSize[0]
    } = params;

    const total = grid.totalCells;
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    const { nx, ny, nz, xMin, yMin, zMin } = grid;
    const dx = grid.boxSize[0] / nx;
    const dy = grid.boxSize[1] / ny;
    const dz = grid.boxSize[2] / nz;

    for (let iz = 0; iz < nz; iz++) {
      const z = zMin + iz * dz;
      const cosKZ = Math.cos(k * z);

      for (let iy = 0; iy < ny; iy++) {
        const y = yMin + iy * dy;
        const sinKY = Math.sin(k * y);
        const cosKY = Math.cos(k * y);

        for (let ix = 0; ix < nx; ix++) {
          const x = xMin + ix * dx;
          const sinKX = Math.sin(k * x);
          const cosKX = Math.cos(k * x);

          const idx = grid.index(ix, iy, iz);
          vx[idx] = U0 * sinKX * cosKY * cosKZ;
          vy[idx] = -U0 * cosKX * sinKY * cosKZ;
          vz[idx] = 0.0;
        }
      }
    }

    return new VectorField3D(grid, vx, vy, vz, { name: 'synthetic_taylor_green', unit: 'km/s' });
  }

  /**
   * Factory: Generates an exact superposition of Potential + Solenoidal ABC + Bulk Velocity modes.
   * 
   * @param {GridIndexer} grid
   * @param {Object} params
   * @param {Array<{ k: Array<number>, amplitude: number }>} [params.potentialModes]
   * @param {Object} [params.solenoidalParams]
   * @param {Array<number>} [params.bulkVelocity]
   * @returns {{ totalField: VectorField3D, purePotential: VectorField3D, pureSolenoidal: VectorField3D, bulkVelocity: VelocityVector }}
   */
  static createSyntheticSuperposition(grid, params = {}) {
    const potModes = params.potentialModes || [
      { k: [(2.0 * Math.PI) / grid.boxSize[0], (2.0 * Math.PI) / grid.boxSize[1], 0], amplitude: 150.0 }
    ];
    const solParams = params.solenoidalParams || { A: 120.0, B: 120.0, C: 120.0, k0: (2.0 * Math.PI) / grid.boxSize[0] };
    const v0Arr = params.bulkVelocity || [300.0, -150.0, 75.0];

    const purePot = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, potModes);
    const pureSol = HelmholtzFFTDecomposer.createSyntheticPureSolenoidalABC(grid, solParams);

    const total = grid.totalCells;
    const totVx = new Float64Array(total);
    const totVy = new Float64Array(total);
    const totVz = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      totVx[i] = purePot.vx[i] + pureSol.vx[i] + v0Arr[0];
      totVy[i] = purePot.vy[i] + pureSol.vy[i] + v0Arr[1];
      totVz[i] = purePot.vz[i] + pureSol.vz[i] + v0Arr[2];
    }

    const totalField = new VectorField3D(grid, totVx, totVy, totVz, { name: 'synthetic_superposition', unit: 'km/s' });
    const bulkVelocity = new VelocityVector(v0Arr[0], v0Arr[1], v0Arr[2]);

    return {
      totalField,
      purePotential: purePot,
      pureSolenoidal: pureSol,
      bulkVelocity
    };
  }

  /**
   * Applies isotropic spatial smoothing filter to a 3D vector field in Fourier space.
   * 
   * @param {VectorField3D} vectorField
   * @param {Object} filterOptions
   * @param {string} [filterOptions.filterType='gaussian'] HelmholtzFilterType
   * @param {number} [filterOptions.smoothingScale=5.0] Smoothing radius R in Mpc/h.
   * @param {number} [filterOptions.cutoffWavenumber] Sharp k-space cutoff wavenumber k_cut in h/Mpc.
   * @returns {VectorField3D} Filtered 3D vector field.
   */
  filterField(vectorField, filterOptions = {}) {
    if (!(vectorField instanceof VectorField3D)) {
      throw new TypeError('filterField: vectorField must be an instance of VectorField3D.');
    }

    const {
      filterType = HelmholtzFilterType.GAUSSIAN,
      smoothingScale = 5.0,
      cutoffWavenumber = (2.0 * Math.PI) / smoothingScale
    } = filterOptions;

    const grid = vectorField.grid;
    const { nx, ny, nz } = grid;
    const fft = this._getFFTEngine(nx, ny, nz);

    const fx = fft.forward(vectorField.vx);
    const fy = fft.forward(vectorField.vy);
    const fz = fft.forward(vectorField.vz);

    const Lx = grid.boxSize[0];
    const Ly = grid.boxSize[1];
    const Lz = grid.boxSize[2];

    const twoPiLx = (2.0 * Math.PI) / Lx;
    const twoPiLy = (2.0 * Math.PI) / Ly;
    const twoPiLz = (2.0 * Math.PI) / Lz;

    const kxArr = new Float64Array(nx);
    const kyArr = new Float64Array(ny);
    const kzArr = new Float64Array(nz);

    for (let ix = 0; ix < nx; ix++) {
      const mx = ix < (nx >> 1) ? ix : ix - nx;
      kxArr[ix] = mx * twoPiLx;
    }
    for (let iy = 0; iy < ny; iy++) {
      const my = iy < (ny >> 1) ? iy : iy - ny;
      kyArr[iy] = my * twoPiLy;
    }
    for (let iz = 0; iz < nz; iz++) {
      const mz = iz < (nz >> 1) ? iz : iz - nz;
      kzArr[iz] = mz * twoPiLz;
    }

    const R = smoothingScale;
    const halfRSq = 0.5 * R * R;

    for (let iz = 0; iz < nz; iz++) {
      const kz = kzArr[iz];
      const kzSq = kz * kz;
      const zOff = iz * nx * ny;

      for (let iy = 0; iy < ny; iy++) {
        const ky = kyArr[iy];
        const kyzSq = ky * ky + kzSq;
        const yOff = zOff + iy * nx;

        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix;
          if (ix === 0 && iy === 0 && iz === 0) continue;

          const kx = kxArr[ix];
          const kMag = Math.sqrt(kx * kx + kyzSq);

          let W = 1.0;
          if (filterType === HelmholtzFilterType.GAUSSIAN) {
            W = Math.exp(-kMag * kMag * halfRSq);
          } else if (filterType === HelmholtzFilterType.TOPHAT_REAL) {
            const kR = kMag * R;
            if (kR < 1e-4) {
              W = 1.0 - 0.1 * kR * kR;
            } else {
              W = 3.0 * (Math.sin(kR) - kR * Math.cos(kR)) / (kR * kR * kR);
            }
          } else if (filterType === HelmholtzFilterType.SHARP_K) {
            W = kMag <= cutoffWavenumber ? 1.0 : 0.0;
          }

          fx.real[idx] *= W;
          fx.imag[idx] *= W;
          fy.real[idx] *= W;
          fy.imag[idx] *= W;
          fz.real[idx] *= W;
          fz.imag[idx] *= W;
        }
      }
    }

    const outVx = fft.inverse(fx.real, fx.imag);
    const outVy = fft.inverse(fy.real, fy.imag);
    const outVz = fft.inverse(fz.real, fz.imag);

    return new VectorField3D(grid, outVx, outVy, outVz, {
      name: `${vectorField.name}_filtered_${filterType}`,
      unit: vectorField.unit
    });
  }

  /**
   * Decomposes vector field into cosmic web multiscale spatial frequency bands.
   * 
   * @param {VectorField3D} vectorField
   * @param {Array<{ band: string, rMin: number, rMax: number }>} [customBands]
   * @returns {Map<string, { bandField: VectorField3D, potField: VectorField3D, solField: VectorField3D, energyFraction: number }>}
   */
  decomposeScaleBands(vectorField, customBands = null) {
    const bands = customBands || [
      { band: ScaleBandType.CLUSTER, rMin: 0.0, rMax: 5.0 },
      { band: ScaleBandType.FILAMENT, rMin: 5.0, rMax: 20.0 },
      { band: ScaleBandType.SUPERCLUSTER, rMin: 20.0, rMax: Infinity }
    ];

    const result = new Map();
    const fullDecomp = this.decompose(vectorField);
    const totalE = fullDecomp.energySpectrum.totalEnergySum;

    for (let b = 0; b < bands.length; b++) {
      const { band, rMin, rMax } = bands[b];
      const kMax = rMin > 0 ? (2.0 * Math.PI) / rMin : Infinity;
      const kMin = Number.isFinite(rMax) ? (2.0 * Math.PI) / rMax : 0.0;

      const grid = vectorField.grid;
      const { nx, ny, nz } = grid;
      const fft = this._getFFTEngine(nx, ny, nz);

      const fx = fft.forward(vectorField.vx);
      const fy = fft.forward(vectorField.vy);
      const fz = fft.forward(vectorField.vz);

      const Lx = grid.boxSize[0];
      const Ly = grid.boxSize[1];
      const Lz = grid.boxSize[2];

      const twoPiLx = (2.0 * Math.PI) / Lx;
      const twoPiLy = (2.0 * Math.PI) / Ly;
      const twoPiLz = (2.0 * Math.PI) / Lz;

      for (let iz = 0; iz < nz; iz++) {
        const mz = iz < (nz >> 1) ? iz : iz - nz;
        const kz = mz * twoPiLz;
        const kzSq = kz * kz;
        const zOff = iz * nx * ny;

        for (let iy = 0; iy < ny; iy++) {
          const my = iy < (ny >> 1) ? iy : iy - ny;
          const ky = my * twoPiLy;
          const kyzSq = ky * ky + kzSq;
          const yOff = zOff + iy * nx;

          for (let ix = 0; ix < nx; ix++) {
            const idx = yOff + ix;
            if (ix === 0 && iy === 0 && iz === 0) {
              if (kMin > 0) {
                fx.real[0] = 0.0;
                fx.imag[0] = 0.0;
                fy.real[0] = 0.0;
                fy.imag[0] = 0.0;
                fz.real[0] = 0.0;
                fz.imag[0] = 0.0;
              }
              continue;
            }

            const mx = ix < (nx >> 1) ? ix : ix - nx;
            const kx = mx * twoPiLx;
            const kMag = Math.sqrt(kx * kx + kyzSq);

            const inBand = kMag >= kMin && kMag < kMax;
            const W = inBand ? 1.0 : 0.0;

            fx.real[idx] *= W;
            fx.imag[idx] *= W;
            fy.real[idx] *= W;
            fy.imag[idx] *= W;
            fz.real[idx] *= W;
            fz.imag[idx] *= W;
          }
        }
      }

      const bandVx = fft.inverse(fx.real, fx.imag);
      const bandVy = fft.inverse(fy.real, fy.imag);
      const bandVz = fft.inverse(fz.real, fz.imag);

      const bandField = new VectorField3D(grid, bandVx, bandVy, bandVz, {
        name: `band_${band}_velocity`,
        unit: vectorField.unit
      });

      const bandDecomp = this.decompose(bandField, { runDiagnostics: false });
      const bandE = bandDecomp.energySpectrum.totalEnergySum;
      const fraction = totalE > 0 ? bandE / totalE : 0.0;

      result.set(band, {
        bandField,
        potField: bandDecomp.vPot,
        solField: bandDecomp.vSol,
        energyFraction: fraction
      });
    }

    return result;
  }
}
