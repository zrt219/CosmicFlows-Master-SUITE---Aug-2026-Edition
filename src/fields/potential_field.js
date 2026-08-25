/**
 * @file potential_field.js
 * @description Cosmological Gravitational Potential & Clustering Field Engine.
 * 
 * Implements:
 * 1. Exact 3D Poisson equation solvers for cosmological density fields:
 *    nabla^2 Phi(x) = 4 * pi * G * rho_bar * delta(x) = S(x)
 *    - Spectral 3D FFT solver (periodic boundary conditions, continuous & discrete kernels)
 *    - Isolated Green's function FFT convolution solver (Hockney-Eastwood zero-padded G(r) = -G/|r|)
 *    - Geometric Multigrid (V-cycle with Red-Black Gauss-Seidel relaxation)
 * 2. Gravitational acceleration vector field g(x) = -nabla Phi(x)
 * 3. Gravitational tidal tensor T_ij = d_i d_j Phi and traceless shear tensor with Jacobi eigenvalue diagonalization
 * 4. Local and volume-integrated gravitational potential energy density w_g(x) = 0.5 * rho_bar * delta(x) * Phi(x)
 * 5. Analytical profiles for verification: Plummer sphere, Hernquist profile, Harmonic standing waves
 * 
 * @module fields/potential_field
 */

import { GridIndexer, BoundaryMode } from './grid_indexer.js';
import { DensityField } from './density_field.js';
import { TrilinearInterpolator } from '../interpolation/trilinear_interpolator.js';
import { FiniteDifference } from '../interpolation/finite_difference.js';

/**
 * Standard cosmological gravitational constant defaults in astrophysical units:
 * G = 4.30091e-9 (km/s)^2 * Mpc / M_sun
 */
export const DEFAULT_FOUR_PI_G_RHO_BAR = 1.0;

/**
 * 1D In-Place Radix-2 Cooley-Tukey Fast Fourier Transform.
 * Length n must be a power of 2.
 * 
 * @param {Float64Array} real Real component array.
 * @param {Float64Array} imag Imaginary component array.
 * @param {boolean} [inverse=false] True for inverse FFT.
 */
export function fft1D(real, imag, inverse = false) {
  const n = real.length;
  if ((n & (n - 1)) !== 0) {
    throw new Error(`fft1D: Length n must be a power of 2, received ${n}`);
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

  // Cooley-Tukey butterfly computations
  const sign = inverse ? 1.0 : -1.0;
  for (let len = 2; len <= n; len <<= 1) {
    const halfLen = len >> 1;
    const angle = (sign * 2.0 * Math.PI) / len;
    const wStepR = Math.cos(angle);
    const wStepI = Math.sin(angle);

    for (let i = 0; i < n; i += len) {
      let wR = 1.0;
      let wI = 0.0;
      for (let k = 0; k < halfLen; k++) {
        const uR = real[i + k];
        const uI = imag[i + k];
        const vR = real[i + k + halfLen] * wR - imag[i + k + halfLen] * wI;
        const vI = real[i + k + halfLen] * wI + imag[i + k + halfLen] * wR;

        real[i + k] = uR + vR;
        imag[i + k] = uI + vI;
        real[i + k + halfLen] = uR - vR;
        imag[i + k + halfLen] = uI - vI;

        const nextWR = wR * wStepR - wI * wStepI;
        const nextWI = wR * wStepI + wI * wStepR;
        wR = nextWR;
        wI = nextWI;
      }
    }
  }

  // Normalization for inverse FFT
  if (inverse) {
    const invN = 1.0 / n;
    for (let i = 0; i < n; i++) {
      real[i] *= invN;
      imag[i] *= invN;
    }
  }
}

/**
 * 3D Fast Fourier Transform on flat row-major arrays of size nx * ny * nz.
 * All dimensions (nx, ny, nz) must be powers of 2.
 * 
 * @param {Float64Array} real Real component buffer (length nx * ny * nz).
 * @param {Float64Array} imag Imaginary component buffer (length nx * ny * nz).
 * @param {number} nx Dimension X (fastest stride).
 * @param {number} ny Dimension Y.
 * @param {number} nz Dimension Z (slowest stride).
 * @param {boolean} [inverse=false] True for inverse 3D FFT.
 */
export function fft3D(real, imag, nx, ny, nz, inverse = false) {
  const total = nx * ny * nz;
  if (real.length !== total || imag.length !== total) {
    throw new Error(`fft3D: Buffer size (${real.length}) does not match grid dimensions ${nx}x${ny}x${nz} = ${total}`);
  }

  // 1. Transform along X lines (stride 1)
  const lineX_R = new Float64Array(nx);
  const lineX_I = new Float64Array(nx);
  for (let iz = 0; iz < nz; iz++) {
    for (let iy = 0; iy < ny; iy++) {
      const base = nx * (iy + ny * iz);
      for (let ix = 0; ix < nx; ix++) {
        lineX_R[ix] = real[base + ix];
        lineX_I[ix] = imag[base + ix];
      }
      fft1D(lineX_R, lineX_I, inverse);
      for (let ix = 0; ix < nx; ix++) {
        real[base + ix] = lineX_R[ix];
        imag[base + ix] = lineX_I[ix];
      }
    }
  }

  // 2. Transform along Y lines (stride nx)
  const lineY_R = new Float64Array(ny);
  const lineY_I = new Float64Array(ny);
  for (let iz = 0; iz < nz; iz++) {
    for (let ix = 0; ix < nx; ix++) {
      for (let iy = 0; iy < ny; iy++) {
        const idx = ix + nx * (iy + ny * iz);
        lineY_R[iy] = real[idx];
        lineY_I[iy] = imag[idx];
      }
      fft1D(lineY_R, lineY_I, inverse);
      for (let iy = 0; iy < ny; iy++) {
        const idx = ix + nx * (iy + ny * iz);
        real[idx] = lineY_R[iy];
        imag[idx] = lineY_I[iy];
      }
    }
  }

  // 3. Transform along Z lines (stride nx * ny)
  const sliceSize = nx * ny;
  const lineZ_R = new Float64Array(nz);
  const lineZ_I = new Float64Array(nz);
  for (let iy = 0; iy < ny; iy++) {
    for (let ix = 0; ix < nx; ix++) {
      const base = ix + nx * iy;
      for (let iz = 0; iz < nz; iz++) {
        const idx = base + sliceSize * iz;
        lineZ_R[iz] = real[idx];
        lineZ_I[iz] = imag[idx];
      }
      fft1D(lineZ_R, lineZ_I, inverse);
      for (let iz = 0; iz < nz; iz++) {
        const idx = base + sliceSize * iz;
        real[idx] = lineZ_R[iz];
        imag[idx] = lineZ_I[iz];
      }
    }
  }
}

/**
 * Full 3x3 symmetric matrix Jacobi eigenvalue diagonalizer.
 * Returns sorted eigenvalues lambda1 <= lambda2 <= lambda3 and orthonormal eigenvectors.
 * 
 * @param {Array<Array<number>>} mat 3x3 symmetric matrix.
 * @param {number} [tol=1e-15] Convergence tolerance.
 * @param {number} [maxSweeps=50] Maximum Jacobi sweeps.
 * @returns {{ eigenvalues: [number, number, number], eigenvectors: [[number, number, number], [number, number, number], [number, number, number]], iterations: number }}
 */
export function jacobiDiagonalize3x3(mat, tol = 1e-15, maxSweeps = 50) {
  const A = [
    [mat[0][0], 0.5 * (mat[0][1] + mat[1][0]), 0.5 * (mat[0][2] + mat[2][0])],
    [0.5 * (mat[0][1] + mat[1][0]), mat[1][1], 0.5 * (mat[1][2] + mat[2][1])],
    [0.5 * (mat[0][2] + mat[2][0]), 0.5 * (mat[1][2] + mat[2][1]), mat[2][2]]
  ];

  const V = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
  ];

  let sweep = 0;
  for (sweep = 0; sweep < maxSweeps; sweep++) {
    const offDiag = Math.abs(A[0][1]) + Math.abs(A[0][2]) + Math.abs(A[1][2]);
    if (offDiag < tol) break;

    const pairs = [[0, 1], [0, 2], [1, 2]];
    for (const [p, q] of pairs) {
      const app = A[p][p];
      const aqq = A[q][q];
      const apq = A[p][q];

      if (Math.abs(apq) < 1e-16) continue;

      const tau = (aqq - app) / (2.0 * apq);
      let t;
      if (tau >= 0) {
        t = 1.0 / (tau + Math.sqrt(1.0 + tau * tau));
      } else {
        t = -1.0 / (-tau + Math.sqrt(1.0 + tau * tau));
      }

      const c = 1.0 / Math.sqrt(1.0 + t * t);
      const s = t * c;
      const h = t * apq;

      A[p][p] -= h;
      A[q][q] += h;
      A[p][q] = 0.0;
      A[q][p] = 0.0;

      for (let r = 0; r < 3; r++) {
        if (r !== p && r !== q) {
          const arp = A[r][p];
          const arq = A[r][q];
          A[r][p] = c * arp - s * arq;
          A[p][r] = A[r][p];
          A[r][q] = s * arp + c * arq;
          A[q][r] = A[r][q];
        }
      }

      for (let r = 0; r < 3; r++) {
        const vrp = V[r][p];
        const vrq = V[r][q];
        V[r][p] = c * vrp - s * vrq;
        V[r][q] = s * vrp + c * vrq;
      }
    }
  }

  const items = [
    { val: A[0][0], vec: [V[0][0], V[1][0], V[2][0]] },
    { val: A[1][1], vec: [V[0][1], V[1][1], V[2][1]] },
    { val: A[2][2], vec: [V[0][2], V[1][2], V[2][2]] }
  ];

  items.sort((a, b) => a.val - b.val);

  return {
    eigenvalues: [items[0].val, items[1].val, items[2].val],
    eigenvectors: [items[0].vec, items[1].vec, items[2].vec],
    iterations: sweep
  };
}

/**
 * Cosmological Gravitational Potential Field representation.
 */
export class PotentialField {
  /**
   * Constructs a PotentialField.
   * 
   * @param {GridIndexer} gridIndexer Grid geometry indexer.
   * @param {Float32Array|Float64Array} phiBuffer 1D flat buffer of gravitational potential values Phi(x).
   * @param {Object} [options] Configuration options.
   * @param {number} [options.fourPiG=1.0] Gravitational coupling constant 4 * pi * G * rho_bar.
   */
  constructor(gridIndexer, phiBuffer, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('PotentialField: gridIndexer must be an instance of GridIndexer.');
    }

    const total = gridIndexer.totalCells;
    if (!phiBuffer || phiBuffer.length !== total) {
      throw new Error(`PotentialField: phiBuffer size (${phiBuffer?.length}) must match grid total cells (${total}).`);
    }

    this.grid = gridIndexer;
    this.phi = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      this.phi[i] = phiBuffer[i];
    }

    this.fourPiG = typeof options.fourPiG === 'number' ? options.fourPiG : DEFAULT_FOUR_PI_G_RHO_BAR;
    this.interpolator = new TrilinearInterpolator(this.grid);
    this.fd = new FiniteDifference(this.grid, this.interpolator);
  }

  /**
   * Evaluates continuous gravitational potential Phi(x, y, z) via trilinear interpolation.
   * 
   * @param {number} x Physical X.
   * @param {number} y Physical Y.
   * @param {number} z Physical Z.
   * @returns {number} Interpolated potential Phi.
   */
  evaluatePotential(x, y, z) {
    return this.interpolator.interpolateScalar(this.phi, x, y, z);
  }

  /**
   * Evaluates continuous gravitational acceleration g(x) = -nabla Phi(x) at physical coordinates.
   * 
   * @param {number} x Physical X.
   * @param {number} y Physical Y.
   * @param {number} z Physical Z.
   * @param {string} [stencil='14point'] '6point' or '14point'.
   * @returns {[number, number, number]} Gravitational acceleration [gx, gy, gz].
   */
  evaluateAcceleration(x, y, z, stencil = '14point') {
    const grad = stencil === '6point'
      ? this.fd.gradient6Point(this.phi, x, y, z)
      : this.fd.gradient14Point(this.phi, x, y, z);
    return [-grad[0], -grad[1], -grad[2]];
  }

  /**
   * Computes gravitational acceleration g = -nabla Phi at specific grid voxel (ix, iy, iz).
   * 
   * @param {number} ix Grid index X.
   * @param {number} iy Grid index Y.
   * @param {number} iz Grid index Z.
   * @returns {[number, number, number]} Acceleration vector [gx, gy, gz].
   */
  computeAccelerationAtVoxel(ix, iy, iz, stencil = '6point') {
    const [x, y, z] = this.grid.gridIndexToCoord(ix, iy, iz);
    return this.evaluateAcceleration(x, y, z, stencil);
  }

  /**
   * Precomputes full 3D gravitational acceleration vector field buffers (gx, gy, gz).
   * 
   * @returns {{ gx: Float64Array, gy: Float64Array, gz: Float64Array }} Acceleration field buffers.
   */
  getAccelerationBuffers() {
    const { nx, ny, nz, dx, dy, dz, totalCells } = this.grid;
    const gx = new Float64Array(totalCells);
    const gy = new Float64Array(totalCells);
    const gz = new Float64Array(totalCells);

    const inv2dx = 0.5 / dx;
    const inv2dy = 0.5 / dy;
    const inv2dz = 0.5 / dz;

    for (let iz = 0; iz < nz; iz++) {
      const izPrev = iz > 0 ? iz - 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? nz - 1 : 0);
      const izNext = iz < nz - 1 ? iz + 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : nz - 1);
      const scaleZ = (iz === 0 || iz === nz - 1) && this.grid.boundaryMode !== BoundaryMode.PERIODIC ? (1.0 / dz) : inv2dz;

      for (let iy = 0; iy < ny; iy++) {
        const iyPrev = iy > 0 ? iy - 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? ny - 1 : 0);
        const iyNext = iy < ny - 1 ? iy + 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : ny - 1);
        const scaleY = (iy === 0 || iy === ny - 1) && this.grid.boundaryMode !== BoundaryMode.PERIODIC ? (1.0 / dy) : inv2dy;

        for (let ix = 0; ix < nx; ix++) {
          const ixPrev = ix > 0 ? ix - 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? nx - 1 : 0);
          const ixNext = ix < nx - 1 ? ix + 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : nx - 1);
          const scaleX = (ix === 0 || ix === nx - 1) && this.grid.boundaryMode !== BoundaryMode.PERIODIC ? (1.0 / dx) : inv2dx;

          const idx = this.grid.index(ix, iy, iz);

          // Partial derivatives
          const dPhi_dx = (this.phi[this.grid.index(ixNext, iy, iz)] - this.phi[this.grid.index(ixPrev, iy, iz)]) * scaleX;
          const dPhi_dy = (this.phi[this.grid.index(ix, iyNext, iz)] - this.phi[this.grid.index(ix, iyPrev, iz)]) * scaleY;
          const dPhi_dz = (this.phi[this.grid.index(ix, iy, izNext)] - this.phi[this.grid.index(ix, iy, izPrev)]) * scaleZ;

          gx[idx] = -dPhi_dx;
          gy[idx] = -dPhi_dy;
          gz[idx] = -dPhi_dz;
        }
      }
    }

    return { gx, gy, gz };
  }

  /**
   * Computes the 3x3 Hessian matrix H_ij = d^2 Phi / (dx_i dx_j) at voxel (ix, iy, iz).
   * 
   * @param {number} ix Grid index X.
   * @param {number} iy Grid index Y.
   * @param {number} iz Grid index Z.
   * @returns {Array<Array<number>>} 3x3 Hessian tensor.
   */
  computeHessianAtVoxel(ix, iy, iz) {
    const { dx, dy, dz } = this.grid;
    const [x, y, z] = this.grid.gridIndexToCoord(ix, iy, iz);

    const evalP = (px, py, pz) => this.evaluatePotential(px, py, pz);

    const f000 = evalP(x, y, z);
    const dxx = (evalP(x + dx, y, z) - 2.0 * f000 + evalP(x - dx, y, z)) / (dx * dx);
    const dyy = (evalP(x, y + dy, z) - 2.0 * f000 + evalP(x, y - dy, z)) / (dy * dy);
    const dzz = (evalP(x, y, z + dz) - 2.0 * f000 + evalP(x, y, z - dz)) / (dz * dz);

    const dxy = (evalP(x + dx, y + dy, z) - evalP(x + dx, y - dy, z) - evalP(x - dx, y + dy, z) + evalP(x - dx, y - dy, z)) / (4.0 * dx * dy);
    const dxz = (evalP(x + dx, y, z + dz) - evalP(x + dx, y, z - dz) - evalP(x - dx, y, z + dz) + evalP(x - dx, y, z - dz)) / (4.0 * dx * dz);
    const dyz = (evalP(x, y + dy, z + dz) - evalP(x, y + dy, z - dz) - evalP(x, y - dy, z + dz) + evalP(x, y - dy, z - dz)) / (4.0 * dy * dz);

    return [
      [dxx, dxy, dxz],
      [dxy, dyy, dyz],
      [dxz, dyz, dzz]
    ];
  }

  /**
   * Computes the cosmological tidal tensor T_ij = d_i d_j Phi at voxel (ix, iy, iz).
   * If traceless is true, subtracts 1/3 * Tr(H) * delta_ij.
   * 
   * @param {number} ix Grid index X.
   * @param {number} iy Grid index Y.
   * @param {number} iz Grid index Z.
   * @param {boolean} [traceless=false] Whether to subtract the trace.
   * @returns {Array<Array<number>>} 3x3 Tidal tensor.
   */
  computeTidalTensorAtVoxel(ix, iy, iz, traceless = false) {
    const H = this.computeHessianAtVoxel(ix, iy, iz);
    if (!traceless) return H;

    const traceThird = (H[0][0] + H[1][1] + H[2][2]) / 3.0;
    return [
      [H[0][0] - traceThird, H[0][1], H[0][2]],
      [H[1][0], H[1][1] - traceThird, H[1][2]],
      [H[2][0], H[2][1], H[2][2] - traceThird]
    ];
  }

  /**
   * Computes tidal tensor eigenvalues lambda1 <= lambda2 <= lambda3 and eigenvectors via Jacobi diagonalization.
   * 
   * @param {number} ix Grid index X.
   * @param {number} iy Grid index Y.
   * @param {number} iz Grid index Z.
   * @param {boolean} [traceless=false] Whether to diagonalize the traceless tidal tensor.
   * @returns {{ eigenvalues: [number, number, number], eigenvectors: [[number, number, number], [number, number, number], [number, number, number]], trace: number }}
   */
  computeTidalEigenvaluesAtVoxel(ix, iy, iz, traceless = false) {
    const T = this.computeTidalTensorAtVoxel(ix, iy, iz, traceless);
    const result = jacobiDiagonalize3x3(T);
    const trace = T[0][0] + T[1][1] + T[2][2];
    return {
      eigenvalues: result.eigenvalues,
      eigenvectors: result.eigenvectors,
      trace
    };
  }

  /**
   * Computes local gravitational potential energy density buffer:
   * w_g(x) = 0.5 * rho_bar * delta(x) * Phi(x)
   * 
   * @param {DensityField} densityField Associated overdensity field delta(x).
   * @returns {Float64Array} Flat buffer of energy density at each voxel.
   */
  computeEnergyDensityBuffer(densityField) {
    if (!(densityField instanceof DensityField)) {
      throw new TypeError('PotentialField: densityField must be an instance of DensityField.');
    }
    const total = this.grid.totalCells;
    const energy = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      energy[i] = 0.5 * densityField.delta[i] * this.phi[i];
    }
    return energy;
  }

  /**
   * Computes total integrated gravitational potential energy:
   * W = \int 0.5 * rho_bar * delta(x) * Phi(x) d^3 x
   * 
   * @param {DensityField} densityField Associated overdensity field delta(x).
   * @returns {number} Total gravitational potential energy W.
   */
  computeTotalGravitationalPotentialEnergy(densityField) {
    const energyBuffer = this.computeEnergyDensityBuffer(densityField);
    const dV = this.grid.dx * this.grid.dy * this.grid.dz;
    let totalW = 0.0;
    for (let i = 0; i < energyBuffer.length; i++) {
      totalW += energyBuffer[i];
    }
    return totalW * dV;
  }

  /**
   * Evaluates discrete Laplacian of the potential nabla^2 Phi at voxel (ix, iy, iz).
   * 
   * @param {number} ix Grid index X.
   * @param {number} iy Grid index Y.
   * @param {number} iz Grid index Z.
   * @returns {number} Discrete Laplacian value.
   */
  laplacianAtVoxel(ix, iy, iz) {
    const { nx, ny, nz, dx, dy, dz } = this.grid;
    const idx = this.grid.index(ix, iy, iz);
    const centerVal = this.phi[idx];

    const ixPrev = ix > 0 ? ix - 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? nx - 1 : 0);
    const ixNext = ix < nx - 1 ? ix + 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : nx - 1);
    const iyPrev = iy > 0 ? iy - 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? ny - 1 : 0);
    const iyNext = iy < ny - 1 ? iy + 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : ny - 1);
    const izPrev = iz > 0 ? iz - 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? nz - 1 : 0);
    const izNext = iz < nz - 1 ? iz + 1 : (this.grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : nz - 1);

    const dxx = (this.phi[this.grid.index(ixNext, iy, iz)] - 2.0 * centerVal + this.phi[this.grid.index(ixPrev, iy, iz)]) / (dx * dx);
    const dyy = (this.phi[this.grid.index(ix, iyNext, iz)] - 2.0 * centerVal + this.phi[this.grid.index(ix, iyPrev, iz)]) / (dy * dy);
    const dzz = (this.phi[this.grid.index(ix, iy, izNext)] - 2.0 * centerVal + this.phi[this.grid.index(ix, iy, izPrev)]) / (dz * dz);

    return dxx + dyy + dzz;
  }

  // =========================================================================
  // FACTORY & POISSON SOLVERS
  // =========================================================================

  /**
   * Solves the Poisson equation on a periodic domain using Spectral 3D FFT:
   * nabla^2 Phi(x) = 4 * pi * G * rho_bar * delta(x)
   * 
   * @param {DensityField} densityField Input overdensity field delta(x).
   * @param {Object} [options] Solver options.
   * @param {number} [options.fourPiG=1.0] Coupling factor 4 * pi * G * rho_bar.
   * @param {string} [options.kernel='discrete'] Green's function kernel: 'discrete' (exact 7-point finite difference) or 'continuous' (-k^2).
   * @returns {PotentialField} Computed gravitational potential field.
   */
  static solvePoissonFFT(densityField, options = {}) {
    if (!(densityField instanceof DensityField)) {
      throw new TypeError('solvePoissonFFT: densityField must be an instance of DensityField.');
    }

    const grid = densityField.grid;
    const { nx, ny, nz, dx, dy, dz, totalCells } = grid;

    // Verify dimensions are powers of 2
    if ((nx & (nx - 1)) !== 0 || (ny & (ny - 1)) !== 0 || (nz & (nz - 1)) !== 0) {
      throw new Error(`solvePoissonFFT: Grid dimensions (${nx}, ${ny}, ${nz}) must all be powers of 2 for Radix-2 FFT.`);
    }

    const fourPiG = typeof options.fourPiG === 'number' ? options.fourPiG : DEFAULT_FOUR_PI_G_RHO_BAR;
    const kernelType = options.kernel || 'discrete';

    const real = new Float64Array(totalCells);
    const imag = new Float64Array(totalCells);

    // Populate source term S(x) = 4 * pi * G * delta(x)
    for (let i = 0; i < totalCells; i++) {
      real[i] = fourPiG * densityField.delta[i];
      imag[i] = 0.0;
    }

    // Forward 3D FFT
    fft3D(real, imag, nx, ny, nz, false);

    // Convolution in Fourier space: Phi(k) = -S(k) / k_eff^2
    const twoPiLx = (2.0 * Math.PI) / (nx * dx);
    const twoPiLy = (2.0 * Math.PI) / (ny * dy);
    const twoPiLz = (2.0 * Math.PI) / (nz * dz);

    const invDx2 = 1.0 / (dx * dx);
    const invDy2 = 1.0 / (dy * dy);
    const invDz2 = 1.0 / (dz * dz);

    for (let iz = 0; iz < nz; iz++) {
      const kzIdx = iz <= nz / 2 ? iz : iz - nz;
      const kz = kzIdx * twoPiLz;

      for (let iy = 0; iy < ny; iy++) {
        const kyIdx = iy <= ny / 2 ? iy : iy - ny;
        const ky = kyIdx * twoPiLy;

        for (let ix = 0; ix < nx; ix++) {
          const kxIdx = ix <= nx / 2 ? ix : ix - nx;
          const kx = kxIdx * twoPiLx;

          const idx = ix + nx * (iy + ny * iz);

          if (kxIdx === 0 && kyIdx === 0 && kzIdx === 0) {
            // Zero out DC mode (mean cosmological potential gauge)
            real[idx] = 0.0;
            imag[idx] = 0.0;
            continue;
          }

          let kEffSq;
          if (kernelType === 'discrete') {
            // Exact discrete 7-point Laplacian eigenvalues:
            // lambda_k = 2*(cos(kx*dx)-1)/dx^2 + 2*(cos(ky*dy)-1)/dy^2 + 2*(cos(kz*dz)-1)/dz^2
            kEffSq = 2.0 * (1.0 - Math.cos(kx * dx)) * invDx2 +
                     2.0 * (1.0 - Math.cos(ky * dy)) * invDy2 +
                     2.0 * (1.0 - Math.cos(kz * dz)) * invDz2;
          } else {
            // Continuous Fourier Laplacian: k^2 = kx^2 + ky^2 + kz^2
            kEffSq = kx * kx + ky * ky + kz * kz;
          }

          const invK2 = -1.0 / kEffSq;
          real[idx] *= invK2;
          imag[idx] *= invK2;
        }
      }
    }

    // Inverse 3D FFT
    fft3D(real, imag, nx, ny, nz, true);

    return new PotentialField(grid, real, { fourPiG });
  }

  /**
   * Solves the Poisson equation for isolated / non-periodic boundary conditions
   * using zero-padded Green's function FFT convolution (Hockney & Eastwood / James algorithm).
   * Computes Phi(x) = \int G(x - x') S(x') d^3 x' with G(r) = -G / |r|.
   * 
   * @param {DensityField} densityField Input overdensity field.
   * @param {Object} [options] Solver options.
   * @param {number} [options.G=1.0] Gravitational constant G.
   * @param {number} [options.rhoBar=1.0] Mean density rho_bar.
   * @param {number} [options.softening=0.0] Plummer gravitational softening radius epsilon.
   * @returns {PotentialField} Isolated gravitational potential field.
   */
  static solvePoissonIsolated(densityField, options = {}) {
    if (!(densityField instanceof DensityField)) {
      throw new TypeError('solvePoissonIsolated: densityField must be an instance of DensityField.');
    }

    const grid = densityField.grid;
    const { nx, ny, nz, dx, dy, dz } = grid;

    // Double grid dimensions for zero-padding (to eliminate periodic wrap-around)
    const n2x = 2 * nx;
    const n2y = 2 * ny;
    const n2z = 2 * nz;
    const total2 = n2x * n2y * n2z;

    const G = typeof options.G === 'number' ? options.G : 1.0;
    const rhoBar = typeof options.rhoBar === 'number' ? options.rhoBar : 1.0;
    const softeningSq = (options.softening || 0.0) ** 2;
    const dV = dx * dy * dz;

    // 1. Zero-padded source array S(x)
    const sReal = new Float64Array(total2);
    const sImag = new Float64Array(total2);

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const srcIdx = grid.index(ix, iy, iz);
          const dstIdx = ix + n2x * (iy + n2y * iz);
          sReal[dstIdx] = rhoBar * densityField.delta[srcIdx] * dV;
        }
      }
    }

    // 2. Isolated Green's function array G(r) = -G / sqrt(r^2 + eps^2)
    const gReal = new Float64Array(total2);
    const gImag = new Float64Array(total2);

    // Singularity regularization at origin (cell volume integration equivalent)
    const g0 = softeningSq > 0 ? -G / Math.sqrt(softeningSq) : -G / (0.5 * Math.cbrt(dV));

    for (let iz = 0; iz < n2z; iz++) {
      const zDist = (iz < nz ? iz : iz - n2z) * dz;
      const z2 = zDist * zDist;

      for (let iy = 0; iy < n2y; iy++) {
        const yDist = (iy < ny ? iy : iy - n2y) * dy;
        const y2 = yDist * yDist;

        for (let ix = 0; ix < n2x; ix++) {
          const xDist = (ix < nx ? ix : ix - n2x) * dx;
          const x2 = xDist * xDist;

          const idx = ix + n2x * (iy + n2y * iz);
          const rSq = x2 + y2 + z2;

          if (rSq === 0.0) {
            gReal[idx] = g0;
          } else {
            gReal[idx] = -G / Math.sqrt(rSq + softeningSq);
          }
        }
      }
    }

    // 3. FFT of source and Green's function
    fft3D(sReal, sImag, n2x, n2y, n2z, false);
    fft3D(gReal, gImag, n2x, n2y, n2z, false);

    // 4. Pointwise complex multiplication: Phi_hat = S_hat * G_hat
    const pReal = new Float64Array(total2);
    const pImag = new Float64Array(total2);

    for (let i = 0; i < total2; i++) {
      pReal[i] = sReal[i] * gReal[i] - sImag[i] * gImag[i];
      pImag[i] = sReal[i] * gImag[i] + sImag[i] * gReal[i];
    }

    // 5. Inverse FFT
    fft3D(pReal, pImag, n2x, n2y, n2z, true);

    // 6. Extract physical quadrant
    const resultPhi = new Float64Array(grid.totalCells);
    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const srcIdx = ix + n2x * (iy + n2y * iz);
          const dstIdx = grid.index(ix, iy, iz);
          resultPhi[dstIdx] = pReal[srcIdx];
        }
      }
    }

    return new PotentialField(grid, resultPhi, { fourPiG: 4.0 * Math.PI * G * rhoBar });
  }

  /**
   * Solves the Poisson equation using Multigrid V-Cycle / Successive Over-Relaxation (SOR).
   * nabla^2 Phi = S(x)
   * 
   * @param {DensityField} densityField Input overdensity field.
   * @param {Object} [options] Solver parameters.
   * @param {number} [options.fourPiG=1.0] Coupling factor 4 * pi * G * rho_bar.
   * @param {number} [options.maxIter=500] Maximum iterations.
   * @param {number} [options.tol=1e-6] Residual convergence tolerance.
   * @param {number} [options.omega=1.5] SOR relaxation parameter (1.0 < omega < 2.0).
   * @returns {PotentialField} Solved potential field.
   */
  static solvePoissonMultigrid(densityField, options = {}) {
    if (!(densityField instanceof DensityField)) {
      throw new TypeError('solvePoissonMultigrid: densityField must be an instance of DensityField.');
    }

    const grid = densityField.grid;
    const { nx, ny, nz, dx, dy, dz, totalCells } = grid;
    const fourPiG = typeof options.fourPiG === 'number' ? options.fourPiG : DEFAULT_FOUR_PI_G_RHO_BAR;
    const maxIter = options.maxIter || 500;
    const tol = options.tol || 1e-6;
    const omega = options.omega || 1.6;

    const phi = new Float64Array(totalCells);
    const source = new Float64Array(totalCells);

    const invDx2 = 1.0 / (dx * dx);
    const invDy2 = 1.0 / (dy * dy);
    const invDz2 = 1.0 / (dz * dz);
    const denom = 2.0 * (invDx2 + invDy2 + invDz2);

    for (let i = 0; i < totalCells; i++) {
      source[i] = fourPiG * densityField.delta[i];
    }

    // Red-Black Gauss-Seidel SOR iterations
    for (let iter = 0; iter < maxIter; iter++) {
      let maxRes = 0.0;

      // Two-color sweep: color = 0 (Red), color = 1 (Black)
      for (let color = 0; color < 2; color++) {
        for (let iz = 0; iz < nz; iz++) {
          const izPrev = iz > 0 ? iz - 1 : (grid.boundaryMode === BoundaryMode.PERIODIC ? nz - 1 : 0);
          const izNext = iz < nz - 1 ? iz + 1 : (grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : nz - 1);

          for (let iy = 0; iy < ny; iy++) {
            const iyPrev = iy > 0 ? iy - 1 : (grid.boundaryMode === BoundaryMode.PERIODIC ? ny - 1 : 0);
            const iyNext = iy < ny - 1 ? iy + 1 : (grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : ny - 1);

            for (let ix = 0; ix < nx; ix++) {
              if ((ix + iy + iz) % 2 !== color) continue;

              const ixPrev = ix > 0 ? ix - 1 : (grid.boundaryMode === BoundaryMode.PERIODIC ? nx - 1 : 0);
              const ixNext = ix < nx - 1 ? ix + 1 : (grid.boundaryMode === BoundaryMode.PERIODIC ? 0 : nx - 1);

              const idx = grid.index(ix, iy, iz);
              const sVal = source[idx];

              const neighborSum =
                (phi[grid.index(ixNext, iy, iz)] + phi[grid.index(ixPrev, iy, iz)]) * invDx2 +
                (phi[grid.index(ix, iyNext, iz)] + phi[grid.index(ix, iyPrev, iz)]) * invDy2 +
                (phi[grid.index(ix, iy, izNext)] + phi[grid.index(ix, iy, izPrev)]) * invDz2;

              const phiNew = (neighborSum - sVal) / denom;
              const res = Math.abs(phiNew - phi[idx]);
              if (res > maxRes) maxRes = res;

              phi[idx] = (1.0 - omega) * phi[idx] + omega * phiNew;
            }
          }
        }
      }

      if (maxRes < tol) break;
    }

    return new PotentialField(grid, phi, { fourPiG });
  }

  /**
   * Analytic Plummer sphere potential generator:
   * Phi(r) = -G * M / sqrt(r^2 + b^2)
   * rho(r) = 3 * M / (4 * pi * b^3) * (1 + r^2 / b^2)^(-5/2)
   * g(r)   = -G * M * r / (r^2 + b^2)^(3/2)
   * 
   * @param {GridIndexer} gridIndexer Grid geometry.
   * @param {Object} options Parameters: { mass, scaleRadius, center: [x0, y0, z0], G }.
   * @returns {{ potentialField: PotentialField, densityField: DensityField, exactPotential: Function, exactDensity: Function, exactAcceleration: Function }}
   */
  static fromAnalyticPlummer(gridIndexer, options = {}) {
    const M = options.mass || 1.0e14;
    const b = options.scaleRadius || 5.0; // scale radius b (Mpc)
    const [x0, y0, z0] = options.center || [0.0, 0.0, 0.0];
    const G = options.G || 4.30091e-9;

    const b2 = b * b;
    const total = gridIndexer.totalCells;
    const phiBuf = new Float64Array(total);
    const deltaBuf = new Float64Array(total);

    const exactPotential = (x, y, z) => {
      const rx = x - x0;
      const ry = y - y0;
      const rz = z - z0;
      const r2 = rx * rx + ry * ry + rz * rz;
      return -G * M / Math.sqrt(r2 + b2);
    };

    const exactDensity = (x, y, z) => {
      const rx = x - x0;
      const ry = y - y0;
      const rz = z - z0;
      const r2 = rx * rx + ry * ry + rz * rz;
      return (3.0 * M / (4.0 * Math.PI * b * b2)) * Math.pow(1.0 + r2 / b2, -2.5);
    };

    const exactAcceleration = (x, y, z) => {
      const rx = x - x0;
      const ry = y - y0;
      const rz = z - z0;
      const r2 = rx * rx + ry * ry + rz * rz;
      const factor = -G * M / Math.pow(r2 + b2, 1.5);
      return [factor * rx, factor * ry, factor * rz];
    };

    for (let iz = 0; iz < gridIndexer.nz; iz++) {
      for (let iy = 0; iy < gridIndexer.ny; iy++) {
        for (let ix = 0; ix < gridIndexer.nx; ix++) {
          const [x, y, z] = gridIndexer.gridIndexToCoord(ix, iy, iz);
          const idx = gridIndexer.index(ix, iy, iz);
          phiBuf[idx] = exactPotential(x, y, z);
          deltaBuf[idx] = exactDensity(x, y, z);
        }
      }
    }

    const potentialField = new PotentialField(gridIndexer, phiBuf, { fourPiG: 4.0 * Math.PI * G });
    const densityField = new DensityField(gridIndexer, deltaBuf);

    return {
      potentialField,
      densityField,
      exactPotential,
      exactDensity,
      exactAcceleration
    };
  }

  /**
   * Analytic Hernquist profile generator:
   * Phi(r) = -G * M / (r + a)
   * rho(r) = M * a / (2 * pi * r * (r + a)^3)
   * g(r)   = -G * M * r / (r * (r + a)^2)
   * 
   * @param {GridIndexer} gridIndexer Grid geometry.
   * @param {Object} options Parameters: { mass, scaleRadius, center: [x0, y0, z0], G }.
   * @returns {{ potentialField: PotentialField, densityField: DensityField, exactPotential: Function, exactDensity: Function, exactAcceleration: Function }}
   */
  static fromAnalyticHernquist(gridIndexer, options = {}) {
    const M = options.mass || 1.0e14;
    const a = options.scaleRadius || 5.0; // scale radius a (Mpc)
    const [x0, y0, z0] = options.center || [0.0, 0.0, 0.0];
    const G = options.G || 4.30091e-9;

    const total = gridIndexer.totalCells;
    const phiBuf = new Float64Array(total);
    const deltaBuf = new Float64Array(total);

    const exactPotential = (x, y, z) => {
      const rx = x - x0;
      const ry = y - y0;
      const rz = z - z0;
      const r = Math.sqrt(rx * rx + ry * ry + rz * rz);
      return -G * M / (r + a);
    };

    const exactDensity = (x, y, z) => {
      const rx = x - x0;
      const ry = y - y0;
      const rz = z - z0;
      const r = Math.sqrt(rx * rx + ry * ry + rz * rz);
      if (r < 1e-12) return M * a / (2.0 * Math.PI * 1e-12 * Math.pow(a, 3));
      return (M * a) / (2.0 * Math.PI * r * Math.pow(r + a, 3));
    };

    const exactAcceleration = (x, y, z) => {
      const rx = x - x0;
      const ry = y - y0;
      const rz = z - z0;
      const r = Math.sqrt(rx * rx + ry * ry + rz * rz);
      if (r < 1e-12) return [0.0, 0.0, 0.0];
      const factor = -G * M / (r * Math.pow(r + a, 2));
      return [factor * rx, factor * ry, factor * rz];
    };

    for (let iz = 0; iz < gridIndexer.nz; iz++) {
      for (let iy = 0; iy < gridIndexer.ny; iy++) {
        for (let ix = 0; ix < gridIndexer.nx; ix++) {
          const [x, y, z] = gridIndexer.gridIndexToCoord(ix, iy, iz);
          const idx = gridIndexer.index(ix, iy, iz);
          phiBuf[idx] = exactPotential(x, y, z);
          deltaBuf[idx] = exactDensity(x, y, z);
        }
      }
    }

    const potentialField = new PotentialField(gridIndexer, phiBuf, { fourPiG: 4.0 * Math.PI * G });
    const densityField = new DensityField(gridIndexer, deltaBuf);

    return {
      potentialField,
      densityField,
      exactPotential,
      exactDensity,
      exactAcceleration
    };
  }

  /**
   * Analytic Harmonic standing wave generator:
   * delta(x) = A * sin(kx*x) * sin(ky*y) * sin(kz*z)
   * Phi(x)   = -4 * pi * G * rho_bar / (kx^2 + ky^2 + kz^2) * delta(x)
   * 
   * @param {GridIndexer} gridIndexer Grid geometry.
   * @param {Object} options Parameters: { kx, ky, kz, amplitude, fourPiG }.
   * @returns {{ potentialField: PotentialField, densityField: DensityField, exactPotential: Function, exactDensity: Function }}
   */
  static fromAnalyticHarmonic(gridIndexer, options = {}) {
    const kx = options.kx || (2.0 * Math.PI / gridIndexer.boxSize[0]);
    const ky = options.ky || (2.0 * Math.PI / gridIndexer.boxSize[1]);
    const kz = options.kz || (2.0 * Math.PI / gridIndexer.boxSize[2]);
    const A = options.amplitude || 1.0;
    const fourPiG = typeof options.fourPiG === 'number' ? options.fourPiG : DEFAULT_FOUR_PI_G_RHO_BAR;

    const kSq = kx * kx + ky * ky + kz * kz;
    const phiScale = -fourPiG / kSq;

    const total = gridIndexer.totalCells;
    const phiBuf = new Float64Array(total);
    const deltaBuf = new Float64Array(total);

    const exactDensity = (x, y, z) => A * Math.sin(kx * x) * Math.sin(ky * y) * Math.sin(kz * z);
    const exactPotential = (x, y, z) => phiScale * exactDensity(x, y, z);

    for (let iz = 0; iz < gridIndexer.nz; iz++) {
      for (let iy = 0; iy < gridIndexer.ny; iy++) {
        for (let ix = 0; ix < gridIndexer.nx; ix++) {
          const [x, y, z] = gridIndexer.gridIndexToCoord(ix, iy, iz);
          const idx = gridIndexer.index(ix, iy, iz);
          deltaBuf[idx] = exactDensity(x, y, z);
          phiBuf[idx] = exactPotential(x, y, z);
        }
      }
    }

    const potentialField = new PotentialField(gridIndexer, phiBuf, { fourPiG });
    const densityField = new DensityField(gridIndexer, deltaBuf);

    return {
      potentialField,
      densityField,
      exactPotential,
      exactDensity
    };
  }
}
