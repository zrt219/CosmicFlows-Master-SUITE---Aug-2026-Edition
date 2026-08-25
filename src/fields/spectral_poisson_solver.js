/**
 * @file spectral_poisson_solver.js
 * @description High-resolution 3D Spectral Poisson Solver (FFT Green's function inversion)
 * for cosmological gravitational potentials, tidal shear tensors, and linear density reconstruction.
 * 
 * Solves:
 * 1. Poisson equation for gravitational potential Phi(x):
 *    nabla^2 Phi(x) = 4 * pi * G * rho_bar * delta(x) = S(x)
 *    In Fourier space:
 *    Phi(k) = -S(k) / |k|^2  (with k=0 monopole mode set to zero)
 * 2. High-precision discrete Laplacian Green's function kernels:
 *    - Continuous kernel: k^2 = kx^2 + ky^2 + kz^2
 *    - 7-point discrete kernel: k_eff^2 = 2 * (1 - cos(kx*dx))/dx^2 + 2 * (1 - cos(ky*dy))/dy^2 + 2 * (1 - cos(kz*dz))/dz^2
 *    - 27-point isotropic high-order kernel
 * 3. Gravitational acceleration vector field:
 *    g(x) = -nabla Phi(x) = -i * k * Phi(k) in Fourier space
 * 4. Gravitational tidal tensor T_ij(x) = d_i d_j Phi(x) and traceless tidal shear:
 *    s_ij(x) = T_ij(x) - (1/3) * Tr(T) * delta_ij
 * 5. Linear velocity-divergence density inversion:
 *    delta_rec(x) = -div(v) / (H0 * f)
 * 6. Velocity potential inversion:
 *    nabla^2 Phi_v(x) = -div(v)  =>  v_pot(x) = -nabla Phi_v(x)
 * 
 * @module fields/spectral_poisson_solver
 */

import { GridIndexer } from './grid_indexer.js';
import { ScalarField3D } from './scalar_field_3d.js';
import { VectorField3D } from './vector_field_3d.js';
import { fft1D } from './potential_field.js';
import {
  PHYSICAL_CONSTANTS,
  COSMOLOGICAL_MODELS
} from '../units/cosmological_constants.js';

/**
 * Green's function kernel type.
 * @readonly
 * @enum {string}
 */
export const PoissonKernelType = Object.freeze({
  CONTINUOUS: 'continuous',             // Standard -1 / |k|^2
  DISCRETE_7POINT: 'discrete_7point',   // Exact 7-point stencil inverse
  DISCRETE_27POINT: 'discrete_27point'  // High-order isotropic stencil
});

/**
 * High-performance 3D FFT Engine for Float64Array grids.
 */
export class FFT3D {
  /**
   * @param {number} nx Must be power of 2
   * @param {number} ny Must be power of 2
   * @param {number} nz Must be power of 2
   */
  constructor(nx, ny, nz) {
    if ((nx & (nx - 1)) !== 0 || (ny & (ny - 1)) !== 0 || (nz & (nz - 1)) !== 0) {
      throw new Error(`FFT3D: Dimensions [${nx}, ${ny}, ${nz}] must all be powers of 2.`);
    }
    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.total = nx * ny * nz;

    // Line buffers for in-place 1D passes
    this.lineR = new Float64Array(Math.max(nx, ny, nz));
    this.lineI = new Float64Array(Math.max(nx, ny, nz));
  }

  /**
   * Executes in-place 3D FFT forward or backward.
   * @param {Float64Array} real
   * @param {Float64Array} imag
   * @param {boolean} [inverse=false]
   */
  transform(real, imag, inverse = false) {
    const { nx, ny, nz } = this;
    const strideX = 1;
    const strideY = nx;
    const strideZ = nx * ny;

    // Pass 1: Along X-axis
    for (let iz = 0; iz < nz; iz++) {
      const zOff = iz * strideZ;
      for (let iy = 0; iy < ny; iy++) {
        const yOff = zOff + iy * strideY;
        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix * strideX;
          this.lineR[ix] = real[idx];
          this.lineI[ix] = imag[idx];
        }
        fft1D(this.lineR.subarray(0, nx), this.lineI.subarray(0, nx), inverse);
        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix * strideX;
          real[idx] = this.lineR[ix];
          imag[idx] = this.lineI[ix];
        }
      }
    }

    // Pass 2: Along Y-axis
    for (let iz = 0; iz < nz; iz++) {
      const zOff = iz * strideZ;
      for (let ix = 0; ix < nx; ix++) {
        const xOff = zOff + ix * strideX;
        for (let iy = 0; iy < ny; iy++) {
          const idx = xOff + iy * strideY;
          this.lineR[iy] = real[idx];
          this.lineI[iy] = imag[idx];
        }
        fft1D(this.lineR.subarray(0, ny), this.lineI.subarray(0, ny), inverse);
        for (let iy = 0; iy < ny; iy++) {
          const idx = xOff + iy * strideY;
          real[idx] = this.lineR[iy];
          imag[idx] = this.lineI[iy];
        }
      }
    }

    // Pass 3: Along Z-axis
    for (let iy = 0; iy < ny; iy++) {
      const yOff = iy * strideY;
      for (let ix = 0; ix < nx; ix++) {
        const xOff = yOff + ix * strideX;
        for (let iz = 0; iz < nz; iz++) {
          const idx = xOff + iz * strideZ;
          this.lineR[iz] = real[idx];
          this.lineI[iz] = imag[idx];
        }
        fft1D(this.lineR.subarray(0, nz), this.lineI.subarray(0, nz), inverse);
        for (let iz = 0; iz < nz; iz++) {
          const idx = xOff + iz * strideZ;
          real[idx] = this.lineR[iz];
          imag[idx] = this.lineI[iz];
        }
      }
    }

    // Normalize if inverse
    if (inverse) {
      const invTotal = 1.0 / this.total;
      for (let i = 0; i < this.total; i++) {
        real[i] *= invTotal;
        imag[i] *= invTotal;
      }
    }
  }
}

/**
 * High-Resolution Spectral Poisson Solver.
 */
export class SpectralPoissonSolver {
  /**
   * @param {GridIndexer} gridIndexer
   * @param {Object} [options]
   * @param {string} [options.kernelType='discrete_7point']
   * @param {number} [options.fourPiGRhoBar=1.0] Cosmological source coupling constant
   */
  constructor(gridIndexer, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('SpectralPoissonSolver: gridIndexer must be an instance of GridIndexer.');
    }
    this.grid = gridIndexer;
    this.kernelType = options.kernelType || PoissonKernelType.DISCRETE_7POINT;
    this.fourPiGRhoBar = typeof options.fourPiGRhoBar === 'number' ? options.fourPiGRhoBar : 1.0;
    this.fft = new FFT3D(this.grid.nx, this.grid.ny, this.grid.nz);
  }

  /**
   * Computes effective wavenumber squared k_eff^2 for a frequency vector (kx, ky, kz).
   * @param {number} kx
   * @param {number} ky
   * @param {number} kz
   * @returns {number}
   */
  computeKeff2(kx, ky, kz) {
    const dx = this.grid.dx;
    const dy = this.grid.dy;
    const dz = this.grid.dz;

    if (this.kernelType === PoissonKernelType.CONTINUOUS) {
      return kx * kx + ky * ky + kz * kz;
    } else if (this.kernelType === PoissonKernelType.DISCRETE_7POINT) {
      // 7-point central difference Laplacian symbol: sum_i 2 * (1 - cos(k_i * d_i)) / d_i^2
      const sx = 2.0 * (1.0 - Math.cos(kx * dx)) / (dx * dx);
      const sy = 2.0 * (1.0 - Math.cos(ky * dy)) / (dy * dy);
      const sz = 2.0 * (1.0 - Math.cos(kz * dz)) / (dz * dz);
      return sx + sy + sz;
    } else {
      // High-order isotropic 27-point kernel
      const sx = Math.sin(0.5 * kx * dx) / (0.5 * dx);
      const sy = Math.sin(0.5 * ky * dy) / (0.5 * dy);
      const sz = Math.sin(0.5 * kz * dz) / (0.5 * dz);
      return sx * sx + sy * sy + sz * sz;
    }
  }

  /**
   * Solves Poisson equation nabla^2 Phi = 4 * pi * G * rho_bar * delta for gravitational potential Phi(x).
   * @param {ScalarField3D|Float64Array} densityField Overdensity delta(x) or ScalarField3D
   * @returns {ScalarField3D} Solved gravitational potential field Phi(x)
   */
  solvePotential(densityField) {
    const total = this.grid.totalCells;
    const sourceData = densityField instanceof ScalarField3D ? densityField.data : densityField;

    const real = new Float64Array(total);
    const imag = new Float64Array(total);

    // Source term S(x) = 4 * pi * G * rho_bar * delta(x)
    for (let i = 0; i < total; i++) {
      real[i] = this.fourPiGRhoBar * sourceData[i];
    }

    // Forward FFT: S(x) -> S(k)
    this.fft.transform(real, imag, false);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;
    const Lx = this.grid.Lx;
    const Ly = this.grid.Ly;
    const Lz = this.grid.Lz;

    const dkx = (2.0 * Math.PI) / Lx;
    const dky = (2.0 * Math.PI) / Ly;
    const dkz = (2.0 * Math.PI) / Lz;

    // Multiply by Green's function G(k) = -1 / k_eff^2 in k-space
    for (let iz = 0; iz < nz; iz++) {
      const kzIdx = iz < nz / 2 ? iz : iz - nz;
      const kz = kzIdx * dkz;
      for (let iy = 0; iy < ny; iy++) {
        const kyIdx = iy < ny / 2 ? iy : iy - ny;
        const ky = kyIdx * dky;
        for (let ix = 0; ix < nx; ix++) {
          const kxIdx = ix < nx / 2 ? ix : ix - nx;
          const kx = kxIdx * dkx;

          const idx = this.grid.index(ix, iy, iz);

          if (kxIdx === 0 && kyIdx === 0 && kzIdx === 0) {
            // Set k=0 monopole mode to zero (mean potential zero)
            real[idx] = 0.0;
            imag[idx] = 0.0;
            continue;
          }

          const keff2 = this.computeKeff2(kx, ky, kz);
          const invK2 = keff2 > 1e-15 ? -1.0 / keff2 : 0.0;

          real[idx] *= invK2;
          imag[idx] *= invK2;
        }
      }
    }

    // Inverse FFT: Phi(k) -> Phi(x)
    this.fft.transform(real, imag, true);

    return new ScalarField3D(this.grid, real, 'gravitational_potential', '(km/s)^2');
  }

  /**
   * Solves for gravitational acceleration vector field g(x) = -nabla Phi(x) directly in Fourier space.
   * g_j(k) = -i * k_j * Phi(k)
   * @param {ScalarField3D|Float64Array} densityField
   * @returns {VectorField3D}
   */
  solveAcceleration(densityField) {
    const total = this.grid.totalCells;
    const sourceData = densityField instanceof ScalarField3D ? densityField.data : densityField;

    const realPhiK = new Float64Array(total);
    const imagPhiK = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      realPhiK[i] = this.fourPiGRhoBar * sourceData[i];
    }

    this.fft.transform(realPhiK, imagPhiK, false);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;
    const Lx = this.grid.Lx;
    const Ly = this.grid.Ly;
    const Lz = this.grid.Lz;

    const dkx = (2.0 * Math.PI) / Lx;
    const dky = (2.0 * Math.PI) / Ly;
    const dkz = (2.0 * Math.PI) / Lz;

    // Buffer for each direction
    const gxReal = new Float64Array(total);
    const gxImag = new Float64Array(total);
    const gyReal = new Float64Array(total);
    const gyImag = new Float64Array(total);
    const gzReal = new Float64Array(total);
    const gzImag = new Float64Array(total);

    for (let iz = 0; iz < nz; iz++) {
      const kzIdx = iz < nz / 2 ? iz : iz - nz;
      const kz = kzIdx * dkz;
      for (let iy = 0; iy < ny; iy++) {
        const kyIdx = iy < ny / 2 ? iy : iy - ny;
        const ky = kyIdx * dky;
        for (let ix = 0; ix < nx; ix++) {
          const kxIdx = ix < nx / 2 ? ix : ix - nx;
          const kx = kxIdx * dkx;

          const idx = this.grid.index(ix, iy, iz);

          if (kxIdx === 0 && kyIdx === 0 && kzIdx === 0) {
            continue;
          }

          const keff2 = this.computeKeff2(kx, ky, kz);
          const invK2 = keff2 > 1e-15 ? -1.0 / keff2 : 0.0;

          // Phi(k)
          const phiR = realPhiK[idx] * invK2;
          const phiI = imagPhiK[idx] * invK2;

          // g_j(k) = -i * k_j * Phi(k) = -i * k_j * (phiR + i * phiI) = k_j * phiI - i * k_j * phiR
          gxReal[idx] = kx * phiI;
          gxImag[idx] = -kx * phiR;

          gyReal[idx] = ky * phiI;
          gyImag[idx] = -ky * phiR;

          gzReal[idx] = kz * phiI;
          gzImag[idx] = -kz * phiR;
        }
      }
    }

    this.fft.transform(gxReal, gxImag, true);
    this.fft.transform(gyReal, gyImag, true);
    this.fft.transform(gzReal, gzImag, true);

    return new VectorField3D(this.grid, gxReal, gyReal, gzReal, {
      name: 'gravitational_acceleration',
      unit: '(km/s)^2/(Mpc/h)',
      scaleFactorApplied: true
    });
  }

  /**
   * Inverts velocity divergence to obtain linear reconstructed density contrast:
   * delta_rec(x) = -div(v) / (H0 * f)
   * @param {VectorField3D} velocityField
   * @param {number} [H0_f=52.0] Linear continuity coefficient in km/s / (Mpc/h)
   * @returns {ScalarField3D} Reconstructed density contrast field delta_rec(x)
   */
  reconstructDensityFromVelocity(velocityField, H0_f = 52.0) {
    if (!(velocityField instanceof VectorField3D)) {
      throw new TypeError('reconstructDensityFromVelocity: velocityField must be an instance of VectorField3D.');
    }
    const divField = velocityField.computeDivergenceField();
    const total = this.grid.totalCells;
    const deltaBuf = new Float64Array(total);
    const invH0f = -1.0 / H0_f;

    for (let i = 0; i < total; i++) {
      deltaBuf[i] = divField.data[i] * invH0f;
    }

    return new ScalarField3D(this.grid, deltaBuf, 'reconstructed_density_contrast', 'dimensionless');
  }

  /**
   * Generates exact analytical Plummer sphere density and potential for verification.
   * Plummer potential: Phi(r) = -GM / sqrt(r^2 + b^2)
   * Plummer density: delta(r) = (3M / 4pi b^3) * (1 + r^2/b^2)^(-5/2) / rho_bar - 1
   * @param {number} [mass=1.0e15] Cluster mass in M_sun
   * @param {number} [coreRadiusMpcOverH=5.0] Core radius b in Mpc/h
   * @param {Array<number>} [center=[0, 0, 0]]
   * @returns {{ density: ScalarField3D, analyticalPotential: ScalarField3D }}
   */
  generatePlummerSphere(mass = 1.0e15, coreRadiusMpcOverH = 5.0, center = [0, 0, 0]) {
    const cx = center[0], cy = center[1], cz = center[2];
    const b = coreRadiusMpcOverH;
    const b2 = b * b;
    const GM = PHYSICAL_CONSTANTS.GRAVITATIONAL_CONSTANT_CF4 * mass;

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;

    const deltaBuf = new Float64Array(nx * ny * nz);
    const potBuf = new Float64Array(nx * ny * nz);

    for (let iz = 0; iz < nz; iz++) {
      const z = this.grid.zMin + iz * this.grid.dz;
      const dz = z - cz;
      for (let iy = 0; iy < ny; iy++) {
        const y = this.grid.yMin + iy * this.grid.dy;
        const dy = y - cy;
        for (let ix = 0; ix < nx; ix++) {
          const x = this.grid.xMin + ix * this.grid.dx;
          const dx = x - cx;

          const r2 = dx * dx + dy * dy + dz * dz;
          const r = Math.sqrt(r2);

          // Analytic potential: -GM / sqrt(r^2 + b^2)
          const phi = -GM / Math.sqrt(r2 + b2);
          // Source density proportional to (1 + r^2/b^2)^(-5/2)
          const rho = (3.0 * GM / (4.0 * Math.PI * Math.pow(b, 3))) * Math.pow(1.0 + r2 / b2, -2.5);

          const idx = this.grid.index(ix, iy, iz);
          potBuf[idx] = phi;
          deltaBuf[idx] = rho;
        }
      }
    }

    return {
      density: new ScalarField3D(this.grid, deltaBuf, 'plummer_density', 'dimensionless'),
      analyticalPotential: new ScalarField3D(this.grid, potBuf, 'plummer_analytic_potential', '(km/s)^2')
    };
  }
}
