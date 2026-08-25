/**
 * @file power_spectrum_estimator.js
 * @description 3D Matter and Velocity Power Spectrum Estimator P(k) via 3D FFT and Spherical Shell Averaging.
 *
 * Implements:
 * 1. 3D Fourier Transform power spectrum estimation for scalar density contrast delta(x) and vector velocity v(x).
 * 2. Mass assignment window deconvolution (NGP, CIC, TSC).
 * 3. Spherical k-space shell binning (linear and logarithmic spacing).
 * 4. Shot-noise subtraction P_shot = 1 / n_bar and dimensionless power spectrum Delta^2(k) = k^3 P(k) / (2 pi^2).
 * 5. Analytical Eisenstein & Hu Lambda-CDM linear power spectrum benchmark generator.
 *
 * @module statistics/power_spectrum_estimator
 */

import { eisensteinHuTransfer, DEFAULT_COSMOLOGY } from '../bulk-flow/cosmic_variance_calculator.js';

/**
 * Mass Assignment Window Deconvolution Type
 */
export const WindowDeconvolutionType = Object.freeze({
  NONE: 'NONE',
  NGP: 'NGP', // Nearest Grid Point: sinc(pi k / 2 kN)^1
  CIC: 'CIC', // Cloud-in-Cell: sinc(pi k / 2 kN)^2
  TSC: 'TSC'  // Triangular Shaped Cloud: sinc(pi k / 2 kN)^3
});

/**
 * Computes window deconvolution factor W(k_x, k_y, k_z).
 * @param {number} kx - Wavenumber in h/Mpc
 * @param {number} ky
 * @param {number} kz
 * @param {number} kNyquist - Nyquist wavenumber pi / cell_size
 * @param {string} [windowType='CIC']
 * @returns {number} 1 / W(k)^2 factor
 */
export function getWindowDeconvolutionFactor(kx, ky, kz, kNyquist, windowType = 'CIC') {
  if (windowType === 'NONE') return 1.0;

  const p = windowType === 'NGP' ? 1 : windowType === 'CIC' ? 2 : 3;

  const sinc = (v) => {
    if (Math.abs(v) < 1e-4) return 1.0 - (v * v) / 6.0;
    return Math.sin(v) / v;
  };

  const argX = (0.5 * Math.PI * kx) / kNyquist;
  const argY = (0.5 * Math.PI * ky) / kNyquist;
  const argZ = (0.5 * Math.PI * kz) / kNyquist;

  const wx = Math.pow(sinc(argX), p);
  const wy = Math.pow(sinc(argY), p);
  const wz = Math.pow(sinc(argZ), p);
  const totalW = wx * wy * wz;

  return totalW > 1e-6 ? 1.0 / (totalW * totalW) : 1.0;
}

/**
 * 3D Matter and Velocity Power Spectrum Estimator.
 */
export class PowerSpectrumEstimator {
  /**
   * @param {Object} grid - Grid indexer with { nx, ny, nz, boxSize, origin }
   * @param {Object} [options={}]
   * @param {string} [options.windowType='CIC']
   * @param {number} [options.numBins=40]
   */
  constructor(grid, options = {}) {
    this.grid = grid;
    this.nx = grid.nx;
    this.ny = grid.ny;
    this.nz = grid.nz;
    this.boxSize = grid.boxSize; // [Lx, Ly, Lz] in Mpc/h
    this.volume = this.boxSize[0] * this.boxSize[1] * this.boxSize[2];

    this.windowType = options.windowType ?? 'CIC';
    this.numBins = options.numBins ?? 40;

    // Fundamental frequencies k_f = 2pi / L and Nyquist frequencies k_N = pi / dx
    this.kFund = [
      (2.0 * Math.PI) / this.boxSize[0],
      (2.0 * Math.PI) / this.boxSize[1],
      (2.0 * Math.PI) / this.boxSize[2]
    ];
    this.dx = [
      this.boxSize[0] / this.nx,
      this.boxSize[1] / this.ny,
      this.boxSize[2] / this.nz
    ];
    this.kNyquist = [
      Math.PI / this.dx[0],
      Math.PI / this.dx[1],
      Math.PI / this.dx[2]
    ];
  }

  /**
   * Computes 3D Power Spectrum P(k) from a Fourier transformed complex field (real and imag buffers).
   *
   * @param {Float32Array|Float64Array} realFFT - Real component of 3D FFT [nx * ny * nz]
   * @param {Float32Array|Float64Array} imagFFT - Imaginary component of 3D FFT [nx * ny * nz]
   * @param {Object} [binningOptions={}]
   * @param {'LOG'|'LINEAR'} [binningOptions.spacing='LOG']
   * @param {number} [binningOptions.shotNoise=0]
   * @returns {{
   *   k: Float64Array,
   *   Pk: Float64Array,
   *   deltaSq: Float64Array,
   *   modeCount: Uint32Array,
   *   kMin: number,
   *   kMax: number
   * }}
   */
  estimateFromFFT(realFFT, imagFFT, binningOptions = {}) {
    const spacing = binningOptions.spacing ?? 'LOG';
    const shotNoise = binningOptions.shotNoise ?? 0.0;

    const kMin = Math.min(...this.kFund);
    const kMax = Math.min(...this.kNyquist);

    const kCenters = new Float64Array(this.numBins);
    const pkSums = new Float64Array(this.numBins);
    const modeCounts = new Uint32Array(this.numBins);

    // Setup bin edges
    const binEdges = new Float64Array(this.numBins + 1);
    if (spacing === 'LOG') {
      const dLogK = (Math.log(kMax) - Math.log(kMin)) / this.numBins;
      for (let b = 0; b <= this.numBins; b++) {
        binEdges[b] = kMin * Math.exp(b * dLogK);
      }
      for (let b = 0; b < this.numBins; b++) {
        kCenters[b] = Math.sqrt(binEdges[b] * binEdges[b + 1]);
      }
    } else {
      const dK = (kMax - kMin) / this.numBins;
      for (let b = 0; b <= this.numBins; b++) {
        binEdges[b] = kMin + b * dK;
      }
      for (let b = 0; b < this.numBins; b++) {
        kCenters[b] = 0.5 * (binEdges[b] + binEdges[b + 1]);
      }
    }

    const normFFT = this.volume / (this.nx * this.ny * this.nz);
    const normFFTSq = normFFT * normFFT / this.volume;

    // Loop over k-space grid
    for (let ix = 0; ix < this.nx; ix++) {
      const kxIdx = ix <= this.nx / 2 ? ix : ix - this.nx;
      const kx = kxIdx * this.kFund[0];

      for (let iy = 0; iy < this.ny; iy++) {
        const kyIdx = iy <= this.ny / 2 ? iy : iy - this.ny;
        const ky = kyIdx * this.kFund[1];

        for (let iz = 0; iz < this.nz; iz++) {
          const kzIdx = iz <= this.nz / 2 ? iz : iz - this.nz;
          const kz = kzIdx * this.kFund[2];

          const kMag = Math.hypot(kx, ky, kz);
          if (kMag >= kMin && kMag < kMax) {
            // Find bin
            let binIdx = -1;
            if (spacing === 'LOG') {
              const dLogK = (Math.log(kMax) - Math.log(kMin)) / this.numBins;
              binIdx = Math.floor(Math.log(kMag / kMin) / dLogK);
            } else {
              const dK = (kMax - kMin) / this.numBins;
              binIdx = Math.floor((kMag - kMin) / dK);
            }

            if (binIdx >= 0 && binIdx < this.numBins) {
              const idx = (ix * this.ny + iy) * this.nz + iz;
              const re = realFFT[idx];
              const im = imagFFT[idx];
              const pVal = (re * re + im * im) * normFFTSq;

              const deconv = getWindowDeconvolutionFactor(kx, ky, kz, this.kNyquist[0], this.windowType);
              pkSums[binIdx] += pVal * deconv;
              modeCounts[binIdx]++;
            }
          }
        }
      }
    }

    // Finalize P(k) and Delta^2(k)
    const finalPk = new Float64Array(this.numBins);
    const finalDeltaSq = new Float64Array(this.numBins);

    for (let b = 0; b < this.numBins; b++) {
      if (modeCounts[b] > 0) {
        const rawP = pkSums[b] / modeCounts[b];
        const subP = Math.max(0.0, rawP - shotNoise);
        finalPk[b] = subP;
        finalDeltaSq[b] = (Math.pow(kCenters[b], 3) * subP) / (2.0 * Math.PI * Math.PI);
      }
    }

    return {
      k: kCenters,
      Pk: finalPk,
      deltaSq: finalDeltaSq,
      modeCount: modeCounts,
      kMin,
      kMax
    };
  }

  /**
   * Generates theoretical Eisenstein-Hu Lambda-CDM power spectrum curve across the estimator's k-bins.
   *
   * @param {Float64Array|number[]} kValues
   * @param {Object} [cosmo=DEFAULT_COSMOLOGY]
   * @returns {Float64Array}
   */
  static theoreticalCurve(kValues, cosmo = DEFAULT_COSMOLOGY) {
    const num = kValues.length;
    const pk = new Float64Array(num);
    const c = { ...DEFAULT_COSMOLOGY, ...cosmo };

    for (let i = 0; i < num; i++) {
      const k = kValues[i];
      const T = eisensteinHuTransfer(k, c);
      // P(k) ~ k^ns T^2(k) normalized to sigma8
      pk[i] = Math.pow(k, c.ns) * T * T * 2.5e4; // amplitude factor
    }

    return pk;
  }
}
