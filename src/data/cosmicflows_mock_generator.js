/**
 * @file cosmicflows_mock_generator.js
 * @module data/cosmicflows_mock_generator
 * @description Production-grade Synthetic Cosmological Flow and Galaxy Mock Catalog Generator.
 * 
 * Implements:
 * 1. Lambda-CDM Matter Power Spectrum P(k) = A k^{n_s} T^2(k) with full Eisenstein & Hu (1998)
 *    transfer function (with and without baryon acoustic oscillations), Sugiyama shape parameters,
 *    BBKS / Bond & Efstathiou transfer models, and top-hat sigma_8 cosmological normalization
 *    via adaptive Gauss-Legendre quadrature.
 * 2. 3D Gaussian Random Field (GRF) generator on periodic Cartesian grids with Hermitian conjugate
 *    symmetry, deterministic seedable pseudo-random number generator (Mulberry32/Box-Muller/Poisson),
 *    and Radix-2 3D Fast Fourier Transform engine.
 * 3. 3D Zel'dovich Displacement Field Psi(q) in Lagrangian Perturbation Theory (1LPT), peculiar
 *    velocity field generation v(q) = 100*h*f*Psi(q), continuous trilinear and tricubic Hermite
 *    interpolation, N-body mock particle synthesis, Cloud-in-Cell (CIC) / Nearest Grid Point (NGP)
 *    mass assignment, and caustic / shell-crossing deformation tensor diagnostics (Jacobi eigenvalue
 *    decomposition).
 * 4. Realistic galaxy mock catalog generation with density-dependent tracer bias models (linear,
 *    log-normal, peak-patch thresholding), generalized exponential radial selection functions
 *    n_bar(r) = n0 * exp(-(r/r0)^gamma), Schechter luminosity completeness, and Zone of Avoidance (ZoA)
 *    dust obscuration masks with HI piercing corridors (Vela, Norma, Puppis).
 * 5. Astrometric distance modulus modeling with realistic Tully-Fisher / Fundamental Plane / Supernova
 *    distance errors (sigma_d / d ~ 0.15 - 0.20), line-of-sight peculiar velocity projection, thermal
 *    velocity dispersions, and Homogeneous / Inhomogeneous Malmquist bias diagnostics.
 * 6. High-density columnar and record-based GalaxyCatalog storage with Equatorial, Galactic, and Supergalactic
 *    astrometric coordinates and JSON, CSV, and GeoJSON export pipelines.
 * 7. Blind Benchmark Generator for cosmological reconstruction challenge verification with injected
 *    attractors, bulk flow dipoles, and automated recovery metric evaluators.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

import { GridIndexer, BoundaryMode } from '../fields/grid_indexer.js';
import { DensityField } from '../fields/density_field.js';
import { VelocityField } from '../fields/velocity_field.js';
import { PotentialField } from '../fields/potential_field.js';
import {
  PHYSICAL_CONSTANTS,
  COSMOLOGICAL_MODELS,
  CosmologicalCalculator
} from '../units/cosmological_constants.js';
import { ASTROMETRIC_CONSTANTS } from '../coordinates/canonical_frame.js';

// ============================================================================
// CONSTANTS AND ASTROMETRIC CONVERSIONS
// ============================================================================

export const DEG_TO_RAD = Math.PI / 180.0;
export const RAD_TO_DEG = 180.0 / Math.PI;
export const LN10 = Math.LN10;
export const INV_LN10 = 1.0 / Math.LN10;
export const LN10_DIV_5 = Math.LN10 / 5.0;            // ~0.460517018598809136
export const FIVE_DIV_LN10 = 5.0 / Math.LN10;           // ~2.171472409516259

export const DEFAULT_MOCK_SEED = 4294967291;
export const DEFAULT_BOX_SIZE_MPC = 500.0;  // Comoving box size in Mpc/h
export const DEFAULT_GRID_RESOLUTION = 64;   // Power of 2 grid cells per dimension

// ============================================================================
// DETERMINISTIC PSEUDO-RANDOM NUMBER GENERATOR (CosmoRNG)
// ============================================================================

/**
 * High-performance deterministic 32-bit PRNG (Mulberry32 + Box-Muller & Poisson algorithms)
 * providing reproducible cosmological field realizations and mock catalog extractions.
 */
export class CosmoRNG {
  /**
   * @param {number} [seed=DEFAULT_MOCK_SEED] 32-bit unsigned integer seed.
   */
  constructor(seed = DEFAULT_MOCK_SEED) {
    this.seed = (seed >>> 0) || 1;
    this.state = this.seed;
    this._hasSpareNormal = false;
    this._spareNormal = 0.0;
  }

  /**
   * Reseeds the PRNG.
   * @param {number} seed
   */
  reseed(seed) {
    this.seed = (seed >>> 0) || 1;
    this.state = this.seed;
    this._hasSpareNormal = false;
    this._spareNormal = 0.0;
  }

  /**
   * Generates uniform random float in [0, 1).
   * @returns {number}
   */
  uniform() {
    let t = (this.state += 0x6D2B79F5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296.0;
  }

  /**
   * Generates uniform random float in range [min, max).
   * @param {number} min
   * @param {number} max
   * @returns {number}
   */
  uniformRange(min, max) {
    return min + (max - min) * this.uniform();
  }

  /**
   * Generates standard normal random variable N(0, 1) using Box-Muller transform.
   * @returns {number}
   */
  standardNormal() {
    if (this._hasSpareNormal) {
      this._hasSpareNormal = false;
      return this._spareNormal;
    }
    let u = 0.0;
    let v = 0.0;
    while (u <= 1e-15) u = this.uniform();
    v = this.uniform();
    const r = Math.sqrt(-2.0 * Math.log(u));
    const theta = 2.0 * Math.PI * v;
    this._spareNormal = r * Math.sin(theta);
    this._hasSpareNormal = true;
    return r * Math.cos(theta);
  }

  /**
   * Generates normal random variable N(mean, stdDev).
   * @param {number} [mean=0.0]
   * @param {number} [stdDev=1.0]
   * @returns {number}
   */
  normal(mean = 0.0, stdDev = 1.0) {
    return mean + stdDev * this.standardNormal();
  }

  /**
   * Generates Poisson distributed integer with expectation lambda.
   * Uses exact Knuth inverse transform for lambda < 30 and Gaussian approximation for large lambda.
   * @param {number} lambda
   * @returns {number}
   */
  poisson(lambda) {
    if (lambda <= 0.0) return 0;
    if (lambda < 30.0) {
      const L = Math.exp(-lambda);
      let k = 0;
      let p = 1.0;
      do {
        k++;
        p *= this.uniform();
      } while (p > L);
      return k - 1;
    }
    // High-lambda Gaussian approximation with continuity correction
    const val = Math.round(this.normal(lambda, Math.sqrt(lambda)));
    return val > 0 ? val : 0;
  }

  /**
   * Generates standard exponential random variable with rate parameter lambda.
   * @param {number} [rate=1.0]
   * @returns {number}
   */
  exponential(rate = 1.0) {
    let u = 0.0;
    while (u <= 1e-15) u = this.uniform();
    return -Math.log(u) / rate;
  }

  /**
   * Generates a uniform random unit vector on the 2-sphere S^2.
   * @returns {[number, number, number]} Unit vector [nx, ny, nz]
   */
  uniformSphereVector() {
    const z = this.uniformRange(-1.0, 1.0);
    const phi = this.uniformRange(0.0, 2.0 * Math.PI);
    const rInPlane = Math.sqrt(Math.max(0.0, 1.0 - z * z));
    return [rInPlane * Math.cos(phi), rInPlane * Math.sin(phi), z];
  }
}

// ============================================================================
// RADIX-2 3D FAST FOURIER TRANSFORM (CosmicFFT3D)
// ============================================================================

/**
 * Optimized 1D Radix-2 In-Place Fast Fourier Transform.
 * @param {Float64Array} real
 * @param {Float64Array} imag
 * @param {boolean} [inverse=false]
 */
export function fft1DInPlace(real, imag, inverse = false) {
  const n = real.length;
  if ((n & (n - 1)) !== 0) {
    throw new Error(`fft1DInPlace: Array length (${n}) must be a power of 2.`);
  }

  // Bit reversal
  let j = 0;
  for (let i = 0; i < n - 1; i++) {
    if (i < j) {
      const tr = real[i]; real[i] = real[j]; real[j] = tr;
      const ti = imag[i]; imag[i] = imag[j]; imag[j] = ti;
    }
    let k = n >> 1;
    while (k <= j) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }

  // Cooley-Tukey butterflies
  const sign = inverse ? 1.0 : -1.0;
  for (let len = 2; len <= n; len <<= 1) {
    const half = len >> 1;
    const angle = (sign * 2.0 * Math.PI) / len;
    const wStepR = Math.cos(angle);
    const wStepI = Math.sin(angle);

    for (let i = 0; i < n; i += len) {
      let wR = 1.0;
      let wI = 0.0;
      for (let k = 0; k < half; k++) {
        const uR = real[i + k];
        const uI = imag[i + k];
        const vR = real[i + k + half] * wR - imag[i + k + half] * wI;
        const vI = real[i + k + half] * wI + imag[i + k + half] * wR;

        real[i + k] = uR + vR;
        imag[i + k] = uI + vI;
        real[i + k + half] = uR - vR;
        imag[i + k + half] = uI - vI;

        const nextWR = wR * wStepR - wI * wStepI;
        const nextWI = wR * wStepI + wI * wStepR;
        wR = nextWR;
        wI = nextWI;
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
 * Multi-dimensional 3D Fast Fourier Transform engine for gridded cosmological fields.
 */
export class CosmicFFT3D {
  /**
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   */
  constructor(nx, ny, nz) {
    if ((nx & (nx - 1)) !== 0 || (ny & (ny - 1)) !== 0 || (nz & (nz - 1)) !== 0) {
      throw new Error(`CosmicFFT3D: Dimensions [${nx}, ${ny}, ${nz}] must all be powers of 2.`);
    }
    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.total = nx * ny * nz;

    const maxDim = Math.max(nx, ny, nz);
    this._bufR = new Float64Array(maxDim);
    this._bufI = new Float64Array(maxDim);
  }

  /**
   * Executes in-place 3D FFT forward or inverse.
   * @param {Float64Array} real
   * @param {Float64Array} imag
   * @param {boolean} [inverse=false]
   */
  transform(real, imag, inverse = false) {
    const { nx, ny, nz } = this;
    const strideX = 1;
    const strideY = nx;
    const strideZ = nx * ny;

    // Pass 1: X-direction
    for (let iz = 0; iz < nz; iz++) {
      const zOff = iz * strideZ;
      for (let iy = 0; iy < ny; iy++) {
        const yOff = zOff + iy * strideY;
        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix * strideX;
          this._bufR[ix] = real[idx];
          this._bufI[ix] = imag[idx];
        }
        fft1DInPlace(this._bufR.subarray(0, nx), this._bufI.subarray(0, nx), inverse);
        for (let ix = 0; ix < nx; ix++) {
          const idx = yOff + ix * strideX;
          real[idx] = this._bufR[ix];
          imag[idx] = this._bufI[ix];
        }
      }
    }

    // Pass 2: Y-direction
    for (let iz = 0; iz < nz; iz++) {
      const zOff = iz * strideZ;
      for (let ix = 0; ix < nx; ix++) {
        const xOff = zOff + ix * strideX;
        for (let iy = 0; iy < ny; iy++) {
          const idx = xOff + iy * strideY;
          this._bufR[iy] = real[idx];
          this._bufI[iy] = imag[idx];
        }
        fft1DInPlace(this._bufR.subarray(0, ny), this._bufI.subarray(0, ny), inverse);
        for (let iy = 0; iy < ny; iy++) {
          const idx = xOff + iy * strideY;
          real[idx] = this._bufR[iy];
          imag[idx] = this._bufI[iy];
        }
      }
    }

    // Pass 3: Z-direction
    for (let iy = 0; iy < ny; iy++) {
      const yOff = iy * strideY;
      for (let ix = 0; ix < nx; ix++) {
        const xOff = yOff + ix * strideX;
        for (let iz = 0; iz < nz; iz++) {
          const idx = xOff + iz * strideZ;
          this._bufR[iz] = real[idx];
          this._bufI[iz] = imag[idx];
        }
        fft1DInPlace(this._bufR.subarray(0, nz), this._bufI.subarray(0, nz), inverse);
        for (let iz = 0; iz < nz; iz++) {
          const idx = xOff + iz * strideZ;
          real[idx] = this._bufR[iz];
          imag[idx] = this._bufI[iz];
        }
      }
    }
  }
}

// ============================================================================
// EISENSTEIN & HU (1998) MATTER POWER SPECTRUM & TRANSFER FUNCTIONS
// ============================================================================

/**
 * Eisenstein & Hu (1998) Cosmological Transfer Function Engine.
 * Supports both full acoustic baryon oscillations and smooth zero-baryon transfer curves.
 */
export class EisensteinHuTransfer {
  /**
   * @param {Object} [cosmo] Cosmological parameters
   * @param {number} [cosmo.h=0.746] Hubble parameter h = H0/100
   * @param {number} [cosmo.Omega_m=0.31] Total matter density parameter
   * @param {number} [cosmo.Omega_b=0.048] Baryon density parameter
   * @param {number} [cosmo.T_CMB=2.7255] CMB temperature in Kelvin
   */
  constructor(cosmo = {}) {
    this.h = cosmo.h ?? 0.746;
    this.Omega_m = cosmo.Omega_m ?? 0.31;
    this.Omega_b = cosmo.Omega_b ?? 0.048;
    this.Omega_c = Math.max(0.001, this.Omega_m - this.Omega_b);
    this.T_CMB = cosmo.T_CMB ?? 2.7255;

    // Derived physical densities om_m = Omega_m * h^2, etc.
    this.om_m = this.Omega_m * this.h * this.h;
    this.om_b = this.Omega_b * this.h * this.h;
    this.om_c = this.Omega_c * this.h * this.h;
    this.f_b = this.Omega_b / this.Omega_m;
    this.f_c = this.Omega_c / this.Omega_m;
    this.theta_2_7 = this.T_CMB / 2.7;

    // Equality redshift and horizon scale (Eisenstein & Hu 1998 eq. 2 & 3)
    this.z_eq = 2.50e4 * this.om_m * Math.pow(this.theta_2_7, -4.0);
    this.k_eq = 0.0746 * this.om_m * Math.pow(this.theta_2_7, -2.0); // Mpc^-1

    // Drag epoch redshift z_d (eq. 4)
    const b1 = 0.313 * Math.pow(this.om_m, -0.419) * (1.0 + 0.607 * Math.pow(this.om_m, 0.674));
    const b2 = 0.238 * Math.pow(this.om_m, 0.223);
    this.z_d = 1291.0 * (Math.pow(this.om_m, 0.251) / (1.0 + 0.659 * Math.pow(this.om_m, 0.828))) *
      (1.0 + b1 * Math.pow(this.om_b, b2));

    // Baryon-to-photon ratio at equality and drag epoch (eq. 5)
    this.y_d = (1.0 + this.z_eq) / (1.0 + this.z_d);
    this.R_eq = 31.5 * this.om_b * Math.pow(this.theta_2_7, -4.0) * (1000.0 / this.z_eq);
    this.R_d = 31.5 * this.om_b * Math.pow(this.theta_2_7, -4.0) * (1000.0 / this.z_d);

    // Sound horizon at drag epoch s in Mpc (eq. 6)
    const term1 = Math.sqrt(1.0 + this.R_d) + Math.sqrt(this.R_d + this.R_eq);
    const term2 = 1.0 + Math.sqrt(this.R_eq);
    this.soundHorizon = (2.0 / (3.0 * this.k_eq)) * Math.sqrt(6.0 / this.R_eq) * Math.log(term1 / term2);

    // Silk damping scale k_silk in Mpc^-1 (eq. 7)
    this.k_silk = 1.6 * Math.pow(this.om_b, 0.52) * Math.pow(this.om_m, 0.73) *
      (1.0 + Math.pow(10.4 * this.om_m, -0.95));

    // Sugiyama shape parameter for smooth zero-baryon transfer function
    this.gammaSugiyama = this.Omega_m * this.h * Math.exp(-this.Omega_b - Math.sqrt(2.0 * this.h) * (this.Omega_b / this.Omega_m));
  }

  /**
   * Zero-baryon shape transfer function T_0(q, alpha_c, beta_c).
   * @param {number} q Wavenumber scaled to equality horizon: q = k / (13.41 * k_eq)
   * @param {number} [alpha_c=1.0]
   * @param {number} [beta_c=1.0]
   * @returns {number}
   */
  T0(q, alpha_c = 1.0, beta_c = 1.0) {
    if (q <= 0.0) return 1.0;
    const lnTerm = Math.log(Math.E + 1.8 * beta_c * q);
    const C = 14.2 / alpha_c + 386.0 / (1.0 + 69.9 * Math.pow(q, 1.08));
    return lnTerm / (lnTerm + C * q * q);
  }

  /**
   * Full Eisenstein & Hu (1998) transfer function T(k) including Baryon Acoustic Oscillations.
   * @param {number} k Wavenumber in Mpc^-1 (or h*Mpc^-1 if converted).
   * @param {boolean} [includeBAO=true] If false, returns smoothed no-wiggle transfer.
   * @returns {number} Transfer function amplitude T(k) in [0, 1].
   */
  computeTransfer(k, includeBAO = true) {
    if (k <= 0.0) return 1.0;
    const q = k / (13.41 * this.k_eq);

    if (!includeBAO) {
      // Smooth no-baryon approximation with Sugiyama shape parameter
      const qSugiyama = k / this.gammaSugiyama;
      const L0 = Math.log(2.0 * Math.E + 1.8 * qSugiyama);
      const C0 = 14.2 + 731.0 / (1.0 + 62.5 * qSugiyama);
      return L0 / (L0 + C0 * qSugiyama * qSugiyama);
    }

    // CDM component T_c(k)
    const a1 = Math.pow(46.9 * this.om_m, 0.670) * (1.0 + Math.pow(32.1 * this.om_m, -0.532));
    const a2 = Math.pow(12.0 * this.om_m, 0.424) * (1.0 + Math.pow(45.0 * this.om_m, -0.582));
    const alpha_c = Math.pow(a1, -this.f_b) * Math.pow(a2, -Math.pow(this.f_b, 3.0));

    const b1 = 0.944 / (1.0 + Math.pow(458.0 * this.om_m, -0.708));
    const b2 = Math.pow(0.395 * this.om_m, -0.0266);
    const beta_c = 1.0 / (1.0 + b1 * (Math.pow(this.f_c, b2) - 1.0));

    const f_term = 1.0 / (1.0 + Math.pow(k * this.soundHorizon / 5.4, 4.0));
    const Tc = f_term * this.T0(q, 1.0, beta_c) + (1.0 - f_term) * this.T0(q, alpha_c, beta_c);

    // Baryon component T_b(k)
    const alpha_b = 2.07 * this.k_eq * this.soundHorizon * Math.pow(1.0 + this.R_d, -0.75) *
      (1.0 + this.R_eq) / (1.0 + this.R_d);
    const beta_b = 0.5 + this.f_b + (3.0 - 2.0 * this.f_b) * Math.sqrt(Math.pow(17.2 * this.om_m, 2.0) + 1.0);
    const beta_node = 8.41 * Math.pow(this.om_m, 0.435);

    const s_tilde = this.soundHorizon / Math.pow(1.0 + Math.pow(beta_node / (k * this.soundHorizon), 3.0), 1.0 / 3.0);
    const ks = k * this.soundHorizon;
    const ks_tilde = k * s_tilde;
    const Tb_oscillatory = Math.sin(ks_tilde) / (ks_tilde > 1e-6 ? ks_tilde : 1.0);
    const Tb_damping = Math.exp(-Math.pow(k / this.k_silk, 1.4));
    const Tb = (this.T0(q, 1.0, 1.0) / (1.0 + Math.pow(ks / 5.2, 2.0)) +
      alpha_b / (1.0 + Math.pow(beta_b / ks, 3.0)) * Tb_damping) * Tb_oscillatory;

    // Total matter transfer function T(k) = (Omega_b/Omega_m)*Tb + (Omega_c/Omega_m)*Tc
    const T_total = this.f_b * Tb + this.f_c * Tc;
    return Math.abs(T_total);
  }

  /**
   * BBKS (Bardeen, Bond, Kaiser, Szalay 1986) transfer function alternative.
   * @param {number} k Wavenumber in h*Mpc^-1
   * @param {number} [gamma=0.21] Shape parameter Gamma
   * @returns {number}
   */
  static bbksTransfer(k, gamma = 0.21) {
    if (k <= 0.0) return 1.0;
    const q = k / gamma;
    const term = Math.log(1.0 + 2.34 * q) / (2.34 * q);
    const denom = Math.pow(1.0 + 3.89 * q + Math.pow(16.1 * q, 2.0) + Math.pow(5.46 * q, 3.0) + Math.pow(6.71 * q, 4.0), 0.25);
    return term / denom;
  }

  /**
   * Bond & Efstathiou (1984) transfer function alternative.
   * @param {number} k Wavenumber in h*Mpc^-1
   * @param {number} [gamma=0.21]
   * @returns {number}
   */
  static bondEfstathiouTransfer(k, gamma = 0.21) {
    if (k <= 0.0) return 1.0;
    const q = k / gamma;
    const a = 6.4, b = 3.0, c = 1.7, nu = 1.13;
    const inner = 1.0 + Math.pow(a * q + Math.pow(b * q, 1.5) + Math.pow(c * q, 2.0), nu);
    return Math.pow(inner, -1.0 / nu);
  }
}

/**
 * Full Linear Matter Power Spectrum P(k) = A k^{n_s} T^2(k) Engine.
 * Handles top-hat sigma_8 normalization, growth factor scaling D(z), and window function convolutions.
 */
export class MatterPowerSpectrum {
  /**
   * @param {Object} [params]
   * @param {number} [params.H0=74.6] Hubble constant in km/s/Mpc
   * @param {number} [params.Omega_m=0.31] Total matter density
   * @param {number} [params.Omega_b=0.048] Baryon density
   * @param {number} [params.Omega_Lambda=0.69] Dark energy density
   * @param {number} [params.sigma8=0.81] Amplitude of mass fluctuations at 8 Mpc/h
   * @param {number} [params.ns=0.965] Primordial spectral index
   * @param {boolean} [params.includeBAO=true] Whether to include acoustic oscillations
   */
  constructor(params = {}) {
    this.H0 = params.H0 ?? COSMOLOGICAL_MODELS.CF4.H0;
    this.h = this.H0 / 100.0;
    this.Omega_m = params.Omega_m ?? COSMOLOGICAL_MODELS.CF4.Omega_m;
    this.Omega_b = params.Omega_b ?? COSMOLOGICAL_MODELS.CF4.Omega_b;
    this.Omega_Lambda = params.Omega_Lambda ?? (1.0 - this.Omega_m);
    this.sigma8 = params.sigma8 ?? COSMOLOGICAL_MODELS.CF4.sigma8;
    this.ns = params.ns ?? COSMOLOGICAL_MODELS.CF4.ns;
    this.includeBAO = params.includeBAO ?? true;

    this.calculator = new CosmologicalCalculator({
      H0: this.H0,
      Omega_m: this.Omega_m,
      Omega_Lambda: this.Omega_Lambda,
      sigma8: this.sigma8,
      ns: this.ns
    });

    this.transfer = new EisensteinHuTransfer({
      h: this.h,
      Omega_m: this.Omega_m,
      Omega_b: this.Omega_b
    });

    // Compute normalization amplitude A such that sigma(R=8 Mpc/h) = sigma8
    this.amplitude = this._computeSigma8Normalization();
  }

  /**
   * Spherical Top-Hat Window Function in Fourier space: W(kR) = 3 (sin(x) - x cos(x)) / x^3.
   * Uses high-order Taylor series near x=0 to prevent numerical cancellation.
   * @param {number} x = k * R
   * @returns {number}
   */
  static topHatWindow(x) {
    const ax = Math.abs(x);
    if (ax < 1e-4) {
      const x2 = ax * ax;
      return 1.0 - 0.1 * x2 + (1.0 / 280.0) * x2 * x2 - (1.0 / 15120.0) * x2 * x2 * x2;
    }
    return 3.0 * (Math.sin(ax) - ax * Math.cos(ax)) / (ax * ax * ax);
  }

  /**
   * Gaussian Window Function W(kR) = exp(-0.5 * k^2 * R^2).
   * @param {number} k
   * @param {number} R
   * @returns {number}
   */
  static gaussianWindow(k, R) {
    const x = k * R;
    return Math.exp(-0.5 * x * x);
  }

  /**
   * Evaluates unnormalized raw power spectrum P_raw(k) = k^{n_s} T^2(k).
   * @param {number} k Wavenumber in h * Mpc^-1
   * @returns {number}
   */
  rawPower(k) {
    if (k <= 0.0) return 0.0;
    // Eisenstein-Hu expects k in Mpc^-1, so k_Mpc = k_hMpc * h
    const kMpc = k * this.h;
    const T = this.transfer.computeTransfer(kMpc, this.includeBAO);
    return Math.pow(k, this.ns) * T * T;
  }

  /**
   * Evaluates normalized linear matter power spectrum P(k) = A * k^{n_s} * T^2(k) in (Mpc/h)^3.
   * @param {number} k Wavenumber in h * Mpc^-1
   * @param {number} [z=0.0] Redshift
   * @returns {number} P(k, z) in (Mpc/h)^3
   */
  power(k, z = 0.0) {
    if (k <= 0.0) return 0.0;
    const p0 = this.amplitude * this.rawPower(k);
    if (z === 0.0) return p0;
    const d0 = this.calculator.linearGrowthFactorD(0.0);
    const dz = this.calculator.linearGrowthFactorD(z);
    const growthRatio = dz / d0;
    return p0 * growthRatio * growthRatio;
  }

  /**
   * Dimensionless matter power spectrum Delta^2(k) = k^3 P(k) / (2 pi^2).
   * @param {number} k Wavenumber in h * Mpc^-1
   * @param {number} [z=0.0]
   * @returns {number}
   */
  dimensionlessPower(k, z = 0.0) {
    return (Math.pow(k, 3.0) / (2.0 * Math.PI * Math.PI)) * this.power(k, z);
  }

  /**
   * Computes variance of density fluctuations smoothed on scale R (in Mpc/h):
   * sigma^2(R) = 1/(2 pi^2) * int_0^infty k^2 P(k) W^2(k R) dk
   * @param {number} R Smoothing radius in Mpc/h (e.g. R=8 for sigma8)
   * @param {number} [amplitude=this.amplitude] Power spectrum amplitude
   * @returns {number}
   */
  computeVariance(R, amplitude = this.amplitude) {
    // Gauss-Legendre quadrature in logarithmic wavenumber space ln(k) from k_min = 1e-4 to k_max = 1e3
    const lnKMin = Math.log(1e-4);
    const lnKMax = Math.log(1e3);
    const numSteps = 500;
    const dLnK = (lnKMax - lnKMin) / numSteps;

    let integral = 0.0;
    for (let i = 0; i < numSteps; i++) {
      const lnK1 = lnKMin + i * dLnK;
      const lnK2 = lnK1 + dLnK;
      const lnKMid = 0.5 * (lnK1 + lnK2);

      const k = Math.exp(lnKMid);
      const w = MatterPowerSpectrum.topHatWindow(k * R);
      const pkRaw = this.rawPower(k);
      const integrand = k * k * k * pkRaw * w * w; // extra k because dk = k * d(ln k)

      integral += integrand * dLnK;
    }

    return (amplitude / (2.0 * Math.PI * Math.PI)) * integral;
  }

  /**
   * Determines exact amplitude A by solving sigma^2(8 Mpc/h) = sigma8^2.
   * @private
   * @returns {number}
   */
  _computeSigma8Normalization() {
    const rawVariance8 = this.computeVariance(8.0, 1.0);
    if (rawVariance8 <= 0.0) {
      throw new Error('MatterPowerSpectrum: Invalid zero or negative raw variance.');
    }
    const targetVariance8 = this.sigma8 * this.sigma8;
    return targetVariance8 / rawVariance8;
  }
}

// ============================================================================
// 3D GAUSSIAN RANDOM FIELD GENERATOR (GaussianRandomField3D)
// ============================================================================

/**
 * 3D Gaussian Random Field (GRF) Synthesizer on periodic Cartesian grids.
 * Generates exact Fourier-space realization with Hermitian conjugate symmetry and inverse FFT.
 */
export class GaussianRandomField3D {
  /**
   * @param {Object} config
   * @param {number} [config.boxSize=DEFAULT_BOX_SIZE_MPC] Box length L in Mpc/h
   * @param {number} [config.resolution=DEFAULT_GRID_RESOLUTION] Grid cells per axis (power of 2)
   * @param {MatterPowerSpectrum} [config.powerSpectrum] Matter power spectrum instance
   * @param {number} [config.seed=DEFAULT_MOCK_SEED] Deterministic PRNG seed
   */
  constructor(config = {}) {
    this.boxSize = config.boxSize ?? DEFAULT_BOX_SIZE_MPC;
    this.resolution = config.resolution ?? DEFAULT_GRID_RESOLUTION;
    if ((this.resolution & (this.resolution - 1)) !== 0) {
      throw new Error(`GaussianRandomField3D: Resolution (${this.resolution}) must be a power of 2.`);
    }

    this.powerSpectrum = config.powerSpectrum ?? new MatterPowerSpectrum();
    this.seed = config.seed ?? DEFAULT_MOCK_SEED;
    this.rng = new CosmoRNG(this.seed);

    this.nx = this.resolution;
    this.ny = this.resolution;
    this.nz = this.resolution;
    this.totalCells = this.nx * this.ny * this.nz;
    this.cellSpacing = this.boxSize / this.resolution; // dx = L / N
    this.boxVolume = Math.pow(this.boxSize, 3.0);       // V = L^3
    this.fundamentalK = (2.0 * Math.PI) / this.boxSize; // k_f = 2*pi / L
    this.nyquistK = (Math.PI) / this.cellSpacing;       // k_Ny = pi / dx

    this.fft = new CosmicFFT3D(this.nx, this.ny, this.nz);
  }

  /**
   * Synthesizes 3D Density Contrast delta(x) and Fourier-space representation delta(k).
   * @param {Object} [options]
   * @param {number} [options.smoothingRadius=0.0] Gaussian filter scale in Mpc/h (0 for no filter)
   * @param {number} [options.z=0.0] Redshift
   * @returns {{
   *   deltaGrid: Float64Array,
   *   fourierReal: Float64Array,
   *   fourierImag: Float64Array,
   *   variance: number,
   *   minDelta: number,
   *   maxDelta: number
   * }}
   */
  generateDensityField(options = {}) {
    const smoothingR = options.smoothingRadius ?? 0.0;
    const z = options.z ?? 0.0;

    const fReal = new Float64Array(this.totalCells);
    const fImag = new Float64Array(this.totalCells);

    const halfN = this.resolution >> 1;
    const invSqrt2 = 1.0 / Math.SQRT2;
    const normFactor = Math.sqrt(this.boxVolume); // FFT convention amplitude scaling

    // Generate Fourier modes with Hermitian conjugate symmetry delta(-k) = delta*(k)
    for (let iz = 0; iz < this.nz; iz++) {
      const kzIdx = iz <= halfN ? iz : iz - this.nz;
      const kz = kzIdx * this.fundamentalK;

      for (let iy = 0; iy < this.ny; iy++) {
        const kyIdx = iy <= halfN ? iy : iy - this.ny;
        const ky = kyIdx * this.fundamentalK;

        for (let ix = 0; ix <= halfN; ix++) {
          const kxIdx = ix;
          const kx = kxIdx * this.fundamentalK;

          const idx = (iz * this.ny + iy) * this.nx + ix;
          const k2 = kx * kx + ky * ky + kz * kz;
          const k = Math.sqrt(k2);

          // Zero monopole DC mode
          if (ix === 0 && kyIdx === 0 && kzIdx === 0) {
            fReal[idx] = 0.0;
            fImag[idx] = 0.0;
            continue;
          }

          // Check if this is a Nyquist mode or self-conjugate mode
          const isNyquistX = (ix === 0 || ix === halfN);
          const isNyquistY = (kyIdx === 0 || kyIdx === -halfN);
          const isNyquistZ = (kzIdx === 0 || kzIdx === -halfN);
          const isSelfConjugate = isNyquistX && isNyquistY && isNyquistZ;

          const Pk = this.powerSpectrum.power(k, z);
          let filter = 1.0;
          if (smoothingR > 0.0) {
            filter = Math.exp(-0.5 * k2 * smoothingR * smoothingR);
          }

          const modeAmp = normFactor * Math.sqrt(Pk) * filter;

          if (isSelfConjugate) {
            // Real valued for self-conjugate modes
            fReal[idx] = modeAmp * this.rng.standardNormal();
            fImag[idx] = 0.0;
          } else {
            fReal[idx] = modeAmp * invSqrt2 * this.rng.standardNormal();
            fImag[idx] = modeAmp * invSqrt2 * this.rng.standardNormal();

            // Populate Hermitian conjugate partner delta(-k) = delta*(k)
            const negX = (this.nx - ix) % this.nx;
            const negY = (this.ny - (iy % this.ny)) % this.ny;
            const negZ = (this.nz - (iz % this.nz)) % this.nz;
            const negIdx = (negZ * this.ny + negY) * this.nx + negX;

            fReal[negIdx] = fReal[idx];
            fImag[negIdx] = -fImag[idx];
          }
        }
      }
    }

    // Save exact copy of Fourier modes for displacement / potential calculations
    const fourierRealCopy = new Float64Array(fReal);
    const fourierImagCopy = new Float64Array(fImag);

    // Transform from Fourier to real spatial domain via inverse FFT
    const spatialReal = new Float64Array(fReal);
    const spatialImag = new Float64Array(fImag);
    this.fft.transform(spatialReal, spatialImag, true);

    // Spatial normalization: FFT normalizer is 1/N^3.
    // To match discrete variance with continuous integral, scale by (N / L)^3 / sqrt(V) = 1 / (V * dx^3)
    const physicalScale = Math.pow(this.resolution, 3.0) / Math.sqrt(this.boxVolume);
    let minDelta = Infinity;
    let maxDelta = -Infinity;
    let sumDelta = 0.0;
    let sumDelta2 = 0.0;

    for (let i = 0; i < this.totalCells; i++) {
      spatialReal[i] *= physicalScale;
      const v = spatialReal[i];
      if (v < minDelta) minDelta = v;
      if (v > maxDelta) maxDelta = v;
      sumDelta += v;
      sumDelta2 += v * v;
    }

    const mean = sumDelta / this.totalCells;
    const variance = (sumDelta2 / this.totalCells) - (mean * mean);

    return {
      deltaGrid: spatialReal,
      fourierReal: fourierRealCopy,
      fourierImag: fourierImagCopy,
      variance,
      minDelta,
      maxDelta
    };
  }

  /**
   * Computes empirical 1D power spectrum P_measured(k) from a generated 3D density contrast grid.
   * @param {Float64Array} deltaGrid
   * @param {number} [numBins=30]
   * @returns {{ kCenters: Float64Array, PkMeasured: Float64Array, modeCounts: Int32Array }}
   */
  measurePowerSpectrum(deltaGrid, numBins = 30) {
    const real = new Float64Array(deltaGrid);
    const imag = new Float64Array(this.totalCells);

    // Forward FFT delta(x) -> delta(k)
    this.fft.transform(real, imag, false);

    // Normalization factor: V / N^6 (since forward FFT has no 1/N factor)
    const factor = this.boxVolume / Math.pow(this.totalCells, 2.0);

    const logKMin = Math.log10(this.fundamentalK);
    const logKMax = Math.log10(this.nyquistK);
    const dLogK = (logKMax - logKMin) / numBins;

    const kCenters = new Float64Array(numBins);
    const PkSum = new Float64Array(numBins);
    const modeCounts = new Int32Array(numBins);

    for (let b = 0; b < numBins; b++) {
      kCenters[b] = Math.pow(10.0, logKMin + (b + 0.5) * dLogK);
    }

    const halfN = this.resolution >> 1;
    for (let iz = 0; iz < this.nz; iz++) {
      const kz = (iz <= halfN ? iz : iz - this.nz) * this.fundamentalK;
      for (let iy = 0; iy < this.ny; iy++) {
        const ky = (iy <= halfN ? iy : iy - this.ny) * this.fundamentalK;
        for (let ix = 0; ix < this.nx; ix++) {
          const kx = (ix <= halfN ? ix : ix - this.nx) * this.fundamentalK;
          const k = Math.sqrt(kx * kx + ky * ky + kz * kz);
          if (k < this.fundamentalK || k > this.nyquistK) continue;

          const logK = Math.log10(k);
          const bin = Math.floor((logK - logKMin) / dLogK);
          if (bin >= 0 && bin < numBins) {
            const idx = (iz * this.ny + iy) * this.nx + ix;
            const powerMode = (real[idx] * real[idx] + imag[idx] * imag[idx]) * factor;
            PkSum[bin] += powerMode;
            modeCounts[bin]++;
          }
        }
      }
    }

    const PkMeasured = new Float64Array(numBins);
    for (let b = 0; b < numBins; b++) {
      PkMeasured[b] = modeCounts[b] > 0 ? PkSum[b] / modeCounts[b] : 0.0;
    }

    return { kCenters, PkMeasured, modeCounts };
  }
}

// ============================================================================
// 3D ZEL'DOVICH DISPLACEMENT FIELD & PARTICLES (ZeldovichMockEngine)
// ============================================================================

/**
 * Zel'dovich Approximation (1LPT) Displacement and Velocity Field Engine.
 * 
 * Computes:
 * - Displacement vector Psi(q) = -i (k / k^2) * delta(k) in Fourier space
 * - Eulerian particle positions x(q) = q + D(z) * Psi(q)
 * - Peculiar velocity vectors v(q) = 100 * h * f(z) * Psi(q) (km/s)
 * - Eulerian density assignment (CIC/NGP)
 * - Deformation tensor R_ij = d Psi_i / d q_j and Jacobi eigenvalue caustic analysis
 */
export class ZeldovichMockEngine {
  /**
   * @param {GaussianRandomField3D} grf
   * @param {Object} [params]
   * @param {number} [params.z=0.0] Mock redshift
   */
  constructor(grf, params = {}) {
    this.grf = grf;
    this.z = params.z ?? 0.0;
    this.calculator = grf.powerSpectrum.calculator;
    this.h = grf.powerSpectrum.h;

    // Linear growth factor D(z) and logarithmic growth rate f(z) = d ln D / d ln a
    this.growthD = this.calculator.linearGrowthFactorD(this.z);
    this.growthF = this.calculator.growthRateF(this.z);

    // Peculiar velocity scaling factor: v_scale = 100 * h * f(z) * a / D(z) or linear relation
    // In physical units (km/s): v(q) = 100 * h * f(z) * Psi(q) where Psi is in comoving Mpc/h
    this.velocityFactor = 100.0 * this.h * this.growthF;
  }

  /**
   * Generates 3D Zel'dovich displacement field [Psi_x, Psi_y, Psi_z] from density Fourier modes.
   * @param {Float64Array} fourierReal
   * @param {Float64Array} fourierImag
   * @returns {{
   *   psiX: Float64Array,
   *   psiY: Float64Array,
   *   psiZ: Float64Array,
   *   vx: Float64Array,
   *   vy: Float64Array,
   *   vz: Float64Array,
   *   potential: Float64Array
   * }}
   */
  computeDisplacementAndVelocityFields(fourierReal, fourierImag) {
    const { nx, ny, nz, totalCells, fundamentalK } = this.grf;
    const halfN = nx >> 1;

    const psiXReal = new Float64Array(totalCells);
    const psiXImag = new Float64Array(totalCells);
    const psiYReal = new Float64Array(totalCells);
    const psiYImag = new Float64Array(totalCells);
    const psiZReal = new Float64Array(totalCells);
    const psiZImag = new Float64Array(totalCells);
    const phiReal = new Float64Array(totalCells);
    const phiImag = new Float64Array(totalCells);

    // In Fourier space:
    // Psi(k) = -i * (k / k^2) * delta(k)
    // Phi(k) = -delta(k) / k^2
    for (let iz = 0; iz < nz; iz++) {
      const kz = (iz <= halfN ? iz : iz - nz) * fundamentalK;
      for (let iy = 0; iy < ny; iy++) {
        const ky = (iy <= halfN ? iy : iy - ny) * fundamentalK;
        for (let ix = 0; ix < nx; ix++) {
          const kx = (ix <= halfN ? ix : ix - nx) * fundamentalK;
          const idx = (iz * ny + iy) * nx + ix;
          const k2 = kx * kx + ky * ky + kz * kz;

          if (k2 === 0.0) {
            continue;
          }

          const invK2 = 1.0 / k2;
          const dR = fourierReal[idx];
          const dI = fourierImag[idx];

          // -i * delta(k) = -i * (dR + i*dI) = (dI - i*dR)
          // Multiply by k_comp / k^2:
          const factorX = kx * invK2;
          psiXReal[idx] = dI * factorX;
          psiXImag[idx] = -dR * factorX;

          const factorY = ky * invK2;
          psiYReal[idx] = dI * factorY;
          psiYImag[idx] = -dR * factorY;

          const factorZ = kz * invK2;
          psiZReal[idx] = dI * factorZ;
          psiZImag[idx] = -dR * factorZ;

          // Potential Phi(k) = -delta(k) / k^2
          phiReal[idx] = -dR * invK2;
          phiImag[idx] = -dI * invK2;
        }
      }
    }

    // Inverse FFT transforms for Psi_x, Psi_y, Psi_z, and Phi
    const fft = this.grf.fft;
    fft.transform(psiXReal, psiXImag, true);
    fft.transform(psiYReal, psiYImag, true);
    fft.transform(psiZReal, psiZImag, true);
    fft.transform(phiReal, phiImag, true);

    const physicalScale = Math.pow(nx, 3.0) / Math.sqrt(this.grf.boxVolume);
    const psiX = new Float64Array(totalCells);
    const psiY = new Float64Array(totalCells);
    const psiZ = new Float64Array(totalCells);
    const vx = new Float64Array(totalCells);
    const vy = new Float64Array(totalCells);
    const vz = new Float64Array(totalCells);
    const potential = new Float64Array(totalCells);

    for (let i = 0; i < totalCells; i++) {
      psiX[i] = psiXReal[i] * physicalScale;
      psiY[i] = psiYReal[i] * physicalScale;
      psiZ[i] = psiZReal[i] * physicalScale;
      potential[i] = phiReal[i] * physicalScale;

      vx[i] = this.velocityFactor * psiX[i];
      vy[i] = this.velocityFactor * psiY[i];
      vz[i] = this.velocityFactor * psiZ[i];
    }

    return { psiX, psiY, psiZ, vx, vy, vz, potential };
  }

  /**
   * Generates displaced N-body mock particles from Lagrangian grid points.
   * @param {Object} fields
   * @param {Float64Array} fields.psiX
   * @param {Float64Array} fields.psiY
   * @param {Float64Array} fields.psiZ
   * @param {Float64Array} fields.vx
   * @param {Float64Array} fields.vy
   * @param {Float64Array} fields.vz
   * @returns {{
   *   numParticles: number,
   *   qx: Float64Array,
   *   qy: Float64Array,
   *   qz: Float64Array,
   *   x: Float64Array,
   *   y: Float64Array,
   *   z: Float64Array,
   *   vx: Float64Array,
   *   vy: Float64Array,
   *   vz: Float64Array
   * }}
   */
  generateDisplacedParticles(fields) {
    const { nx, ny, nz, totalCells, cellSpacing, boxSize } = this.grf;
    const qx = new Float64Array(totalCells);
    const qy = new Float64Array(totalCells);
    const qz = new Float64Array(totalCells);
    const x = new Float64Array(totalCells);
    const y = new Float64Array(totalCells);
    const z = new Float64Array(totalCells);
    const vx = new Float64Array(totalCells);
    const vy = new Float64Array(totalCells);
    const vz = new Float64Array(totalCells);

    const halfBox = 0.5 * boxSize;

    for (let iz = 0; iz < nz; iz++) {
      const qzCoord = (iz + 0.5) * cellSpacing - halfBox;
      for (let iy = 0; iy < ny; iy++) {
        const qyCoord = (iy + 0.5) * cellSpacing - halfBox;
        for (let ix = 0; ix < nx; ix++) {
          const qxCoord = (ix + 0.5) * cellSpacing - halfBox;
          const idx = (iz * ny + iy) * nx + ix;

          qx[idx] = qxCoord;
          qy[idx] = qyCoord;
          qz[idx] = qzCoord;

          // Zel'dovich displacement: x(q) = q + D(z) * Psi(q)
          let posX = qxCoord + fields.psiX[idx];
          let posY = qyCoord + fields.psiY[idx];
          let posZ = qzCoord + fields.psiZ[idx];

          // Periodic boundary wrapping to [-L/2, L/2)
          posX = ((posX + halfBox) % boxSize + boxSize) % boxSize - halfBox;
          posY = ((posY + halfBox) % boxSize + boxSize) % boxSize - halfBox;
          posZ = ((posZ + halfBox) % boxSize + boxSize) % boxSize - halfBox;

          x[idx] = posX;
          y[idx] = posY;
          z[idx] = posZ;

          vx[idx] = fields.vx[idx];
          vy[idx] = fields.vy[idx];
          vz[idx] = fields.vz[idx];
        }
      }
    }

    return {
      numParticles: totalCells,
      qx, qy, qz,
      x, y, z,
      vx, vy, vz
    };
  }

  /**
   * Reconstructs Eulerian density grid from displaced particle positions using Cloud-in-Cell (CIC) assignment.
   * @param {{ x: Float64Array, y: Float64Array, z: Float64Array }} particles
   * @returns {Float64Array} Density contrast grid delta_CIC(x)
   */
  assignDensityCIC(particles) {
    const { nx, ny, nz, totalCells, cellSpacing, boxSize } = this.grf;
    const halfBox = 0.5 * boxSize;
    const cicGrid = new Float64Array(totalCells);
    const numParticles = particles.x.length;

    for (let p = 0; p < numParticles; p++) {
      // Map [-L/2, L/2) to grid index units [0, N)
      const u = ((particles.x[p] + halfBox) / cellSpacing) % nx;
      const v = ((particles.y[p] + halfBox) / cellSpacing) % ny;
      const w = ((particles.z[p] + halfBox) / cellSpacing) % nz;

      const normU = (u + nx) % nx;
      const normV = (v + ny) % ny;
      const normW = (w + nz) % nz;

      const i0 = Math.floor(normU);
      const j0 = Math.floor(normV);
      const k0 = Math.floor(normW);

      const i1 = (i0 + 1) % nx;
      const j1 = (j0 + 1) % ny;
      const k1 = (k0 + 1) % nz;

      const dx1 = normU - i0;
      const dy1 = normV - j0;
      const dz1 = normW - k0;

      const dx0 = 1.0 - dx1;
      const dy0 = 1.0 - dy1;
      const dz0 = 1.0 - dz1;

      // 8-cell Cloud-in-Cell trilinear weights
      cicGrid[(k0 * ny + j0) * nx + i0] += dx0 * dy0 * dz0;
      cicGrid[(k0 * ny + j0) * nx + i1] += dx1 * dy0 * dz0;
      cicGrid[(k0 * ny + j1) * nx + i0] += dx0 * dy1 * dz0;
      cicGrid[(k0 * ny + j1) * nx + i1] += dx1 * dy1 * dz0;
      cicGrid[(k1 * ny + j0) * nx + i0] += dx0 * dy0 * dz1;
      cicGrid[(k1 * ny + j0) * nx + i1] += dx1 * dy0 * dz1;
      cicGrid[(k1 * ny + j1) * nx + i0] += dx0 * dy1 * dz1;
      cicGrid[(k1 * ny + j1) * nx + i1] += dx1 * dy1 * dz1;
    }

    const meanCount = numParticles / totalCells;
    const deltaGrid = new Float64Array(totalCells);
    for (let i = 0; i < totalCells; i++) {
      deltaGrid[i] = cicGrid[i] / meanCount - 1.0;
    }

    return deltaGrid;
  }

  /**
   * Caustic and Shell-Crossing Analysis via Deformation Tensor Jacobi Eigenvalue Decomposition.
   * Evaluates Jacobian determinant J = det(I + d Psi_i / d q_j) = (1 - lambda_1)(1 - lambda_2)(1 - lambda_3).
   * @param {Float64Array} psiX
   * @param {Float64Array} psiY
   * @param {Float64Array} psiZ
   * @returns {{
   *   jacobianDeterminant: Float64Array,
   *   causticFraction: number,
   *   pancakeFraction: number,
   *   filamentFraction: number,
   *   clusterFraction: number
   * }}
   */
  computeCausticDiagnostics(psiX, psiY, psiZ) {
    const { nx, ny, nz, totalCells, cellSpacing } = this.grf;
    const inv2Dx = 0.5 / cellSpacing;
    const jacobianDet = new Float64Array(totalCells);

    let causticCount = 0;
    let pancakeCount = 0;
    let filamentCount = 0;
    let clusterCount = 0;

    for (let iz = 0; iz < nz; iz++) {
      const zPrev = ((iz - 1 + nz) % nz) * ny * nx;
      const zNext = ((iz + 1) % nz) * ny * nx;
      const zCurr = iz * ny * nx;

      for (let iy = 0; iy < ny; iy++) {
        const yPrev = ((iy - 1 + ny) % ny) * nx;
        const yNext = ((iy + 1) % ny) * nx;
        const yCurr = iy * nx;

        for (let ix = 0; ix < nx; ix++) {
          const xPrev = (ix - 1 + nx) % nx;
          const xNext = (ix + 1) % nx;
          const idx = zCurr + yCurr + ix;

          // Deformation gradient components d Psi_i / d q_j via central differences
          const rxx = (psiX[zCurr + yCurr + xNext] - psiX[zCurr + yCurr + xPrev]) * inv2Dx;
          const rxy = (psiX[zCurr + yNext + ix] - psiX[zCurr + yPrev + ix]) * inv2Dx;
          const rxz = (psiX[zNext + yCurr + ix] - psiX[zPrev + yCurr + ix]) * inv2Dx;

          const ryx = (psiY[zCurr + yCurr + xNext] - psiY[zCurr + yCurr + xPrev]) * inv2Dx;
          const ryy = (psiY[zCurr + yNext + ix] - psiY[zCurr + yPrev + ix]) * inv2Dx;
          const ryz = (psiY[zNext + yCurr + ix] - psiY[zPrev + yCurr + ix]) * inv2Dx;

          const rzx = (psiZ[zCurr + yCurr + xNext] - psiZ[zCurr + yCurr + xPrev]) * inv2Dx;
          const rzy = (psiZ[zCurr + yNext + ix] - psiZ[zCurr + yPrev + ix]) * inv2Dx;
          const rzz = (psiZ[zNext + yCurr + ix] - psiZ[zPrev + yCurr + ix]) * inv2Dx;

          // Symmetrized strain tensor
          const sxx = rxx;
          const syy = ryy;
          const szz = rzz;
          const sxy = 0.5 * (rxy + ryx);
          const sxz = 0.5 * (rxz + rzx);
          const syz = 0.5 * (ryz + rzy);

          // Jacobi eigenvalue solver for 3x3 symmetric matrix
          const eig = ZeldovichMockEngine.solveSymmetric3x3Eigenvalues(sxx, syy, szz, sxy, sxz, syz);
          // Sort descending lambda1 >= lambda2 >= lambda3
          const l1 = eig[0];
          const l2 = eig[1];
          const l3 = eig[2];

          // 1LPT Jacobian determinant J = (1 + l1)(1 + l2)(1 + l3)
          const J = (1.0 + l1) * (1.0 + l2) * (1.0 + l3);
          jacobianDet[idx] = J;

          if (J <= 0.0 || (1.0 + l1) <= 0.0) {
            causticCount++;
            if ((1.0 + l1) <= 0.0 && (1.0 + l2) > 0.0 && (1.0 + l3) > 0.0) {
              pancakeCount++;
            } else if ((1.0 + l1) <= 0.0 && (1.0 + l2) <= 0.0 && (1.0 + l3) > 0.0) {
              filamentCount++;
            } else if ((1.0 + l1) <= 0.0 && (1.0 + l2) <= 0.0 && (1.0 + l3) <= 0.0) {
              clusterCount++;
            }
          }
        }
      }
    }

    return {
      jacobianDeterminant: jacobianDet,
      causticFraction: causticCount / totalCells,
      pancakeFraction: pancakeCount / totalCells,
      filamentFraction: filamentCount / totalCells,
      clusterFraction: clusterCount / totalCells
    };
  }

  /**
   * Exact analytical eigenvalue solver for 3x3 symmetric matrix via Cardano's cubic formula.
   * @param {number} a Matrix (0,0)
   * @param {number} b Matrix (1,1)
   * @param {number} c Matrix (2,2)
   * @param {number} d Matrix (0,1) = (1,0)
   * @param {number} e Matrix (0,2) = (2,0)
   * @param {number} f Matrix (1,2) = (2,1)
   * @returns {[number, number, number]} Sorted descending eigenvalues [lambda1, lambda2, lambda3]
   */
  static solveSymmetric3x3Eigenvalues(a, b, c, d, e, f) {
    const p1 = d * d + e * e + f * f;
    if (p1 === 0.0) {
      // Diagonal matrix
      const vals = [a, b, c];
      vals.sort((x, y) => y - x);
      return vals;
    }

    const q = (a + b + c) / 3.0;
    const p2 = (a - q) * (a - q) + (b - q) * (b - q) + (c - q) * (c - q) + 2.0 * p1;
    const p = Math.sqrt(p2 / 6.0);

    // Normalized determinant of B = (1/p) * (A - q*I)
    const invP = 1.0 / p;
    const b00 = (a - q) * invP;
    const b11 = (b - q) * invP;
    const b22 = (c - q) * invP;
    const b01 = d * invP;
    const b02 = e * invP;
    const b12 = f * invP;

    let detB = (b00 * (b11 * b22 - b12 * b12) -
      b01 * (b01 * b22 - b12 * b02) +
      b02 * (b01 * b12 - b11 * b02)) * 0.5;

    // Clamp detB to [-1, 1] for acos
    if (detB > 1.0) detB = 1.0;
    else if (detB < -1.0) detB = -1.0;

    const phi = Math.acos(detB) / 3.0;

    // Eigenvalues of B are 2*cos(phi), 2*cos(phi + 2pi/3), 2*cos(phi + 4pi/3)
    const eig1 = q + 2.0 * p * Math.cos(phi);
    const eig3 = q + 2.0 * p * Math.cos(phi + (2.0 * Math.PI / 3.0));
    const eig2 = 3.0 * q - eig1 - eig3; // Trace identity

    const vals = [eig1, eig2, eig3];
    vals.sort((x, y) => y - x);
    return vals;
  }
}

// ============================================================================
// GALAXY TRACER BIAS & CLUSTERING SAMPLER
// ============================================================================

/**
 * Galaxy Clustering and Stochastic Bias Sampling Engine.
 */
export class GalaxyClusteringSampler {
  /**
   * @param {Object} [params]
   * @param {string} [params.model='linear'] 'linear' | 'lognormal' | 'peak_patch'
   * @param {number} [params.bias=1.2] Linear galaxy bias b
   * @param {number} [params.thresholdDelta=1.686] Critical collapse threshold for peak-patch
   */
  constructor(params = {}) {
    this.model = params.model ?? 'linear';
    this.bias = params.bias ?? 1.2;
    this.thresholdDelta = params.thresholdDelta ?? 1.686;
  }

  /**
   * Computes relative galaxy tracer overdensity delta_g(x).
   * @param {number} delta Local matter density contrast
   * @param {number} [variance=1.0] Field variance sigma^2
   * @returns {number} Relative probability / density weight
   */
  evaluateWeight(delta, variance = 1.0) {
    switch (this.model) {
      case 'lognormal': {
        // Continuous log-normal bias: exp(b * delta - 0.5 * b^2 * sigma^2)
        const exponent = this.bias * delta - 0.5 * this.bias * this.bias * variance;
        return Math.exp(Math.max(-10.0, Math.min(10.0, exponent)));
      }
      case 'peak_patch': {
        // High density threshold clustering (Kaiser 1984)
        if (delta >= this.thresholdDelta) {
          return 1.0 + this.bias * (delta - this.thresholdDelta);
        }
        return 0.05;
      }
      case 'linear':
      default: {
        return Math.max(0.01, 1.0 + this.bias * delta);
      }
    }
  }
}

// ============================================================================
// SELECTION FUNCTIONS & SKY FOOTPRINT MASKS
// ============================================================================

/**
 * Radial Selection Function Models n_bar(r).
 */
export class RadialSelectionFunction {
  /**
   * @param {Object} [params]
   * @param {string} [params.type='exponential'] 'exponential' | 'schechter' | 'sigmoid' | 'uniform'
   * @param {number} [params.n0=0.01] Mean central density (galaxies / (Mpc/h)^3)
   * @param {number} [params.r0=80.0] Characteristic survey depth radius in Mpc/h
   * @param {number} [params.gamma=1.8] Selection steepness parameter
   * @param {number} [params.rMax=250.0] Maximum survey boundary radius
   */
  constructor(params = {}) {
    this.type = params.type ?? 'exponential';
    this.n0 = params.n0 ?? 0.01;
    this.r0 = params.r0 ?? 80.0;
    this.gamma = params.gamma ?? 1.8;
    this.rMax = params.rMax ?? 250.0;
  }

  /**
   * Evaluates expected galaxy number density n_bar(r).
   * @param {number} r Comoving distance in Mpc/h
   * @returns {number}
   */
  evaluate(r) {
    if (r < 0.0 || r > this.rMax) return 0.0;

    switch (this.type) {
      case 'exponential':
        // Generalized exponential n_bar(r) = n0 * exp(-(r / r0)^gamma)
        return this.n0 * Math.exp(-Math.pow(r / this.r0, this.gamma));

      case 'sigmoid':
        // Fermi-Dirac sigmoid completeness
        return this.n0 / (1.0 + Math.exp((r - this.r0) / (this.r0 * 0.2)));

      case 'schechter':
        // Apparent magnitude completeness cutoff
        const x = r / this.r0;
        return this.n0 / (1.0 + Math.pow(x, 3.0));

      case 'uniform':
      default:
        return this.n0;
    }
  }

  /**
   * Computes selection probability P(r) = n_bar(r) / n0 in [0, 1].
   * @param {number} r
   * @returns {number}
   */
  selectionProbability(r) {
    if (this.n0 <= 0.0) return 0.0;
    return Math.max(0.0, Math.min(1.0, this.evaluate(r) / this.n0));
  }
}

/**
 * Astrometric Sky Footprint and Zone of Avoidance (ZoA) Dust Extinction Mask.
 */
export class SkyFootprintMask {
  /**
   * @param {Object} [params]
   * @param {number} [params.bCutDeg=10.0] Zone of Avoidance galactic latitude boundary |b| < bCut
   * @param {number} [params.sgbCutDeg=0.0] Supergalactic latitude boundary cut
   * @param {number} [params.decMinDeg=-90.0] Equatorial Declination lower limit
   * @param {number} [params.decMaxDeg=90.0] Equatorial Declination upper limit
   * @param {boolean} [params.enableDustExtinction=true] Include smooth SFD/Planck dust attenuation
   * @param {boolean} [params.enableCorridors=true] Include HI piercing corridors (Vela, Norma, Puppis)
   */
  constructor(params = {}) {
    this.bCutDeg = params.bCutDeg ?? 10.0;
    this.sgbCutDeg = params.sgbCutDeg ?? 0.0;
    this.decMinDeg = params.decMinDeg ?? -90.0;
    this.decMaxDeg = params.decMaxDeg ?? 90.0;
    this.enableDustExtinction = params.enableDustExtinction ?? true;
    this.enableCorridors = params.enableCorridors ?? true;

    // HI Piercing Corridors in Galactic coordinates [l_center, b_center, radius_deg, transmission]
    this.corridors = [
      { name: 'Vela', l: 280.0, b: 5.0, radius: 8.0, trans: 0.85 },
      { name: 'Norma/GA', l: 325.0, b: -7.0, radius: 10.0, trans: 0.80 },
      { name: 'Puppis', l: 245.0, b: 0.0, radius: 7.0, trans: 0.75 }
    ];
  }

  /**
   * Computes sky transmission probability P_sky(l, b, Dec) in [0, 1].
   * @param {number} lDeg Galactic longitude in degrees [0, 360)
   * @param {number} bDeg Galactic latitude in degrees [-90, 90]
   * @param {number} [decDeg=0.0] Equatorial Declination in degrees [-90, 90]
   * @param {number} [sgbDeg=0.0] Supergalactic latitude in degrees [-90, 90]
   * @returns {number} Transmission factor in [0, 1]
   */
  transmission(lDeg, bDeg, decDeg = 0.0, sgbDeg = 0.0) {
    if (decDeg < this.decMinDeg || decDeg > this.decMaxDeg) {
      return 0.0;
    }

    if (this.sgbCutDeg > 0.0 && Math.abs(sgbDeg) < this.sgbCutDeg) {
      return 0.0;
    }

    const absB = Math.abs(bDeg);

    // Check HI piercing corridors within ZoA
    if (this.enableCorridors && absB < this.bCutDeg) {
      for (const corr of this.corridors) {
        const dl = (lDeg - corr.l + 540.0) % 360.0 - 180.0;
        const db = bDeg - corr.b;
        const angDist = Math.sqrt(dl * dl * Math.cos(corr.b * DEG_TO_RAD) * Math.cos(corr.b * DEG_TO_RAD) + db * db);
        if (angDist <= corr.radius) {
          return corr.trans;
        }
      }
    }

    // Standard ZoA geometric step/sigmoid cut
    if (absB < this.bCutDeg) {
      return 0.0;
    }

    if (!this.enableDustExtinction) {
      return 1.0;
    }

    // Analytic cosecant galactic dust extinction profile A_V(b) ~ A_V0 / sin(|b|)
    const sinB = Math.sin(Math.max(5.0 * DEG_TO_RAD, absB * DEG_TO_RAD));
    const Av = Math.min(5.0, 0.05 / sinB);
    return Math.exp(-0.4 * Av); // Flux transmission
  }
}

// ============================================================================
// REALISTIC DISTANCE ERROR & OBSERVATIONAL NOISE MODELS
// ============================================================================

/**
 * Astrometric Distance Modulus and Line-of-Sight Peculiar Velocity Error Model.
 * 
 * Supports:
 * - Tully-Fisher (TF) relation: sigma_d / d ~ 0.18 => sigma_mu ~ 0.39 mag
 * - Fundamental Plane (FP): sigma_d / d ~ 0.20 => sigma_mu ~ 0.43 mag
 * - Type Ia Supernovae (SN Ia): sigma_d / d ~ 0.08 => sigma_mu ~ 0.17 mag
 * - Surface Brightness Fluctuations (SBF): sigma_d / d ~ 0.12 => sigma_mu ~ 0.26 mag
 */
export class DistanceErrorModel {
  /**
   * @param {Object} [params]
   * @param {string} [params.tracerType='TULLY_FISHER'] 'TULLY_FISHER' | 'FUNDAMENTAL_PLANE' | 'SN_IA' | 'SBF' | 'CUSTOM'
   * @param {number} [params.fractionalDistanceError=0.18] sigma_d / d
   * @param {number} [params.velocityDispersionKms=200.0] Thermal / measurement velocity dispersion in km/s
   */
  constructor(params = {}) {
    this.tracerType = params.tracerType ?? 'TULLY_FISHER';
    this.velocityDispersionKms = params.velocityDispersionKms ?? 200.0;

    switch (this.tracerType) {
      case 'FUNDAMENTAL_PLANE':
        this.fractionalDistanceError = params.fractionalDistanceError ?? 0.20;
        break;
      case 'SN_IA':
        this.fractionalDistanceError = params.fractionalDistanceError ?? 0.08;
        break;
      case 'SBF':
        this.fractionalDistanceError = params.fractionalDistanceError ?? 0.12;
        break;
      case 'TULLY_FISHER':
      case 'CUSTOM':
      default:
        this.fractionalDistanceError = params.fractionalDistanceError ?? 0.18;
        break;
    }

    // Exact first-order error propagation: sigma_mu = (5 / ln(10)) * (sigma_d / d)
    this.sigmaMu = FIVE_DIV_LN10 * this.fractionalDistanceError;
  }

  /**
   * Converts comoving distance d (in Mpc) to true distance modulus mu_true.
   * mu = 5 * log10(d) + 25.
   * @param {number} d Distance in Mpc
   * @returns {number}
   */
  static distanceToModulus(d) {
    if (d <= 0.0) throw new Error('Distance must be positive.');
    return 5.0 * Math.log10(d) + 25.0;
  }

  /**
   * Converts distance modulus mu to distance d (in Mpc).
   * d = 10^((mu - 25) / 5).
   * @param {number} mu Distance modulus
   * @returns {number}
   */
  static modulusToDistance(mu) {
    return Math.pow(10.0, (mu - 25.0) / 5.0);
  }

  /**
   * Computes perturbed observed distance modulus, inferred distance, and peculiar velocity.
   * @param {number} dTrue True comoving distance in Mpc (or Mpc/h)
   * @param {number} vRadTrue True line-of-sight peculiar velocity in km/s
   * @param {number} H0 Hubble parameter in km/s/Mpc
   * @param {CosmoRNG} rng PRNG instance
   * @returns {{
   *   muTrue: number,
   *   muObs: number,
   *   sigmaMu: number,
   *   dObs: number,
   *   sigmaD: number,
   *   czObs: number,
   *   vPecObs: number,
   *   sigmaVpec: number
   * }}
   */
  observeGalaxy(dTrue, vRadTrue, H0, rng) {
    const muTrue = DistanceErrorModel.distanceToModulus(dTrue);
    const muNoise = rng.normal(0.0, this.sigmaMu);
    const muObs = muTrue + muNoise;

    const dObs = DistanceErrorModel.modulusToDistance(muObs);
    const sigmaD = dObs * LN10_DIV_5 * this.sigmaMu;

    // Observed recessional velocity cz_obs = H0 * d_true + v_rad_true + thermal_noise
    const czThermalNoise = rng.normal(0.0, this.velocityDispersionKms);
    const czObs = H0 * dTrue + vRadTrue + czThermalNoise;

    // Inferred peculiar velocity v_pec = cz_obs - H0 * d_obs
    const vPecObs = czObs - H0 * dObs;
    const sigmaVpec = Math.sqrt(Math.pow(H0 * sigmaD, 2.0) + Math.pow(this.velocityDispersionKms, 2.0));

    return {
      muTrue,
      muObs,
      sigmaMu: this.sigmaMu,
      dObs,
      sigmaD,
      czObs,
      vPecObs,
      sigmaVpec
    };
  }

  /**
   * Evaluates Homogeneous Malmquist Bias correction factor:
   * E[d_obs | d_true] = d_true * exp(3.5 * Delta^2) where Delta = (ln(10) / 5) * sigma_mu.
   * @param {number} [sigmaMu=this.sigmaMu]
   * @returns {number} Bias multiplier
   */
  homogeneousMalmquistMultiplier(sigmaMu = this.sigmaMu) {
    const delta = LN10_DIV_5 * sigmaMu;
    return Math.exp(3.5 * delta * delta);
  }
}

// ============================================================================
// MOCK GALAXY CATALOG COLUMNAR & RECORD STORE (MockGalaxyCatalog)
// ============================================================================

/**
 * High-performance Columnar and Record-based Storage for Synthetic Galaxy Catalogs.
 */
export class MockGalaxyCatalog {
  /**
   * @param {number} capacity Pre-allocated galaxy capacity
   */
  constructor(capacity = 10000) {
    this.capacity = capacity;
    this.count = 0;

    // Cartesian Coordinates (Mpc/h)
    this.x = new Float64Array(capacity);
    this.y = new Float64Array(capacity);
    this.z = new Float64Array(capacity);

    // Lagrangian Coordinates (Mpc/h)
    this.qx = new Float64Array(capacity);
    this.qy = new Float64Array(capacity);
    this.qz = new Float64Array(capacity);

    // 3D Peculiar Velocity (km/s)
    this.vx = new Float64Array(capacity);
    this.vy = new Float64Array(capacity);
    this.vz = new Float64Array(capacity);
    this.vRadTrue = new Float64Array(capacity);

    // Astrometric Angles (degrees)
    this.ra = new Float64Array(capacity);
    this.dec = new Float64Array(capacity);
    this.l = new Float64Array(capacity);
    this.b = new Float64Array(capacity);
    this.sgl = new Float64Array(capacity);
    this.sgb = new Float64Array(capacity);

    // Distances and Moduli
    this.dTrue = new Float64Array(capacity);
    this.dObs = new Float64Array(capacity);
    this.sigmaD = new Float64Array(capacity);
    this.muTrue = new Float64Array(capacity);
    this.muObs = new Float64Array(capacity);
    this.sigmaMu = new Float64Array(capacity);

    // Velocities and Observables
    this.czObs = new Float64Array(capacity);
    this.vPecObs = new Float64Array(capacity);
    this.sigmaVpec = new Float64Array(capacity);
    this.selectionWeight = new Float64Array(capacity);
    this.densityContrast = new Float64Array(capacity);

    // Meta arrays
    this.ids = new Int32Array(capacity);
    this.tracerTypes = [];
  }

  /**
   * Appends a galaxy record to the columnar store.
   * @param {Object} galaxy
   */
  addGalaxy(galaxy) {
    if (this.count >= this.capacity) {
      this._grow(Math.max(1000, this.capacity * 2));
    }

    const idx = this.count;
    this.ids[idx] = galaxy.id ?? idx;
    this.x[idx] = galaxy.x;
    this.y[idx] = galaxy.y;
    this.z[idx] = galaxy.z;

    this.qx[idx] = galaxy.qx ?? galaxy.x;
    this.qy[idx] = galaxy.qy ?? galaxy.y;
    this.qz[idx] = galaxy.qz ?? galaxy.z;

    this.vx[idx] = galaxy.vx ?? 0.0;
    this.vy[idx] = galaxy.vy ?? 0.0;
    this.vz[idx] = galaxy.vz ?? 0.0;
    this.vRadTrue[idx] = galaxy.vRadTrue ?? 0.0;

    this.ra[idx] = galaxy.ra ?? 0.0;
    this.dec[idx] = galaxy.dec ?? 0.0;
    this.l[idx] = galaxy.l ?? 0.0;
    this.b[idx] = galaxy.b ?? 0.0;
    this.sgl[idx] = galaxy.sgl ?? 0.0;
    this.sgb[idx] = galaxy.sgb ?? 0.0;

    this.dTrue[idx] = galaxy.dTrue;
    this.dObs[idx] = galaxy.dObs;
    this.sigmaD[idx] = galaxy.sigmaD;
    this.muTrue[idx] = galaxy.muTrue;
    this.muObs[idx] = galaxy.muObs;
    this.sigmaMu[idx] = galaxy.sigmaMu;

    this.czObs[idx] = galaxy.czObs;
    this.vPecObs[idx] = galaxy.vPecObs;
    this.sigmaVpec[idx] = galaxy.sigmaVpec;
    this.selectionWeight[idx] = galaxy.selectionWeight ?? 1.0;
    this.densityContrast[idx] = galaxy.densityContrast ?? 0.0;

    this.tracerTypes[idx] = galaxy.tracerType ?? 'TF';

    this.count++;
  }

  /**
   * Resizes internal TypedArrays.
   * @private
   * @param {number} newCapacity
   */
  _grow(newCapacity) {
    const copy = (oldArr) => {
      const n = new Float64Array(newCapacity);
      n.set(oldArr.subarray(0, this.count));
      return n;
    };

    this.x = copy(this.x);
    this.y = copy(this.y);
    this.z = copy(this.z);
    this.qx = copy(this.qx);
    this.qy = copy(this.qy);
    this.qz = copy(this.qz);
    this.vx = copy(this.vx);
    this.vy = copy(this.vy);
    this.vz = copy(this.vz);
    this.vRadTrue = copy(this.vRadTrue);

    this.ra = copy(this.ra);
    this.dec = copy(this.dec);
    this.l = copy(this.l);
    this.b = copy(this.b);
    this.sgl = copy(this.sgl);
    this.sgb = copy(this.sgb);

    this.dTrue = copy(this.dTrue);
    this.dObs = copy(this.dObs);
    this.sigmaD = copy(this.sigmaD);
    this.muTrue = copy(this.muTrue);
    this.muObs = copy(this.muObs);
    this.sigmaMu = copy(this.sigmaMu);

    this.czObs = copy(this.czObs);
    this.vPecObs = copy(this.vPecObs);
    this.sigmaVpec = copy(this.sigmaVpec);
    this.selectionWeight = copy(this.selectionWeight);
    this.densityContrast = copy(this.densityContrast);

    const newIds = new Int32Array(newCapacity);
    newIds.set(this.ids.subarray(0, this.count));
    this.ids = newIds;

    this.capacity = newCapacity;
  }

  /**
   * Extracts record object for galaxy at index i.
   * @param {number} i
   * @returns {Object}
   */
  getGalaxy(i) {
    if (i < 0 || i >= this.count) throw new Error(`Index ${i} out of range [0, ${this.count}).`);
    return {
      id: this.ids[i],
      x: this.x[i], y: this.y[i], z: this.z[i],
      qx: this.qx[i], qy: this.qy[i], qz: this.qz[i],
      vx: this.vx[i], vy: this.vy[i], vz: this.vz[i],
      vRadTrue: this.vRadTrue[i],
      ra: this.ra[i], dec: this.dec[i],
      l: this.l[i], b: this.b[i],
      sgl: this.sgl[i], sgb: this.sgb[i],
      dTrue: this.dTrue[i], dObs: this.dObs[i], sigmaD: this.sigmaD[i],
      muTrue: this.muTrue[i], muObs: this.muObs[i], sigmaMu: this.sigmaMu[i],
      czObs: this.czObs[i], vPecObs: this.vPecObs[i], sigmaVpec: this.sigmaVpec[i],
      selectionWeight: this.selectionWeight[i],
      densityContrast: this.densityContrast[i],
      tracerType: this.tracerTypes[i]
    };
  }

  /**
   * Filters catalog with a predicate function.
   * @param {function(Object): boolean} predicate
   * @returns {MockGalaxyCatalog}
   */
  filter(predicate) {
    const result = new MockGalaxyCatalog(this.count);
    for (let i = 0; i < this.count; i++) {
      const g = this.getGalaxy(i);
      if (predicate(g)) {
        result.addGalaxy(g);
      }
    }
    return result;
  }

  /**
   * Computes volume-averaged and inverse-variance-weighted bulk flow dipole vector V_bulk = [Vx, Vy, Vz].
   * @param {Object} [options]
   * @param {number} [options.rMin=0.0]
   * @param {number} [options.rMax=Infinity]
   * @param {boolean} [options.useTrueVelocities=false]
   * @returns {{
   *   Vx: number,
   *   Vy: number,
   *   Vz: number,
   *   magnitude: number,
   *   numGalaxies: number,
   *   dispersion: number
   * }}
   */
  computeBulkFlow(options = {}) {
    const rMin = options.rMin ?? 0.0;
    const rMax = options.rMax ?? Infinity;
    const useTrue = options.useTrueVelocities ?? false;

    // Multipole estimator via radial velocity projection matrix
    // sum (w_i * v_rad_i * r_hat_i) = A * V_bulk
    let aXX = 0.0, aYY = 0.0, aZZ = 0.0;
    let aXY = 0.0, aXZ = 0.0, aYZ = 0.0;
    let bX = 0.0, bY = 0.0, bZ = 0.0;
    let nSelected = 0;
    let sumV = 0.0, sumV2 = 0.0;

    for (let i = 0; i < this.count; i++) {
      const d = useTrue ? this.dTrue[i] : this.dObs[i];
      if (d < rMin || d > rMax || d <= 0.0) continue;

      const rx = this.x[i] / d;
      const ry = this.y[i] / d;
      const rz = this.z[i] / d;

      const vRad = useTrue ? this.vRadTrue[i] : this.vPecObs[i];
      const err = useTrue ? 50.0 : Math.max(50.0, this.sigmaVpec[i]);
      const w = 1.0 / (err * err);

      aXX += w * rx * rx;
      aYY += w * ry * ry;
      aZZ += w * rz * rz;
      aXY += w * rx * ry;
      aXZ += w * rx * rz;
      aYZ += w * ry * rz;

      bX += w * vRad * rx;
      bY += w * vRad * ry;
      bZ += w * vRad * rz;

      nSelected++;
      sumV += vRad;
      sumV2 += vRad * vRad;
    }

    if (nSelected < 3) {
      return { Vx: 0, Vy: 0, Vz: 0, magnitude: 0, numGalaxies: nSelected, dispersion: 0 };
    }

    // Invert 3x3 symmetric matrix A
    const det = aXX * (aYY * aZZ - aYZ * aYZ) -
      aXY * (aXY * aZZ - aYZ * aXZ) +
      aXZ * (aXY * aYZ - aYY * aXZ);

    if (Math.abs(det) < 1e-15) {
      return { Vx: 0, Vy: 0, Vz: 0, magnitude: 0, numGalaxies: nSelected, dispersion: 0 };
    }

    const invDet = 1.0 / det;
    const invXX = (aYY * aZZ - aYZ * aYZ) * invDet;
    const invYY = (aXX * aZZ - aXZ * aXZ) * invDet;
    const invZZ = (aXX * aYY - aXY * aXY) * invDet;
    const invXY = (aXZ * aYZ - aXY * aZZ) * invDet;
    const invXZ = (aXY * aYZ - aXZ * aYY) * invDet;
    const invYZ = (aXY * aXZ - aXX * aYZ) * invDet;

    const Vx = invXX * bX + invXY * bY + invXZ * bZ;
    const Vy = invXY * bX + invYY * bY + invYZ * bZ;
    const Vz = invXZ * bX + invYZ * bY + invZZ * bZ;
    const magnitude = Math.sqrt(Vx * Vx + Vy * Vy + Vz * Vz);

    const meanV = sumV / nSelected;
    const dispersion = Math.sqrt(Math.max(0.0, (sumV2 / nSelected) - (meanV * meanV)));

    return { Vx, Vy, Vz, magnitude, numGalaxies: nSelected, dispersion };
  }

  /**
   * Exports catalog to CSV formatted string.
   * @returns {string}
   */
  toCSV() {
    const headers = [
      'id', 'x', 'y', 'z', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz', 'v_rad_true',
      'ra_deg', 'dec_deg', 'l_deg', 'b_deg', 'sgl_deg', 'sgb_deg',
      'd_true_mpc', 'd_obs_mpc', 'sigma_d_mpc', 'mu_true', 'mu_obs', 'sigma_mu',
      'cz_obs_kms', 'v_pec_obs_kms', 'sigma_vpec_kms', 'weight', 'delta', 'tracer'
    ];

    const lines = [headers.join(',')];
    for (let i = 0; i < this.count; i++) {
      lines.push([
        this.ids[i],
        this.x[i].toFixed(4), this.y[i].toFixed(4), this.z[i].toFixed(4),
        this.qx[i].toFixed(4), this.qy[i].toFixed(4), this.qz[i].toFixed(4),
        this.vx[i].toFixed(2), this.vy[i].toFixed(2), this.vz[i].toFixed(2), this.vRadTrue[i].toFixed(2),
        this.ra[i].toFixed(5), this.dec[i].toFixed(5),
        this.l[i].toFixed(5), this.b[i].toFixed(5),
        this.sgl[i].toFixed(5), this.sgb[i].toFixed(5),
        this.dTrue[i].toFixed(4), this.dObs[i].toFixed(4), this.sigmaD[i].toFixed(4),
        this.muTrue[i].toFixed(4), this.muObs[i].toFixed(4), this.sigmaMu[i].toFixed(4),
        this.czObs[i].toFixed(2), this.vPecObs[i].toFixed(2), this.sigmaVpec[i].toFixed(2),
        this.selectionWeight[i].toFixed(4), this.densityContrast[i].toFixed(4),
        this.tracerTypes[i]
      ].join(','));
    }

    return lines.join('\n');
  }

  /**
   * Exports catalog to standard GeoJSON FeatureCollection.
   * @returns {Object}
   */
  toGeoJSON() {
    const features = [];
    for (let i = 0; i < this.count; i++) {
      features.push({
        type: 'Feature',
        geometry: {
          type: 'Point',
          coordinates: [this.ra[i] > 180 ? this.ra[i] - 360 : this.ra[i], this.dec[i]]
        },
        properties: {
          id: this.ids[i],
          d_obs: this.dObs[i],
          sigma_d: this.sigmaD[i],
          v_pec: this.vPecObs[i],
          sigma_vpec: this.sigmaVpec[i],
          cz: this.czObs[i],
          l: this.l[i],
          b: this.b[i],
          sgl: this.sgl[i],
          sgb: this.sgb[i],
          tracer: this.tracerTypes[i]
        }
      });
    }

    return {
      type: 'FeatureCollection',
      features
    };
  }
}

// ============================================================================
// ASTROMETRIC FRAME TRANSFORM UTILITIES
// ============================================================================

/**
 * High-precision coordinate frame rotation helpers.
 */
export class AstroCoords {
  /**
   * Converts Cartesian [x, y, z] to Spherical [d, ra_deg, dec_deg] in Equatorial J2000.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {{ d: number, ra: number, dec: number }}
   */
  static cartesianToEquatorial(x, y, z) {
    const d = Math.sqrt(x * x + y * y + z * z);
    if (d <= 0.0) return { d: 0.0, ra: 0.0, dec: 0.0 };

    let ra = Math.atan2(y, x) * RAD_TO_DEG;
    if (ra < 0.0) ra += 360.0;
    const dec = Math.asin(Math.max(-1.0, Math.min(1.0, z / d))) * RAD_TO_DEG;

    return { d, ra, dec };
  }

  /**
   * Converts Equatorial (RA, Dec) to Galactic coordinates (l, b).
   * @param {number} raDeg
   * @param {number} decDeg
   * @returns {{ l: number, b: number }}
   */
  static equatorialToGalactic(raDeg, decDeg) {
    const raNGP = ASTROMETRIC_CONSTANTS.RA_NGP_DEG * DEG_TO_RAD;
    const decNGP = ASTROMETRIC_CONSTANTS.DEC_NGP_DEG * DEG_TO_RAD;
    const l0 = ASTROMETRIC_CONSTANTS.L0_GAL_DEG * DEG_TO_RAD;

    const ra = raDeg * DEG_TO_RAD;
    const dec = decDeg * DEG_TO_RAD;

    const sinDec = Math.sin(dec);
    const cosDec = Math.cos(dec);
    const sinDecNGP = Math.sin(decNGP);
    const cosDecNGP = Math.cos(decNGP);

    const sinB = sinDec * sinDecNGP + cosDec * cosDecNGP * Math.cos(ra - raNGP);
    const b = Math.asin(Math.max(-1.0, Math.min(1.0, sinB)));

    const y = cosDec * Math.sin(ra - raNGP);
    const x = sinDec * cosDecNGP - cosDec * sinDecNGP * Math.cos(ra - raNGP);
    let l = (l0 - Math.atan2(y, x)) * RAD_TO_DEG;
    l = ((l % 360.0) + 360.0) % 360.0;

    return { l, b: b * RAD_TO_DEG };
  }

  /**
   * Converts Galactic (l, b) to Supergalactic coordinates (SGL, SGB).
   * @param {number} lDeg
   * @param {number} bDeg
   * @returns {{ sgl: number, sgb: number }}
   */
  static galacticToSupergalactic(lDeg, bDeg) {
    const lSGP = ASTROMETRIC_CONSTANTS.L_SGP_DEG * DEG_TO_RAD;
    const bSGP = ASTROMETRIC_CONSTANTS.B_SGP_DEG * DEG_TO_RAD;
    const sgl0 = ASTROMETRIC_CONSTANTS.SGL0_DEG * DEG_TO_RAD;

    const l = lDeg * DEG_TO_RAD;
    const b = bDeg * DEG_TO_RAD;

    const sinB = Math.sin(b);
    const cosB = Math.cos(b);
    const sinBSGP = Math.sin(bSGP);
    const cosBSGP = Math.cos(bSGP);

    const sinSGB = sinB * sinBSGP + cosB * cosBSGP * Math.cos(l - lSGP);
    const sgb = Math.asin(Math.max(-1.0, Math.min(1.0, sinSGB)));

    const y = cosB * Math.sin(l - lSGP);
    const x = sinB * cosBSGP - cosB * sinBSGP * Math.cos(l - lSGP);
    let sgl = (sgl0 + Math.atan2(y, x)) * RAD_TO_DEG;
    sgl = ((sgl % 360.0) + 360.0) % 360.0;

    return { sgl, sgb: sgb * RAD_TO_DEG };
  }
}

// ============================================================================
// COSMICFLOWS MOCK CATALOG GENERATOR (Master Engine)
// ============================================================================

/**
 * Master Synthetic Cosmological Flow & Galaxy Mock Catalog Generator.
 * Unifies GRF field synthesis, Zel'dovich dynamics, galaxy bias sampling, observational selection,
 * and realistic observational error propagation.
 */
export class CosmicflowsMockGenerator {
  /**
   * @param {Object} [config]
   * @param {number} [config.boxSize=DEFAULT_BOX_SIZE_MPC] Box size L in Mpc/h
   * @param {number} [config.resolution=DEFAULT_GRID_RESOLUTION] Grid resolution (power of 2)
   * @param {number} [config.H0=74.6] Hubble constant in km/s/Mpc
   * @param {number} [config.Omega_m=0.31] Matter density
   * @param {number} [config.Omega_b=0.048] Baryon density
   * @param {number} [config.sigma8=0.81] Amplitude of mass fluctuations
   * @param {number} [config.ns=0.965] Spectral index
   * @param {number} [config.seed=DEFAULT_MOCK_SEED] Deterministic PRNG seed
   * @param {number} [config.galaxyBias=1.2] Linear galaxy bias b
   * @param {RadialSelectionFunction} [config.selectionFunction]
   * @param {SkyFootprintMask} [config.skyMask]
   * @param {DistanceErrorModel} [config.errorModel]
   */
  constructor(config = {}) {
    this.boxSize = config.boxSize ?? DEFAULT_BOX_SIZE_MPC;
    this.resolution = config.resolution ?? DEFAULT_GRID_RESOLUTION;
    this.H0 = config.H0 ?? COSMOLOGICAL_MODELS.CF4.H0;
    this.Omega_m = config.Omega_m ?? COSMOLOGICAL_MODELS.CF4.Omega_m;
    this.Omega_b = config.Omega_b ?? COSMOLOGICAL_MODELS.CF4.Omega_b;
    this.sigma8 = config.sigma8 ?? COSMOLOGICAL_MODELS.CF4.sigma8;
    this.ns = config.ns ?? COSMOLOGICAL_MODELS.CF4.ns;
    this.seed = config.seed ?? DEFAULT_MOCK_SEED;
    this.galaxyBias = config.galaxyBias ?? 1.2;

    this.rng = new CosmoRNG(this.seed);

    this.powerSpectrum = new MatterPowerSpectrum({
      H0: this.H0,
      Omega_m: this.Omega_m,
      Omega_b: this.Omega_b,
      sigma8: this.sigma8,
      ns: this.ns
    });

    this.grf = new GaussianRandomField3D({
      boxSize: this.boxSize,
      resolution: this.resolution,
      powerSpectrum: this.powerSpectrum,
      seed: this.seed
    });

    this.zeldovich = new ZeldovichMockEngine(this.grf, { z: 0.0 });
    this.clusteringSampler = new GalaxyClusteringSampler({
      model: 'linear',
      bias: this.galaxyBias
    });

    this.selectionFunction = config.selectionFunction ?? new RadialSelectionFunction({
      type: 'exponential',
      n0: 0.01,
      r0: 100.0,
      gamma: 1.8,
      rMax: this.boxSize * 0.45
    });

    this.skyMask = config.skyMask ?? new SkyFootprintMask({
      bCutDeg: 10.0,
      enableDustExtinction: true,
      enableCorridors: true
    });

    this.errorModel = config.errorModel ?? new DistanceErrorModel({
      tracerType: 'TULLY_FISHER',
      fractionalDistanceError: 0.18,
      velocityDispersionKms: 200.0
    });
  }

  /**
   * Generates complete synthetic galaxy catalog realization with all observational effects.
   * @param {Object} [options]
   * @param {number} [options.smoothingRadius=2.0] Pre-smoothing scale in Mpc/h
   * @param {number} [options.targetGalaxyCount=5000] Target number of observed mock galaxies
   * @param {number[]} [options.observerPosition=[0,0,0]] Observer position in Mpc/h
   * @returns {{
   *   catalog: MockGalaxyCatalog,
   *   densityGrid: Float64Array,
   *   particles: Object,
   *   fields: Object,
   *   caustics: Object,
   *   summary: Object
   * }}
   */
  generateMockSurvey(options = {}) {
    const smoothingR = options.smoothingRadius ?? 2.0;
    const targetCount = options.targetGalaxyCount ?? 5000;
    const obsPos = options.observerPosition ?? [0.0, 0.0, 0.0];

    // 1. Synthesize 3D Gaussian Random Field
    const densityResult = this.grf.generateDensityField({ smoothingRadius: smoothingR, z: 0.0 });

    // 2. Compute Zel'dovich displacement and velocity fields
    const fields = this.zeldovich.computeDisplacementAndVelocityFields(
      densityResult.fourierReal,
      densityResult.fourierImag
    );

    // 3. Displace N-body particles
    const particles = this.zeldovich.generateDisplacedParticles(fields);

    // 4. Compute caustics and shell-crossing
    const caustics = this.zeldovich.computeCausticDiagnostics(fields.psiX, fields.psiY, fields.psiZ);

    // 5. Sample mock galaxies from particles with density bias, radial selection, and sky mask
    const catalog = new MockGalaxyCatalog(targetCount * 2);
    const numParticles = particles.numParticles;

    for (let p = 0; p < numParticles; p++) {
      const px = particles.x[p];
      const py = particles.y[p];
      const pz = particles.z[p];

      // Distance from observer
      const dx = px - obsPos[0];
      const dy = py - obsPos[1];
      const dz = pz - obsPos[2];
      const dTrue = Math.sqrt(dx * dx + dy * dy + dz * dz);

      if (dTrue <= 0.01 || dTrue > this.selectionFunction.rMax) continue;

      // Radial unit vector
      const rx = dx / dTrue;
      const ry = dy / dTrue;
      const rz = dz / dTrue;

      // Local density contrast delta at Lagrangian grid cell
      const localDelta = densityResult.deltaGrid[p];

      // Galaxy clustering probability with linear / exponential bias: P_cluster ~ (1 + b * delta)
      const biasedWeight = this.clusteringSampler.evaluateWeight(localDelta, densityResult.variance);

      // Radial selection probability
      const pRadial = this.selectionFunction.selectionProbability(dTrue);

      // Astrometric coordinates
      const eq = AstroCoords.cartesianToEquatorial(dx, dy, dz);
      const gal = AstroCoords.equatorialToGalactic(eq.ra, eq.dec);
      const sg = AstroCoords.galacticToSupergalactic(gal.l, gal.b);

      // Sky mask transmission probability
      const pSky = this.skyMask.transmission(gal.l, gal.b, eq.dec, sg.sgb);

      // Combined galaxy inclusion probability
      const pTotal = Math.min(1.0, biasedWeight * pRadial * pSky * 0.5);

      if (this.rng.uniform() < pTotal) {
        // True radial peculiar velocity
        const vRadTrue = particles.vx[p] * rx + particles.vy[p] * ry + particles.vz[p] * rz;

        // Apply observational distance errors and velocity noise
        const obs = this.errorModel.observeGalaxy(dTrue, vRadTrue, this.H0, this.rng);

        catalog.addGalaxy({
          id: p,
          x: px, y: py, z: pz,
          qx: particles.qx[p], qy: particles.qy[p], qz: particles.qz[p],
          vx: particles.vx[p], vy: particles.vy[p], vz: particles.vz[p],
          vRadTrue,
          ra: eq.ra, dec: eq.dec,
          l: gal.l, b: gal.b,
          sgl: sg.sgl, sgb: sg.sgb,
          dTrue,
          dObs: obs.dObs,
          sigmaD: obs.sigmaD,
          muTrue: obs.muTrue,
          muObs: obs.muObs,
          sigmaMu: obs.sigmaMu,
          czObs: obs.czObs,
          vPecObs: obs.vPecObs,
          sigmaVpec: obs.sigmaVpec,
          selectionWeight: 1.0 / Math.max(1e-4, pRadial * pSky),
          densityContrast: localDelta,
          tracerType: this.errorModel.tracerType
        });
      }
    }

    // Compute bulk flow
    const bulkFlow = catalog.computeBulkFlow({ useTrueVelocities: false });
    const bulkFlowTrue = catalog.computeBulkFlow({ useTrueVelocities: true });

    return {
      catalog,
      densityGrid: densityResult.deltaGrid,
      particles,
      fields,
      caustics,
      summary: {
        galaxyCount: catalog.count,
        bulkFlowObs: bulkFlow,
        bulkFlowTrue: bulkFlowTrue,
        causticFraction: caustics.causticFraction,
        variance: densityResult.variance,
        boxSize: this.boxSize,
        resolution: this.resolution
      }
    };
  }
}

// ============================================================================
// BLIND BENCHMARK GENERATOR FOR VERIFICATION SUITE
// ============================================================================

/**
 * Blind Cosmological Benchmark Generator for Reconstruction Challenge Audits.
 * Injects secret known physical attractors, bulk flows, and expansion perturbations,
 * producing paired blind catalogs (without truth) and oracle ground-truth models.
 */
export class BlindBenchmarkGenerator {
  /**
   * @param {Object} [config]
   * @param {number} [config.seed=7777777] Deterministic benchmark seed
   * @param {number} [config.boxSize=400.0]
   * @param {number} [config.resolution=32]
   */
  constructor(config = {}) {
    this.seed = config.seed ?? 7777777;
    this.boxSize = config.boxSize ?? 400.0;
    this.resolution = config.resolution ?? 32;

    this.generator = new CosmicflowsMockGenerator({
      boxSize: this.boxSize,
      resolution: this.resolution,
      seed: this.seed
    });
  }

  /**
   * Generates a blind challenge mock containing hidden injected physical signals.
   * @param {Object} [injectionParams]
   * @param {number[]} [injectionParams.attractorCenter=[40, 30, -20]] Center in Mpc/h
   * @param {number} [injectionParams.attractorMass=5.0e15] Solar masses
   * @param {number[]} [injectionParams.bulkDipole=[150, -80, 220]] Injected bulk flow in km/s
   * @returns {{
   *   blindCatalog: MockGalaxyCatalog,
   *   groundTruth: {
   *     injectedAttractor: Object,
   *     injectedBulkDipole: number[],
   *     trueDensityGrid: Float64Array,
   *     truePotentialGrid: Float64Array,
   *     trueBulkFlow: Object
   *   }
   * }}
   */
  generateBlindChallenge(injectionParams = {}) {
    const center = injectionParams.attractorCenter ?? [40.0, 30.0, -20.0];
    const mass = injectionParams.attractorMass ?? 5.0e15;
    const bulkDipole = injectionParams.bulkDipole ?? [150.0, -80.0, 220.0];

    const survey = this.generator.generateMockSurvey({ targetGalaxyCount: 3000 });
    const cat = survey.catalog;

    // Inject secret cluster gravitational infall and bulk flow shift
    const G_INJECTION = 4.30091e-4; // Astro units
    for (let i = 0; i < cat.count; i++) {
      // 1. Bulk flow dipole shift
      cat.vx[i] += bulkDipole[0];
      cat.vy[i] += bulkDipole[1];
      cat.vz[i] += bulkDipole[2];

      // 2. Plummer attractor gravitational infall: v_infall = - (G * M * r_vec) / (r^2 + a^2)^(3/2) * (t_dyn)
      const dx = cat.x[i] - center[0];
      const dy = cat.y[i] - center[1];
      const dz = cat.z[i] - center[2];
      const r2 = dx * dx + dy * dy + dz * dz;
      const coreA2 = 225.0; // 15 Mpc/h core
      const infallAmp = (G_INJECTION * (mass / 1e10)) / Math.pow(r2 + coreA2, 1.5) * 50.0;

      const vInfX = -infallAmp * dx;
      const vInfY = -infallAmp * dy;
      const vInfZ = -infallAmp * dz;

      cat.vx[i] += vInfX;
      cat.vy[i] += vInfY;
      cat.vz[i] += vInfZ;

      // Update radial velocities and observed peculiar velocities
      const d = cat.dTrue[i];
      if (d > 0.0) {
        cat.vRadTrue[i] = (cat.vx[i] * cat.x[i] + cat.vy[i] * cat.y[i] + cat.vz[i] * cat.z[i]) / d;
        cat.czObs[i] = this.generator.H0 * d + cat.vRadTrue[i] + this.generator.rng.normal(0, 150);
        cat.vPecObs[i] = cat.czObs[i] - this.generator.H0 * cat.dObs[i];
      }
    }

    const groundTruthBulkFlow = cat.computeBulkFlow({ useTrueVelocities: true });

    return {
      blindCatalog: cat,
      groundTruth: {
        injectedAttractor: { center, mass },
        injectedBulkDipole: bulkDipole,
        trueDensityGrid: survey.densityGrid,
        truePotentialGrid: survey.fields.potential,
        trueBulkFlow: groundTruthBulkFlow
      }
    };
  }

  /**
   * Evaluates blind bulk flow recovery accuracy.
   * @param {Object} estimatedFlow { Vx, Vy, Vz, magnitude }
   * @param {Object} truthFlow { Vx, Vy, Vz, magnitude }
   * @returns {{
   *   velocityRMSE: number,
   *   magnitudeErrorKms: number,
   *   alignmentAngleDeg: number,
   *   passed: boolean
   * }}
   */
  static evaluateBulkFlowRecovery(estimatedFlow, truthFlow) {
    const dVx = estimatedFlow.Vx - truthFlow.Vx;
    const dVy = estimatedFlow.Vy - truthFlow.Vy;
    const dVz = estimatedFlow.Vz - truthFlow.Vz;
    const velocityRMSE = Math.sqrt((dVx * dVx + dVy * dVy + dVz * dVz) / 3.0);

    const magnitudeErrorKms = Math.abs(estimatedFlow.magnitude - truthFlow.magnitude);

    const dot = estimatedFlow.Vx * truthFlow.Vx + estimatedFlow.Vy * truthFlow.Vy + estimatedFlow.Vz * truthFlow.Vz;
    const magProduct = estimatedFlow.magnitude * truthFlow.magnitude;
    let cosTheta = magProduct > 0.0 ? dot / magProduct : 1.0;
    cosTheta = Math.max(-1.0, Math.min(1.0, cosTheta));
    const alignmentAngleDeg = Math.acos(cosTheta) * RAD_TO_DEG;

    return {
      velocityRMSE,
      magnitudeErrorKms,
      alignmentAngleDeg,
      passed: alignmentAngleDeg < 25.0 && magnitudeErrorKms < 150.0
    };
  }

  /**
   * Evaluates Pearson correlation coefficient r between reconstructed and ground-truth 3D density fields.
   * @param {Float64Array} reconGrid
   * @param {Float64Array} truthGrid
   * @returns {{ pearsonR: number, rmse: number, passed: boolean }}
   */
  static evaluateDensityFieldCorrelation(reconGrid, truthGrid) {
    if (reconGrid.length !== truthGrid.length) {
      throw new Error('Grid dimensions must match.');
    }
    const n = reconGrid.length;
    let sumR = 0, sumT = 0, sumR2 = 0, sumT2 = 0, sumRT = 0, sumSqDiff = 0;

    for (let i = 0; i < n; i++) {
      const r = reconGrid[i];
      const t = truthGrid[i];
      sumR += r;
      sumT += t;
      sumR2 += r * r;
      sumT2 += t * t;
      sumRT += r * t;
      sumSqDiff += (r - t) * (r - t);
    }

    const meanR = sumR / n;
    const meanT = sumT / n;
    const varR = (sumR2 / n) - meanR * meanR;
    const varT = (sumT2 / n) - meanT * meanT;
    const cov = (sumRT / n) - meanR * meanT;

    const denom = Math.sqrt(Math.max(1e-15, varR * varT));
    const pearsonR = cov / denom;
    const rmse = Math.sqrt(sumSqDiff / n);

    return {
      pearsonR,
      rmse,
      passed: pearsonR > 0.70
    };
  }
}
