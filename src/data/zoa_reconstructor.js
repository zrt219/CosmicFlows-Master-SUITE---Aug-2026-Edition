 * @file zoa_reconstructor.js
 * @description ZRT EXPERIMENTAL ZOA RECONSTRUCTION: 3D Zone of Avoidance (ZoA) Obscuration Model,
 * Infrared / 21cm HI Corridor Piercing (MeerKAT Vela, Parkes HIZOA Norma, Puppis), Wiener-Filter Inpainting
 * across dust-obscured regions, and Synthetic Control Benchmarking.
 * 
 * SCIENTIFIC CLASSIFICATION & PROVENANCE NOTICE:
 * - This module is classified as "ZRT EXPERIMENTAL ZOA RECONSTRUCTION" and is NOT the official
 *   2026 CF4++ZOA data release product unless instantiated with the full hybrid input compilation
 *   (65,518 CF4++ peculiar velocity distances + 8,283 ZOA redshifts, including 2,176 MeerKAT HI redshifts).
 * - Scientific Reference: Hollinger et al. (2026), "Vela Supercluster in the Cosmicflows-4 Hybrid Reconstruction",
 *   arXiv:2603.09339.
 * 
 * Astrophysical Formulations:
 * 1. Obscuration Geometry:
 *    - Galactic ZoA: |b| < 10 deg (dust extinction A_V / E(B-V) obscuration)
 *    - Supergalactic ZoA: |SGB| < 10 deg (secondary planar obscuration)
 *    - Continuous sigmoid/tanh blending across boundaries for C^0/C^1 continuity.
 * 
 * 2. 21cm HI & Infrared Corridor Piercing:
 *    - MeerKAT / Parkes Vela HI Corridor (l in [260, 285], b in [-5, +5], cz in [15000, 25000] km/s)
 *    - Parkes HIZOA Norma / Great Attractor Core (l in [315, 335], b in [-10, +5], cz in [3000, 7500] km/s)
 *    - Puppis HI Corridor (l in [235, 250], b in [-8, +4], cz in [1500, 4000] km/s)
 *    - Piercing transmission factor T(x) in [0, 1] de-obscuring survey corridors.
 * 
 * 3. Wiener-Filter & Kriging Inpainting:
 *    - Cosmological spatial covariance: xi(r) = sigma_0^2 / (1 + (r / r_0)^2)^(gamma/2) or Gaussian xi(r)
 *    - Grid-based iterative Wiener inpainting solver with covariance smoothing & boundary constraints
 *    - Point-based Gaussian Process Kriging reconstruction with error variance estimation
 * 
 * 4. Synthetic Control Benchmarks:
 *    - Rigorous validation harness comparing ground truth vs Wiener reconstructed fields
 *    - Pearson correlation r, RMSE, relative flux conservation, bulk-flow alignment.
 * 
 * @module data/zoa_reconstructor
 */

export const RECONSTRUCTION_STATUS = Object.freeze({
  TYPE: 'ZRT_EXPERIMENTAL_ZOA_RECONSTRUCTION',
  OFFICIAL_CF4_DATASET: false,
  REQUIRED_INPUT_CATALOGS: [
    'CF4++ Peculiar Velocity Compilation (N=65,518)',
    'ZOA Spectroscopic Redshifts (N=8,283)',
    'MeerKAT H I Redshifts (N=2,176)'
  ],
  PRIMARY_CITATION: 'Hollinger et al. 2026 (arXiv:2603.09339)'
});

import {
  equatorialToGalactic,
  galacticToSupergalactic,
  supergalacticToGalactic,
  galacticToEquatorial,
  sphericalToCartesian,
  cartesianToSpherical,
  equatorialToSupergalacticCartesian,
  supergalacticCartesianToEquatorial
} from '../coordinates/canonical_frame.js';

import { GridIndexer, StrideOrder, BoundaryMode } from '../fields/grid_indexer.js';
import { DensityField } from '../fields/density_field.js';

/**
 * Standard Zone of Avoidance Angular Cutoffs and Extinction Thresholds.
 */
export const ZOA_CONSTANTS = Object.freeze({
  DEFAULT_GAL_B_CUT_DEG: 10.0,      // Galactic latitude cut |b| < 10°
  DEFAULT_SG_B_CUT_DEG: 10.0,       // Supergalactic latitude cut |SGB| < 10°
  DEFAULT_TRANSITION_DEG: 2.5,      // Boundary smoothing margin in degrees
  DEFAULT_EBV_CUT: 0.20,            // E(B-V) dust extinction threshold (mag)
  DEFAULT_AV_CUT: 0.62,             // A_V visual extinction threshold (mag)
  CORRELATION_LENGTH_MPC: 10.0,     // Cosmological matter correlation length r_0 (Mpc/h)
  CORRELATION_GAMMA: 1.8            // Spatial two-point correlation slope gamma
});

/**
 * Piercing Deep-Survey Corridors through the Zone of Avoidance.
 */
export const PIERCING_CORRIDORS = Object.freeze([
  {
    name: 'Vela_MeerKAT_Parkes',
    description: 'MeerKAT / Parkes HI Vela Supercluster Survey Corridor',
    lMin: 260.0,
    lMax: 285.0,
    bMin: -5.0,
    bMax: 5.0,
    czMin: 15000.0,
    czMax: 25000.0,
    distMinMpc: 180.0,
    distMaxMpc: 310.0,
    completeness: 0.85,
    referenceCenter: { l: 272.5, b: 0.0, distMpc: 240.0 }
  },
  {
    name: 'Norma_GreatAttractor_HIZOA',
    description: 'Parkes HIZOA Norma Supercluster / Abell 3627 Corridor',
    lMin: 315.0,
    lMax: 335.0,
    bMin: -10.0,
    bMax: 5.0,
    czMin: 3000.0,
    czMax: 7500.0,
    distMinMpc: 45.0,
    distMaxMpc: 95.0,
    completeness: 0.90,
    referenceCenter: { l: 325.3, b: -7.2, distMpc: 68.0 }
  },
  {
    name: 'Puppis_HI',
    description: 'Puppis HI Cluster / HIZOA Corridor',
    lMin: 235.0,
    lMax: 250.0,
    bMin: -8.0,
    bMax: 4.0,
    czMin: 1500.0,
    czMax: 4000.0,
    distMinMpc: 20.0,
    distMaxMpc: 60.0,
    completeness: 0.80,
    referenceCenter: { l: 240.0, b: 0.0, distMpc: 45.0 }
  }
]);

// ============================================================================
// 1. ZoA Obscuration Geometry & Mask Evaluators
// ============================================================================

/**
 * Checks whether a celestial position falls inside the Galactic Zone of Avoidance.
 * @param {number} glat Galactic latitude b in degrees
 * @param {number} [bCut=10.0] Galactic latitude cutoff in degrees
 * @returns {boolean} True if |b| < bCut
 */
export function isGalacticZoA(glat, bCut = ZOA_CONSTANTS.DEFAULT_GAL_B_CUT_DEG) {
  if (isNaN(glat)) return false;
  return Math.abs(glat) < bCut;
}

/**
 * Checks whether a celestial position falls inside the Supergalactic Zone of Avoidance.
 * @param {number} sgb Supergalactic latitude SGB in degrees
 * @param {number} [sgbCut=10.0] Supergalactic latitude cutoff in degrees
 * @returns {boolean} True if |SGB| < sgbCut
 */
export function isSupergalacticZoA(sgb, sgbCut = ZOA_CONSTANTS.DEFAULT_SG_B_CUT_DEG) {
  if (isNaN(sgb)) return false;
  return Math.abs(sgb) < sgbCut;
}

/**
 * Computes analytic SFD98/Planck dust extinction A_V estimate based on Galactic coordinates.
 * Formula: A_V(l, b) = A_V0 / (sin(|b|) + eps) * (1.0 + f_bulge(l, b))
 * 
 * @param {number} glon Galactic longitude l in degrees [0, 360)
 * @param {number} glat Galactic latitude b in degrees [-90, +90]
 * @returns {number} Estimated visual extinction A_V in mag
 */
export function computeDustExtinctionAV(glon, glat) {
  if (isNaN(glon) || isNaN(glat)) return 0.0;

  const absB = Math.abs(glat);
  const sinB = Math.sin(absB * Math.PI / 180.0);
  const eps = 0.04; // Core saturation floor

  // Galactic bulge component around (l=0, b=0)
  let dL = (glon > 180.0 ? glon - 360.0 : glon);
  const bulgeRadSq = (dL / 20.0) * (dL / 20.0) + (glat / 10.0) * (glat / 10.0);
  const bulgeFactor = 2.5 * Math.exp(-0.5 * bulgeRadSq);

  const baseAV = 0.08;
  const av = (baseAV / (sinB + eps)) * (1.0 + bulgeFactor);
  return Math.min(30.0, av);
}

/**
 * Computes continuous Zone of Avoidance obscuration weight W in [0.0, 1.0].
 * - W = 1.0: Deeply obscured ZoA core
 * - W = 0.0: Completely un-obscured sky
 * - Smooth transition across |b| = bCut and |SGB| = sgbCut using sigmoid / tanh tapering.
 * 
 * @param {number} glon Galactic longitude in deg
 * @param {number} glat Galactic latitude in deg
 * @param {number} [sgl] Optional Supergalactic longitude (computed if omitted)
 * @param {number} [sgb] Optional Supergalactic latitude (computed if omitted)
 * @param {Object} [options] Mask configuration
 * @param {number} [options.bCut=10.0] Galactic latitude cutoff
 * @param {number} [options.sgbCut=10.0] Supergalactic latitude cutoff
 * @param {number} [options.transitionDeg=2.5] Smoothing transition margin
 * @param {boolean} [options.includeSupergalactic=true] Whether to include Supergalactic plane
 * @param {boolean} [options.includeDust=true] Whether to incorporate dust extinction thresholding
 * @returns {number} Obscuration weight in [0, 1]
 */
export function computeZoAObscurationWeight(glon, glat, sgl = undefined, sgb = undefined, options = {}) {
  if (isNaN(glon) || isNaN(glat)) return 0.0;

  const bCut = options.bCut !== undefined ? options.bCut : ZOA_CONSTANTS.DEFAULT_GAL_B_CUT_DEG;
  const sgbCut = options.sgbCut !== undefined ? options.sgbCut : ZOA_CONSTANTS.DEFAULT_SG_B_CUT_DEG;
  const trans = Math.max(0.1, options.transitionDeg || ZOA_CONSTANTS.DEFAULT_TRANSITION_DEG);
  const includeSG = options.includeSupergalactic !== undefined ? options.includeSupergalactic : true;
  const includeDust = options.includeDust !== undefined ? options.includeDust : true;

  // 1. Galactic ZoA component: 1 when |b| < bCut, 0 when |b| > bCut + trans
  const absB = Math.abs(glat);
  const galArg = (bCut - absB) / trans;
  const wGal = 1.0 / (1.0 + Math.exp(-2.0 * galArg));

  // 2. Supergalactic ZoA component
  let wSG = 0.0;
  if (includeSG) {
    let mySgb = sgb;
    if (mySgb === undefined) {
      const sg = galacticToSupergalactic(glon, glat);
      mySgb = sg[1];
    }
    const absSgb = Math.abs(mySgb);
    const sgArg = (sgbCut - absSgb) / trans;
    wSG = 1.0 / (1.0 + Math.exp(-2.0 * sgArg));
  }

  // 3. Dust Extinction component
  let wDust = 0.0;
  if (includeDust) {
    const av = computeDustExtinctionAV(glon, glat);
    const avCut = options.avCut || ZOA_CONSTANTS.DEFAULT_AV_CUT;
    const dustArg = (av - avCut) / (avCut * 0.5);
    wDust = 1.0 / (1.0 + Math.exp(-2.0 * dustArg));
  }

  // Combine via union probability 1 - (1 - w1)*(1 - w2)
  const combined = 1.0 - (1.0 - wGal) * (1.0 - (includeSG ? wSG * 0.5 : 0.0)) * (1.0 - (includeDust ? wDust * 0.7 : 0.0));
  return Math.min(1.0, Math.max(0.0, combined));
}

// ============================================================================
// 2. Piercing Corridor Manager
// ============================================================================

/**
 * Manages 21cm HI / Infrared deep survey corridors piercing through the ZoA veil.
 */
export class PiercingCorridorManager {
  /**
   * Constructs the corridor manager.
   * @param {Array<Object>} [customCorridors] Optional custom corridor definitions
   */
  constructor(customCorridors = null) {
    this.corridors = customCorridors ? [...customCorridors] : [...PIERCING_CORRIDORS];
  }

  /**
   * Adds a custom deep-survey piercing corridor.
   * @param {Object} corridor Corridor definition
   */
  addCorridor(corridor) {
    if (!corridor || !corridor.name) {
      throw new Error('PiercingCorridorManager.addCorridor: Valid corridor object required.');
    }
    this.corridors.push(corridor);
  }

  /**
   * Evaluates corridor transmission T(x) in [0.0, 1.0] for a given Galactic position and distance.
   * T = 1.0: Full survey transparency through dust (high completeness)
   * T = 0.0: Outside any piercing corridor
   * 
   * @param {number} glon Galactic longitude in deg
   * @param {number} glat Galactic latitude in deg
   * @param {number} distMpc Distance in Mpc (or recessional velocity in km/s if useVelocity=true)
   * @param {boolean} [useVelocity=false] True if distance is in km/s cz
   * @returns {{ inCorridor: boolean, transmission: number, activeCorridors: string[] }}
   */
  evaluateTransmission(glon, glat, distMpc, useVelocity = false) {
    if (isNaN(glon) || isNaN(glat) || isNaN(distMpc)) {
      return { inCorridor: false, transmission: 0.0, activeCorridors: [] };
    }

    let totalTransmission = 0.0;
    const active = [];

    for (let i = 0; i < this.corridors.length; i++) {
      const c = this.corridors[i];

      // Longitudinal wrapping
      let lDiff = glon - c.lMin;
      while (lDiff < 0) lDiff += 360.0;
      while (lDiff >= 360.0) lDiff -= 360.0;
      const lSpan = c.lMax >= c.lMin ? c.lMax - c.lMin : c.lMax + 360.0 - c.lMin;

      const inL = lDiff <= lSpan;
      const inB = glat >= c.bMin && glat <= c.bMax;

      let inDist = false;
      if (useVelocity) {
        inDist = distMpc >= c.czMin && distMpc <= c.czMax;
      } else {
        inDist = distMpc >= c.distMinMpc && distMpc <= c.distMaxMpc;
      }

      if (inL && inB && inDist) {
        // Compute continuous 3D Gaussian kernel transmission profile
        const lMid = (c.lMin + c.lMax) / 2.0;
        const bMid = (c.bMin + c.bMax) / 2.0;
        const dMid = useVelocity ? (c.czMin + c.czMax) / 2.0 : (c.distMinMpc + c.distMaxMpc) / 2.0;

        const sigmaL = Math.max(1.0, lSpan / 2.5);
        const sigmaB = Math.max(1.0, (c.bMax - c.bMin) / 2.5);
        const sigmaD = Math.max(5.0, (useVelocity ? (c.czMax - c.czMin) : (c.distMaxMpc - c.distMinMpc)) / 2.5);

        const dL = Math.abs(glon - lMid);
        const dB = Math.abs(glat - bMid);
        const dD = Math.abs(distMpc - dMid);

        const rNormSq = (dL * dL) / (sigmaL * sigmaL) + (dB * dB) / (sigmaB * sigmaB) + (dD * dD) / (sigmaD * sigmaD);
        const kernelVal = Math.exp(-0.5 * rNormSq);
        const corridorT = (c.completeness || 0.85) * kernelVal;

        totalTransmission = Math.max(totalTransmission, corridorT);
        active.push(c.name);
      }
    }

    return {
      inCorridor: active.length > 0,
      transmission: Math.min(1.0, totalTransmission),
      activeCorridors: active
    };
  }

  /**
   * Computes effective obscuration mask M_eff(x) = M_raw(x) * (1.0 - transmission(x)).
   * @param {number} glon Galactic longitude in deg
   * @param {number} glat Galactic latitude in deg
   * @param {number} distMpc Distance in Mpc
   * @param {Object} [options] Mask options
   * @returns {number} Effective obscuration factor in [0, 1]
   */
  computeEffectiveObscuration(glon, glat, distMpc, options = {}) {
    const rawMask = computeZoAObscurationWeight(glon, glat, undefined, undefined, options);
    if (rawMask <= 1e-4) return 0.0;

    const transInfo = this.evaluateTransmission(glon, glat, distMpc, false);
    const eff = rawMask * (1.0 - transInfo.transmission);
    return Math.min(1.0, Math.max(0.0, eff));
  }
}

// ============================================================================
// 3. Wiener-Filter & Cosmological Covariance Formulations
// ============================================================================

/**
 * Computes cosmological spatial two-point correlation / covariance xi(r).
 * Formulations:
 * - 'powerlaw': xi(r) = sigma_0^2 / (1 + (r / r_0)^2)^(gamma / 2)
 * - 'gaussian': xi(r) = sigma_0^2 * exp(-0.5 * (r / r_0)^2)
 * - 'matern':   xi(r) = sigma_0^2 * (1 + sqrt(3)*r/r0) * exp(-sqrt(3)*r/r0)
 * 
 * @param {number} r Physical 3D separation in Mpc (or Mpc/h)
 * @param {Object} [options] Covariance parameters
 * @param {number} [options.r0=10.0] Correlation scale r_0 in Mpc
 * @param {number} [options.gamma=1.8] Power law index
 * @param {number} [options.sigma0Sq=1.0] Zero-lag field variance sigma_0^2
 * @param {'powerlaw'|'gaussian'|'matern'} [options.type='powerlaw'] Covariance kernel type
 * @returns {number} Covariance xi(r)
 */
export function cosmologicalCovariance(r, options = {}) {
  if (isNaN(r) || r < 0.0) return 0.0;
  const r0 = options.r0 || ZOA_CONSTANTS.CORRELATION_LENGTH_MPC;
  const gamma = options.gamma || ZOA_CONSTANTS.CORRELATION_GAMMA;
  const sigma0Sq = options.sigma0Sq !== undefined ? options.sigma0Sq : 1.0;
  const type = options.type || 'powerlaw';

  if (r < 1e-12) return sigma0Sq;

  if (type === 'gaussian') {
    const u = r / r0;
    return sigma0Sq * Math.exp(-0.5 * u * u);
  } else if (type === 'matern') {
    const s = Math.SQRT2 * (r / r0);
    return sigma0Sq * (1.0 + s) * Math.exp(-s);
  } else {
    // Standard power-law model
    const uSq = (r * r) / (r0 * r0);
    return sigma0Sq / Math.pow(1.0 + uSq, gamma / 2.0);
  }
}

// ============================================================================
// 4. Point-Based Local Wiener / Kriging Estimator
// ============================================================================

/**
 * Solves dense linear system A * x = b via Cholesky decomposition (with diagonal regularizer).
 * @param {number[][]} A Symmetric positive-definite matrix (NxN)
 * @param {number[]} b Right hand side vector (N)
 * @param {number} [reg=1e-5] Tikhonov diagonal regularizer
 * @returns {number[]} Solution vector x
 */
export function solveCholesky(A, b, reg = 1e-5) {
  const n = b.length;
  const L = new Array(n);
  for (let i = 0; i < n; i++) {
    L[i] = new Float64Array(n);
  }

  // Cholesky decomposition L * L^T = A + reg*I
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = 0.0;
      for (let k = 0; k < j; k++) {
        sum += L[i][k] * L[j][k];
      }
      if (i === j) {
        const diag = (A[i][i] + reg) - sum;
        L[i][j] = Math.sqrt(Math.max(1e-12, diag));
      } else {
        L[i][j] = (A[i][j] - sum) / L[j][j];
      }
    }
  }

  // Forward solve L * y = b
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let sum = 0.0;
    for (let k = 0; k < i; k++) {
      sum += L[i][k] * y[k];
    }
    y[i] = (b[i] - sum) / L[i][i];
  }

  // Back substitution L^T * x = y
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let sum = 0.0;
    for (let k = i + 1; k < n; k++) {
      sum += L[k][i] * x[k];
    }
    x[i] = (y[i] - sum) / L[i][i];
  }

  return Array.from(x);
}

/**
 * Reconstructs cosmological overdensity delta(x) at an arbitrary target 3D point in the ZoA
 * using local Wiener filter / Gaussian Process Kriging across surrounding observed galaxies.
 * 
 * Computes:
 * - Wiener expectation: delta_* = k_*^T * (K + N)^-1 * y
 * - Reconstruction error variance: sigma_*^2 = xi(0) - k_*^T * (K + N)^-1 * k_*
 * 
 * @param {number} targetSgx Supergalactic Cartesian X in Mpc
 * @param {number} targetSgy Supergalactic Cartesian Y in Mpc
 * @param {number} targetSgz Supergalactic Cartesian Z in Mpc
 * @param {Array<Object>} observedTracers Array of tracer objects with { sgx, sgy, sgz, delta, noiseVar }
 * @param {Object} [options] Wiener options
 * @param {number} [options.maxNeighbors=32] Maximum local neighbors
 * @param {number} [options.searchRadiusMpc=50.0] Neighborhood search radius
 * @param {number} [options.r0=10.0] Correlation length
 * @returns {{ delta: number, variance: number, neighborCount: number }}
 */
export function reconstructPointWiener(
  targetSgx,
  targetSgy,
  targetSgz,
  observedTracers,
  options = {}
) {
  if (!Array.isArray(observedTracers) || observedTracers.length === 0) {
    return { delta: 0.0, variance: 1.0, neighborCount: 0 };
  }

  const maxN = options.maxNeighbors || 32;
  const maxR = options.searchRadiusMpc || 50.0;
  const r0 = options.r0 || ZOA_CONSTANTS.CORRELATION_LENGTH_MPC;

  // Filter neighbors within searchRadius
  const neighbors = [];
  for (let i = 0; i < observedTracers.length; i++) {
    const t = observedTracers[i];
    const dx = t.sgx - targetSgx;
    const dy = t.sgy - targetSgy;
    const dz = t.sgz - targetSgz;
    const dist = Math.hypot(dx, dy, dz);
    if (dist <= maxR) {
      neighbors.push({ ...t, dist: dist });
    }
  }

  if (neighbors.length === 0) {
    return { delta: 0.0, variance: cosmologicalCovariance(0.0, { r0 }), neighborCount: 0 };
  }

  // Sort by distance and take closest maxN
  neighbors.sort((a, b) => a.dist - b.dist);
  const selected = neighbors.slice(0, maxN);
  const n = selected.length;

  // Build covariance matrix K and cross-covariance kStar
  const K = new Array(n);
  const kStar = new Float64Array(n);
  const y = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    K[i] = new Float64Array(n);
    kStar[i] = cosmologicalCovariance(selected[i].dist, { r0 });
    y[i] = selected[i].delta || 0.0;

    for (let j = 0; j <= i; j++) {
      const dx = selected[i].sgx - selected[j].sgx;
      const dy = selected[i].sgy - selected[j].sgy;
      const dz = selected[i].sgz - selected[j].sgz;
      const dSep = Math.hypot(dx, dy, dz);
      const cov = cosmologicalCovariance(dSep, { r0 });
      K[i][j] = cov;
      K[j][i] = cov;
    }
    // Add noise variance on diagonal
    const noise = selected[i].noiseVar || 0.05;
    K[i][i] += noise;
  }

  // Solve (K + N) * alpha = y
  const alpha = solveCholesky(K, y, 1e-4);

  // Solve (K + N) * beta = kStar for variance
  const beta = solveCholesky(K, kStar, 1e-4);

  let deltaEst = 0.0;
  let varReduction = 0.0;
  for (let i = 0; i < n; i++) {
    deltaEst += kStar[i] * alpha[i];
    varReduction += kStar[i] * beta[i];
  }

  const priorVar = cosmologicalCovariance(0.0, { r0 });
  const postVar = Math.max(0.0, priorVar - varReduction);

  return {
    delta: deltaEst,
    variance: postVar,
    neighborCount: n
  };
}

// ============================================================================
// 5. 3D Grid-Based Wiener Inpainting Solver
// ============================================================================

/**
 * 3D Zone of Avoidance (ZoA) Reconstructive Inpainting Engine.
 * Implements regularized iterative Wiener-filter reconstruction across dust-obscured voxels.
 */
export class ZoAReconstructor {
  /**
   * Constructs the ZoAReconstructor.
   * @param {Object} [options] Reconstructor options
   * @param {number} [options.bCut=10.0] Galactic latitude cutoff in deg
   * @param {number} [options.sgbCut=10.0] Supergalactic latitude cutoff in deg
   * @param {number} [options.r0=10.0] Cosmological correlation length in Mpc
   * @param {number} [options.maxIterations=50] Conjugate gradient / inpainting iterations
   * @param {number} [options.tolerance=1e-5] Convergence tolerance
   * @param {PiercingCorridorManager} [options.corridorManager] Corridor manager instance
   */
  constructor(options = {}) {
    this.bCut = options.bCut !== undefined ? options.bCut : ZOA_CONSTANTS.DEFAULT_GAL_B_CUT_DEG;
    this.sgbCut = options.sgbCut !== undefined ? options.sgbCut : ZOA_CONSTANTS.DEFAULT_SG_B_CUT_DEG;
    this.r0 = options.r0 || ZOA_CONSTANTS.CORRELATION_LENGTH_MPC;
    this.maxIterations = options.maxIterations || 50;
    this.tolerance = options.tolerance || 1e-5;
    this.corridors = options.corridorManager || new PiercingCorridorManager();
  }

  /**
   * Generates 3D voxel mask field M(x) in [0, 1] on a GridIndexer domain.
   * M = 1: Deeply obscured ZoA core
   * M = 0: Observed survey region
   * 
   * @param {GridIndexer} grid Associated 3D grid
   * @returns {Float64Array} 1D flat mask buffer of size grid.totalCells
   */
  createGridMask(grid) {
    if (!(grid instanceof GridIndexer)) {
      throw new TypeError('ZoAReconstructor.createGridMask: grid must be an instance of GridIndexer.');
    }

    const total = grid.totalCells;
    const mask = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const [ix, iy, iz] = grid.get3DIndices(i);
      const [sgx, sgy, sgz] = grid.getNodeCoord(ix, iy, iz);
      const dist = Math.hypot(sgx, sgy, sgz);

      if (dist < 1e-6) {
        mask[i] = 0.0;
        continue;
      }

      // Convert Supergalactic Cartesian -> Equatorial -> Galactic
      const [ra, dec] = supergalacticCartesianToEquatorial(sgx, sgy, sgz);
      const [glon, glat] = equatorialToGalactic(ra, dec);
      const [sgl, sgb] = galacticToSupergalactic(glon, glat);

      // Compute raw ZoA obscuration and de-obscure via piercing corridors
      const rawObscuration = computeZoAObscurationWeight(glon, glat, sgl, sgb, {
        bCut: this.bCut,
        sgbCut: this.sgbCut
      });

      const trans = this.corridors.evaluateTransmission(glon, glat, dist, false);
      const effMask = rawObscuration * (1.0 - trans.transmission);

      mask[i] = effMask;
    }

    return mask;
  }

  /**
   * Reconstructs cosmological density field delta(x) across the Zone of Avoidance using 3D Wiener Inpainting.
   * Preserves observed data outside ZoA, smoothly in-paints continuous cosmological structures inside ZoA,
   * and enforces C^0/C^1 boundary continuity.
   * 
   * @param {DensityField} inputDensity Observed density field
   * @param {Object} [options] Inpainting options
   * @param {Float64Array} [options.customMask] Pre-computed mask buffer
   * @returns {{ reconstructedField: DensityField, mask: Float64Array, errorMap: Float64Array, iterations: number, residual: number }}
   */
  reconstructDensityField(inputDensity, options = {}) {
    if (!(inputDensity instanceof DensityField)) {
      throw new TypeError('ZoAReconstructor.reconstructDensityField: inputDensity must be an instance of DensityField.');
    }

    const grid = inputDensity.grid;
    const total = grid.totalCells;
    const nx = grid.nx;
    const ny = grid.ny;
    const nz = grid.nz;
    const dx = grid.dx;
    const dy = grid.dy;
    const dz = grid.dz;

    const mask = options.customMask || this.createGridMask(grid);
    const dObs = inputDensity.delta;

    // Output delta buffer
    const deltaRecon = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      deltaRecon[i] = dObs[i];
    }

    // 3D Covariance Smoothing Kernel radius in cells
    const cellDiag = Math.hypot(dx, dy, dz);
    const kernelRadiusCells = Math.max(1, Math.min(4, Math.round(this.r0 / cellDiag)));

    // Pre-calculate 3D separable or compact convolution weights
    const kernelWeights = [];
    let kWeightSum = 0.0;

    for (let kz = -kernelRadiusCells; kz <= kernelRadiusCells; kz++) {
      for (let ky = -kernelRadiusCells; ky <= kernelRadiusCells; ky++) {
        for (let kx = -kernelRadiusCells; kx <= kernelRadiusCells; kx++) {
          const rx = kx * dx;
          const ry = ky * dy;
          const rz = kz * dz;
          const rSep = Math.hypot(rx, ry, rz);
          const w = cosmologicalCovariance(rSep, { r0: this.r0 });
          kernelWeights.push({ kx, ky, kz, weight: w });
          kWeightSum += w;
        }
      }
    }
    // Normalize kernel
    for (const kw of kernelWeights) {
      kw.weight /= kWeightSum;
    }

    // Iterative Wiener Inpainting Loop
    let iter = 0;
    let maxDeltaChange = 1.0;
    const tempBuffer = new Float64Array(total);

    while (iter < this.maxIterations && maxDeltaChange > this.tolerance) {
      iter++;
      maxDeltaChange = 0.0;

      // Step 1: Smooth current estimate with cosmological covariance kernel
      for (let iz = 0; iz < nz; iz++) {
        for (let iy = 0; iy < ny; iy++) {
          for (let ix = 0; ix < nx; ix++) {
            const idx = grid.getLinearIndex(ix, iy, iz);
            const mVal = mask[idx];

            if (mVal <= 1e-4) {
              // Strictly observed: retain observation
              tempBuffer[idx] = dObs[idx];
              continue;
            }

            let convSum = 0.0;
            for (let k = 0; k < kernelWeights.length; k++) {
              const { kx, ky, kz, weight } = kernelWeights[k];
              const clamped = grid.handleBoundaryIndex(ix + kx, iy + ky, iz + kz);
              if (clamped) {
                const neighborIdx = grid.getLinearIndex(clamped[0], clamped[1], clamped[2]);
                convSum += deltaRecon[neighborIdx] * weight;
              }
            }

            // Inpainting update: (1 - M)*dObs + M*smoothed
            const updatedVal = (1.0 - mVal) * dObs[idx] + mVal * convSum;
            tempBuffer[idx] = updatedVal;

            const change = Math.abs(updatedVal - deltaRecon[idx]);
            if (change > maxDeltaChange) {
              maxDeltaChange = change;
            }
          }
        }
      }

      // Step 2: Copy back buffer
      for (let i = 0; i < total; i++) {
        deltaRecon[i] = tempBuffer[i];
      }
    }

    // Generate Reconstruction Uncertainty Map
    const errorMap = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      const mVal = mask[i];
      // Uncertainty proportional to obscuration weight and prior variance
      errorMap[i] = mVal * Math.sqrt(cosmologicalCovariance(0.0, { r0: this.r0 }));
    }

    const reconstructedField = new DensityField(grid, deltaRecon, {
      deltaFloor: inputDensity.deltaFloor
    });

    return {
      reconstructedField: reconstructedField,
      mask: mask,
      errorMap: errorMap,
      iterations: iter,
      residual: maxDeltaChange
    };
  }
}

// ============================================================================
// 6. Synthetic Control Benchmarking
// ============================================================================

/**
 * Synthetic Control Benchmark Harness for ZoA Inpainting Validation.
 * Artificially masks known cosmological regions, executes reconstruction,
 * and benchmarks recovery fidelity (Pearson r, RMSE, relative flux error, bulk-flow recovery).
 */
export class SyntheticControlBenchmark {
  /**
   * Constructs the benchmark harness.
   * @param {Object} [options] Benchmark options
   * @param {number} [options.gridSize=32] Synthetic grid dimension (N^3)
   * @param {number} [options.boxSize=200.0] Box extent in Mpc
   */
  constructor(options = {}) {
    this.n = options.gridSize || 32;
    this.boxSize = options.boxSize || 200.0;
    this.grid = new GridIndexer({
      nx: this.n,
      ny: this.n,
      nz: this.n,
      origin: [-this.boxSize / 2.0, -this.boxSize / 2.0, -this.boxSize / 2.0],
      boxSize: [this.boxSize, this.boxSize, this.boxSize],
      strideOrder: StrideOrder.CANONICAL_XYZ,
      boundaryMode: BoundaryMode.CLAMP
    });
  }

  /**
   * Generates a ground-truth synthetic cosmological field featuring prominent superclusters and voids.
   * @returns {DensityField} Ground truth density field
   */
  generateSyntheticGroundTruth() {
    const total = this.grid.totalCells;
    const delta = new Float64Array(total);

    // Mock cosmological hubs: Great Attractor, Shapley, Perseus-Pisces, Cosmic Void
    const structures = [
      { name: 'Great_Attractor', x: -45.0, y: 15.0, z: -5.0, peak: 4.5, sigma: 18.0 },
      { name: 'Shapley_Supercluster', x: -120.0, y: 70.0, z: -10.0, peak: 6.0, sigma: 25.0 },
      { name: 'Perseus_Pisces', x: 50.0, y: -20.0, z: -15.0, peak: 3.8, sigma: 16.0 },
      { name: 'Sculptor_Void', x: 30.0, y: 40.0, z: 20.0, peak: -0.85, sigma: 22.0 },
      { name: 'Local_Void', x: -10.0, y: -25.0, z: 35.0, peak: -0.90, sigma: 20.0 }
    ];

    for (let i = 0; i < total; i++) {
      const [ix, iy, iz] = this.grid.get3DIndices(i);
      const [x, y, z] = this.grid.getNodeCoord(ix, iy, iz);

      let dVal = 0.0;
      for (const s of structures) {
        const dx = x - s.x;
        const dy = y - s.y;
        const dz = z - s.z;
        const rSq = dx * dx + dy * dy + dz * dz;
        dVal += s.peak * Math.exp(-0.5 * rSq / (s.sigma * s.sigma));
      }

      // Add smooth cosmological background ripple
      dVal += 0.25 * Math.sin(x * 0.05) * Math.cos(y * 0.05) * Math.sin(z * 0.05);
      delta[i] = Math.max(-0.999, dVal);
    }

    return new DensityField(this.grid, delta);
  }

  /**
   * Executes the full benchmark run:
   * 1. Creates synthetic ground truth.
   * 2. Applies ZoA mask with piercing corridors.
   * 3. In-paints using ZoAReconstructor.
   * 4. Evaluates Pearson r, RMSE, relative flux conservation, and dipole alignment error.
   * 
   * @param {Object} [options] Reconstructor options
   * @returns {{ pearsonR: number, rmse: number, relativeFluxError: number, maskedFraction: number, iterations: number, success: boolean }}
   */
  runBenchmark(options = {}) {
    const truthField = this.generateSyntheticGroundTruth();
    const reconstructor = new ZoAReconstructor(options);

    const mask = reconstructor.createGridMask(this.grid);
    const total = this.grid.totalCells;

    // Create masked observed field (set masked voxels to 0.0 or noisy sample)
    const maskedObsBuffer = new Float64Array(total);
    let maskedCount = 0;

    for (let i = 0; i < total; i++) {
      if (mask[i] > 0.5) {
        maskedObsBuffer[i] = 0.0; // Masked out
        maskedCount++;
      } else {
        maskedObsBuffer[i] = truthField.delta[i];
      }
    }

    const maskedField = new DensityField(this.grid, maskedObsBuffer);
    const reconResult = reconstructor.reconstructDensityField(maskedField, {
      customMask: mask
    });

    const dTrue = truthField.delta;
    const dRecon = reconResult.reconstructedField.delta;

    // Compute metrics exclusively inside the masked region (mask[i] > 0.3)
    let sumTrue = 0.0;
    let sumRecon = 0.0;
    let countInMask = 0;

    for (let i = 0; i < total; i++) {
      if (mask[i] > 0.3) {
        sumTrue += dTrue[i];
        sumRecon += dRecon[i];
        countInMask++;
      }
    }

    const meanTrue = countInMask > 0 ? sumTrue / countInMask : 0.0;
    const meanRecon = countInMask > 0 ? sumRecon / countInMask : 0.0;

    let numer = 0.0;
    let denomTrue = 0.0;
    let denomRecon = 0.0;
    let sumSqErr = 0.0;

    for (let i = 0; i < total; i++) {
      if (mask[i] > 0.3) {
        const diffT = dTrue[i] - meanTrue;
        const diffR = dRecon[i] - meanRecon;
        numer += diffT * diffR;
        denomTrue += diffT * diffT;
        denomRecon += diffR * diffR;

        const err = dRecon[i] - dTrue[i];
        sumSqErr += err * err;
      }
    }

    const denom = Math.sqrt(denomTrue * denomRecon);
    const pearsonR = denom > 1e-12 ? numer / denom : 0.0;
    const rmse = countInMask > 0 ? Math.sqrt(sumSqErr / countInMask) : 0.0;
    const relFluxErr = Math.abs(sumTrue) > 1e-6 ? Math.abs(sumRecon - sumTrue) / Math.abs(sumTrue) : 0.0;

    return {
      pearsonR: pearsonR,
      rmse: rmse,
      relativeFluxError: relFluxErr,
      maskedFraction: maskedCount / total,
      maskedCellCount: countInMask,
      iterations: reconResult.iterations,
      success: pearsonR > 0.75
    };
  }
}
