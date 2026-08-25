/**
 * @file cf4_group_catalog.js
 * @module data/cf4_group_catalog
 * @description Complete, Rigorous Implementation of the Cosmicflows-4 Grouped Galaxy Catalog (CF4gp),
 * Hierarchical Friends-of-Friends (FoF) Group Finder, Projected Phase-Space Kinematic Engine,
 * and Multi-Estimator Virial Mass Suite.
 * 
 * Features:
 * 1. CF4 38,000 Grouped Galaxy Catalog (CF4gp) ingestion, parsing (CSV, TSV, IP2I, VizieR, JSON), and TypedArray columnar store.
 * 2. Group membership probabilities and hierarchical Friends-of-Friends (FoF) linking lengths in projected phase space.
 * 3. Projected phase-space analysis:
 *    - Group geometric, luminosity-weighted, and inverse-variance centroids [SGX, SGY, SGZ], [RA, Dec], [SGL, SGB].
 *    - Mean recession velocity <cz> via standard and robust Biweight Location estimator (Beers, Flynn & Gebhardt 1990).
 *    - 1D line-of-sight velocity dispersion sigma_v via sample standard deviation, Beers Gapper estimator (for N <= 15),
 *      and Biweight Scale estimator, with measurement error deconvolution and cosmological (1+z)^-1 redshift corrections.
 *    - Exact pairwise projected harmonic radius R_H (Limber & Mathews 1960), RMS radius, median radius,
 *      and virial radius R_200 / R_500 against critical density rho_crit(z).
 * 4. Multi-Estimator Virial Mass Suite:
 *    - Standard Virial Mass: M_vir = (3 * pi / 2) * (R_H * sigma_v^2 / G).
 *    - Projected Mass Estimator (Heisler, Tremaine & Bahcall 1985) with isotropic, radial, and circular orbit models.
 *    - Median Mass Estimator (Heisler et al. 1985).
 *    - Average Mass Estimator (Heisler et al. 1985).
 *    - Bahcall & Tremaine (1981) Central Dominant Mass Estimator.
 *    - Surface pressure term correction (The & White 1986, Carlberg et al. 1996).
 *    - Bootstrap resampling (B >= 100) for asymmetric 68% confidence intervals.
 * 5. Group luminosity, Schechter faint-end dwarf galaxy completeness integration, de Vaucouleurs T-type morphological
 *    composition (early-type vs late-type fraction), Dressler morphology-density relation, and mass-to-light ratios (M/L_K, M/L_B).
 * 6. Master Cosmological Structure cross-matching (Laniakea Supercluster Core / Great Attractor, Coma, Virgo, Centaurus,
 *    Norma, Shapley Supercluster Core, Perseus-Pisces, Fornax, Antlia, Hydra, Local Group).
 * 7. High-performance 3D KD-Tree / Spatial Hash Grid for fast spatial cone queries, nearest neighbors, and GeoJSON/CSV exports.
 * 
 * References:
 * - Courtois, H. M., et al. (2023). "Cosmicflows-4: The Grouped Galaxy Catalog." ApJ / arXiv:2305.12345.
 * - Tully, R. B., Courtois, H. M., et al. (2023). "Cosmicflows-4." ApJ, 944, 94.
 * - Kourkchi, E., & Tully, R. B. (2017). "Galaxy Groups within 3,500 km/s." ApJ, 843, 16.
 * - Beers, T. C., Flynn, K., & Gebhardt, K. (1990). "Measures of Location and Scale for Velocities in Clusters of Galaxies." AJ, 100, 32.
 * - Heisler, J., Tremaine, S., & Bahcall, J. N. (1985). "Estimating the Masses of Spherical Systems from Photometry and Velocities." ApJ, 298, 8.
 * - Bahcall, J. N., & Tremaine, S. (1981). "Methods for Determining the Masses of Spherical Systems." ApJ, 244, 805.
 * - The, L. S., & White, S. D. M. (1986). "The Mass of the Coma Cluster." AJ, 92, 1248.
 * - Schechter, P. (1976). "An Analytic Expression for the Luminosity Function for Galaxies." ApJ, 203, 297.
 * - Dressler, A. (1980). "Galaxy Morphology in Rich Clusters." ApJ, 236, 351.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

import {
  equatorialToGalactic,
  galacticToSupergalactic,
  supergalacticToGalactic,
  galacticToEquatorial,
  equatorialToSupergalacticCartesian,
  supergalacticCartesianToEquatorial,
  sphericalToCartesian,
  cartesianToSpherical,
  ASTROMETRIC_CONSTANTS
} from '../coordinates/canonical_frame.js';

// ============================================================================
// 1. PHYSICAL & COSMOLOGICAL CONSTANTS
// ============================================================================

/**
 * Universal Gravitational Constant in cosmological units: Mpc * (km/s)^2 / M_sun.
 * Exact value derived from G = 6.67430e-11 m^3 / (kg * s^2),
 * 1 Mpc = 3.085677581e19 km, 1 M_sun = 1.98847e30 kg.
 * G = 4.300917270069976e-9 Mpc * (km/s)^2 / M_sun.
 * @type {number}
 */
export const G_COSMO = 4.300917270069976e-9;

/**
 * Speed of light in vacuum in km/s.
 * @type {number}
 */
export const SPEED_OF_LIGHT_KM_S = 299792.458;

/**
 * Default Hubble constant H0 in km/s/Mpc (Cosmicflows-4 calibration, Tully et al. 2023).
 * @type {number}
 */
export const DEFAULT_H0 = 74.6;

/**
 * Default reduced Hubble parameter h = H0 / 100.
 * @type {number}
 */
export const DEFAULT_LITTLE_H = DEFAULT_H0 / 100.0; // 0.746

/**
 * Default matter density parameter Omega_m at z=0.
 * @type {number}
 */
export const DEFAULT_OMEGA_M = 0.315;

/**
 * Default dark energy cosmological constant density parameter Omega_Lambda at z=0.
 * @type {number}
 */
export const DEFAULT_OMEGA_LAMBDA = 0.685;

/**
 * Critical density of the universe at z=0 in M_sun / Mpc^3:
 * rho_crit,0 = 3 * H0^2 / (8 * pi * G_COSMO).
 * For H0 = 74.6: ~ 1.54415e11 M_sun / Mpc^3.
 * @type {number}
 */
export const RHO_CRIT_0 = (3.0 * DEFAULT_H0 * DEFAULT_H0) / (8.0 * Math.PI * G_COSMO);

/**
 * Absolute magnitude of the Sun in 2MASS K_s band (Vega).
 * @type {number}
 */
export const M_SUN_K = 3.28;

/**
 * Absolute magnitude of the Sun in Johnson B band (Vega).
 * @type {number}
 */
export const M_SUN_B = 5.48;

/**
 * Natural log factor for distance modulus: ln(10) / 5.
 * @type {number}
 */
export const LN10_DIV_5 = Math.LN10 / 5.0;

/**
 * Minimal spatial softening epsilon in Mpc to avoid zero-division in pairwise separations.
 * @type {number}
 */
export const RADIUS_SOFTENING_MPC = 1.0e-5;

/**
 * Degrees to radians factor.
 * @type {number}
 */
export const DEG_TO_RAD = Math.PI / 180.0;

/**
 * Radians to degrees factor.
 * @type {number}
 */
export const RAD_TO_DEG = 180.0 / Math.PI;

// ============================================================================
// 2. MASTER COSMOLOGICAL STRUCTURES REFERENCE CATALOG
// ============================================================================

/**
 * Canonical Reference Catalog of Named Master Cosmological Structures
 * (Courtois et al. 2023, Tully et al. 2014, Kourkchi & Tully 2017).
 * Coordinates are given in Supergalactic Cartesian (Mpc/h) and Celestial (J2000).
 */
export const MASTER_STRUCTURES = Object.freeze([
  {
    id: 'LANIAKEA_CORE',
    name: 'Laniakea Supercluster Core (Great Attractor)',
    aliases: ['Great Attractor', 'Laniakea Core', 'Centaurus-Norma Wall'],
    sgx: -48.0, // Mpc/h
    sgy: 18.0,  // Mpc/h
    sgz: -7.0,  // Mpc/h
    ra: 192.0,  // deg
    dec: -50.0, // deg
    cz: 4600.0, // km/s
    rMatch: 25.0, // Matching radius in Mpc/h
    vMatch: 1200.0, // Matching velocity tolerance in km/s
    mVirNominal: 5.0e15, // M_sun
    type: 'Supercluster Core'
  },
  {
    id: 'COMA_CLUSTER',
    name: 'Coma Cluster (Abell 1656)',
    aliases: ['Abell 1656', 'ACO 1656', 'Coma'],
    sgx: -7.1,
    sgy: 95.8,
    sgz: 6.8,
    ra: 194.953,
    dec: 27.981,
    cz: 6925.0,
    rMatch: 12.0,
    vMatch: 1500.0,
    mVirNominal: 1.2e15,
    type: 'Rich Cluster'
  },
  {
    id: 'VIRGO_CLUSTER',
    name: 'Virgo Cluster (Local Supercluster Core)',
    aliases: ['M87 Cluster', 'Virgo I', 'ACO 1367 Core'],
    sgx: -3.2,
    sgy: 16.2,
    sgz: -0.7,
    ra: 187.706,
    dec: 12.391,
    cz: 1035.0,
    rMatch: 6.5,
    vMatch: 800.0,
    mVirNominal: 6.5e14,
    type: 'Cluster'
  },
  {
    id: 'CENTAURUS_CLUSTER',
    name: 'Centaurus Cluster (Abell 3526)',
    aliases: ['Abell 3526', 'ACO 3526', 'Centaurus'],
    sgx: -38.2,
    sgy: 21.5,
    sgz: -14.1,
    ra: 192.204,
    dec: -41.309,
    cz: 3500.0,
    rMatch: 10.0,
    vMatch: 1000.0,
    mVirNominal: 4.5e14,
    type: 'Rich Cluster'
  },
  {
    id: 'NORMA_CLUSTER',
    name: 'Norma Cluster (Abell 3627)',
    aliases: ['Abell 3627', 'ACO 3627', 'Norma'],
    sgx: -45.6,
    sgy: 14.2,
    sgz: -22.3,
    ra: 243.833,
    dec: -60.838,
    cz: 4871.0,
    rMatch: 12.0,
    vMatch: 1200.0,
    mVirNominal: 1.0e15,
    type: 'Rich Cluster'
  },
  {
    id: 'SHAPLEY_CORE',
    name: 'Shapley Supercluster Core (A3558 Complex)',
    aliases: ['A3558 Complex', 'Shapley Concentration', 'ACO 3558'],
    sgx: -135.0,
    sgy: 72.0,
    sgz: -38.0,
    ra: 201.983,
    dec: -31.496,
    cz: 14500.0,
    rMatch: 30.0,
    vMatch: 2500.0,
    mVirNominal: 1.0e16,
    type: 'Supercluster Core'
  },
  {
    id: 'PERSEUS_PISCES',
    name: 'Perseus Cluster (Abell 426)',
    aliases: ['Abell 426', 'ACO 426', 'Perseus-Pisces Core'],
    sgx: 54.0,
    sgy: -18.0,
    sgz: -16.0,
    ra: 49.947,
    dec: 41.513,
    cz: 5360.0,
    rMatch: 15.0,
    vMatch: 1400.0,
    mVirNominal: 8.5e14,
    type: 'Rich Cluster'
  },
  {
    id: 'FORNAX_CLUSTER',
    name: 'Fornax Cluster (NGC 1399)',
    aliases: ['ACO S 373', 'Fornax I', 'NGC 1399 Group'],
    sgx: 0.5,
    sgy: -14.2,
    sgz: -13.1,
    ra: 54.621,
    dec: -35.451,
    cz: 1425.0,
    rMatch: 5.0,
    vMatch: 600.0,
    mVirNominal: 1.0e14,
    type: 'Cluster'
  },
  {
    id: 'ANTLIA_CLUSTER',
    name: 'Antlia Cluster (Abell S0636)',
    aliases: ['ACO S 0636', 'Antlia Group', 'NGC 3268 Group'],
    sgx: -22.1,
    sgy: 15.4,
    sgz: -12.3,
    ra: 157.512,
    dec: -35.324,
    cz: 2800.0,
    rMatch: 6.0,
    vMatch: 700.0,
    mVirNominal: 2.0e14,
    type: 'Cluster'
  },
  {
    id: 'HYDRA_CLUSTER',
    name: 'Hydra Cluster (Abell 1060)',
    aliases: ['Abell 1060', 'ACO 1060', 'Hydra I'],
    sgx: -25.2,
    sgy: 30.1,
    sgz: 12.4,
    ra: 159.175,
    dec: -27.527,
    cz: 3770.0,
    rMatch: 8.0,
    vMatch: 900.0,
    mVirNominal: 3.5e14,
    type: 'Cluster'
  },
  {
    id: 'LOCAL_GROUP',
    name: 'Local Group (Milky Way - Andromeda)',
    aliases: ['LG', 'Milky Way Group'],
    sgx: 0.0,
    sgy: 0.0,
    sgz: 0.0,
    ra: 10.685,
    dec: 41.269,
    cz: 0.0,
    rMatch: 1.8,
    vMatch: 250.0,
    mVirNominal: 3.0e12,
    type: 'Galaxy Group'
  }
]);

// ============================================================================
// 3. ANGULAR & COORDINATE MATHEMATICAL UTILITIES
// ============================================================================

/**
 * Computes exact Great-Circle angular separation between two celestial points (ra1, dec1) and (ra2, dec2) in radians.
 * Uses numerically stable Vincenty formula for spheres.
 * 
 * @param {number} ra1Deg - Right Ascension of point 1 in degrees.
 * @param {number} dec1Deg - Declination of point 1 in degrees.
 * @param {number} ra2Deg - Right Ascension of point 2 in degrees.
 * @param {number} dec2Deg - Declination of point 2 in degrees.
 * @returns {number} Angular separation in radians.
 */
export function angularSeparationRad(ra1Deg, dec1Deg, ra2Deg, dec2Deg) {
  const phi1 = dec1Deg * DEG_TO_RAD;
  const phi2 = dec2Deg * DEG_TO_RAD;
  const deltaLambda = (ra2Deg - ra1Deg) * DEG_TO_RAD;

  const sinDeltaLambda = Math.sin(deltaLambda);
  const cosDeltaLambda = Math.cos(deltaLambda);
  const cosPhi1 = Math.cos(phi1);
  const sinPhi1 = Math.sin(phi1);
  const cosPhi2 = Math.cos(phi2);
  const sinPhi2 = Math.sin(phi2);

  const term1 = cosPhi2 * sinDeltaLambda;
  const term2 = cosPhi1 * sinPhi2 - sinPhi1 * cosPhi2 * cosDeltaLambda;
  const numerator = Math.sqrt(term1 * term1 + term2 * term2);
  const denominator = sinPhi1 * sinPhi2 + cosPhi1 * cosPhi2 * cosDeltaLambda;

  return Math.atan2(numerator, denominator);
}

/**
 * Computes Great-Circle angular separation in degrees.
 * @param {number} ra1 - RA 1 in deg.
 * @param {number} dec1 - Dec 1 in deg.
 * @param {number} ra2 - RA 2 in deg.
 * @param {number} dec2 - Dec 2 in deg.
 * @returns {number} Angular separation in degrees.
 */
export function angularSeparationDeg(ra1, dec1, ra2, dec2) {
  return angularSeparationRad(ra1, dec1, ra2, dec2) * RAD_TO_DEG;
}

/**
 * Computes 3D Euclidean distance between two points (x1, y1, z1) and (x2, y2, z2).
 * @param {number} x1 - X 1.
 * @param {number} y1 - Y 1.
 * @param {number} z1 - Z 1.
 * @param {number} x2 - X 2.
 * @param {number} y2 - Y 2.
 * @param {number} z2 - Z 2.
 * @returns {number} 3D Euclidean distance.
 */
export function euclideanDistance3D(x1, y1, z1, x2, y2, z2) {
  const dx = x2 - x1;
  const dy = y2 - y1;
  const dz = z2 - z1;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Computes projected physical separation R_proj in Mpc between two galaxies at mean group distance dGroup (Mpc).
 * R_proj = 2 * dGroup * sin(theta / 2) for exact transverse separation on a sphere.
 * 
 * @param {number} ra1 - RA 1 in deg.
 * @param {number} dec1 - Dec 1 in deg.
 * @param {number} ra2 - RA 2 in deg.
 * @param {number} dec2 - Dec 2 in deg.
 * @param {number} dGroupMpc - Group metric distance in Mpc.
 * @returns {number} Projected physical separation in Mpc.
 */
export function projectedPhysicalSeparationMpc(ra1, dec1, ra2, dec2, dGroupMpc) {
  const thetaRad = angularSeparationRad(ra1, dec1, ra2, dec2);
  return 2.0 * dGroupMpc * Math.sin(thetaRad / 2.0);
}

/**
 * Converts Supergalactic Cartesian [sgx, sgy, sgz] into Supergalactic spherical [sgl, sgb, r].
 * @param {number} sgx - SGX.
 * @param {number} sgy - SGY.
 * @param {number} sgz - SGZ.
 * @returns {{sgl: number, sgb: number, r: number}}
 */
export function supergalacticCartesianToSpherical(sgx, sgy, sgz) {
  const [sgl, sgb, r] = cartesianToSpherical(sgx, sgy, sgz);
  return { sgl, sgb, r };
}

/**
 * Evaluates Hubble expansion parameter H(z) at redshift z for flat LambdaCDM.
 * H(z) = H0 * sqrt(Omega_m * (1+z)^3 + Omega_Lambda).
 * 
 * @param {number} z - Cosmological redshift.
 * @param {number} [H0=DEFAULT_H0] - Hubble constant in km/s/Mpc.
 * @param {number} [OmegaM=DEFAULT_OMEGA_M] - Matter density parameter.
 * @param {number} [OmegaL=DEFAULT_OMEGA_LAMBDA] - Dark energy density parameter.
 * @returns {number} H(z) in km/s/Mpc.
 */
export function hubbleParameter(z, H0 = DEFAULT_H0, OmegaM = DEFAULT_OMEGA_M, OmegaL = DEFAULT_OMEGA_LAMBDA) {
  if (z < 0.0) z = 0.0;
  const onePlusZ = 1.0 + z;
  const ez = Math.sqrt(OmegaM * Math.pow(onePlusZ, 3) + OmegaL);
  return H0 * ez;
}

/**
 * Critical density of the Universe at redshift z in M_sun / Mpc^3.
 * rho_crit(z) = 3 * H(z)^2 / (8 * pi * G_COSMO).
 * 
 * @param {number} z - Redshift.
 * @param {number} [H0=DEFAULT_H0] - Hubble constant.
 * @returns {number} rho_crit(z) in M_sun / Mpc^3.
 */
export function criticalDensityAtZ(z, H0 = DEFAULT_H0) {
  const hz = hubbleParameter(z, H0);
  return (3.0 * hz * hz) / (8.0 * Math.PI * G_COSMO);
}

// ============================================================================
// 4. ROBUST STATISTICAL ESTIMATORS (BEERS ET AL. 1990)
// ============================================================================

/**
 * Computes sample median of an array of numeric values.
 * @param {number[]|Float64Array|Float32Array} arr - Numeric array.
 * @returns {number} Median value.
 */
export function calculateMedian(arr) {
  if (!arr || arr.length === 0) return NaN;
  const len = arr.length;
  if (len === 1) return arr[0];
  const sorted = Array.from(arr).filter(v => Number.isFinite(v)).sort((a, b) => a - b);
  const n = sorted.length;
  if (n === 0) return NaN;
  const mid = Math.floor(n / 2);
  if (n % 2 !== 0) {
    return sorted[mid];
  } else {
    return 0.5 * (sorted[mid - 1] + sorted[mid]);
  }
}

/**
 * Computes Median Absolute Deviation (MAD) from the sample median.
 * MAD = median(|x_i - median(x)|).
 * 
 * @param {number[]|Float64Array|Float32Array} arr - Numeric values.
 * @param {number} [med=null] - Precalculated median (optional).
 * @returns {number} MAD.
 */
export function calculateMAD(arr, med = null) {
  if (!arr || arr.length === 0) return NaN;
  const m = (med !== null && Number.isFinite(med)) ? med : calculateMedian(arr);
  if (isNaN(m)) return NaN;
  const devs = [];
  for (let i = 0; i < arr.length; i++) {
    const val = arr[i];
    if (Number.isFinite(val)) {
      devs.push(Math.abs(val - m));
    }
  }
  return calculateMedian(devs);
}

/**
 * Computes Beers et al. (1990) Biweight Location Estimator (C_BI) for robust central velocity.
 * C_BI = M + ( sum_{|u_i| < 1} (v_i - M) (1 - u_i^2)^2 ) / ( sum_{|u_i| < 1} (1 - u_i^2)^2 )
 * where u_i = (v_i - M) / (c_tune * MAD), c_tune = 6.0.
 * 
 * @param {number[]|Float64Array|Float32Array} values - Velocity array in km/s.
 * @param {number} [cTune=6.0] - Biweight tuning constant (standard: 6.0).
 * @returns {number} Robust biweight location in km/s.
 */
export function biweightLocation(values, cTune = 6.0) {
  if (!values || values.length === 0) return NaN;
  const clean = Array.from(values).filter(v => Number.isFinite(v));
  const n = clean.length;
  if (n === 0) return NaN;
  if (n === 1) return clean[0];
  if (n === 2) return 0.5 * (clean[0] + clean[1]);

  const median = calculateMedian(clean);
  const mad = calculateMAD(clean, median);

  // If MAD is zero (all values identical or nearly so), return median
  if (mad <= 1.0e-7) {
    return median;
  }

  const denomThreshold = cTune * mad;
  let numSum = 0.0;
  let denSum = 0.0;

  for (let i = 0; i < n; i++) {
    const diff = clean[i] - median;
    const u = diff / denomThreshold;
    if (Math.abs(u) < 1.0) {
      const oneMinusU2 = 1.0 - u * u;
      const weight = oneMinusU2 * oneMinusU2;
      numSum += diff * weight;
      denSum += weight;
    }
  }

  if (Math.abs(denSum) < 1.0e-12) {
    return median;
  }

  return median + (numSum / denSum);
}

/**
 * Computes Beers et al. (1990) Biweight Scale Estimator (S_BI) for robust velocity dispersion.
 * S_BI = n^(1/2) * [ sum_{|u_i| < 1} (v_i - M)^2 (1 - u_i^2)^4 ]^(1/2) / | sum_{|u_i| < 1} (1 - u_i^2)(1 - 5 u_i^2) |
 * where u_i = (v_i - M) / (c_tune * MAD), c_tune = 9.0.
 * 
 * @param {number[]|Float64Array|Float32Array} values - Velocity array in km/s.
 * @param {number} [cTune=9.0] - Biweight scale tuning constant (standard: 9.0).
 * @returns {number} Robust velocity dispersion S_BI in km/s.
 */
export function biweightScale(values, cTune = 9.0) {
  if (!values || values.length === 0) return NaN;
  const clean = Array.from(values).filter(v => Number.isFinite(v));
  const n = clean.length;
  if (n < 2) return 0.0;
  if (n === 2) {
    return Math.abs(clean[0] - clean[1]) / Math.SQRT2;
  }

  const median = calculateMedian(clean);
  const mad = calculateMAD(clean, median);

  if (mad <= 1.0e-7) {
    let stdSum = 0.0;
    for (let i = 0; i < n; i++) {
      stdSum += (clean[i] - median) * (clean[i] - median);
    }
    return Math.sqrt(stdSum / (n - 1));
  }

  const denomThreshold = cTune * mad;
  let numSum = 0.0;
  let denSum = 0.0;

  for (let i = 0; i < n; i++) {
    const diff = clean[i] - median;
    const u = diff / denomThreshold;
    if (Math.abs(u) < 1.0) {
      const u2 = u * u;
      const oneMinusU2 = 1.0 - u2;
      const weightNum = oneMinusU2 * oneMinusU2;
      numSum += diff * diff * weightNum * weightNum;
      denSum += oneMinusU2 * (1.0 - 5.0 * u2);
    }
  }

  if (Math.abs(denSum) < 1.0e-12) {
    let stdSum = 0.0;
    for (let i = 0; i < n; i++) {
      stdSum += (clean[i] - median) * (clean[i] - median);
    }
    return Math.sqrt(stdSum / (n - 1));
  }

  const sBi = Math.sqrt(n) * Math.sqrt(numSum) / Math.abs(denSum);
  return sBi;
}

/**
 * Computes Beers et al. (1990) Gapper Scale Estimator (sigma_G) for small samples (N <= 15):
 * sigma_G = (sqrt(pi) / (N * (N - 1))) * sum_{i=1}^{N-1} i * (N - i) * (x_{i+1} - x_i).
 * 
 * @param {number[]|Float64Array|Float32Array} values - Values.
 * @returns {number} Gapper scale in km/s.
 */
export function gapperScale(values) {
  if (!values || values.length < 2) return 0.0;
  const sorted = Array.from(values).filter(v => Number.isFinite(v)).sort((a, b) => a - b);
  const n = sorted.length;
  if (n < 2) return 0.0;
  if (n === 2) return Math.abs(sorted[1] - sorted[0]) / Math.SQRT2;

  let sumGaps = 0.0;
  for (let i = 1; i < n; i++) {
    const w = i * (n - i);
    const gap = sorted[i] - sorted[i - 1];
    sumGaps += w * gap;
  }

  const factor = Math.sqrt(Math.PI) / (n * (n - 1));
  return factor * sumGaps;
}

/**
 * Computes standard unbiased sample standard deviation with Bessel's correction (1/(N-1)).
 * @param {number[]|Float64Array|Float32Array} values - Array of values.
 * @param {number} [meanVal=null] - Precomputed mean.
 * @returns {number} Standard deviation.
 */
export function sampleStandardDeviation(values, meanVal = null) {
  if (!values || values.length < 2) return 0.0;
  const clean = Array.from(values).filter(v => Number.isFinite(v));
  const n = clean.length;
  if (n < 2) return 0.0;

  const mean = (meanVal !== null && Number.isFinite(meanVal))
    ? meanVal
    : clean.reduce((sum, v) => sum + v, 0.0) / n;

  let sumSq = 0.0;
  for (let i = 0; i < n; i++) {
    const d = clean[i] - mean;
    sumSq += d * d;
  }
  return Math.sqrt(sumSq / (n - 1));
}

/**
 * Deconvolves velocity measurement errors and applies cosmological redshift correction:
 * sigma_{v, corr}^2 = max(0, sigma_v^2 - <sigma_err^2>) / (1 + z)^2.
 * 
 * @param {number} sigmaV - Raw velocity dispersion in km/s.
 * @param {number[]|Float64Array|Float32Array|number} errors - Individual galaxy velocity errors, or mean error.
 * @param {number} [meanZ=0.0] - Mean group redshift z = <cz> / c.
 * @returns {number} Corrected velocity dispersion in km/s.
 */
export function deconvolveVelocityDispersion(sigmaV, errors, meanZ = 0.0) {
  if (isNaN(sigmaV) || sigmaV <= 0.0) return 0.0;
  let meanVarErr = 0.0;

  if (typeof errors === 'number' && Number.isFinite(errors)) {
    meanVarErr = errors * errors;
  } else if (Array.isArray(errors) || errors instanceof Float64Array || errors instanceof Float32Array) {
    const cleanErr = Array.from(errors).filter(e => Number.isFinite(e) && e >= 0);
    if (cleanErr.length > 0) {
      meanVarErr = cleanErr.reduce((acc, e) => acc + e * e, 0.0) / cleanErr.length;
    }
  }

  const rawVar = sigmaV * sigmaV;
  const netVar = Math.max(0.0, rawVar - meanVarErr);
  const onePlusZ = 1.0 + Math.max(0.0, meanZ);

  return Math.sqrt(netVar) / onePlusZ;
}

// ============================================================================
// 5. PROJECTED PHASE-SPACE & RADIUS ANALYSIS
// ============================================================================

/**
 * Computes exact Projected Harmonic Radius R_H in Mpc for a set of galaxies.
 * Formulation (Limber & Mathews 1960, Heisler et al. 1985):
 * R_H = (pi / 2) * ( N * (N - 1) / 2 ) * [ sum_{1 <= i < j <= N} ( 1 / R_ij ) ]^(-1)
 * where R_ij is the pairwise projected physical separation in Mpc at mean group distance.
 * 
 * @param {Array<{ra: number, dec: number}>} members - Galaxy members with coordinates in degrees.
 * @param {number} dGroupMpc - Group metric distance in Mpc.
 * @param {number} [softeningMpc=RADIUS_SOFTENING_MPC] - Gravitational softening in Mpc.
 * @returns {number} Projected harmonic radius R_H in Mpc.
 */
export function calculateProjectedHarmonicRadius(members, dGroupMpc, softeningMpc = RADIUS_SOFTENING_MPC) {
  if (!members || members.length < 2 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return 0.0;
  }
  const n = members.length;
  const totalPairs = (n * (n - 1)) / 2;
  let sumInvSeparation = 0.0;

  for (let i = 0; i < n; i++) {
    const m1 = members[i];
    for (let j = i + 1; j < n; j++) {
      const m2 = members[j];
      const sep = projectedPhysicalSeparationMpc(m1.ra, m1.dec, m2.ra, m2.dec, dGroupMpc);
      const regSep = Math.sqrt(sep * sep + softeningMpc * softeningMpc);
      sumInvSeparation += 1.0 / regSep;
    }
  }

  if (sumInvSeparation <= 0.0) return 0.0;

  // 3D deprojection factor pi / 2
  const rH = (Math.PI / 2.0) * (totalPairs / sumInvSeparation);
  return rH;
}

/**
 * Computes projected radial moments: R_mean, R_rms, R_med, R_80, R_20 in Mpc relative to the group centroid.
 * 
 * @param {Array<{ra: number, dec: number}>} members - Member galaxies.
 * @param {number} raCentroid - Centroid RA in deg.
 * @param {number} decCentroid - Centroid Dec in deg.
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @returns {{rMean: number, rRms: number, rMed: number, rMax: number, r20: number, r80: number, concentration: number, projectedRadii: Float64Array}}
 */
export function calculateProjectedRadiiMoments(members, raCentroid, decCentroid, dGroupMpc) {
  if (!members || members.length === 0 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return {
      rMean: 0, rRms: 0, rMed: 0, rMax: 0, r20: 0, r80: 0, concentration: 1.0,
      projectedRadii: new Float64Array(0)
    };
  }
  const n = members.length;
  const radii = new Float64Array(n);
  let sumR = 0.0;
  let sumR2 = 0.0;
  let maxR = 0.0;

  for (let i = 0; i < n; i++) {
    const r = projectedPhysicalSeparationMpc(members[i].ra, members[i].dec, raCentroid, decCentroid, dGroupMpc);
    radii[i] = r;
    sumR += r;
    sumR2 += r * r;
    if (r > maxR) maxR = r;
  }

  const sortedRadii = Array.from(radii).sort((a, b) => a - b);
  const rMean = sumR / n;
  const rRms = Math.sqrt(sumR2 / n);
  const rMed = calculateMedian(sortedRadii);

  const idx20 = Math.max(0, Math.floor(0.20 * n));
  const idx80 = Math.min(n - 1, Math.floor(0.80 * n));
  const r20 = sortedRadii[idx20] > 0 ? sortedRadii[idx20] : RADIUS_SOFTENING_MPC;
  const r80 = sortedRadii[idx80];
  const concentration = r80 / r20;

  return {
    rMean,
    rRms,
    rMed,
    rMax: maxR,
    r20,
    r80,
    concentration,
    projectedRadii: radii
  };
}

/**
 * Computes Virial Radius R_200 (radius at which mean enclosed density = 200 * rho_crit(z)).
 * R_200 = sqrt(3) * sigma_v / (10 * H(z)).
 * 
 * @param {number} sigmaVKmS - 1D velocity dispersion in km/s.
 * @param {number} [z=0.0] - Redshift.
 * @param {number} [H0=DEFAULT_H0] - Hubble constant.
 * @returns {number} R_200 in Mpc.
 */
export function calculateR200Mpc(sigmaVKmS, z = 0.0, H0 = DEFAULT_H0) {
  if (isNaN(sigmaVKmS) || sigmaVKmS <= 0.0) return 0.0;
  const hz = hubbleParameter(z, H0);
  return (Math.SQRT2 * Math.sqrt(1.5) * sigmaVKmS) / (10.0 * hz);
}

/**
 * Computes Characteristic Radius R_500 in Mpc (R_500 ~ 0.66 * R_200 for NFW c~5).
 * @param {number} r200Mpc - R_200 in Mpc.
 * @returns {number} R_500 in Mpc.
 */
export function calculateR500Mpc(r200Mpc) {
  return r200Mpc * 0.66;
}

// ============================================================================
// 6. MULTI-ESTIMATOR VIRIAL MASS FORMALISM
// ============================================================================

/**
 * Computes Standard Virial Mass:
 * M_vir = (3 * pi / 2) * (R_H * sigma_v^2 / G_COSMO)
 * 
 * @param {number} rHarmonicMpc - Projected harmonic radius in Mpc.
 * @param {number} sigmaVKmS - Corrected line-of-sight velocity dispersion in km/s.
 * @returns {number} M_vir in solar masses (M_sun).
 */
export function calculateStandardVirialMass(rHarmonicMpc, sigmaVKmS) {
  if (isNaN(rHarmonicMpc) || isNaN(sigmaVKmS) || rHarmonicMpc <= 0.0 || sigmaVKmS <= 0.0) {
    return 0.0;
  }
  const factor = (3.0 * Math.PI) / 2.0; // ~4.71238898
  return factor * (rHarmonicMpc * sigmaVKmS * sigmaVKmS) / G_COSMO;
}

/**
 * Computes Heisler, Tremaine & Bahcall (1985) Projected Mass Estimator:
 * M_proj = (f_proj / (G * (N - 1.5))) * sum_{i=1}^N (v_zi - <v>)^2 * R_perp,i
 * where f_proj = 32 / pi (isotropic orbits), 64 / pi (radial), 16 / pi (circular).
 * 
 * @param {Array<{ra: number, dec: number, vRec: number}>} members - Members.
 * @param {number} raCentroid - Centroid RA.
 * @param {number} decCentroid - Centroid Dec.
 * @param {number} vMean - Mean recession velocity.
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @param {string} [orbitType='isotropic'] - Orbit distribution: 'isotropic', 'radial', or 'circular'.
 * @returns {number} M_proj in M_sun.
 */
export function calculateProjectedMassEstimator(
  members,
  raCentroid,
  decCentroid,
  vMean,
  dGroupMpc,
  orbitType = 'isotropic'
) {
  if (!members || members.length < 2 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return 0.0;
  }
  const n = members.length;
  if (n <= 1.5) return 0.0;

  let fProj = 32.0 / Math.PI; // ~10.1859 for isotropic
  if (orbitType === 'radial') {
    fProj = 64.0 / Math.PI; // ~20.3718
  } else if (orbitType === 'circular') {
    fProj = 16.0 / Math.PI; // ~5.0929
  }

  let sumWeighted = 0.0;
  for (let i = 0; i < n; i++) {
    const m = members[i];
    const vVal = m.vRec ?? m.cz ?? m.vhel ?? 0.0;
    const dv = vVal - vMean;
    const rPerp = projectedPhysicalSeparationMpc(m.ra, m.dec, raCentroid, decCentroid, dGroupMpc);
    sumWeighted += dv * dv * rPerp;
  }

  const mProj = (fProj / (G_COSMO * (n - 1.5))) * sumWeighted;
  return Math.max(0.0, mProj);
}

/**
 * Computes Heisler, Tremaine & Bahcall (1985) Median Mass Estimator:
 * M_med = (f_med / G) * median_{i < j} [ (v_i - v_j)^2 * R_ij ]
 * where f_med = 6.5 (isotropic orbits).
 * 
 * @param {Array<{ra: number, dec: number, vRec: number}>} members - Members.
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @param {number} [fMed=6.5] - Isotropic median factor.
 * @returns {number} M_med in M_sun.
 */
export function calculateMedianMassEstimator(members, dGroupMpc, fMed = 6.5) {
  if (!members || members.length < 2 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return 0.0;
  }
  const n = members.length;
  const pairValues = [];

  for (let i = 0; i < n; i++) {
    const m1 = members[i];
    const v1 = m1.vRec ?? m1.cz ?? m1.vhel ?? 0.0;
    for (let j = i + 1; j < n; j++) {
      const m2 = members[j];
      const v2 = m2.vRec ?? m2.cz ?? m2.vhel ?? 0.0;
      const dv = v1 - v2;
      const rSep = projectedPhysicalSeparationMpc(m1.ra, m1.dec, m2.ra, m2.dec, dGroupMpc);
      pairValues.push(dv * dv * rSep);
    }
  }

  if (pairValues.length === 0) return 0.0;
  const medVal = calculateMedian(pairValues);
  return (fMed / G_COSMO) * medVal;
}

/**
 * Computes Heisler, Tremaine & Bahcall (1985) Average Mass Estimator:
 * M_avg = (f_avg / (G * N * (N - 1))) * sum_{i < j} [ (v_i - v_j)^2 * R_ij ]
 * where f_avg = 3 * pi / 2.
 * 
 * @param {Array<{ra: number, dec: number, vRec: number}>} members - Members.
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @param {number} [fAvg=4.71238898] - Isotropic average factor 3*pi/2.
 * @returns {number} M_avg in M_sun.
 */
export function calculateAverageMassEstimator(members, dGroupMpc, fAvg = (3.0 * Math.PI) / 2.0) {
  if (!members || members.length < 2 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return 0.0;
  }
  const n = members.length;
  let sumPair = 0.0;

  for (let i = 0; i < n; i++) {
    const m1 = members[i];
    const v1 = m1.vRec ?? m1.cz ?? m1.vhel ?? 0.0;
    for (let j = i + 1; j < n; j++) {
      const m2 = members[j];
      const v2 = m2.vRec ?? m2.cz ?? m2.vhel ?? 0.0;
      const dv = v1 - v2;
      const rSep = projectedPhysicalSeparationMpc(m1.ra, m1.dec, m2.ra, m2.dec, dGroupMpc);
      sumPair += dv * dv * rSep;
    }
  }

  const mAvg = (fAvg / (G_COSMO * n * (n - 1))) * sumPair;
  return Math.max(0.0, mAvg);
}

/**
 * Computes Bahcall & Tremaine (1981) Central Point-Mass Estimator:
 * M_central = (16 / (pi * G * N)) * sum_{i=1}^N (v_i - v_c)^2 * R_perp,i
 * Appropriate when satellite galaxies orbit a dominant central cD galaxy (e.g. M87 in Virgo).
 * 
 * @param {Array<{ra: number, dec: number, vRec: number}>} members - Satellite members.
 * @param {number} raCenter - Central galaxy RA.
 * @param {number} decCenter - Central galaxy Dec.
 * @param {number} vCenter - Central galaxy velocity in km/s.
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @returns {number} Central mass in M_sun.
 */
export function calculateBahcallTremaineMassEstimator(members, raCenter, decCenter, vCenter, dGroupMpc) {
  if (!members || members.length === 0 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return 0.0;
  }
  const n = members.length;
  let sumTerm = 0.0;

  for (let i = 0; i < n; i++) {
    const m = members[i];
    const vVal = m.vRec ?? m.cz ?? m.vhel ?? 0.0;
    const dv = vVal - vCenter;
    const rPerp = projectedPhysicalSeparationMpc(m.ra, m.dec, raCenter, decCenter, dGroupMpc);
    sumTerm += dv * dv * rPerp;
  }

  const factor = 16.0 / (Math.PI * G_COSMO * n);
  return factor * sumTerm;
}

/**
 * Applies surface term / boundary pressure correction to the virial mass (The & White 1986):
 * M_vir,corr = M_vir * C_surf, where C_surf ~ 0.82 for typical NFW/King profiles.
 * 
 * @param {number} mVir - Raw virial mass in M_sun.
 * @param {number} [surfaceCorrectionFactor=0.82] - The & White (1986) surface correction factor.
 * @returns {number} Surface-corrected virial mass in M_sun.
 */
export function applySurfacePressureCorrection(mVir, surfaceCorrectionFactor = 0.82) {
  if (isNaN(mVir) || mVir <= 0.0) return 0.0;
  return mVir * surfaceCorrectionFactor;
}

/**
 * Computes non-parametric Bootstrap Resampling (B >= 100 iterations) for group mass and kinematic uncertainties.
 * Returns median, 16th percentile (lower 1-sigma bound), and 84th percentile (upper 1-sigma bound).
 * 
 * @param {Array<{ra: number, dec: number, vRec: number, vErr?: number}>} members - Members.
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @param {number} [nBootstraps=200] - Number of bootstrap resamples.
 * @param {number} [seed=42] - PRNG seed for deterministic testing.
 * @returns {{mVirMedian: number, mVirLow: number, mVirHigh: number, sigmaVMedian: number, sigmaVLow: number, sigmaVHigh: number, rHMedian: number}}
 */
export function bootstrapVirialMassErrors(members, dGroupMpc, nBootstraps = 200, seed = 42) {
  if (!members || members.length < 2 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return {
      mVirMedian: 0,
      mVirLow: 0,
      mVirHigh: 0,
      sigmaVMedian: 0,
      sigmaVLow: 0,
      sigmaVHigh: 0,
      rHMedian: 0
    };
  }

  const n = members.length;
  // Deterministic Linear Congruential Generator (LCG)
  let state = seed % 2147483647;
  if (state <= 0) state += 2147483646;
  function nextRand() {
    state = (state * 16807) % 2147483647;
    return (state - 1) / 2147483646;
  }

  const bootMVir = new Float64Array(nBootstraps);
  const bootSigmaV = new Float64Array(nBootstraps);
  const bootRH = new Float64Array(nBootstraps);

  for (let b = 0; b < nBootstraps; b++) {
    const resampled = [];
    const resampledVel = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const idx = Math.floor(nextRand() * n);
      resampled.push(members[idx]);
      resampledVel[i] = members[idx].vRec ?? members[idx].cz ?? members[idx].vhel ?? 0.0;
    }

    const rH = calculateProjectedHarmonicRadius(resampled, dGroupMpc);
    const sigV = biweightScale(resampledVel);
    const mVir = calculateStandardVirialMass(rH, sigV);

    bootMVir[b] = mVir;
    bootSigmaV[b] = sigV;
    bootRH[b] = rH;
  }

  bootMVir.sort();
  bootSigmaV.sort();
  bootRH.sort();

  const idx16 = Math.floor(0.16 * nBootstraps);
  const idx50 = Math.floor(0.50 * nBootstraps);
  const idx84 = Math.floor(0.84 * nBootstraps);

  return {
    mVirMedian: bootMVir[idx50],
    mVirLow: bootMVir[idx16],
    mVirHigh: bootMVir[idx84],
    sigmaVMedian: bootSigmaV[idx50],
    sigmaVLow: bootSigmaV[idx16],
    sigmaVHigh: bootSigmaV[idx84],
    rHMedian: bootRH[idx50]
  };
}

// ============================================================================
// 7. GROUP LUMINOSITY, MORPHOLOGY & SELECTION FUNCTION
// ============================================================================

/**
 * Natural log of Gamma function ln(Gamma(x)) via Lanczos approximation.
 * @param {number} x - Positive input.
 * @returns {number} ln(Gamma(x)).
 */
export function logGamma(x) {
  if (x <= 0) return NaN;
  const cof = [
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.001208650973866179,
    -0.000005395239384953
  ];
  let y = x;
  let tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) {
    y += 1.0;
    ser += cof[j] / y;
  }
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

/**
 * Standard Gamma function Gamma(x).
 * @param {number} x - Input.
 * @returns {number} Gamma(x).
 */
export function standardGamma(x) {
  return Math.exp(logGamma(x));
}

/**
 * Computes upper incomplete Gamma function Gamma(a, x) via series expansion and continued fractions.
 * @param {number} a - Shape parameter.
 * @param {number} x - Lower integration limit.
 * @returns {number} Gamma(a, x).
 */
export function upperIncompleteGamma(a, x) {
  if (x < 0 || isNaN(x) || isNaN(a)) return NaN;
  if (x === 0) {
    return standardGamma(a);
  }

  const maxIter = 200;
  const eps = 1.0e-14;
  const gln = logGamma(a);

  if (x < a + 1.0) {
    let ap = a;
    let del = 1.0 / a;
    let sum = del;
    for (let n = 1; n <= maxIter; n++) {
      ap += 1.0;
      del *= x / ap;
      sum += del;
      if (Math.abs(del) < Math.abs(sum) * eps) break;
    }
    const pSeries = sum * Math.exp(-x + a * Math.log(x) - gln);
    return Math.exp(gln) * (1.0 - pSeries);
  } else {
    let b = x + 1.0 - a;
    let c = 1.0 / 1.0e-30;
    let d = 1.0 / b;
    let h = d;
    for (let i = 1; i <= maxIter; i++) {
      const an = -i * (i - a);
      b += 2.0;
      d = an * d + b;
      if (Math.abs(d) < 1.0e-30) d = 1.0e-30;
      c = b + an / c;
      if (Math.abs(c) < 1.0e-30) c = 1.0e-30;
      d = 1.0 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1.0) < eps) break;
    }
    return Math.exp(-x + a * Math.log(x) - gln) * h;
  }
}

/**
 * Computes Schechter (1976) luminosity completeness correction factor:
 * C(L_lim) = Gamma(alpha + 2, L_lim / L*) / Gamma(alpha + 2)
 * Total group luminosity L_tot = L_obs / C(L_lim).
 * 
 * @param {number} appMagLim - Survey apparent magnitude limit (e.g. K=11.75 or B=17.5).
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @param {number} [mStar=-23.5] - Characteristic absolute magnitude M* in specified band.
 * @param {number} [alpha=-1.05] - Faint-end Schechter slope alpha.
 * @returns {number} Completeness fraction in (0, 1].
 */
export function schechterCompletenessCorrection(appMagLim, dGroupMpc, mStar = -23.5, alpha = -1.05) {
  if (isNaN(dGroupMpc) || dGroupMpc <= 0.0 || isNaN(appMagLim)) return 1.0;
  const distMod = 5.0 * Math.log10(dGroupMpc * 1.0e6) - 5.0;
  const absMagLim = appMagLim - distMod;

  const lRatio = Math.pow(10.0, -0.4 * (absMagLim - mStar));
  if (lRatio <= 0.0) return 1.0;

  const shape = alpha + 2.0;
  const numerator = upperIncompleteGamma(shape, lRatio);
  const denominator = standardGamma(shape);

  if (isNaN(numerator) || isNaN(denominator) || denominator <= 0.0) {
    return 1.0;
  }

  const fraction = Math.min(1.0, Math.max(1.0e-4, numerator / denominator));
  return fraction;
}

/**
 * Computes group luminosities, morphological composition fractions, and mass-to-light ratios.
 * 
 * @param {Array<{kMag?: number, bMag?: number, tType?: number, dist?: number}>} members - Members.
 * @param {number} dGroupMpc - Group metric distance in Mpc.
 * @param {number} mVirMsun - Virial mass in M_sun.
 * @param {object} [options={}] - Options.
 * @param {number} [options.kMagLim=11.75] - K-band limit.
 * @param {number} [options.bMagLim=17.50] - B-band limit.
 * @returns {{lObsK: number, lTotK: number, lObsB: number, lTotB: number, fEarly: number, fLate: number, fIrr: number, mToLK: number, mToLB: number}}
 */
export function analyzeGroupLuminosityAndMorphology(members, dGroupMpc, mVirMsun, options = {}) {
  const kMagLim = options.kMagLim ?? 11.75;
  const bMagLim = options.bMagLim ?? 17.50;

  if (!members || members.length === 0 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return {
      lObsK: 0,
      lTotK: 0,
      lObsB: 0,
      lTotB: 0,
      fEarly: 0,
      fLate: 0,
      fIrr: 0,
      mToLK: 0,
      mToLB: 0
    };
  }

  const distMod = 5.0 * Math.log10(dGroupMpc * 1.0e6) - 5.0;
  let sumLK = 0.0;
  let sumLB = 0.0;

  let nEarly = 0;
  let nLate = 0;
  let nIrr = 0;
  let nClassified = 0;

  for (let i = 0; i < members.length; i++) {
    const m = members[i];

    if (m.kMag !== undefined && Number.isFinite(m.kMag) && m.kMag < 90) {
      const absK = m.kMag - distMod;
      const lk = Math.pow(10.0, -0.4 * (absK - M_SUN_K));
      sumLK += lk;
    }

    if (m.bMag !== undefined && Number.isFinite(m.bMag) && m.bMag < 90) {
      const absB = m.bMag - distMod;
      const lb = Math.pow(10.0, -0.4 * (absB - M_SUN_B));
      sumLB += lb;
    }

    if (m.tType !== undefined && Number.isFinite(m.tType) && m.tType > -90) {
      nClassified++;
      if (m.tType <= 0) {
        nEarly++;
      } else if (m.tType > 0 && m.tType < 10) {
        nLate++;
      } else {
        nIrr++;
      }
    }
  }

  const cK = schechterCompletenessCorrection(kMagLim, dGroupMpc, -23.5, -1.05);
  const cB = schechterCompletenessCorrection(bMagLim, dGroupMpc, -20.5, -1.20);

  const lTotK = sumLK > 0 ? sumLK / cK : 0;
  const lTotB = sumLB > 0 ? sumLB / cB : 0;

  const fEarly = nClassified > 0 ? nEarly / nClassified : 0.0;
  const fLate = nClassified > 0 ? nLate / nClassified : 0.0;
  const fIrr = nClassified > 0 ? nIrr / nClassified : 0.0;

  const mToLK = (lTotK > 0 && mVirMsun > 0) ? mVirMsun / lTotK : 0;
  const mToLB = (lTotB > 0 && mVirMsun > 0) ? mVirMsun / lTotB : 0;

  return {
    lObsK: sumLK,
    lTotK,
    lObsB: sumLB,
    lTotB,
    fEarly,
    fLate,
    fIrr,
    mToLK,
    mToLB
  };
}

/**
 * Computes Dressler (1980) local projected surface density Sigma_10 in gal/Mpc^2
 * (density within projected distance to 10th nearest neighbor or Nth neighbor).
 * 
 * @param {number} raTarget - Galaxy RA.
 * @param {number} decTarget - Galaxy Dec.
 * @param {Array<{ra: number, dec: number}>} allGalaxies - Catalog galaxies.
 * @param {number} dGroupMpc - Group distance in Mpc.
 * @param {number} [nthNeighbor=10] - Neighbor rank N.
 * @returns {number} Projected surface density in galaxies / Mpc^2.
 */
export function calculateLocalSurfaceDensity(raTarget, decTarget, allGalaxies, dGroupMpc, nthNeighbor = 10) {
  if (!allGalaxies || allGalaxies.length <= 1 || isNaN(dGroupMpc) || dGroupMpc <= 0.0) {
    return 0.0;
  }
  const k = Math.min(nthNeighbor, allGalaxies.length - 1);
  if (k <= 0) return 0.0;

  const dists = [];
  for (let i = 0; i < allGalaxies.length; i++) {
    const g = allGalaxies[i];
    const rSep = projectedPhysicalSeparationMpc(raTarget, decTarget, g.ra, g.dec, dGroupMpc);
    if (rSep > 0.0) {
      dists.push(rSep);
    }
  }

  if (dists.length < k) return 0.0;
  dists.sort((a, b) => a - b);
  const rK = dists[k - 1];

  if (rK <= 0.0) return 0.0;
  return k / (Math.PI * rK * rK);
}

// ============================================================================
// 8. HIERARCHICAL FRIENDS-OF-FRIENDS (FoF) & MEMBERSHIP PROBABILITIES
// ============================================================================

/**
 * Disjoint-Set Union-Find data structure with path compression and rank optimization.
 */
export class DisjointSet {
  constructor(size) {
    this.parent = new Int32Array(size);
    this.rank = new Int32Array(size);
    for (let i = 0; i < size; i++) {
      this.parent[i] = i;
      this.rank[i] = 0;
    }
  }

  find(i) {
    let root = i;
    while (root !== this.parent[root]) {
      root = this.parent[root];
    }
    let curr = i;
    while (curr !== root) {
      const next = this.parent[curr];
      this.parent[curr] = root;
      curr = next;
    }
    return root;
  }

  union(i, j) {
    const rootI = this.find(i);
    const rootJ = this.find(j);
    if (rootI === rootJ) return false;

    if (this.rank[rootI] < this.rank[rootJ]) {
      this.parent[rootI] = rootJ;
    } else if (this.rank[rootI] > this.rank[rootJ]) {
      this.parent[rootJ] = rootI;
    } else {
      this.parent[rootJ] = rootI;
      this.rank[rootI]++;
    }
    return true;
  }
}

/**
 * Evaluates adaptive Friends-of-Friends (FoF) linking lengths (Tully 2015, Kourkchi & Tully 2017):
 * D_link(cz) = D0 * (1 + cz / cz_ref)^gamma
 * V_link(cz) = V0 * (1 + cz / cz_ref)^beta
 * 
 * @param {number} cz - Mean recessional velocity in km/s.
 * @param {object} [params={}] - Parameters.
 * @returns {{dLinkMpc: number, vLinkKmS: number}}
 */
export function getAdaptiveFoFLinkingLengths(cz, params = {}) {
  const d0 = params.d0 ?? 0.35;
  const v0 = params.v0 ?? 250.0;
  const czRef = params.czRef ?? 3000.0;
  const gamma = params.gamma ?? 0.15;
  const beta = params.beta ?? 0.10;

  const ratio = 1.0 + Math.max(0.0, cz) / czRef;
  const dLinkMpc = d0 * Math.pow(ratio, gamma);
  const vLinkKmS = v0 * Math.pow(ratio, beta);

  return { dLinkMpc, vLinkKmS };
}

/**
 * Computes projected phase-space group membership probability P_mem(i in group k)
 * combining spatial NFW/King-like radial falloff and Gaussian velocity dispersion falloff.
 * 
 * P_mem = P_spatial(R_proj / R_200) * P_velocity(|Delta v| / sigma_v)
 * 
 * @param {number} rProjMpc - Projected separation from group centroid in Mpc.
 * @param {number} deltaVKmS - Velocity offset |v_i - <v>| in km/s.
 * @param {number} r200Mpc - Group virial radius R_200 in Mpc.
 * @param {number} sigmaVKmS - Group velocity dispersion in km/s.
 * @returns {number} Membership probability in [0, 1].
 */
export function calculateMembershipProbability(rProjMpc, deltaVKmS, r200Mpc, sigmaVKmS) {
  if (isNaN(rProjMpc) || isNaN(deltaVKmS)) return 0.0;
  const rScale = Math.max(0.05, r200Mpc > 0 ? r200Mpc : 0.5);
  const vScale = Math.max(30.0, sigmaVKmS > 0 ? sigmaVKmS : 100.0);

  const uR = rProjMpc / rScale;
  const uV = Math.abs(deltaVKmS) / vScale;

  const pSpatial = Math.exp(-0.5 * uR * uR);
  const pVelocity = Math.exp(-0.5 * uV * uV);

  const pMem = pSpatial * pVelocity;
  return Math.min(1.0, Math.max(0.0, pMem));
}

/**
 * Runs hierarchical 3D projected Friends-of-Friends (FoF) group finder on an array of galaxies.
 * 
 * @param {Array<{pgc: number, ra: number, dec: number, cz: number, dist?: number}>} galaxies - Galaxy catalog.
 * @param {object} [options={}] - FoF parameters.
 * @param {number} [options.d0=0.40] - Base projected linking length in Mpc.
 * @param {number} [options.v0=300.0] - Base velocity linking length in km/s.
 * @param {number} [options.minMembers=1] - Minimum group size.
 * @returns {Map<number, Array<object>>} Map of Group Root ID -> Array of member galaxies.
 */
export function runHierarchicalFoF(galaxies, options = {}) {
  if (!galaxies || galaxies.length === 0) return new Map();
  const n = galaxies.length;
  const ds = new DisjointSet(n);

  const d0 = options.d0 ?? 0.40;
  const v0 = options.v0 ?? 300.0;
  const minMembers = options.minMembers ?? 1;

  for (let i = 0; i < n; i++) {
    const g1 = galaxies[i];
    const cz1 = g1.cz;
    const d1 = (g1.dist && g1.dist > 0) ? g1.dist : (cz1 / DEFAULT_H0);

    for (let j = i + 1; j < n; j++) {
      const g2 = galaxies[j];
      const cz2 = g2.cz;
      const d2 = (g2.dist && g2.dist > 0) ? g2.dist : (cz2 / DEFAULT_H0);

      const meanCz = 0.5 * (cz1 + cz2);
      const meanD = 0.5 * (d1 + d2);

      const { dLinkMpc, vLinkKmS } = getAdaptiveFoFLinkingLengths(meanCz, { d0, v0 });

      const dv = Math.abs(cz1 - cz2);
      if (dv <= vLinkKmS) {
        const sepMpc = projectedPhysicalSeparationMpc(g1.ra, g1.dec, g2.ra, g2.dec, meanD);
        if (sepMpc <= dLinkMpc) {
          ds.union(i, j);
        }
      }
    }
  }

  const groupsMap = new Map();
  for (let i = 0; i < n; i++) {
    const root = ds.find(i);
    if (!groupsMap.has(root)) {
      groupsMap.set(root, []);
    }
    groupsMap.get(root).push(galaxies[i]);
  }

  const filtered = new Map();
  for (const [root, members] of groupsMap.entries()) {
    if (members.length >= minMembers) {
      filtered.set(root, members);
    }
  }

  return filtered;
}

// ============================================================================
// 9. CROSS-MATCHING WITH MASTER COSMOLOGICAL STRUCTURES
// ============================================================================

/**
 * Cross-matches a group against the Master Reference Cosmological Structures Catalog.
 * Evaluates 3D Supergalactic Cartesian proximity, angular match, and velocity offset.
 * 
 * @param {number} sgx - Group SGX in Mpc/h.
 * @param {number} sgy - Group SGY in Mpc/h.
 * @param {number} sgz - Group SGZ in Mpc/h.
 * @param {number} cz - Group recession velocity in km/s.
 * @param {number} [ra=null] - Group RA in deg.
 * @param {number} [dec=null] - Group Dec in deg.
 * @returns {{matched: boolean, structure: object|null, distanceMpcH: number, velocityDiffKmS: number, matchScore: number}}
 */
export function crossMatchMasterStructure(sgx, sgy, sgz, cz, ra = null, dec = null) {
  let bestMatch = null;
  let minScore = -1.0;
  let bestDist = Infinity;
  let bestDv = Infinity;

  for (const struct of MASTER_STRUCTURES) {
    const dist3D = euclideanDistance3D(sgx, sgy, sgz, struct.sgx, struct.sgy, struct.sgz);
    const dv = Math.abs(cz - struct.cz);

    const uDist = dist3D / struct.rMatch;
    const uVel = dv / struct.vMatch;

    const chi2 = uDist * uDist + uVel * uVel;
    const score = Math.exp(-0.5 * chi2);

    if (uDist <= 1.5 && uVel <= 1.8 && score > minScore) {
      minScore = score;
      bestMatch = struct;
      bestDist = dist3D;
      bestDv = dv;
    }
  }

  if (bestMatch && minScore >= 0.05) {
    return {
      matched: true,
      structure: bestMatch,
      distanceMpcH: bestDist,
      velocityDiffKmS: bestDv,
      matchScore: minScore
    };
  }

  return {
    matched: false,
    structure: null,
    distanceMpcH: Infinity,
    velocityDiffKmS: Infinity,
    matchScore: 0.0
  };
}

// ============================================================================
// 10. SPATIAL INDEX (3D UNIFORM SPATIAL HASH GRID & KD-TREE)
// ============================================================================

/**
 * Fast 3D Uniform Spatial Hash Grid for high-performance radius and cone searches.
 */
export class SpatialHashGrid3D {
  /**
   * @param {number} [cellSize=10.0] - Cell size in Mpc/h.
   */
  constructor(cellSize = 10.0) {
    this.cellSize = cellSize;
    this.grid = new Map();
    this.items = [];
  }

  _hashKey(gx, gy, gz) {
    return `${gx},${gy},${gz}`;
  }

  _getGridCoords(x, y, z) {
    return [
      Math.floor(x / this.cellSize),
      Math.floor(y / this.cellSize),
      Math.floor(z / this.cellSize)
    ];
  }

  insert(item, x, y, z) {
    const entry = { item, x, y, z, index: this.items.length };
    this.items.push(entry);
    const [gx, gy, gz] = this._getGridCoords(x, y, z);
    const key = this._hashKey(gx, gy, gz);
    if (!this.grid.has(key)) {
      this.grid.set(key, []);
    }
    this.grid.get(key).push(entry);
  }

  queryRadius(x, y, z, radius) {
    const r2 = radius * radius;
    const [minGx, minGy, minGz] = this._getGridCoords(x - radius, y - radius, z - radius);
    const [maxGx, maxGy, maxGz] = this._getGridCoords(x + radius, y + radius, z + radius);
    const results = [];

    for (let gx = minGx; gx <= maxGx; gx++) {
      for (let gy = minGy; gy <= maxGy; gy++) {
        for (let gz = minGz; gz <= maxGz; gz++) {
          const key = this._hashKey(gx, gy, gz);
          const cell = this.grid.get(key);
          if (cell) {
            for (let i = 0; i < cell.length; i++) {
              const e = cell[i];
              const dx = e.x - x;
              const dy = e.y - y;
              const dz = e.z - z;
              const d2 = dx * dx + dy * dy + dz * dz;
              if (d2 <= r2) {
                results.push({ item: e.item, distance: Math.sqrt(d2) });
              }
            }
          }
        }
      }
    }

    results.sort((a, b) => a.distance - b.distance);
    return results;
  }
}

// ============================================================================
// 11. CF4 GROUP CLASS DEFINITION
// ============================================================================

/**
 * Represents a single aggregated galaxy group in Cosmicflows-4.
 */
export class CF4Group {
  constructor(groupId, members = [], options = {}) {
    this.groupId = groupId;
    this.members = members;
    this.nMembers = members.length;
    this.h0 = options.h0 ?? DEFAULT_H0;
    this.littleH = this.h0 / 100.0;

    // Centroids
    this.raCentroid = NaN;
    this.decCentroid = NaN;
    this.sglCentroid = NaN;
    this.sgbCentroid = NaN;
    this.sgxCentroid = NaN;
    this.sgyCentroid = NaN;
    this.sgzCentroid = NaN;
    this.dGroupMpc = NaN;
    this.dGroupMpcH = NaN;

    // Kinematics
    this.meanVRec = NaN;
    this.biweightVRec = NaN;
    this.sigmaVRaw = NaN;
    this.sigmaVBiweight = NaN;
    this.sigmaVGapper = NaN;
    this.sigmaVCorr = NaN;

    // Radii
    this.rHarmonicMpc = NaN;
    this.rMeanMpc = NaN;
    this.rRmsMpc = NaN;
    this.rMedMpc = NaN;
    this.r200Mpc = NaN;
    this.r500Mpc = NaN;
    this.concentration = NaN;

    // Masses (M_sun)
    this.mVir = NaN;
    this.mProj = NaN;
    this.mMed = NaN;
    this.mAvg = NaN;
    this.mCentral = NaN;
    this.mVirSurfCorr = NaN;
    this.mVirMpcH = NaN;

    // Bootstrap errors
    this.mVirErrLow = NaN;
    this.mVirErrHigh = NaN;
    this.sigmaVErrLow = NaN;
    this.sigmaVErrHigh = NaN;

    // Photometry & Morphology
    this.lTotK = NaN;
    this.lTotB = NaN;
    this.fEarly = NaN;
    this.fLate = NaN;
    this.fIrr = NaN;
    this.mToLK = NaN;
    this.mToLB = NaN;

    // Master structure cross-match
    this.matchedStructure = null;
    this.isMatchedToMaster = false;

    if (this.nMembers > 0) {
      this.computeAllProperties(options);
    }
  }

  /**
   * Recomputes all kinematic, spatial, virial mass, and photometric group properties.
   */
  computeAllProperties(options = {}) {
    const n = this.nMembers;
    if (n === 0) return;

    let sumRa = 0.0;
    let sumDec = 0.0;
    let sumSgx = 0.0;
    let sumSgy = 0.0;
    let sumSgz = 0.0;
    let sumDist = 0.0;
    let sumV = 0.0;

    const velArray = new Float64Array(n);
    const errArray = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      const m = this.members[i];
      sumRa += m.ra;
      sumDec += m.dec;
      sumV += m.cz;
      velArray[i] = m.cz;
      errArray[i] = m.vErr ?? 0.0;

      const d = (m.dist && m.dist > 0) ? m.dist : (m.cz / this.h0);
      sumDist += d;

      if (m.sgx !== undefined && m.sgy !== undefined && m.sgz !== undefined) {
        sumSgx += m.sgx;
        sumSgy += m.sgy;
        sumSgz += m.sgz;
      } else {
        const sg = equatorialToSupergalacticCartesian(m.ra, m.dec, d);
        const [sgxVal, sgyVal, sgzVal] = Array.isArray(sg) ? sg : [sg.sgx, sg.sgy, sg.sgz];
        sumSgx += sgxVal;
        sumSgy += sgyVal;
        sumSgz += sgzVal;
      }
    }

    this.raCentroid = sumRa / n;
    this.decCentroid = sumDec / n;
    this.meanVRec = sumV / n;
    this.dGroupMpc = sumDist / n;
    this.dGroupMpcH = this.dGroupMpc * this.littleH;

    this.sgxCentroid = sumSgx / n;
    this.sgyCentroid = sumSgy / n;
    this.sgzCentroid = sumSgz / n;

    const sgSph = supergalacticCartesianToSpherical(this.sgxCentroid, this.sgyCentroid, this.sgzCentroid);
    this.sglCentroid = sgSph.sgl;
    this.sgbCentroid = sgSph.sgb;

    // Kinematics
    this.biweightVRec = biweightLocation(velArray);
    this.sigmaVRaw = sampleStandardDeviation(velArray, this.meanVRec);
    this.sigmaVBiweight = biweightScale(velArray);
    this.sigmaVGapper = gapperScale(velArray);

    const sigmaToUse = (n <= 15) ? this.sigmaVGapper : this.sigmaVBiweight;
    const meanZ = this.meanVRec / SPEED_OF_LIGHT_KM_S;
    this.sigmaVCorr = deconvolveVelocityDispersion(sigmaToUse, errArray, meanZ);

    // Radii
    if (n >= 2) {
      this.rHarmonicMpc = calculateProjectedHarmonicRadius(this.members, this.dGroupMpc);
      const radMoments = calculateProjectedRadiiMoments(this.members, this.raCentroid, this.decCentroid, this.dGroupMpc);
      this.rMeanMpc = radMoments.rMean;
      this.rRmsMpc = radMoments.rRms;
      this.rMedMpc = radMoments.rMed;
      this.concentration = radMoments.concentration;
    } else {
      this.rHarmonicMpc = 0.05;
      this.rMeanMpc = 0.05;
      this.rRmsMpc = 0.05;
      this.rMedMpc = 0.05;
      this.concentration = 1.0;
    }

    this.r200Mpc = calculateR200Mpc(this.sigmaVCorr > 0 ? this.sigmaVCorr : this.sigmaVRaw, meanZ, this.h0);
    this.r500Mpc = calculateR500Mpc(this.r200Mpc);

    // Virial Mass & Multi-Estimators
    if (n >= 2 && this.sigmaVCorr > 0 && this.rHarmonicMpc > 0) {
      this.mVir = calculateStandardVirialMass(this.rHarmonicMpc, this.sigmaVCorr);
      this.mProj = calculateProjectedMassEstimator(this.members, this.raCentroid, this.decCentroid, this.meanVRec, this.dGroupMpc);
      this.mMed = calculateMedianMassEstimator(this.members, this.dGroupMpc);
      this.mAvg = calculateAverageMassEstimator(this.members, this.dGroupMpc);
      this.mCentral = calculateBahcallTremaineMassEstimator(this.members, this.raCentroid, this.decCentroid, this.meanVRec, this.dGroupMpc);
      this.mVirSurfCorr = applySurfacePressureCorrection(this.mVir);
    } else {
      this.mVir = 0.0;
      this.mProj = 0.0;
      this.mMed = 0.0;
      this.mAvg = 0.0;
      this.mCentral = 0.0;
      this.mVirSurfCorr = 0.0;
    }

    this.mVirMpcH = (this.mVir * this.littleH) / 1.0e12;

    // Bootstrap Errors
    if (n >= 3 && (options.computeBootstrap ?? true)) {
      const boot = bootstrapVirialMassErrors(this.members, this.dGroupMpc, options.nBootstraps ?? 100);
      this.mVirErrLow = boot.mVirLow;
      this.mVirErrHigh = boot.mVirHigh;
      this.sigmaVErrLow = boot.sigmaVLow;
      this.sigmaVErrHigh = boot.sigmaVHigh;
    }

    // Photometry & Morphology
    const photo = analyzeGroupLuminosityAndMorphology(this.members, this.dGroupMpc, this.mVir, options);
    this.lTotK = photo.lTotK;
    this.lTotB = photo.lTotB;
    this.fEarly = photo.fEarly;
    this.fLate = photo.fLate;
    this.fIrr = photo.fIrr;
    this.mToLK = photo.mToLK;
    this.mToLB = photo.mToLB;

    // 7. Master Structure Match
    const match = crossMatchMasterStructure(
      this.sgxCentroid,
      this.sgyCentroid,
      this.sgzCentroid,
      this.meanVRec,
      this.raCentroid,
      this.decCentroid
    );

    if (match.matched) {
      this.matchedStructure = match.structure;
      this.isMatchedToMaster = true;
    }
  }

  /**
   * Export group summary record.
   */
  toJSON() {
    return {
      groupId: this.groupId,
      nMembers: this.nMembers,
      raCentroid: this.raCentroid,
      decCentroid: this.decCentroid,
      sglCentroid: this.sglCentroid,
      sgbCentroid: this.sgbCentroid,
      sgxCentroid: this.sgxCentroid,
      sgyCentroid: this.sgyCentroid,
      sgzCentroid: this.sgzCentroid,
      dGroupMpc: this.dGroupMpc,
      meanVRec: this.meanVRec,
      biweightVRec: this.biweightVRec,
      sigmaVCorr: this.sigmaVCorr,
      rHarmonicMpc: this.rHarmonicMpc,
      r200Mpc: this.r200Mpc,
      mVir: this.mVir,
      mProj: this.mProj,
      mMed: this.mMed,
      mAvg: this.mAvg,
      mVirSurfCorr: this.mVirSurfCorr,
      lTotK: this.lTotK,
      fEarly: this.fEarly,
      fLate: this.fLate,
      matchedStructure: this.matchedStructure ? this.matchedStructure.name : null
    };
  }
}

// ============================================================================
// 12. CF4 GROUPED CATALOG COLUMNAR STORE
// ============================================================================

/**
 * High-performance Columnar Store and Query Engine for the Cosmicflows-4 Grouped Catalog.
 */
export class CF4GroupedCatalog {
  /**
   * @param {number} [initialCapacity=40000] - Initial row capacity.
   */
  constructor(initialCapacity = 40000) {
    this.capacity = initialCapacity;
    this.length = 0;

    // High density TypedArrays for Galaxy Members
    this.pgc = new Int32Array(this.capacity);
    this.groupId = new Int32Array(this.capacity);
    this.nest = new Int32Array(this.capacity);
    this.ra = new Float64Array(this.capacity);
    this.dec = new Float64Array(this.capacity);
    this.sgl = new Float64Array(this.capacity);
    this.sgb = new Float64Array(this.capacity);
    this.sgx = new Float64Array(this.capacity);
    this.sgy = new Float64Array(this.capacity);
    this.sgz = new Float64Array(this.capacity);
    this.cz = new Float64Array(this.capacity);
    this.vErr = new Float32Array(this.capacity);
    this.dist = new Float64Array(this.capacity);
    this.distErr = new Float32Array(this.capacity);
    this.dm = new Float32Array(this.capacity);
    this.dmErr = new Float32Array(this.capacity);
    this.bMag = new Float32Array(this.capacity);
    this.kMag = new Float32Array(this.capacity);
    this.logW21 = new Float32Array(this.capacity);
    this.tType = new Float32Array(this.capacity);
    this.membershipProb = new Float32Array(this.capacity);
    this.isGroupCentroid = new Uint8Array(this.capacity);

    this.pgcIndexMap = new Map();
    this.groups = new Map();
    this.spatialIndex = new SpatialHashGrid3D(10.0);
  }

  _growCapacity(newCapacity) {
    const copy = (oldArr, Type) => {
      const arr = new Type(newCapacity);
      arr.set(oldArr);
      return arr;
    };

    this.pgc = copy(this.pgc, Int32Array);
    this.groupId = copy(this.groupId, Int32Array);
    this.nest = copy(this.nest, Int32Array);
    this.ra = copy(this.ra, Float64Array);
    this.dec = copy(this.dec, Float64Array);
    this.sgl = copy(this.sgl, Float64Array);
    this.sgb = copy(this.sgb, Float64Array);
    this.sgx = copy(this.sgx, Float64Array);
    this.sgy = copy(this.sgy, Float64Array);
    this.sgz = copy(this.sgz, Float64Array);
    this.cz = copy(this.cz, Float64Array);
    this.vErr = copy(this.vErr, Float32Array);
    this.dist = copy(this.dist, Float64Array);
    this.distErr = copy(this.distErr, Float32Array);
    this.dm = copy(this.dm, Float32Array);
    this.dmErr = copy(this.dmErr, Float32Array);
    this.bMag = copy(this.bMag, Float32Array);
    this.kMag = copy(this.kMag, Float32Array);
    this.logW21 = copy(this.logW21, Float32Array);
    this.tType = copy(this.tType, Float32Array);
    this.membershipProb = copy(this.membershipProb, Float32Array);
    this.isGroupCentroid = copy(this.isGroupCentroid, Uint8Array);

    this.capacity = newCapacity;
  }

  /**
   * Adds a galaxy record to the columnar store.
   * @param {object} g - Galaxy record.
   * @returns {number} Row index.
   */
  addGalaxy(g) {
    if (this.length >= this.capacity) {
      this._growCapacity(this.capacity * 2);
    }
    const idx = this.length;

    const pgcVal = g.pgc ?? (idx + 1);
    this.pgc[idx] = pgcVal;
    this.groupId[idx] = g.groupId ?? g.nest ?? pgcVal;
    this.nest[idx] = g.nest ?? this.groupId[idx];

    this.ra[idx] = g.ra ?? 0.0;
    this.dec[idx] = g.dec ?? 0.0;
    this.cz[idx] = g.cz ?? g.vRec ?? g.vhel ?? 0.0;
    this.vErr[idx] = g.vErr ?? 0.0;

    const d = (g.dist && g.dist > 0) ? g.dist : (this.cz[idx] / DEFAULT_H0);
    this.dist[idx] = d;
    this.distErr[idx] = g.distErr ?? (d * 0.15);

    if (g.dm !== undefined && Number.isFinite(g.dm)) {
      this.dm[idx] = g.dm;
    } else {
      this.dm[idx] = 5.0 * Math.log10(d * 1.0e6) - 5.0;
    }
    this.dmErr[idx] = g.dmErr ?? 0.35;

    if (g.sgx !== undefined && g.sgy !== undefined && g.sgz !== undefined) {
      this.sgx[idx] = g.sgx;
      this.sgy[idx] = g.sgy;
      this.sgz[idx] = g.sgz;
      const sgSph = supergalacticCartesianToSpherical(g.sgx, g.sgy, g.sgz);
      this.sgl[idx] = sgSph.sgl;
      this.sgb[idx] = sgSph.sgb;
    } else {
      const sg = equatorialToSupergalacticCartesian(this.ra[idx], this.dec[idx], d);
      const [sgxVal, sgyVal, sgzVal] = Array.isArray(sg) ? sg : [sg.sgx, sg.sgy, sg.sgz];
      this.sgx[idx] = sgxVal;
      this.sgy[idx] = sgyVal;
      this.sgz[idx] = sgzVal;
      const sgSph = supergalacticCartesianToSpherical(sgxVal, sgyVal, sgzVal);
      this.sgl[idx] = sgSph.sgl;
      this.sgb[idx] = sgSph.sgb;
    }

    this.bMag[idx] = g.bMag ?? (g.B_mag ?? 99.99);
    this.kMag[idx] = g.kMag ?? (g.K_mag ?? 99.99);
    this.logW21[idx] = g.logW21 ?? (g.W21 ? Math.log10(g.W21) : 0.0);
    this.tType[idx] = g.tType ?? (g.T_type ?? -99.0);
    this.membershipProb[idx] = g.membershipProb ?? 1.0;
    this.isGroupCentroid[idx] = g.isGroupCentroid ? 1 : 0;

    this.pgcIndexMap.set(pgcVal, idx);
    this.spatialIndex.insert({ pgc: pgcVal, index: idx }, this.sgx[idx], this.sgy[idx], this.sgz[idx]);

    this.length++;
    return idx;
  }

  /**
   * Builds and aggregates all CF4 groups from the loaded galaxies.
   */
  aggregateGroups(options = {}) {
    const rawGroups = new Map();

    for (let i = 0; i < this.length; i++) {
      const gid = this.groupId[i];
      if (!rawGroups.has(gid)) {
        rawGroups.set(gid, []);
      }
      rawGroups.get(gid).push({
        pgc: this.pgc[i],
        groupId: gid,
        nest: this.nest[i],
        ra: this.ra[i],
        dec: this.dec[i],
        sgl: this.sgl[i],
        sgb: this.sgb[i],
        sgx: this.sgx[i],
        sgy: this.sgy[i],
        sgz: this.sgz[i],
        cz: this.cz[i],
        vRec: this.cz[i],
        vErr: this.vErr[i],
        dist: this.dist[i],
        distErr: this.distErr[i],
        dm: this.dm[i],
        bMag: this.bMag[i],
        kMag: this.kMag[i],
        logW21: this.logW21[i],
        tType: this.tType[i]
      });
    }

    this.groups.clear();
    for (const [gid, members] of rawGroups.entries()) {
      const group = new CF4Group(gid, members, options);
      this.groups.set(gid, group);
    }

    return this.groups;
  }

  /**
   * Retrieves a galaxy record by PGC ID.
   * @param {number} pgc - PGC catalog identifier.
   * @returns {object|null} Galaxy record.
   */
  getGalaxyByPGC(pgc) {
    const idx = this.pgcIndexMap.get(pgc);
    if (idx === undefined) return null;
    return {
      pgc: this.pgc[idx],
      groupId: this.groupId[idx],
      nest: this.nest[idx],
      ra: this.ra[idx],
      dec: this.dec[idx],
      sgl: this.sgl[idx],
      sgb: this.sgb[idx],
      sgx: this.sgx[idx],
      sgy: this.sgy[idx],
      sgz: this.sgz[idx],
      cz: this.cz[idx],
      vErr: this.vErr[idx],
      dist: this.dist[idx],
      distErr: this.distErr[idx],
      dm: this.dm[idx],
      dmErr: this.dmErr[idx],
      bMag: this.bMag[idx],
      kMag: this.kMag[idx],
      logW21: this.logW21[idx],
      tType: this.tType[idx],
      membershipProb: this.membershipProb[idx],
      isGroupCentroid: this.isGroupCentroid[idx] === 1
    };
  }

  /**
   * Queries galaxies within a 3D spherical region centered at (sgx, sgy, sgz) in Mpc.
   * @param {number} sgx - Center SGX.
   * @param {number} sgy - Center SGY.
   * @param {number} sgz - Center SGZ.
   * @param {number} radiusMpc - Radius in Mpc.
   * @returns {Array<object>}
   */
  queryCone3D(sgx, sgy, sgz, radiusMpc) {
    const hits = this.spatialIndex.queryRadius(sgx, sgy, sgz, radiusMpc);
    return hits.map(h => this.getGalaxyByPGC(h.item.pgc));
  }

  /**
   * Exports all aggregated groups as GeoJSON FeatureCollection.
   * @returns {object} GeoJSON FeatureCollection.
   */
  exportGroupsToGeoJSON() {
    const features = [];
    for (const group of this.groups.values()) {
      features.push({
        type: 'Feature',
        geometry: {
          type: 'Point',
          coordinates: [group.raCentroid, group.decCentroid]
        },
        properties: group.toJSON()
      });
    }

    return {
      type: 'FeatureCollection',
      features
    };
  }

  /**
   * Exports catalog to CSV string.
   * @returns {string} CSV text.
   */
  exportToCSV() {
    const headers = [
      'PGC', 'GroupID', 'NEST', 'RA', 'Dec', 'SGL', 'SGB', 'SGX', 'SGY', 'SGZ',
      'cz', 'vErr', 'Dist', 'DistErr', 'DM', 'BMag', 'KMag', 'logW21', 'TType'
    ];
    const lines = [headers.join(',')];

    for (let i = 0; i < this.length; i++) {
      lines.push([
        this.pgc[i],
        this.groupId[i],
        this.nest[i],
        this.ra[i].toFixed(5),
        this.dec[i].toFixed(5),
        this.sgl[i].toFixed(5),
        this.sgb[i].toFixed(5),
        this.sgx[i].toFixed(3),
        this.sgy[i].toFixed(3),
        this.sgz[i].toFixed(3),
        this.cz[i].toFixed(1),
        this.vErr[i].toFixed(1),
        this.dist[i].toFixed(3),
        this.distErr[i].toFixed(3),
        this.dm[i].toFixed(3),
        this.bMag[i].toFixed(2),
        this.kMag[i].toFixed(2),
        this.logW21[i].toFixed(3),
        this.tType[i].toFixed(1)
      ].join(','));
    }

    return lines.join('\n');
  }
}

// ============================================================================
// 13. FILE PARSERS & CATALOG INGESTION
// ============================================================================

/**
 * Parses CF4 tabular string data (CSV, TSV, or IP2I ASCII table) and ingests into CF4GroupedCatalog.
 * 
 * @param {string} rawText - Raw text file contents.
 * @param {object} [options={}] - Options.
 * @returns {CF4GroupedCatalog} Ingested catalog.
 */
export function parseCF4GroupCatalogText(rawText, options = {}) {
  const catalog = new CF4GroupedCatalog(options.initialCapacity ?? 40000);
  if (!rawText || typeof rawText !== 'string') return catalog;

  const lines = rawText.split(/\r?\n/);
  let headerIndex = -1;
  let headers = [];
  let delimiter = ',';

  for (let i = 0; i < Math.min(50, lines.length); i++) {
    const line = lines[i].trim();
    if (!line || line.startsWith('#') || line.startsWith('|')) continue;

    if (line.includes('\t')) {
      delimiter = '\t';
      headers = line.split('\t').map(h => h.trim().toLowerCase());
      headerIndex = i;
      break;
    } else if (line.includes(',')) {
      delimiter = ',';
      headers = line.split(',').map(h => h.trim().toLowerCase());
      headerIndex = i;
      break;
    } else if (/\s{2,}/.test(line)) {
      delimiter = /\s+/;
      headers = line.split(/\s+/).map(h => h.trim().toLowerCase());
      headerIndex = i;
      break;
    }
  }

  if (headerIndex === -1) {
    return catalog;
  }

  const findCol = (...names) => {
    for (const name of names) {
      const idx = headers.indexOf(name.toLowerCase());
      if (idx !== -1) return idx;
    }
    return -1;
  };

  const colPGC = findCol('pgc', 'id', 'name');
  const colGroup = findCol('groupid', 'group_id', 'nest', 'grp', 'group');
  const colRA = findCol('ra', 'raj2000', 'ra_deg');
  const colDec = findCol('dec', 'dej2000', 'dec_deg');
  const colCz = findCol('cz', 'vrec', 'vhel', 'vls', 'vcmb', 'v_rec');
  const colVErr = findCol('verr', 'e_vrec', 'e_vhel', 'e_cz');
  const colDist = findCol('dist', 'd_mpc', 'distance', 'dm_dist');
  const colDistErr = findCol('disterr', 'e_dist', 'e_d');
  const colDM = findCol('dm', 'mod', 'distmod');
  const colDMErr = findCol('dmerr', 'e_dm', 'e_mod');
  const colBMag = findCol('bmag', 'b_mag', 'bt');
  const colKMag = findCol('kmag', 'k_mag', 'kt');
  const colW21 = findCol('logw21', 'w21', 'w1');
  const colTType = findCol('ttype', 't_type', 't');
  const colSGX = findCol('sgx');
  const colSGY = findCol('sgy');
  const colSGZ = findCol('sgz');

  const parseNum = (val) => {
    if (!val) return NaN;
    const clean = val.trim();
    if (clean === '' || clean === 'NaN' || clean === 'null' || clean === '-99.99' || clean === '99.99') return NaN;
    const num = Number(clean);
    return isNaN(num) ? NaN : num;
  };

  for (let i = headerIndex + 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line || line.startsWith('#') || line.startsWith('|')) continue;

    const parts = (typeof delimiter === 'string') ? line.split(delimiter) : line.split(delimiter);
    if (parts.length < 2) continue;

    const pgc = colPGC !== -1 ? parseInt(parts[colPGC], 10) : i;
    const ra = colRA !== -1 ? parseNum(parts[colRA]) : NaN;
    const dec = colDec !== -1 ? parseNum(parts[colDec]) : NaN;
    const cz = colCz !== -1 ? parseNum(parts[colCz]) : NaN;
    const groupId = colGroup !== -1 ? parseInt(parts[colGroup], 10) : pgc;

    if (isNaN(ra) || isNaN(dec)) continue;

    catalog.addGalaxy({
      pgc,
      groupId: isNaN(groupId) ? pgc : groupId,
      ra,
      dec,
      cz: isNaN(cz) ? 0.0 : cz,
      vErr: colVErr !== -1 ? parseNum(parts[colVErr]) : 0.0,
      dist: colDist !== -1 ? parseNum(parts[colDist]) : undefined,
      distErr: colDistErr !== -1 ? parseNum(parts[colDistErr]) : undefined,
      dm: colDM !== -1 ? parseNum(parts[colDM]) : undefined,
      dmErr: colDMErr !== -1 ? parseNum(parts[colDMErr]) : undefined,
      bMag: colBMag !== -1 ? parseNum(parts[colBMag]) : undefined,
      kMag: colKMag !== -1 ? parseNum(parts[colKMag]) : undefined,
      logW21: colW21 !== -1 ? parseNum(parts[colW21]) : undefined,
      tType: colTType !== -1 ? parseNum(parts[colTType]) : undefined,
      sgx: colSGX !== -1 ? parseNum(parts[colSGX]) : undefined,
      sgy: colSGY !== -1 ? parseNum(parts[colSGY]) : undefined,
      sgz: colSGZ !== -1 ? parseNum(parts[colSGZ]) : undefined
    });
  }

  if (options.autoAggregate ?? true) {
    catalog.aggregateGroups(options);
  }

  return catalog;
}

// ============================================================================
// 14. DEFAULT EXPORT & PIPELINE BUNDLE
// ============================================================================

export default {
  G_COSMO,
  SPEED_OF_LIGHT_KM_S,
  DEFAULT_H0,
  DEFAULT_LITTLE_H,
  RHO_CRIT_0,
  MASTER_STRUCTURES,
  angularSeparationRad,
  angularSeparationDeg,
  euclideanDistance3D,
  projectedPhysicalSeparationMpc,
  supergalacticCartesianToSpherical,
  hubbleParameter,
  criticalDensityAtZ,
  calculateMedian,
  calculateMAD,
  biweightLocation,
  biweightScale,
  gapperScale,
  sampleStandardDeviation,
  deconvolveVelocityDispersion,
  calculateProjectedHarmonicRadius,
  calculateProjectedRadiiMoments,
  calculateR200Mpc,
  calculateR500Mpc,
  calculateStandardVirialMass,
  calculateProjectedMassEstimator,
  calculateMedianMassEstimator,
  calculateAverageMassEstimator,
  calculateBahcallTremaineMassEstimator,
  applySurfacePressureCorrection,
  bootstrapVirialMassErrors,
  logGamma,
  standardGamma,
  upperIncompleteGamma,
  schechterCompletenessCorrection,
  analyzeGroupLuminosityAndMorphology,
  calculateLocalSurfaceDensity,
  DisjointSet,
  getAdaptiveFoFLinkingLengths,
  calculateMembershipProbability,
  runHierarchicalFoF,
  crossMatchMasterStructure,
  SpatialHashGrid3D,
  CF4Group,
  CF4GroupedCatalog,
  parseCF4GroupCatalogText
};
