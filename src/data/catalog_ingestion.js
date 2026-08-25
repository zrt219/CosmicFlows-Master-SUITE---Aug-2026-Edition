/**
 * @file catalog_ingestion.js
 * @description Production-grade Astrometric Catalog Ingestion Pipeline for the full 56,000+ galaxy
 * Cosmicflows-4 (CF4) distance compilation and 31,561 grouped catalog.
 * 
 * Implements:
 * 1. Binary & ASCII Parsers:
 *    - FITS Binary Table parser (XTENSION='BINTABLE', big-endian decoding for 'E', 'D', 'J', 'I', 'B', 'K', 'A', 'L').
 *    - IP2I custom ASCII table parser with header metadata detection and sentinel handling.
 *    - CDS VizieR catalogue parser (TSV / byte-column formats with metadata extraction).
 * 2. Cosmographic & Distance Modulus Relations:
 *    - Distance modulus inversion: mu = 5*log10(d) + 25 <-> d = 10^((mu - 25)/5) Mpc.
 *    - Rigorous first-order error propagation: sigma_d = d * (ln(10)/5) * sigma_mu.
 *    - Recessional velocity conversions (cz_cmb, cz_helio, cz_lg) and Hubble relations.
 * 3. Peculiar Velocity & Tully-Fisher Inversion:
 *    - Linear peculiar velocity: v_pec = v_rec - H0 * d.
 *    - Relativistic/cosmographic expansion correction: v_pec = (v_rec - H0*d) / (1 + H0*d/c).
 *    - Tully-Fisher absolute magnitude & distance modulus estimation: M = a*(logW - 2.5) + b.
 * 4. Malmquist Bias Corrections:
 *    - Homogeneous Malmquist correction: d_corr = d_obs * exp(-1.5 * (ln(10)/5)^2 * sigma_mu^2).
 *    - Inhomogeneous Malmquist bias (IMB): Line-of-sight density integration with posterior expectation.
 * 5. Selection Functions & Effective Volume Weights:
 *    - Sigmoid/Fermi-Dirac apparent magnitude selection function.
 *    - Schechter luminosity function distance completeness integration.
 * 6. Grouped Catalog Aggregation:
 *    - Inverse-variance weighted group distance modulus, velocity dispersion, centroids, multiplicity.
 * 7. GalaxyCatalog Columnar Store:
 *    - High-performance TypedArray storage with fast slicing, filtering, and GeoJSON export.
 * 
 * @module data/catalog_ingestion
 */

import {
  equatorialToGalactic,
  galacticToSupergalactic,
  equatorialToSupergalacticCartesian,
  sphericalToCartesian,
  cartesianToSpherical,
  ASTROMETRIC_CONSTANTS
} from '../coordinates/canonical_frame.js';

// Physical and Cosmological Constants
export const SPEED_OF_LIGHT_KM_S = 299792.458;
export const DEFAULT_H0 = ASTROMETRIC_CONSTANTS.DEFAULT_H0 || 74.6;
export const LN10_DIV_5 = Math.LN10 / 5.0; // ~0.4605170185988092
export const HOMOGENEOUS_MALMQUIST_FACTOR = 1.5 * LN10_DIV_5 * LN10_DIV_5; // ~0.3181138243
export const HOMOGENEOUS_MALMQUIST_DM_FACTOR = 3.0 * Math.LN10 / 25.0; // ~0.2763102112

/**
 * Standard CF4 Distance Measurement Methodology Codes.
 * @readonly
 * @enum {string}
 */
export const DistanceMethod = Object.freeze({
  TF: 'TF',           // Tully-Fisher relation (HI / Optical / Infrared)
  FP: 'FP',           // Fundamental Plane (Early-type galaxies)
  SNIA: 'SNIa',       // Type Ia Supernovae
  SBF: 'SBF',         // Surface Brightness Fluctuations
  TRGB: 'TRGB',       // Tip of the Red Giant Branch
  CEPHEID: 'Cepheid', // Classical Cepheid Period-Luminosity
  OTHER: 'Other'
});

// ============================================================================
// 1. Distance Modulus & Peculiar Velocity Mathematical Formulations
// ============================================================================

/**
 * Converts distance modulus mu to physical distance in Mpc (or Mpc/h).
 * Formula: d = 10^((mu - 25.0) / 5.0)
 * 
 * @param {number} mu Distance modulus (mag), e.g. 31.0
 * @param {number} [H0=74.6] Hubble constant in km/s/Mpc
 * @param {boolean} [inMpch=false] If true, returns distance in h^-1 Mpc (where h = H0 / 100)
 * @returns {number} Distance in Mpc (or h^-1 Mpc)
 */
export function distanceModulusToDistance(mu, H0 = DEFAULT_H0, inMpch = false) {
  if (typeof mu !== 'number' || isNaN(mu)) {
    return NaN;
  }
  const dMpc = Math.pow(10.0, (mu - 25.0) / 5.0);
  if (!inMpch) {
    return dMpc;
  }
  const h = H0 / 100.0;
  return dMpc * h;
}

/**
 * Converts physical distance in Mpc to distance modulus mu.
 * Formula: mu = 5.0 * log10(d_Mpc) + 25.0
 * 
 * @param {number} distMpc Distance in Mpc (> 0)
 * @returns {number} Distance modulus in mag
 */
export function distanceToDistanceModulus(distMpc) {
  if (typeof distMpc !== 'number' || isNaN(distMpc) || distMpc <= 0.0) {
    return NaN;
  }
  return 5.0 * Math.log10(distMpc) + 25.0;
}

/**
 * Computes Gaussian distance error sigma_d from distance modulus error sigma_mu.
 * Analytic derivative: sigma_d = d * (ln(10) / 5) * sigma_mu
 * 
 * @param {number} distMpc Distance in Mpc
 * @param {number} sigmaMu Distance modulus uncertainty in mag
 * @returns {number} Distance uncertainty in Mpc
 */
export function distanceModulusUncertaintyToDistance(distMpc, sigmaMu) {
  if (isNaN(distMpc) || isNaN(sigmaMu) || distMpc <= 0.0 || sigmaMu < 0.0) {
    return NaN;
  }
  return distMpc * LN10_DIV_5 * sigmaMu;
}

/**
 * Converts distance uncertainty sigma_d to distance modulus uncertainty sigma_mu.
 * Formula: sigma_mu = sigma_d / (d * (ln(10) / 5))
 * 
 * @param {number} distMpc Distance in Mpc
 * @param {number} sigmaDist Distance uncertainty in Mpc
 * @returns {number} Distance modulus uncertainty in mag
 */
export function distanceUncertaintyToDistanceModulus(distMpc, sigmaDist) {
  if (isNaN(distMpc) || isNaN(sigmaDist) || distMpc <= 0.0 || sigmaDist < 0.0) {
    return NaN;
  }
  return sigmaDist / (distMpc * LN10_DIV_5);
}

/**
 * Calculates peculiar velocity v_pec along the line of sight.
 * Linear form: v_pec = v_rec - H0 * d
 * Relativistic form: v_pec = (v_rec - H0 * d) / (1 + H0 * d / c)
 * 
 * @param {number} vRec Recessional velocity in km/s (e.g. CMB frame)
 * @param {number} distMpc Physical distance in Mpc
 * @param {number} [H0=74.6] Hubble constant in km/s/Mpc
 * @param {Object} [options] Calculation options
 * @param {boolean} [options.relativistic=false] Whether to apply first-order relativistic expansion correction
 * @returns {number} Peculiar velocity in km/s
 */
export function calculatePeculiarVelocity(vRec, distMpc, H0 = DEFAULT_H0, options = {}) {
  if (isNaN(vRec) || isNaN(distMpc) || isNaN(H0) || distMpc < 0.0) {
    return NaN;
  }
  const vHubble = H0 * distMpc;
  const vLinear = vRec - vHubble;
  if (!options.relativistic) {
    return vLinear;
  }
  const denominator = 1.0 + vHubble / SPEED_OF_LIGHT_KM_S;
  return vLinear / denominator;
}

/**
 * Computes formal peculiar velocity uncertainty including distance error,
 * redshift measurement error, and intrinsic cosmic thermal dispersion.
 * Formula: sigma_v = sqrt(sigma_vrec^2 + (H0 * sigma_d)^2 + sigma_thermal^2)
 * 
 * @param {number} distMpc Distance in Mpc
 * @param {number} sigmaMu Distance modulus uncertainty in mag
 * @param {number} [sigmaVrec=0.0] Recessional velocity uncertainty in km/s
 * @param {number} [H0=74.6] Hubble constant in km/s/Mpc
 * @param {number} [thermalDispersion=150.0] Thermal cosmic velocity dispersion in km/s
 * @returns {number} Peculiar velocity uncertainty in km/s
 */
export function peculiarVelocityError(
  distMpc,
  sigmaMu,
  sigmaVrec = 0.0,
  H0 = DEFAULT_H0,
  thermalDispersion = 150.0
) {
  if (isNaN(distMpc) || isNaN(sigmaMu) || distMpc <= 0.0) {
    return NaN;
  }
  const sigmaD = distanceModulusUncertaintyToDistance(distMpc, sigmaMu);
  const sigmaHubble = H0 * sigmaD;
  const varTotal = (sigmaVrec * sigmaVrec) + (sigmaHubble * sigmaHubble) + (thermalDispersion * thermalDispersion);
  return Math.sqrt(varTotal);
}

/**
 * Tully-Fisher relation absolute magnitude estimator.
 * Formula: M = M0 + a * (log10(W21) - logW0)
 * 
 * @param {number} logW21 Logarithm of HI 21cm line width log10(W_21 in km/s)
 * @param {number} [slope=-8.0] TF slope (e.g. -8.0 in W1 band, -7.5 in I-band)
 * @param {number} [zeroPoint=-20.5] Absolute magnitude at logW0
 * @param {number} [logW0=2.5] Reference log linewidth (log10(316.2 km/s))
 * @returns {number} Inferred absolute magnitude M
 */
export function tullyFisherAbsoluteMagnitude(logW21, slope = -8.0, zeroPoint = -20.5, logW0 = 2.5) {
  if (isNaN(logW21)) {
    return NaN;
  }
  return zeroPoint + slope * (logW21 - logW0);
}

// ============================================================================
// 2. Malmquist Bias Correction Models
// ============================================================================

/**
 * Applies Homogeneous Malmquist Bias correction for a uniform spatial distribution.
 * 
 * For distance:
 *   d_true = d_obs * exp( -1.5 * (ln(10)/5)^2 * sigma_mu^2 )
 *   d_true = d_obs * exp( -0.318114 * sigma_mu^2 )
 * 
 * For distance modulus:
 *   mu_true = mu_obs - 3.0 * ln(10) / 25.0 * sigma_mu^2
 *   mu_true = mu_obs - 0.276310 * sigma_mu^2
 * 
 * @param {number} value Distance in Mpc (if mode='distance') or distance modulus in mag (if mode='modulus')
 * @param {number} sigmaMu Distance modulus uncertainty in mag
 * @param {Object} [options] Configuration options
 * @param {'distance'|'modulus'|'forward'} [options.mode='distance'] Calculation mode
 * @returns {number} Bias-corrected quantity
 */
export function correctHomogeneousMalmquist(value, sigmaMu, options = {}) {
  const mode = options.mode || 'distance';
  if (isNaN(value) || isNaN(sigmaMu) || sigmaMu < 0.0) {
    return NaN;
  }
  const sigmaMuSq = sigmaMu * sigmaMu;

  if (mode === 'distance') {
    if (value <= 0.0) return NaN;
    const factor = Math.exp(-HOMOGENEOUS_MALMQUIST_FACTOR * sigmaMuSq);
    return value * factor;
  } else if (mode === 'modulus') {
    const deltaMu = HOMOGENEOUS_MALMQUIST_DM_FACTOR * sigmaMuSq;
    return value - deltaMu;
  } else if (mode === 'forward') {
    // Expected observed modulus given true modulus
    const deltaMu = HOMOGENEOUS_MALMQUIST_DM_FACTOR * sigmaMuSq;
    return value + deltaMu;
  }
  throw new Error(`correctHomogeneousMalmquist: Unknown mode '${mode}'. Use 'distance', 'modulus', or 'forward'.`);
}

/**
 * Applies Inhomogeneous Malmquist Bias (IMB) correction incorporating local cosmological
 * density fluctuations delta(x) along the galaxy line of sight.
 * 
 * Evaluates posterior expectation:
 * E[d | mu_obs, hat{n}] = \frac{\int d * d^2 (1 + delta(d hat{n})) * N(mu_obs | 5*log10(d)+25, sigma_mu^2) dd}{\int d^2 (1 + delta(d hat{n})) * N(mu_obs | 5*log10(d)+25, sigma_mu^2) dd}
 * 
 * @param {number} sgx Supergalactic Cartesian X unit vector / position
 * @param {number} sgy Supergalactic Cartesian Y unit vector / position
 * @param {number} sgz Supergalactic Cartesian Z unit vector / position
 * @param {number} muObs Observed distance modulus (mag)
 * @param {number} sigmaMu Distance modulus uncertainty (mag)
 * @param {Function} [densityInterpolator=null] Density field sampling function delta(x, y, z)
 * @param {Object} [options] Numerical quadrature options
 * @param {number} [options.numSamples=64] Number of Simpson's quadrature integration steps
 * @param {number} [options.sigmaSpan=4.0] Number of standard deviations to integrate across
 * @returns {{ correctedDistance: number, correctedModulus: number, biasCorrection: number, posteriorStd: number }}
 */
export function correctInhomogeneousMalmquist(
  sgx,
  sgy,
  sgz,
  muObs,
  sigmaMu,
  densityInterpolator = null,
  options = {}
) {
  if (isNaN(muObs) || isNaN(sigmaMu) || sigmaMu <= 0.0) {
    return {
      correctedDistance: NaN,
      correctedModulus: NaN,
      biasCorrection: NaN,
      posteriorStd: NaN
    };
  }

  const dObs = distanceModulusToDistance(muObs);
  const sigmaD = distanceModulusUncertaintyToDistance(dObs, sigmaMu);

  // Direction unit vector
  const norm = Math.hypot(sgx, sgy, sgz);
  const nx = norm > 1e-12 ? sgx / norm : 1.0;
  const ny = norm > 1e-12 ? sgy / norm : 0.0;
  const nz = norm > 1e-12 ? sgz / norm : 0.0;

  const nSteps = Math.max(16, (options.numSamples || 64) | 1); // Ensure odd count for Simpson's rule
  const nSigma = options.sigmaSpan || 4.0;
  const dMin = Math.max(0.1, dObs - nSigma * sigmaD);
  const dMax = dObs + nSigma * sigmaD;
  const deltaD = (dMax - dMin) / (nSteps - 1);

  let numIntegral = 0.0;
  let denIntegral = 0.0;
  let varIntegral = 0.0;

  const twoSigmaMuSq = 2.0 * sigmaMu * sigmaMu;

  for (let i = 0; i < nSteps; i++) {
    const d = dMin + i * deltaD;
    const simpsonWeight = (i === 0 || i === nSteps - 1) ? 1.0 : (i % 2 === 1 ? 4.0 : 2.0);

    const px = d * nx;
    const py = d * ny;
    const pz = d * nz;

    let delta = 0.0;
    if (typeof densityInterpolator === 'function') {
      try {
        delta = Math.max(-0.999, densityInterpolator(px, py, pz));
      } catch {
        delta = 0.0;
      }
    }

    const muModel = 5.0 * Math.log10(d) + 25.0;
    const diffMu = muObs - muModel;
    const likelihood = Math.exp(-(diffMu * diffMu) / twoSigmaMuSq);

    // Volume element * (1 + delta) * Gaussian likelihood
    const integrand = (d * d) * (1.0 + delta) * likelihood * simpsonWeight;

    denIntegral += integrand;
    numIntegral += d * integrand;
    varIntegral += (d * d) * integrand;
  }

  if (denIntegral <= 0.0 || isNaN(denIntegral)) {
    const fallbackD = correctHomogeneousMalmquist(dObs, sigmaMu, { mode: 'distance' });
    return {
      correctedDistance: fallbackD,
      correctedModulus: distanceToDistanceModulus(fallbackD),
      biasCorrection: fallbackD - dObs,
      posteriorStd: sigmaD
    };
  }

  const dMean = numIntegral / denIntegral;
  const dVar = Math.max(0.0, (varIntegral / denIntegral) - (dMean * dMean));
  const postStd = Math.sqrt(dVar);
  const muMean = distanceToDistanceModulus(dMean);

  return {
    correctedDistance: dMean,
    correctedModulus: muMean,
    biasCorrection: dMean - dObs,
    posteriorStd: postStd
  };
}

// ============================================================================
// 3. Selection Functions & Survey Weights
// ============================================================================

/**
 * Evaluates smooth Fermi-Dirac / sigmoid apparent magnitude selection probability.
 * S(m) = 1 / (1 + exp((m - m_lim) / delta_m))
 * 
 * @param {number} appMag Apparent magnitude m
 * @param {number} [magLim=17.5] 50% completeness magnitude threshold
 * @param {number} [deltaM=0.4] Soft transition width in mag
 * @returns {number} Selection probability in [0, 1]
 */
export function evaluateMagnitudeSelection(appMag, magLim = 17.5, deltaM = 0.4) {
  if (isNaN(appMag)) {
    return 0.0;
  }
  const arg = (appMag - magLim) / Math.max(1e-4, deltaM);
  if (arg > 50.0) return 0.0;
  if (arg < -50.0) return 1.0;
  return 1.0 / (1.0 + Math.exp(arg));
}

/**
 * Incomplete upper Gamma function Gamma(a, x) = int_x^inf t^(a-1) e^-t dt
 * via continued fractions / series expansion for Schechter luminosity integration.
 * 
 * @param {number} a Shape parameter (alpha + 1)
 * @param {number} x Lower integration bound (x >= 0)
 * @returns {number} Upper incomplete gamma value
 */
export function upperIncompleteGamma(a, x) {
  if (x <= 0.0) {
    return gammaLanczos(a);
  }
  if (x > 100.0) {
    return 0.0;
  }

  const maxIter = 100;
  const eps = 1e-12;
  let b = x + 1.0 - a;
  let c = 1.0 / 1e-30;
  let d = 1.0 / b;
  let h = d;

  for (let i = 1; i <= maxIter; i++) {
    const an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (Math.abs(d) < 1e-30) d = 1e-30;
    c = b + an / c;
    if (Math.abs(c) < 1e-30) c = 1e-30;
    d = 1.0 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1.0) < eps) break;
  }

  const factor = Math.exp(-x + a * Math.log(x));
  return factor * h;
}

/**
 * Lanczos approximation for standard Gamma(z).
 * @param {number} z Real number (z > 0)
 * @returns {number} Gamma(z)
 */
function gammaLanczos(z) {
  if (z < 0.5) {
    return Math.PI / (Math.sin(Math.PI * z) * gammaLanczos(1.0 - z));
  }
  z -= 1.0;
  const g = 7;
  const C = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109583111912,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  ];
  let x = C[0];
  for (let i = 1; i < g + 2; i++) {
    x += C[i] / (z + i);
  }
  const t = z + g + 0.5;
  return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
}

/**
 * Evaluates distance selection completeness Phi(d) for a flux-limited survey
 * following a Schechter luminosity function.
 * 
 * @param {number} distMpc Distance in Mpc
 * @param {number} [magLim=17.5] Apparent magnitude limit
 * @param {number} [mStar=-20.5] Characteristic absolute magnitude M*
 * @param {number} [alpha=-1.0] Schechter faint-end slope
 * @returns {number} Completeness fraction in [0, 1]
 */
export function evaluateSchechterDistanceSelection(
  distMpc,
  magLim = 17.5,
  mStar = -20.5,
  alpha = -1.0
) {
  if (isNaN(distMpc) || distMpc <= 0.0) {
    return 0.0;
  }
  const mLimAbs = magLim - 5.0 * Math.log10(distMpc) - 25.0;
  const lRatio = Math.pow(10.0, -0.4 * (mLimAbs - mStar));
  const a = alpha + 1.0;

  if (a <= 0.0) {
    return Math.exp(-lRatio);
  }

  const gammaUpper = upperIncompleteGamma(a, lRatio);
  const gammaTotal = gammaLanczos(a);
  const fraction = gammaUpper / gammaTotal;
  return Math.min(1.0, Math.max(0.0, fraction));
}

// ============================================================================
// 4. Grouped Catalog Aggregation Engine
// ============================================================================

/**
 * Aggregates individual galaxy measurements into a grouped CF4 cluster catalog.
 * 
 * @param {Array<Object>} galaxies Array of galaxy objects
 * @param {Object} [options] Aggregation configuration
 * @param {string} [options.groupKey='groupId'] Object property identifying group membership
 * @param {number} [options.H0=74.6] Hubble constant in km/s/Mpc
 * @returns {Array<Object>} Array of grouped cluster objects
 */
export function aggregateGroupedCatalog(galaxies, options = {}) {
  if (!Array.isArray(galaxies) || galaxies.length === 0) {
    return [];
  }

  const groupKey = options.groupKey || 'groupId';
  const H0 = options.H0 || DEFAULT_H0;

  // Bucket galaxies by group ID
  const groupMap = new Map();
  let ungroupedCounter = 1;

  for (let i = 0; i < galaxies.length; i++) {
    const gal = galaxies[i];
    let gid = gal[groupKey];
    if (gid === undefined || gid === null || gid === '' || gid <= 0) {
      gid = `single_${ungroupedCounter++}`;
    } else {
      gid = String(gid);
    }

    if (!groupMap.has(gid)) {
      groupMap.set(gid, []);
    }
    groupMap.get(gid).push(gal);
  }

  const aggregated = [];

  for (const [gid, members] of groupMap.entries()) {
    const nMem = members.length;

    let sumWeightMu = 0.0;
    let sumWeight = 0.0;
    let sumVrec = 0.0;
    let sumVrecWeight = 0.0;

    let sumX = 0.0;
    let sumY = 0.0;
    let sumZ = 0.0;
    let sumCoordWeight = 0.0;

    const methodsSet = new Set();

    for (let j = 0; j < nMem; j++) {
      const g = members[j];
      const mu = typeof g.mu === 'number' && !isNaN(g.mu) ? g.mu : g.mod;
      const dmu = typeof g.dmu === 'number' && !isNaN(g.dmu) ? g.dmu : (g.dmod || 0.35);
      const vrec = typeof g.cz === 'number' && !isNaN(g.cz) ? g.cz : (g.vrec || g.vcmb || 0.0);
      if (g.method) methodsSet.add(g.method);

      const wMu = 1.0 / Math.max(1e-4, dmu * dmu);
      if (!isNaN(mu)) {
        sumWeightMu += mu * wMu;
        sumWeight += wMu;
      }

      const wV = 1.0;
      sumVrec += vrec * wV;
      sumVrecWeight += wV;

      // Coordinate aggregation
      const ra = g.ra || 0.0;
      const dec = g.dec || 0.0;
      const [x, y, z] = sphericalToCartesian(ra, dec, 1.0);
      sumX += x;
      sumY += y;
      sumZ += z;
      sumCoordWeight += 1.0;
    }

    const meanMu = sumWeight > 0.0 ? sumWeightMu / sumWeight : NaN;
    const errMu = sumWeight > 0.0 ? 1.0 / Math.sqrt(sumWeight) : 0.4;
    const meanVrec = sumVrecWeight > 0.0 ? sumVrec / sumVrecWeight : 0.0;

    // Velocity dispersion
    let velDisp = 0.0;
    if (nMem > 1) {
      let sumSqDiff = 0.0;
      for (let j = 0; j < nMem; j++) {
        const g = members[j];
        const vrec = typeof g.cz === 'number' && !isNaN(g.cz) ? g.cz : (g.vrec || g.vcmb || 0.0);
        const diff = vrec - meanVrec;
        sumSqDiff += diff * diff;
      }
      velDisp = Math.sqrt(sumSqDiff / (nMem - 1));
    }

    // Mean Sky coordinates
    const meanX = sumX / sumCoordWeight;
    const meanY = sumY / sumCoordWeight;
    const meanZ = sumZ / sumCoordWeight;
    const [meanRa, meanDec] = cartesianToSpherical(meanX, meanY, meanZ);

    const [meanL, meanB] = equatorialToGalactic(meanRa, meanDec);
    const [meanSgl, meanSgb] = galacticToSupergalactic(meanL, meanB);

    const distMpc = distanceModulusToDistance(meanMu, H0);
    const distErrMpc = distanceModulusUncertaintyToDistance(distMpc, errMu);
    const [sgx, sgy, sgz] = equatorialToSupergalacticCartesian(meanRa, meanDec, distMpc);
    const vPec = calculatePeculiarVelocity(meanVrec, distMpc, H0);
    const vPecErr = peculiarVelocityError(distMpc, errMu, velDisp / Math.sqrt(nMem), H0);

    aggregated.push({
      groupId: gid.startsWith('single_') ? null : (isNaN(Number(gid)) ? gid : Number(gid)),
      groupLabel: gid,
      multiplicity: nMem,
      ra: meanRa,
      dec: meanDec,
      glon: meanL,
      glat: meanB,
      sgl: meanSgl,
      sgb: meanSgb,
      sgx: sgx,
      sgy: sgy,
      sgz: sgz,
      cz: meanVrec,
      velDispersion: velDisp,
      mu: meanMu,
      dmu: errMu,
      distMpc: distMpc,
      distErrMpc: distErrMpc,
      vpec: vPec,
      vpecErr: vPecErr,
      methods: Array.from(methodsSet),
      members: members
    });
  }

  return aggregated;
}

// ============================================================================
// 5. Binary & ASCII File Format Parsers
// ============================================================================

/**
 * Decodes 80-character FITS header cards into key-value map.
 * @param {Uint8Array} headerBytes Raw header bytes (multiple of 2880)
 * @returns {Map<string, string|number|boolean>} Map of FITS keywords
 */
export function decodeFITSHeaderCards(headerBytes) {
  const cards = new Map();
  const nCards = Math.floor(headerBytes.length / 80);
  const decoder = new TextDecoder('ascii');

  for (let i = 0; i < nCards; i++) {
    const cardStr = decoder.decode(headerBytes.subarray(i * 80, (i + 1) * 80));
    const eqIdx = cardStr.indexOf('=');
    if (eqIdx === -1) {
      const key = cardStr.substring(0, 8).trim();
      if (key === 'END') {
        cards.set('END', true);
        break;
      }
      continue;
    }

    const key = cardStr.substring(0, eqIdx).trim();
    let valComment = cardStr.substring(eqIdx + 1);
    const slashIdx = valComment.indexOf('/');
    let valStr = (slashIdx !== -1 ? valComment.substring(0, slashIdx) : valComment).trim();

    if (valStr.startsWith("'")) {
      const endQuote = valStr.indexOf("'", 1);
      const strVal = endQuote !== -1 ? valStr.substring(1, endQuote).trim() : valStr.substring(1).trim();
      cards.set(key, strVal);
    } else if (valStr === 'T') {
      cards.set(key, true);
    } else if (valStr === 'F') {
      cards.set(key, false);
    } else {
      const num = Number(valStr);
      cards.set(key, isNaN(num) ? valStr : num);
    }
  }

  return cards;
}

/**
 * Parses FITS Binary Tables (BINTABLE extension) adhering to IAU FITS Standard 4.0.
 * 
 * @param {ArrayBuffer|Uint8Array} buffer Input FITS file buffer
 * @param {Object} [options] Parser options
 * @param {number} [options.extensionIndex=1] Target BINTABLE extension index (1 = first extension)
 * @returns {{ headers: Map<string, any>, columnNames: string[], rowCount: number, columns: Object, rows: Object[] }}
 */
export function parseFITSBinaryTable(buffer, options = {}) {
  const bytes = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  const dataView = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const decoder = new TextDecoder('ascii');

  let offset = 0;
  let extCount = 0;
  const targetExt = options.extensionIndex !== undefined ? options.extensionIndex : 1;

  while (offset < bytes.byteLength) {
    const headerStart = offset;
    let headerEnd = headerStart;
    let foundEnd = false;

    // Scan for 2880-byte block containing 'END     '
    while (headerEnd < bytes.byteLength) {
      const block = bytes.subarray(headerEnd, headerEnd + 2880);
      headerEnd += 2880;
      for (let i = 0; i < 2880; i += 80) {
        const key = decoder.decode(block.subarray(i, i + 8)).trim();
        if (key === 'END') {
          foundEnd = true;
          break;
        }
      }
      if (foundEnd) break;
    }

    const headerBytes = bytes.subarray(headerStart, headerEnd);
    const header = decodeFITSHeaderCards(headerBytes);
    offset = headerEnd;

    const xtension = header.get('XTENSION');
    const naxis1 = Number(header.get('NAXIS1') || 0);
    const naxis2 = Number(header.get('NAXIS2') || 0);
    const pcount = Number(header.get('PCOUNT') || 0);
    const gcount = Number(header.get('GCOUNT') || 1);

    const dataSize = Math.ceil((naxis1 * naxis2 * gcount + pcount) / 2880) * 2880;

    if (xtension === 'BINTABLE' || xtension === 'TABLE') {
      extCount++;
      if (extCount === targetExt || targetExt === 0) {
        const tfields = Number(header.get('TFIELDS') || 0);
        const colDefs = [];

        let currentByteOffset = 0;
        for (let f = 1; f <= tfields; f++) {
          const ttype = String(header.get(`TTYPE${f}`) || `COL${f}`).trim();
          const tform = String(header.get(`TFORM${f}`) || '1E').trim();
          const tscal = Number(header.get(`TSCAL${f}`) !== undefined ? header.get(`TSCAL${f}`) : 1.0);
          const tzero = Number(header.get(`TZERO${f}`) !== undefined ? header.get(`TZERO${f}`) : 0.0);
          const tnull = header.get(`TNULL${f}`);

          const match = tform.match(/^(\d*)([A-Z])$/);
          const repeat = match && match[1] ? parseInt(match[1], 10) : 1;
          const typeCode = match ? match[2] : 'E';

          let byteLength = 0;
          switch (typeCode) {
            case 'A': byteLength = repeat; break;
            case 'L': byteLength = repeat; break;
            case 'B': byteLength = repeat; break;
            case 'I': byteLength = repeat * 2; break;
            case 'J': byteLength = repeat * 4; break;
            case 'K': byteLength = repeat * 8; break;
            case 'E': byteLength = repeat * 4; break;
            case 'D': byteLength = repeat * 8; break;
            default: byteLength = repeat * 4; break;
          }

          colDefs.push({
            name: ttype,
            form: tform,
            type: typeCode,
            repeat: repeat,
            scale: tscal,
            zero: tzero,
            nullVal: tnull,
            colOffset: currentByteOffset,
            byteLength: byteLength
          });

          currentByteOffset += byteLength;
        }

        const dataStartOffset = offset;
        const rowCount = naxis2;
        const columns = {};
        for (const col of colDefs) {
          if (col.type === 'E' || col.type === 'D') {
            columns[col.name] = new Float64Array(rowCount);
          } else if (col.type === 'J' || col.type === 'I' || col.type === 'B') {
            columns[col.name] = new Int32Array(rowCount);
          } else if (col.type === 'K') {
            columns[col.name] = new Float64Array(rowCount);
          } else {
            columns[col.name] = new Array(rowCount);
          }
        }

        const rows = new Array(rowCount);

        for (let r = 0; r < rowCount; r++) {
          const rowBytePos = dataStartOffset + r * naxis1;
          const rowObj = {};

          for (const col of colDefs) {
            const fieldPos = rowBytePos + col.colOffset;
            let val;

            if (fieldPos + col.byteLength > bytes.byteLength) {
              val = NaN;
            } else {
              switch (col.type) {
                case 'E': {
                  const raw = dataView.getFloat32(fieldPos, false);
                  val = isNaN(raw) ? NaN : raw * col.scale + col.zero;
                  break;
                }
                case 'D': {
                  const raw = dataView.getFloat64(fieldPos, false);
                  val = isNaN(raw) ? NaN : raw * col.scale + col.zero;
                  break;
                }
                case 'J': {
                  const raw = dataView.getInt32(fieldPos, false);
                  val = raw === col.nullVal ? NaN : raw * col.scale + col.zero;
                  break;
                }
                case 'I': {
                  const raw = dataView.getInt16(fieldPos, false);
                  val = raw === col.nullVal ? NaN : raw * col.scale + col.zero;
                  break;
                }
                case 'B': {
                  const raw = dataView.getUint8(fieldPos);
                  val = raw === col.nullVal ? NaN : raw * col.scale + col.zero;
                  break;
                }
                case 'K': {
                  const rawBig = dataView.getBigInt64(fieldPos, false);
                  val = Number(rawBig) * col.scale + col.zero;
                  break;
                }
                case 'A': {
                  const strBytes = bytes.subarray(fieldPos, fieldPos + col.byteLength);
                  val = decoder.decode(strBytes).trim();
                  break;
                }
                case 'L': {
                  const charCode = bytes[fieldPos];
                  val = (charCode === 84 || charCode === 116);
                  break;
                }
                default:
                  val = dataView.getFloat32(fieldPos, false);
              }
            }

            columns[col.name][r] = val;
            rowObj[col.name] = val;
          }
          rows[r] = rowObj;
        }

        return {
          headers: header,
          columnNames: colDefs.map(c => c.name),
          rowCount: rowCount,
          columns: columns,
          rows: rows
        };
      }
    }

    offset += dataSize;
  }

  throw new Error(`parseFITSBinaryTable: No valid BINTABLE extension found at index ${targetExt}.`);
}

/**
 * Parses IP2I Lyon / Cosmicflows custom ASCII tables.
 * 
 * @param {string} textData Raw ASCII string content
 * @param {Object} [options] Parser configuration
 * @returns {{ metadata: Object, columnNames: string[], rows: Object[], rowCount: number }}
 */
export function parseIP2IASCIITable(textData, options = {}) {
  if (typeof textData !== 'string') {
    throw new TypeError('parseIP2IASCIITable: textData must be a string.');
  }

  const lines = textData.split(/\r?\n/);
  const metadata = {};
  let headerCols = options.columnNames ? [...options.columnNames] : null;
  const rows = [];

  const SENTINELS = new Set(['-99.99', '-99.9', '-999', '99.99', '99.9', '9999.0', '9999', 'NaN', 'nan', 'null', '...']);

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line) continue;

    if (line.startsWith('#')) {
      const commentContent = line.substring(1).trim();
      const metaMatch = commentContent.match(/^([A-Za-z0-9_]+)\s*[:=]\s*(.+)$/);
      if (metaMatch) {
        const k = metaMatch[1];
        const v = metaMatch[2].trim();
        const num = Number(v);
        metadata[k] = isNaN(num) ? v : num;
      }

      if (!headerCols && (commentContent.toLowerCase().startsWith('columns:') || commentContent.toLowerCase().startsWith('cols:'))) {
        const colListStr = commentContent.substring(commentContent.indexOf(':') + 1).trim();
        headerCols = colListStr.split(/[\s,|]+/).filter(c => c.length > 0);
      }
      continue;
    }

    const tokens = line.split(/[\s,|]+/).filter(t => t.length > 0);
    if (tokens.length === 0) continue;

    if (!headerCols) {
      const isHeader = tokens.some(t => isNaN(Number(t)));
      if (isHeader) {
        headerCols = tokens.map(t => t.trim());
        continue;
      } else {
        headerCols = tokens.map((_, idx) => `col_${idx}`);
      }
    }

    const rowObj = {};
    for (let c = 0; c < headerCols.length; c++) {
      const colName = headerCols[c];
      const token = c < tokens.length ? tokens[c] : '';

      if (SENTINELS.has(token) || token === '') {
        rowObj[colName] = NaN;
      } else {
        const num = Number(token);
        rowObj[colName] = isNaN(num) ? token : num;
      }
    }
    rows.push(rowObj);
  }

  return {
    metadata: metadata,
    columnNames: headerCols || [],
    rows: rows,
    rowCount: rows.length
  };
}

/**
 * Parses CDS VizieR astronomical catalogue formats.
 * 
 * @param {string} textData Raw CDS VizieR text string
 * @param {Object} [options] Parser options
 * @returns {{ metadata: Object, columnNames: string[], rows: Object[], rowCount: number }}
 */
export function parseCDSVizieRCatalog(textData, options = {}) {
  if (typeof textData !== 'string') {
    throw new TypeError('parseCDSVizieRCatalog: textData must be a string.');
  }

  const lines = textData.split(/\r?\n/);
  const metadata = {};
  const rows = [];
  let columnNames = [];
  let dataStarted = false;

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const trimmed = line.trim();
    if (!trimmed) continue;

    if (trimmed.startsWith('#')) {
      const meta = trimmed.substring(1).trim();
      const match = meta.match(/^([A-Za-z0-9_]+)\s*[:=]\s*(.+)$/);
      if (match) {
        metadata[match[1]] = match[2].trim();
      }

      // Check if comment line defines column names
      if (!dataStarted && (meta.includes('|') || meta.toLowerCase().startsWith('columns:') || meta.toLowerCase().startsWith('cols:'))) {
        const colStr = meta.includes(':') ? meta.substring(meta.indexOf(':') + 1) : meta;
        columnNames = colStr.split(/[\t|]/).map(c => c.trim()).filter(c => c.length > 0);
      }
      continue;
    }

    if (trimmed.startsWith('---') || trimmed.includes('---|---')) {
      dataStarted = true;
      continue;
    }

    if (!dataStarted) {
      if (trimmed.includes('\t') || trimmed.includes('|')) {
        columnNames = trimmed.split(/[\t|]/).map(c => c.trim()).filter(c => c.length > 0);
      }
      continue;
    }

    const tokens = (trimmed.includes('\t') ? trimmed.split('\t') : (trimmed.includes('|') ? trimmed.split('|') : trimmed.split(/\s+/)))
      .map(t => t.trim())
      .filter(t => t.length > 0);

    if (tokens.length === 0) continue;

    const rowObj = {};
    for (let c = 0; c < columnNames.length; c++) {
      const col = columnNames[c];
      const tok = c < tokens.length ? tokens[c] : '';
      if (tok === '' || tok === '---' || tok === 'NaN' || tok === 'null') {
        rowObj[col] = NaN;
      } else {
        const num = Number(tok);
        rowObj[col] = isNaN(num) ? tok : num;
      }
    }
    rows.push(rowObj);
  }

  return {
    metadata: metadata,
    columnNames: columnNames,
    rows: rows,
    rowCount: rows.length
  };
}

// ============================================================================
// 6. GalaxyCatalog Class (High-Performance Columnar Storage)
// ============================================================================

/**
 * GalaxyCatalog represents a high-density, columnar-backed astronomical catalog.
 */
export class GalaxyCatalog {
  /**
   * Constructs a GalaxyCatalog instance.
   * @param {number} count Total number of galaxy records
   */
  constructor(count) {
    this.count = count;
    this.pgc = new Int32Array(count);
    this.ra = new Float64Array(count);
    this.dec = new Float64Array(count);
    this.glon = new Float64Array(count);
    this.glat = new Float64Array(count);
    this.sgl = new Float64Array(count);
    this.sgb = new Float64Array(count);
    this.sgx = new Float64Array(count);
    this.sgy = new Float64Array(count);
    this.sgz = new Float64Array(count);
    this.cz = new Float64Array(count);
    this.mu = new Float64Array(count);
    this.dmu = new Float64Array(count);
    this.dist = new Float64Array(count);
    this.distErr = new Float64Array(count);
    this.vpec = new Float64Array(count);
    this.vpecErr = new Float64Array(count);
    this.weight = new Float64Array(count);
    this.groupId = new Int32Array(count);
    this.method = new Array(count);
  }

  /**
   * Retrieves single galaxy record by index as a JavaScript object.
   * @param {number} index Galaxy row index
   * @returns {Object} Galaxy data object
   */
  getGalaxy(index) {
    if (index < 0 || index >= this.count) {
      throw new RangeError(`GalaxyCatalog.getGalaxy: Index ${index} out of range [0, ${this.count - 1}].`);
    }
    return {
      index: index,
      pgc: this.pgc[index],
      ra: this.ra[index],
      dec: this.dec[index],
      glon: this.glon[index],
      glat: this.glat[index],
      sgl: this.sgl[index],
      sgb: this.sgb[index],
      sgx: this.sgx[index],
      sgy: this.sgy[index],
      sgz: this.sgz[index],
      cz: this.cz[index],
      mu: this.mu[index],
      dmu: this.dmu[index],
      dist: this.dist[index],
      distErr: this.distErr[index],
      vpec: this.vpec[index],
      vpecErr: this.vpecErr[index],
      weight: this.weight[index],
      groupId: this.groupId[index],
      method: this.method[index]
    };
  }

  /**
   * Filters the catalog by a predicate function, returning a new sliced GalaxyCatalog.
   * @param {Function} predicate Filter function
   * @returns {GalaxyCatalog} Sliced sub-catalog
   */
  filter(predicate) {
    const matchingIndices = [];
    for (let i = 0; i < this.count; i++) {
      const g = this.getGalaxy(i);
      if (predicate(g, i)) {
        matchingIndices.push(i);
      }
    }

    const sub = new GalaxyCatalog(matchingIndices.length);
    for (let j = 0; j < matchingIndices.length; j++) {
      const src = matchingIndices[j];
      sub.pgc[j] = this.pgc[src];
      sub.ra[j] = this.ra[src];
      sub.dec[j] = this.dec[src];
      sub.glon[j] = this.glon[src];
      sub.glat[j] = this.glat[src];
      sub.sgl[j] = this.sgl[src];
      sub.sgb[j] = this.sgb[src];
      sub.sgx[j] = this.sgx[src];
      sub.sgy[j] = this.sgy[src];
      sub.sgz[j] = this.sgz[src];
      sub.cz[j] = this.cz[src];
      sub.mu[j] = this.mu[src];
      sub.dmu[j] = this.dmu[src];
      sub.dist[j] = this.dist[src];
      sub.distErr[j] = this.distErr[src];
      sub.vpec[j] = this.vpec[src];
      sub.vpecErr[j] = this.vpecErr[src];
      sub.weight[j] = this.weight[src];
      sub.groupId[j] = this.groupId[src];
      sub.method[j] = this.method[src];
    }
    return sub;
  }

  /**
   * Computes bounding box in Supergalactic Cartesian coordinates.
   * @returns {{ xMin: number, xMax: number, yMin: number, yMax: number, zMin: number, zMax: number }}
   */
  getBoundingBox() {
    let xMin = Infinity, xMax = -Infinity;
    let yMin = Infinity, yMax = -Infinity;
    let zMin = Infinity, zMax = -Infinity;

    for (let i = 0; i < this.count; i++) {
      const x = this.sgx[i];
      const y = this.sgy[i];
      const z = this.sgz[i];
      if (!isNaN(x) && x < xMin) xMin = x;
      if (!isNaN(x) && x > xMax) xMax = x;
      if (!isNaN(y) && y < yMin) yMin = y;
      if (!isNaN(y) && y > yMax) yMax = y;
      if (!isNaN(z) && z < zMin) zMin = z;
      if (!isNaN(z) && z > zMax) zMax = z;
    }

    return { xMin, xMax, yMin, yMax, zMin, zMax };
  }

  /**
   * Exports catalog to standard GeoJSON FeatureCollection format.
   * @returns {Object} GeoJSON FeatureCollection
   */
  toGeoJSON() {
    const features = new Array(this.count);
    for (let i = 0; i < this.count; i++) {
      features[i] = {
        type: 'Feature',
        geometry: {
          type: 'Point',
          coordinates: [this.ra[i], this.dec[i], this.dist[i]]
        },
        properties: {
          pgc: this.pgc[i],
          sgl: this.sgl[i],
          sgb: this.sgb[i],
          sgx: this.sgx[i],
          sgy: this.sgy[i],
          sgz: this.sgz[i],
          cz: this.cz[i],
          vpec: this.vpec[i],
          mu: this.mu[i],
          weight: this.weight[i],
          groupId: this.groupId[i],
          method: this.method[i]
        }
      };
    }
    return {
      type: 'FeatureCollection',
      features: features
    };
  }
}

// ============================================================================
// 7. Full Catalog Ingestion Pipeline
// ============================================================================

function resolveRowField(row, aliases, fallback = NaN) {
  for (const alias of aliases) {
    if (row[alias] !== undefined && row[alias] !== null && row[alias] !== '') {
      const val = row[alias];
      const num = Number(val);
      return isNaN(num) ? val : num;
    }
  }
  return fallback;
}

/**
 * CatalogIngestionPipeline provides end-to-end ingestion, validation, coordinate transformation,
 * Malmquist bias correction, and selection function weighting.
 */
export class CatalogIngestionPipeline {
  /**
   * Constructs the ingestion pipeline.
   * @param {Object} [options] Pipeline configuration
   */
  constructor(options = {}) {
    this.H0 = options.H0 || DEFAULT_H0;
    this.applyMalmquist = options.applyMalmquist !== undefined ? options.applyMalmquist : true;
    this.malmquistModel = options.malmquistModel || 'homogeneous';
    this.densityInterpolator = options.densityInterpolator || null;
    this.applySelectionWeights = options.applySelectionWeights !== undefined ? options.applySelectionWeights : true;
    this.magLim = options.magLim || 17.5;
  }

  /**
   * Ingests parsed raw galaxy records into a structured GalaxyCatalog.
   * @param {Array<Object>} rawRows Array of record objects
   * @returns {GalaxyCatalog} Calibrated GalaxyCatalog
   */
  processRecords(rawRows) {
    if (!Array.isArray(rawRows)) {
      throw new TypeError('CatalogIngestionPipeline.processRecords: rawRows must be an array.');
    }

    const n = rawRows.length;
    const cat = new GalaxyCatalog(n);

    for (let i = 0; i < n; i++) {
      const row = rawRows[i];

      const pgc = resolveRowField(row, ['PGC', 'pgc', 'ID', 'id', 'Name', 'name', 'objid'], i + 1);
      cat.pgc[i] = typeof pgc === 'number' ? pgc : parseInt(pgc, 10) || (i + 1);

      const ra = resolveRowField(row, ['RA', 'ra', 'RAJ2000', 'radeg', 'ALPHA', 'alpha'], 0.0);
      const dec = resolveRowField(row, ['DEC', 'dec', 'DEJ2000', 'dedeg', 'DELTA', 'delta'], 0.0);
      cat.ra[i] = ra;
      cat.dec[i] = dec;

      const cz = resolveRowField(row, ['cz', 'CZ', 'Vcmb', 'vcmb', 'V_cmb', 'vrec', 'Vrec', 'Vhel', 'vhel'], 0.0);
      cat.cz[i] = cz;

      let mu = resolveRowField(row, ['DM', 'dm', 'mu', 'MOD', 'mod', 'distmod', 'm-M', 'dist_mod'], NaN);
      let dmu = resolveRowField(row, ['eDM', 'edm', 'dmu', 'e_DM', 'e_mod', 'sig_mu', 'dmod'], 0.35);

      let distDirect = resolveRowField(row, ['dist', 'Dist', 'D', 'd_mpc', 'distance'], NaN);
      if (isNaN(mu) && !isNaN(distDirect) && distDirect > 0.0) {
        mu = distanceToDistanceModulus(distDirect);
      }

      let distMpc = NaN;
      if (!isNaN(mu)) {
        if (this.applyMalmquist) {
          if (this.malmquistModel === 'homogeneous') {
            const dObs = distanceModulusToDistance(mu, this.H0);
            distMpc = correctHomogeneousMalmquist(dObs, dmu, { mode: 'distance' });
            mu = distanceToDistanceModulus(distMpc);
          } else if (this.malmquistModel === 'inhomogeneous' && typeof this.densityInterpolator === 'function') {
            const [sgx0, sgy0, sgz0] = equatorialToSupergalacticCartesian(ra, dec, 1.0);
            const imb = correctInhomogeneousMalmquist(sgx0, sgy0, sgz0, mu, dmu, this.densityInterpolator);
            distMpc = imb.correctedDistance;
            mu = imb.correctedModulus;
          } else {
            distMpc = distanceModulusToDistance(mu, this.H0);
          }
        } else {
          distMpc = distanceModulusToDistance(mu, this.H0);
        }
      } else if (!isNaN(distDirect)) {
        distMpc = distDirect;
      }

      cat.mu[i] = mu;
      cat.dmu[i] = dmu;
      cat.dist[i] = distMpc;
      cat.distErr[i] = distanceModulusUncertaintyToDistance(distMpc, dmu);

      const [glon, glat] = equatorialToGalactic(ra, dec);
      const [sgl, sgb] = galacticToSupergalactic(glon, glat);
      const [sgx, sgy, sgz] = equatorialToSupergalacticCartesian(ra, dec, isNaN(distMpc) ? 1.0 : distMpc);

      cat.glon[i] = glon;
      cat.glat[i] = glat;
      cat.sgl[i] = sgl;
      cat.sgb[i] = sgb;
      cat.sgx[i] = sgx;
      cat.sgy[i] = sgy;
      cat.sgz[i] = sgz;

      const vpec = calculatePeculiarVelocity(cz, distMpc, this.H0);
      const vpecErr = peculiarVelocityError(distMpc, dmu, 0.0, this.H0);
      cat.vpec[i] = vpec;
      cat.vpecErr[i] = vpecErr;

      let w = 1.0;
      if (this.applySelectionWeights && !isNaN(distMpc) && distMpc > 0.0) {
        const phi = evaluateSchechterDistanceSelection(distMpc, this.magLim);
        w = 1.0 / Math.max(0.01, phi);
      }
      cat.weight[i] = w;

      const gid = resolveRowField(row, ['GroupId', 'group_id', 'NestId', 'nest', 'grp', 'group'], -1);
      cat.groupId[i] = typeof gid === 'number' ? gid : (parseInt(gid, 10) || -1);

      const meth = resolveRowField(row, ['Method', 'method', 'METH', 'type'], DistanceMethod.OTHER);
      cat.method[i] = String(meth);
    }

    return cat;
  }

  /**
   * Ingests FITS file ArrayBuffer directly.
   * @param {ArrayBuffer|Uint8Array} buffer FITS binary table buffer
   * @param {Object} [options] FITS options
   * @returns {GalaxyCatalog} Calibrated GalaxyCatalog
   */
  ingestFITS(buffer, options = {}) {
    const parsed = parseFITSBinaryTable(buffer, options);
    return this.processRecords(parsed.rows);
  }

  /**
   * Ingests IP2I ASCII table string directly.
   * @param {string} asciiString Raw IP2I ASCII text
   * @param {Object} [options] Parser options
   * @returns {GalaxyCatalog} Calibrated GalaxyCatalog
   */
  ingestIP2IASCII(asciiString, options = {}) {
    const parsed = parseIP2IASCIITable(asciiString, options);
    return this.processRecords(parsed.rows);
  }

  /**
   * Ingests CDS VizieR catalogue string directly.
   * @param {string} cdsString Raw CDS VizieR text
   * @param {Object} [options] Parser options
   * @returns {GalaxyCatalog} Calibrated GalaxyCatalog
   */
  ingestCDS(cdsString, options = {}) {
    const parsed = parseCDSVizieRCatalog(cdsString, options);
    return this.processRecords(parsed.rows);
  }
}
