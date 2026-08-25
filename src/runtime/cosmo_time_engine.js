/**
 * @file cosmo_time_engine.js
 * @module runtime/cosmo_time_engine
 * @description Cosmological Time Evolution Physics Engine for the ZRT Cosmicflows Workbench.
 * Implements physical cosmic time evolution across t in [-13.8 Gyr, +10.0 Gyr], encompassing:
 *  - Exact Friedmann scale factor a(t) for flat Lambda-CDM cosmology
 *  - Redshift z(a) and inverse relations t(a), t(z)
 *  - Hubble parameter H(z) and normalized expansion rate E(z)
 *  - Density parameter evolutions Omega_m(z), Omega_Lambda(z), Omega_r(z), Omega_k(z)
 *  - Carroll, Press & Turner (1992) linear growth factor D+(z) and exact quadrature integral D+(a)
 *  - Logarithmic growth rate f(z) = dln(D+)/dln(a) ~ Omega_m(z)^0.55
 *  - Cosmological distance-redshift integrals (Comoving, Luminosity, Angular Diameter, Lookback Time)
 *  - Backward Lagrangian Zel'dovich de-advection (z -> infty) for primordial un-clustering
 *  - Forward non-linear Eulerian sink accretion (t -> +10 Gyr) for Shapley Core and Great Attractor basins
 *  - High-performance in-place zero-allocation buffer mutation for 60+ FPS real-time scrubbing
 *  - Simulation transport controller (play/pause/speed/presets) and telemetry generation
 *
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Speed of light in vacuum (km/s).
 * @constant {number}
 */
export const SPEED_OF_LIGHT_KM_S = 299792.458;

/**
 * Standard Megaparsec to kilometer conversion constant.
 * @constant {number}
 */
export const MPC_TO_KM = 3.08567758149137e19;

/**
 * Gigayear to seconds conversion constant.
 * @constant {number}
 */
export const GYR_TO_SECONDS = 3.15576e16;

/**
 * Standard Wiener Filter displacement-to-velocity scale factor (km/s / (h^-1 Mpc)).
 * H0 * f(Omega_m) ~= 74.6 * 0.315^0.55 ~= 52.0 km/s / (h^-1 Mpc).
 * @constant {number}
 */
export const DEFAULT_VELOCITY_SCALE_FACTOR = 52.0;

/**
 * Default fiducial Lambda-CDM Cosmological Parameters (Cosmicflows-4 / Planck fiducial baseline).
 * @type {Readonly<{H0: number, Omega_m: number, Omega_L: number, Omega_r: number, Omega_k: number, t0: number, gamma: number}>}
 */
export const DEFAULT_COSMO_PARAMS = Object.freeze({
  H0: 74.6,           // Hubble constant in km/s/Mpc
  Omega_m: 0.315,     // Total matter density parameter at z=0
  Omega_L: 0.685,     // Dark energy (cosmological constant) density parameter at z=0
  Omega_r: 0.0,       // Radiation density parameter at z=0 (negligible in low-z cosmography)
  Omega_k: 0.0,       // Spatial curvature density parameter (flat universe = 0)
  t0: 13.787,         // Age of the universe in Gyr
  gamma: 0.55         // Growth index gamma (Peebles / Wang & Steinhardt)
});

/**
 * Backwards-compatible alias for default cosmological parameters.
 */
export const COSMO_PARAMS = DEFAULT_COSMO_PARAMS;

/**
 * Standard Cosmological Time Presets (in Delta t Gyr relative to present t=0).
 * @enum {number}
 */
export const TIME_PRESETS = Object.freeze({
  BIG_BANG: -13.78,       // Near primordial recombination / initial singularity limit
  COSMIC_NOON: -10.4,     // z ~ 2 peak star formation epoch
  REIONIZATION: -12.8,    // z ~ 6 reionization boundary
  GALAXY_FORM: -11.5,     // z ~ 3 early proto-cluster collapse
  PRESENT: 0.0,           // Present day epoch (z = 0, a = 1, D+ = 1)
  FUTURE_5: 5.0,          // t = +5.0 Gyr future dark-energy-dominated expansion
  FUTURE_10: 10.0         // t = +10.0 Gyr asymptotic attractor condensation
});

/**
 * Default Major Cosmological Attractor Sinks (Supergalactic Cartesian coordinates in km/s).
 * @type {ReadonlyArray<{id: string, name: string, x: number, y: number, z: number, massWeight: number, radiusScale: number, infallRate: number}>}
 */
export const DEFAULT_SINKS = Object.freeze([
  {
    id: 'shapley_core',
    name: 'Shapley Supercluster Core',
    x: 7200.0,
    y: -2400.0,
    z: 8600.0,
    massWeight: 1.0,
    radiusScale: 7500.0,
    infallRate: 0.35
  },
  {
    id: 'great_attractor',
    name: 'Great Attractor / Norma-Centaurus',
    x: -4800.0,
    y: 3900.0,
    z: 850.0,
    massWeight: 0.85,
    radiusScale: 6500.0,
    infallRate: 0.30
  },
  {
    id: 'perseus_pisces',
    name: 'Perseus-Pisces Supercluster',
    x: 4800.0,
    y: -1800.0,
    z: -3200.0,
    massWeight: 0.70,
    radiusScale: 6000.0,
    infallRate: 0.25
  },
  {
    id: 'coma_cluster',
    name: 'Coma Supercluster Core',
    x: 0.0,
    y: 7000.0,
    z: 1200.0,
    massWeight: 0.60,
    radiusScale: 5000.0,
    infallRate: 0.20
  }
]);

// ---------------------------------------------------------------------------
// Helper: Parameter normalization
// ---------------------------------------------------------------------------
function getParam(params, key, fallback) {
  if (params && typeof params[key] === 'number' && !isNaN(params[key])) {
    return params[key];
  }
  return fallback;
}

// ---------------------------------------------------------------------------
// Pure Cosmological Mathematical Functions
// ---------------------------------------------------------------------------

/**
 * Calculates exact cosmic time t0 (age of universe in Gyr) from cosmological parameters for flat Lambda-CDM.
 * Formula: t0 = (2 / (3 * H0 * sqrt(Omega_L))) * asinh(sqrt(Omega_L / Omega_m))
 * 
 * @param {number} [H0=74.6] Hubble constant in km/s/Mpc
 * @param {number} [Omega_m=0.315] Matter density parameter
 * @param {number} [Omega_L=0.685] Dark energy density parameter
 * @returns {number} Age of the universe in Gyr
 */
export function calculateCosmicAge(H0 = DEFAULT_COSMO_PARAMS.H0, Omega_m = DEFAULT_COSMO_PARAMS.Omega_m, Omega_L = DEFAULT_COSMO_PARAMS.Omega_L) {
  if (H0 <= 0 || Omega_m <= 0) {
    throw new RangeError('calculateCosmicAge: H0 and Omega_m must be positive numbers.');
  }
  const H0_Gyr = H0 * (GYR_TO_SECONDS / MPC_TO_KM);
  if (Omega_L <= 0) {
    // Einstein-de Sitter limit (Omega_m = 1, Omega_L = 0): t0 = 2 / (3 * H0)
    return 2.0 / (3.0 * H0_Gyr);
  }
  const factor = 2.0 / (3.0 * H0_Gyr * Math.sqrt(Omega_L));
  const arg = Math.sqrt(Omega_L / Omega_m);
  return factor * Math.asinh(arg);
}

/**
 * Evaluates the exact Friedmann scale factor a(t) for flat Lambda-CDM cosmology.
 *
 * Analytic Derivation:
 * In a spatially flat universe (Omega_k = 0) with pressureless matter and cosmological constant:
 *   a(t) = (Omega_m / Omega_L)^(1/3) * [ sinh( (t_cosmic / t0) * asinh(sqrt(Omega_L / Omega_m)) ) ]^(2/3)
 * For Einstein-de Sitter (Omega_L = 0, Omega_m = 1):
 *   a(t) = (t_cosmic / t0)^(2/3)
 *
 * @param {number} tGyr Offset from present epoch in Gyr, t in [-13.8, +10.0]
 * @param {Object} [params=DEFAULT_COSMO_PARAMS] Cosmological parameters
 * @returns {number} Scale factor a(t), with a(0) = 1.0 identically
 */
export function scaleFactorAtTime(tGyr, params = DEFAULT_COSMO_PARAMS) {
  const t0 = getParam(params, 't0', DEFAULT_COSMO_PARAMS.t0);
  const Om = getParam(params, 'Omega_m', DEFAULT_COSMO_PARAMS.Omega_m);
  const OL = getParam(params, 'Omega_L', DEFAULT_COSMO_PARAMS.Omega_L);

  if (typeof tGyr !== 'number' || isNaN(tGyr)) {
    return 1.0;
  }

  // Exact present day invariant
  if (Math.abs(tGyr) < 1e-12) {
    return 1.0;
  }

  const t_cosmic = Math.max(1e-7, t0 + tGyr);

  if (OL <= 1e-9) {
    // Einstein-de Sitter limit
    return Math.pow(t_cosmic / t0, 2.0 / 3.0);
  }

  const ratio = Om / OL;
  const factor = Math.cbrt(ratio); // (Omega_m / Omega_L)^(1/3)
  const alpha = Math.asinh(Math.sqrt(OL / Om));
  const sinhArg = (t_cosmic / t0) * alpha;

  const s = Math.sinh(sinhArg);
  if (s <= 0) {
    return 1e-6;
  }

  const a = factor * Math.pow(s, 2.0 / 3.0);
  return Math.max(1e-6, a);
}

/**
 * Inverts the scale factor a to find the exact offset cosmic time Delta t in Gyr.
 *
 * Inversion Formula:
 *   t_cosmic(a) = t0 * asinh( sqrt(Omega_L / Omega_m) * a^(3/2) ) / asinh( sqrt(Omega_L / Omega_m) )
 *   Delta t(a) = t_cosmic(a) - t0
 *
 * @param {number} a Scale factor (must be > 0)
 * @param {Object} [params=DEFAULT_COSMO_PARAMS] Cosmological parameters
 * @returns {number} Delta t in Gyr (t=0 at a=1)
 */
export function timeAtScaleFactor(a, params = DEFAULT_COSMO_PARAMS) {
  if (typeof a !== 'number' || isNaN(a) || a <= 0) {
    throw new RangeError(`timeAtScaleFactor: Scale factor 'a' must be a positive number, got ${a}`);
  }
  if (Math.abs(a - 1.0) < 1e-12) {
    return 0.0;
  }

  const t0 = getParam(params, 't0', DEFAULT_COSMO_PARAMS.t0);
  const Om = getParam(params, 'Omega_m', DEFAULT_COSMO_PARAMS.Omega_m);
  const OL = getParam(params, 'Omega_L', DEFAULT_COSMO_PARAMS.Omega_L);

  if (OL <= 1e-9) {
    // EdS limit: a = (t_cosmic / t0)^(2/3) -> t_cosmic = t0 * a^(3/2)
    return t0 * Math.pow(a, 1.5) - t0;
  }

  const alpha = Math.asinh(Math.sqrt(OL / Om));
  const arg = Math.sqrt(OL / Om) * Math.pow(a, 1.5);
  const t_cosmic = (t0 * Math.asinh(arg)) / alpha;

  return t_cosmic - t0;
}

/**
 * Converts scale factor a to cosmological redshift z.
 * z = 1/a - 1
 * 
 * @param {number} a Scale factor
 * @returns {number} Redshift z
 */
export function redshiftAtScaleFactor(a) {
  if (typeof a !== 'number' || isNaN(a) || a <= 0) {
    throw new RangeError(`redshiftAtScaleFactor: Scale factor 'a' must be positive, got ${a}`);
  }
  return 1.0 / a - 1.0;
}

/**
 * Converts cosmological redshift z to scale factor a.
 * a = 1 / (1 + z)
 * 
 * @param {number} z Redshift (z > -1)
 * @returns {number} Scale factor a
 */
export function scaleFactorAtRedshift(z) {
  if (typeof z !== 'number' || isNaN(z) || z <= -1.0) {
    throw new RangeError(`scaleFactorAtRedshift: Redshift 'z' must be > -1.0, got ${z}`);
  }
  return 1.0 / (1.0 + z);
}

/**
 * Dimensionless Hubble parameter E(a) = H(a) / H0.
 * E(a) = sqrt( Omega_m * a^-3 + Omega_r * a^-4 + Omega_k * a^-2 + Omega_L )
 * 
 * @param {number} a Scale factor
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} Dimensionless expansion factor E(a)
 */
export function expansionRateE(a, params = DEFAULT_COSMO_PARAMS) {
  if (typeof a !== 'number' || isNaN(a) || a <= 0) {
    throw new RangeError(`expansionRateE: Scale factor 'a' must be positive, got ${a}`);
  }
  const Om = getParam(params, 'Omega_m', DEFAULT_COSMO_PARAMS.Omega_m);
  const OL = getParam(params, 'Omega_L', DEFAULT_COSMO_PARAMS.Omega_L);
  const Or = getParam(params, 'Omega_r', DEFAULT_COSMO_PARAMS.Omega_r);
  const Ok = getParam(params, 'Omega_k', DEFAULT_COSMO_PARAMS.Omega_k);

  const a2 = a * a;
  const a3 = a2 * a;
  const a4 = a3 * a;

  const val = Om / a3 + Or / a4 + Ok / a2 + OL;
  return Math.sqrt(Math.max(1e-12, val));
}

/**
 * Hubble parameter H(z) in km/s/Mpc.
 * 
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} H(z) in km/s/Mpc
 */
export function hubbleParameter(z, params = DEFAULT_COSMO_PARAMS) {
  const H0 = getParam(params, 'H0', DEFAULT_COSMO_PARAMS.H0);
  const a = scaleFactorAtRedshift(z);
  return H0 * expansionRateE(a, params);
}

/**
 * Time-dependent matter density parameter Omega_m(z).
 * Omega_m(z) = (Omega_m,0 * (1+z)^3) / E(z)^2
 * 
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} Omega_m(z)
 */
export function matterDensityParameter(z, params = DEFAULT_COSMO_PARAMS) {
  const Om = getParam(params, 'Omega_m', DEFAULT_COSMO_PARAMS.Omega_m);
  const a = scaleFactorAtRedshift(z);
  const E = expansionRateE(a, params);
  return (Om / (a * a * a)) / (E * E);
}

/**
 * Time-dependent dark energy density parameter Omega_Lambda(z).
 * Omega_Lambda(z) = Omega_L,0 / E(z)^2
 * 
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} Omega_Lambda(z)
 */
export function darkEnergyDensityParameter(z, params = DEFAULT_COSMO_PARAMS) {
  const OL = getParam(params, 'Omega_L', DEFAULT_COSMO_PARAMS.Omega_L);
  const a = scaleFactorAtRedshift(z);
  const E = expansionRateE(a, params);
  return OL / (E * E);
}

/**
 * Carroll, Press & Turner (1992) / Eisenstein & Hu (1999) linear growth factor D+(z).
 * Normalized such that D+(z=0) = 1.0 identically.
 *
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} Linear growth factor D+(z)
 */
export function linearGrowthFactorCPT(z, params = DEFAULT_COSMO_PARAMS) {
  if (Math.abs(z) < 1e-12) {
    return 1.0;
  }

  const Om0 = getParam(params, 'Omega_m', DEFAULT_COSMO_PARAMS.Omega_m);
  const OL0 = getParam(params, 'Omega_L', DEFAULT_COSMO_PARAMS.Omega_L);

  const a = scaleFactorAtRedshift(z);
  const Om_z = matterDensityParameter(z, params);
  const OL_z = darkEnergyDensityParameter(z, params);

  function carrollG(om, ol) {
    const denom = Math.pow(om, 4.0 / 7.0) - ol + (1.0 + 0.5 * om) * (1.0 + ol / 70.0);
    return (2.5 * om) / Math.max(1e-6, denom);
  }

  const g_z = carrollG(Om_z, OL_z);
  const g_0 = carrollG(Om0, OL0);

  const D_plus = (a * g_z) / g_0;
  return Math.max(0.0, D_plus);
}

/**
 * Exact Numerical Quadrature Integration of Linear Growth Factor D+(a).
 *
 * @param {number} a Scale factor
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @param {number} [nSteps=256] Number of quadrature subdivisions
 * @returns {number} Exact normalized linear growth factor D+(a) (D+(1) = 1.0)
 */
export function linearGrowthFactorExact(a, params = DEFAULT_COSMO_PARAMS, nSteps = 256) {
  if (typeof a !== 'number' || isNaN(a) || a <= 0) {
    return 0.0;
  }
  if (Math.abs(a - 1.0) < 1e-12) {
    return 1.0;
  }

  const Om0 = getParam(params, 'Omega_m', DEFAULT_COSMO_PARAMS.Omega_m);
  const OL0 = getParam(params, 'Omega_L', DEFAULT_COSMO_PARAMS.Omega_L);

  const uMax = Math.sqrt(a);
  const h = uMax / nSteps;

  let sum = 0.0;
  for (let i = 0; i <= nSteps; i++) {
    const u = i * h;
    const u2 = u * u;
    const u4 = u2 * u2;
    const u6 = u4 * u2;
    const denom = Math.pow(Om0 + OL0 * u6, 1.5);
    const f = denom > 0 ? (2.0 * u4) / denom : 0.0;

    let weight = 2.0;
    if (i === 0 || i === nSteps) weight = 1.0;
    else if (i % 2 === 1) weight = 4.0;
    sum += weight * f;
  }
  const integral_a = (h / 3.0) * sum;
  const unnorm_D_a = expansionRateE(a, params) * integral_a;

  const uMax1 = 1.0;
  const h1 = uMax1 / nSteps;
  let sum1 = 0.0;
  for (let i = 0; i <= nSteps; i++) {
    const u = i * h1;
    const u2 = u * u;
    const u4 = u2 * u2;
    const u6 = u4 * u2;
    const denom = Math.pow(Om0 + OL0 * u6, 1.5);
    const f = denom > 0 ? (2.0 * u4) / denom : 0.0;

    let weight = 2.0;
    if (i === 0 || i === nSteps) weight = 1.0;
    else if (i % 2 === 1) weight = 4.0;
    sum1 += weight * f;
  }
  const integral_1 = (h1 / 3.0) * sum1;
  const unnorm_D_1 = expansionRateE(1.0, params) * integral_1;

  if (unnorm_D_1 <= 0) return a;
  return unnorm_D_a / unnorm_D_1;
}

/**
 * Logarithmic growth rate f(z) = dln(D+) / dln(a) ~= Omega_m(z)^gamma.
 *
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} Logarithmic growth rate f(z)
 */
export function growthRate(z, params = DEFAULT_COSMO_PARAMS) {
  const gamma = getParam(params, 'gamma', DEFAULT_COSMO_PARAMS.gamma);
  const Om_z = matterDensityParameter(z, params);
  return Math.pow(Math.max(1e-6, Om_z), gamma);
}

// ---------------------------------------------------------------------------
// Cosmological Distance Integrals
// ---------------------------------------------------------------------------

/**
 * Computes line-of-sight Comoving Distance d_C(z) in Mpc.
 *
 * @param {number} z Redshift (z >= 0)
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @param {number} [nSteps=128]
 * @returns {number} Comoving distance d_C(z) in Mpc
 */
export function comovingDistance(z, params = DEFAULT_COSMO_PARAMS, nSteps = 128) {
  if (z <= 0) return 0.0;
  const H0 = getParam(params, 'H0', DEFAULT_COSMO_PARAMS.H0);
  const c = SPEED_OF_LIGHT_KM_S;

  const h = z / nSteps;
  let sum = 0.0;
  for (let i = 0; i <= nSteps; i++) {
    const zi = i * h;
    const ai = scaleFactorAtRedshift(zi);
    const E = expansionRateE(ai, params);
    const invE = 1.0 / E;

    let weight = 2.0;
    if (i === 0 || i === nSteps) weight = 1.0;
    else if (i % 2 === 1) weight = 4.0;
    sum += weight * invE;
  }
  const integral = (h / 3.0) * sum;
  return (c / H0) * integral;
}

/**
 * Computes Luminosity Distance d_L(z) in Mpc for flat universe.
 *
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} Luminosity distance in Mpc
 */
export function luminosityDistance(z, params = DEFAULT_COSMO_PARAMS) {
  if (z <= 0) return 0.0;
  return (1.0 + z) * comovingDistance(z, params);
}

/**
 * Computes Angular Diameter Distance d_A(z) in Mpc for flat universe.
 *
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {number} Angular diameter distance in Mpc
 */
export function angularDiameterDistance(z, params = DEFAULT_COSMO_PARAMS) {
  if (z <= 0) return 0.0;
  return comovingDistance(z, params) / (1.0 + z);
}

/**
 * Computes Lookback Time t_L(z) in Gyr.
 *
 * @param {number} z Redshift
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @param {number} [nSteps=128]
 * @returns {number} Lookback time in Gyr
 */
export function lookbackTime(z, params = DEFAULT_COSMO_PARAMS, nSteps = 128) {
  if (z <= 0) return 0.0;
  const H0 = getParam(params, 'H0', DEFAULT_COSMO_PARAMS.H0);
  const H0_Gyr = H0 * (GYR_TO_SECONDS / MPC_TO_KM);

  const h = z / nSteps;
  let sum = 0.0;
  for (let i = 0; i <= nSteps; i++) {
    const zi = i * h;
    const ai = scaleFactorAtRedshift(zi);
    const E = expansionRateE(ai, params);
    const f = 1.0 / ((1.0 + zi) * E);

    let weight = 2.0;
    if (i === 0 || i === nSteps) weight = 1.0;
    else if (i % 2 === 1) weight = 4.0;
    sum += weight * f;
  }
  const integral = (h / 3.0) * sum;
  return integral / H0_Gyr;
}

// ---------------------------------------------------------------------------
// Unified Cosmological State Computation
// ---------------------------------------------------------------------------

/**
 * Primary Cosmological State Evaluation Function.
 * Computes full analytical cosmological state at cosmic time offset tGyr in [-13.8, +10.0] Gyr.
 *
 * @param {number} tGyr Offset from present epoch in Gyr
 * @param {Object} [params=DEFAULT_COSMO_PARAMS]
 * @returns {Object} Comprehensive cosmological telemetry state
 */
export function computeCosmology(tGyr, params = DEFAULT_COSMO_PARAMS) {
  const t0 = getParam(params, 't0', DEFAULT_COSMO_PARAMS.t0);
  const Om0 = getParam(params, 'Omega_m', DEFAULT_COSMO_PARAMS.Omega_m);
  const OL0 = getParam(params, 'Omega_L', DEFAULT_COSMO_PARAMS.Omega_L);
  const H0 = getParam(params, 'H0', DEFAULT_COSMO_PARAMS.H0);

  const dtClamped = Math.max(-13.8, Math.min(10.0, typeof tGyr === 'number' && !isNaN(tGyr) ? tGyr : 0.0));
  const t_cosmic = Math.max(1e-6, t0 + dtClamped);

  const a = scaleFactorAtTime(dtClamped, params);
  const z = redshiftAtScaleFactor(a);
  const E_a = expansionRateE(a, params);
  const H_z = H0 * E_a;
  const Om_z = matterDensityParameter(z, params);
  const OL_z = darkEnergyDensityParameter(z, params);

  const D_plus = linearGrowthFactorCPT(z, params);
  const f_z = growthRate(z, params);

  return {
    t: dtClamped,
    t_cosmic: t_cosmic,
    a: a,
    scaleFactor: a,
    z: z,
    redshift: z,
    D_plus: D_plus,
    dPlus: D_plus,
    D: D_plus,
    growthFactor: D_plus,
    H_z: H_z,
    Hz: H_z,
    H: H_z,
    hubble: H_z,
    Omega_m: Om_z,
    omegaM: Om_z,
    Omega_L: OL_z,
    Omega_Lambda: OL_z,
    omegaL: OL_z,
    f: f_z,
    f_z: f_z,
    growthRate: f_z,
    E: E_a,
    params: { ...params }
  };
}

// ---------------------------------------------------------------------------
// Particle De-Advection & Eulerian Sink Dynamics
// ---------------------------------------------------------------------------

/**
 * Backward Lagrangian Zel'dovich De-Advection (z -> infty, t < 0).
 * In linear perturbation theory: x(t) = q + D+(t) * Psi(q).
 *
 * @param {Float32Array|Float64Array|number[]} positions Flat [x, y, z, x, y, z, ...] array
 * @param {Float32Array|Float64Array|number[]} velocities Flat [vx, vy, vz, ...] array (km/s)
 * @param {number} tGyr Offset time in Gyr (tGyr <= 0 for past)
 * @param {Object} [options]
 * @param {number} [options.velocityScale=52.0]
 * @param {Float32Array|Float64Array} [options.outBuffer] Target buffer for zero-alloc output
 * @returns {Float32Array|Float64Array} De-advected positions array
 */
export function deAdvectLagrangianParticles(positions, velocities, tGyr, options = {}) {
  if (!positions || !velocities) {
    throw new TypeError('deAdvectLagrangianParticles: positions and velocities must be defined.');
  }
  const count = Math.floor(positions.length / 3);
  const out = options.outBuffer || new Float32Array(positions.length);

  const cosmo = computeCosmology(tGyr);
  const D_plus = cosmo.D_plus;
  const dispFactor = Math.max(0.0, 1.0 - D_plus);

  for (let i = 0; i < count; i++) {
    const idx = i * 3;
    const bx = positions[idx];
    const by = positions[idx + 1];
    const bz = positions[idx + 2];

    const vx = velocities[idx] || 0.0;
    const vy = velocities[idx + 1] || 0.0;
    const vz = velocities[idx + 2] || 0.0;

    out[idx] = bx - vx * dispFactor;
    out[idx + 1] = by - vy * dispFactor;
    out[idx + 2] = bz - vz * dispFactor;
  }

  return out;
}

/**
 * Forward Non-Linear Eulerian Sink Accretion (t > 0, t -> +10 Gyr).
 *
 * @param {Float32Array|Float64Array|number[]} positions Base position array [x, y, z, ...]
 * @param {Float32Array|Float64Array|number[]} velocities Base velocity array [vx, vy, vz, ...]
 * @param {number} tGyr Future offset time in Gyr (tGyr >= 0)
 * @param {Array<Object>} [sinks=DEFAULT_SINKS] Array of sink definitions
 * @param {Object} [options]
 * @param {Float32Array|Float64Array} [options.outBuffer]
 * @returns {Float32Array|Float64Array} Evolved future coordinates
 */
export function evolveEulerianSinks(positions, velocities, tGyr, sinks = DEFAULT_SINKS, options = {}) {
  if (!positions) {
    throw new TypeError('evolveEulerianSinks: positions must be defined.');
  }
  const count = Math.floor(positions.length / 3);
  const out = options.outBuffer || new Float32Array(positions.length);

  const dtClamped = Math.max(0.0, Math.min(10.0, tGyr));
  const futS = dtClamped / 10.0;
  const flowCoeff = 0.2 * dtClamped;

  for (let i = 0; i < count; i++) {
    const idx = i * 3;
    const bx = positions[idx];
    const by = positions[idx + 1];
    const bz = positions[idx + 2];

    let vx = 0.0, vy = 0.0, vz = 0.0;
    if (velocities) {
      vx = velocities[idx] || 0.0;
      vy = velocities[idx + 1] || 0.0;
      vz = velocities[idx + 2] || 0.0;
    }

    let totalInfX = 0.0;
    let totalInfY = 0.0;
    let totalInfZ = 0.0;

    for (let s = 0; s < sinks.length; s++) {
      const sink = sinks[s];
      const dx = sink.x - bx;
      const dy = sink.y - by;
      const dz = sink.z - bz;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const maxRadius = sink.radiusScale * 2.0;

      if (dist < maxRadius && dist > 1.0) {
        const infFrac = sink.infallRate * futS * Math.exp(-dist / sink.radiusScale);
        totalInfX += dx * infFrac;
        totalInfY += dy * infFrac;
        totalInfZ += dz * infFrac;
      }
    }

    out[idx] = bx + vx * flowCoeff + totalInfX;
    out[idx + 1] = by + vy * flowCoeff + totalInfY;
    out[idx + 2] = bz + vz * flowCoeff + totalInfZ;
  }

  return out;
}

/**
 * Unified Particle Mesh Cosmic Evolution (Seamless t in [-13.8, +10.0] Gyr).
 *
 * @param {Float32Array|Float64Array} basePositions Original z=0 galaxy coordinates
 * @param {Float32Array|Float64Array} baseVelocities Original peculiar velocity vectors
 * @param {Float32Array|Float64Array} targetBuffer Active rendering coordinate buffer to mutate
 * @param {number} tGyr Cosmic time offset in Gyr
 * @param {Object} [options]
 * @returns {Float32Array|Float64Array} Mutated targetBuffer
 */
export function evolveParticleMesh(basePositions, baseVelocities, targetBuffer, tGyr, options = {}) {
  if (!basePositions || !targetBuffer) {
    throw new TypeError('evolveParticleMesh: basePositions and targetBuffer must be defined.');
  }

  const dt = Math.max(-13.8, Math.min(10.0, typeof tGyr === 'number' && !isNaN(tGyr) ? tGyr : 0.0));
  const sinks = options.sinks || DEFAULT_SINKS;

  if (dt < 0) {
    return deAdvectLagrangianParticles(basePositions, baseVelocities, dt, {
      outBuffer: targetBuffer,
      velocityScale: options.velocityScale || DEFAULT_VELOCITY_SCALE_FACTOR
    });
  } else {
    return evolveEulerianSinks(basePositions, baseVelocities, dt, sinks, {
      outBuffer: targetBuffer
    });
  }
}

/**
 * Evolves a 3D scalar density contrast field delta(x, t) across cosmic time.
 *
 * @param {Float32Array|Float64Array} densityBuffer 3D scalar density grid (delta = rho/rho_bar - 1)
 * @param {Float32Array|Float64Array} targetBuffer Target buffer to populate
 * @param {number} tGyr Cosmic time offset in Gyr
 * @returns {Float32Array|Float64Array} Evolved density field
 */
export function evolveDensityField(densityBuffer, targetBuffer, tGyr) {
  if (!densityBuffer || !targetBuffer || densityBuffer.length !== targetBuffer.length) {
    throw new Error('evolveDensityField: Buffer sizes must match.');
  }

  const cosmo = computeCosmology(tGyr);
  const D_plus = cosmo.D_plus;
  const dt = cosmo.t;
  const len = densityBuffer.length;

  if (dt <= 0) {
    for (let i = 0; i < len; i++) {
      targetBuffer[i] = densityBuffer[i] * D_plus;
    }
  } else {
    const growth = 1.0 + 0.2 * dt;
    for (let i = 0; i < len; i++) {
      const delta = densityBuffer[i];
      if (delta > 0) {
        targetBuffer[i] = delta * Math.pow(growth, 1.2);
      } else {
        targetBuffer[i] = Math.max(-1.0, delta * growth);
      }
    }
  }

  return targetBuffer;
}

// ---------------------------------------------------------------------------
// Cosmological Time Engine Class
// ---------------------------------------------------------------------------

/**
 * Cosmological Time Evolution Engine Class.
 */
export class CosmoTimeEngine {
  /**
   * @param {Object} [customParams] Optional cosmological parameter overrides
   */
  constructor(customParams = {}) {
    this.params = { ...DEFAULT_COSMO_PARAMS, ...customParams };
    this.currentTime = 0.0; // Gyr offset from present
    this.isPlayingState = false;
    this.playbackSpeed = 1.0; // Gyr per real-time second
    this.PRESETS = TIME_PRESETS;
    this.PARAMS = this.params;
    this.DEFAULT_SINKS = DEFAULT_SINKS;
    this.lastFrameTimestamp = null;
  }

  /**
   * Evaluates cosmological state at specified time offset or current simulation time.
   * @param {number} [tGyr]
   * @returns {Object}
   */
  computeCosmology(tGyr) {
    const t = typeof tGyr === 'number' ? tGyr : this.currentTime;
    return computeCosmology(t, this.params);
  }

  /**
   * Sets simulation time in Gyr relative to present day.
   * @param {number} tGyr Time in [-13.8, +10.0] Gyr
   * @returns {number} Clamped time applied
   */
  setTime(tGyr) {
    this.currentTime = Math.max(-13.8, Math.min(10.0, typeof tGyr === 'number' && !isNaN(tGyr) ? tGyr : 0.0));
    return this.currentTime;
  }

  /**
   * Gets current simulation cosmic time offset in Gyr.
   * @returns {number}
   */
  getTime() {
    return this.currentTime;
  }

  /**
   * Starts transport playback loop.
   */
  play() {
    this.isPlayingState = true;
    this.lastFrameTimestamp = null;
  }

  /**
   * Pauses transport playback loop.
   */
  pause() {
    this.isPlayingState = false;
    this.lastFrameTimestamp = null;
  }

  /**
   * Returns current transport playing status.
   * @returns {boolean}
   */
  isPlaying() {
    return this.isPlayingState;
  }

  /**
   * Sets playback animation rate in Gyr per second.
   * @param {number} speed
   */
  setPlaybackSpeed(speed) {
    if (typeof speed === 'number' && !isNaN(speed) && speed > 0) {
      this.playbackSpeed = speed;
    }
  }

  /**
   * Alias for setPlaybackSpeed.
   * @param {number} speed
   */
  setSpeed(speed) {
    this.setPlaybackSpeed(speed);
  }

  /**
   * Advances simulation time by dtGyr.
   * @param {number} dtGyr
   * @returns {number} Updated time
   */
  step(dtGyr) {
    let nextTime = this.currentTime + dtGyr;
    if (nextTime > 10.0) {
      nextTime = -13.78;
    } else if (nextTime < -13.8) {
      nextTime = 10.0;
    }
    this.currentTime = nextTime;
    return this.currentTime;
  }

  /**
   * Performs an animation frame tick using elapsed wall-clock milliseconds.
   * @param {number} nowMs Current timestamp from performance.now()
   * @returns {Object} Telemetry of updated state
   */
  tick(nowMs) {
    if (!this.isPlayingState) {
      return this.computeCosmology(this.currentTime);
    }

    if (this.lastFrameTimestamp !== null) {
      const elapsedSec = (nowMs - this.lastFrameTimestamp) / 1000.0;
      const dtGyr = Math.min(0.5, elapsedSec * this.playbackSpeed);
      this.step(dtGyr);
    }
    this.lastFrameTimestamp = nowMs;
    return this.computeCosmology(this.currentTime);
  }

  /**
   * Generates formatted HUD Telemetry object for UI rendering.
   * @param {number} [tGyr]
   * @returns {Object} Formatted display strings and raw values
   */
  getHUDTelemetry(tGyr) {
    const cosmo = this.computeCosmology(tGyr);
    return {
      raw: cosmo,
      formattedTime: (cosmo.t >= 0 ? '+' : '') + cosmo.t.toFixed(2) + ' Gyr',
      formattedCosmicAge: cosmo.t_cosmic.toFixed(3) + ' Gyr',
      formattedRedshift: cosmo.z >= 100 ? Math.round(cosmo.z).toLocaleString() : cosmo.z.toFixed(2),
      formattedScaleFactor: cosmo.a.toFixed(3),
      formattedGrowthFactor: cosmo.D_plus.toFixed(3),
      formattedHubble: Math.round(cosmo.H_z).toLocaleString() + ' km/s/Mpc',
      formattedOmegaM: (cosmo.Omega_m * 100).toFixed(1) + '%',
      formattedOmegaL: (cosmo.Omega_L * 100).toFixed(1) + '%'
    };
  }

  /**
   * De-advects or collapses a mesh buffer in-place.
   * @param {Float32Array|Float64Array} basePositions
   * @param {Float32Array|Float64Array} baseVelocities
   * @param {Float32Array|Float64Array} targetPositions
   * @param {number} [tGyr]
   * @returns {Float32Array|Float64Array}
   */
  evolveMesh(basePositions, baseVelocities, targetPositions, tGyr) {
    const t = typeof tGyr === 'number' ? tGyr : this.currentTime;
    return evolveParticleMesh(basePositions, baseVelocities, targetPositions, t, {
      sinks: this.DEFAULT_SINKS
    });
  }
}

/**
 * Singleton Default Export of Cosmological Time Engine.
 */
export const cosmoTimeEngine = new CosmoTimeEngine();
