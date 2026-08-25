/**
 * @file cosmological_constants.js
 * @description Cosmological parameters, Lambda-CDM background evolution, growth factors,
 * and astrophysical constants for the CosmicFlows-4 Research Workbench.
 * 
 * Implements:
 * 1. Physical constants in SI, CGS, and astrophysical units (CODATA 2022 / IAU 2015).
 * 2. Cosmological parameter sets:
 *    - CF4 Standard (H0 = 74.6 km/s/Mpc, Omega_m = 0.31, Omega_Lambda = 0.69)
 *    - Planck 2018 (H0 = 67.4 km/s/Mpc, Omega_m = 0.315, Omega_Lambda = 0.685)
 *    - WMAP 9-Year (H0 = 70.0 km/s/Mpc, Omega_m = 0.279, Omega_Lambda = 0.721)
 *    - Einstein-de Sitter fiducial benchmark (Omega_m = 1.0, Omega_Lambda = 0.0)
 * 3. Distance measures:
 *    - Expansion rate E(z) = H(z) / H0
 *    - Comoving radial distance D_C(z) = c / H0 * int_0^z dz' / E(z')
 *    - Transverse comoving distance D_M(z)
 *    - Angular diameter distance D_A(z) = D_M(z) / (1 + z)
 *    - Luminosity distance D_L(z) = (1 + z) * D_M(z)
 *    - Lookback time t_L(z) and cosmic age t(z)
 * 4. Linear perturbation theory:
 *    - Matter density parameter evolution Omega_m(a) = Omega_m0 * a^-3 / E^2(a)
 *    - Linear growth factor D(a) via exact numerical quadrature of Heath (1977) integral:
 *      D(a) = 5 * Omega_m0 / 2 * E(a) * int_0^a da' / (a' * E(a'))^3
 *    - Peebles growth rate parameter: f(a) = d ln D / d ln a ~= Omega_m(a)^0.55
 *    - Linear velocity-density continuity scaling factor: H(a) * f(a) (km/s / (Mpc/h))
 * 
 * @module units/cosmological_constants
 */

/**
 * Standard Astrometric and Physical Constants.
 * @readonly
 */
export const PHYSICAL_CONSTANTS = Object.freeze({
  // Speed of light in vacuum
  SPEED_OF_LIGHT_SI: 299792458.0,             // m / s
  SPEED_OF_LIGHT_KMS: 299792.458,             // km / s
  SPEED_OF_LIGHT_CGS: 2.99792458e10,          // cm / s

  // Newtonian Gravitational Constant
  GRAVITATIONAL_CONSTANT_SI: 6.67430e-11,     // m^3 / (kg * s^2)
  GRAVITATIONAL_CONSTANT_CGS: 6.67430e-8,     // cm^3 / (g * s^2)
  // G in (km/s)^2 * (Mpc/h) / (10^10 M_sun/h)
  GRAVITATIONAL_CONSTANT_ASTRO: 4.30091e-4,   // (km/s)^2 * Mpc / (10^10 M_sun)
  GRAVITATIONAL_CONSTANT_CF4: 4.30091e-9,     // (km/s)^2 * Mpc / M_sun

  // Solar Mass and Luminosity
  SOLAR_MASS_KG: 1.98847e30,                  // kg
  SOLAR_MASS_G: 1.98847e33,                   // g
  SOLAR_LUMINOSITY_WATTS: 3.828e26,           // W
  SOLAR_LUMINOSITY_ERG_S: 3.828e33,           // erg / s

  // Distances
  METERS_PER_MPC: 3.085677581491367e22,       // m
  METERS_PER_KPC: 3.085677581491367e19,       // m
  METERS_PER_PC: 3.085677581491367e16,        // m
  METERS_PER_AU: 1.495978707e11,              // m
  METERS_PER_LY: 9.4607304725808e15,          // m
  KM_PER_MPC: 3.085677581491367e19,           // km
  KM_PER_KPC: 3.085677581491367e16,           // km

  // Times
  SECONDS_PER_JULIAN_YEAR: 31557600.0,        // s (365.25 days)
  SECONDS_PER_MYR: 3.15576e13,                // s
  SECONDS_PER_GYR: 3.15576e16,                // s

  // Critical density conversion constant: 3 * H100^2 / (8 * pi * G) in M_sun / Mpc^3
  // H100 = 100 km/s/Mpc = 3.240779289e-18 s^-1
  CRITICAL_DENSITY_COEFF_MSUN_MPC3: 2.77536627e11, // (M_sun / Mpc^3) / h^2
  CRITICAL_DENSITY_COEFF_KG_M3: 1.87847e-26        // (kg / m^3) / h^2
});

/**
 * Standard Cosmological Models.
 * @readonly
 */
export const COSMOLOGICAL_MODELS = Object.freeze({
  CF4: Object.freeze({
    name: 'CosmicFlows-4 Standard',
    H0: 74.6,             // km/s/Mpc
    h: 0.746,
    Omega_m: 0.31,
    Omega_Lambda: 0.69,
    Omega_b: 0.048,
    Omega_r: 8.48e-5,
    Omega_k: 0.0,
    w0: -1.0,
    wa: 0.0,
    sigma8: 0.81,
    ns: 0.965,
    growthIndexGamma: 0.55
  }),
  PLANCK18: Object.freeze({
    name: 'Planck 2018 (TT,TE,EE+lowE+lensing)',
    H0: 67.4,             // km/s/Mpc
    h: 0.674,
    Omega_m: 0.315,
    Omega_Lambda: 0.685,
    Omega_b: 0.0493,
    Omega_r: 9.16e-5,
    Omega_k: 0.0,
    w0: -1.0,
    wa: 0.0,
    sigma8: 0.811,
    ns: 0.9649,
    growthIndexGamma: 0.55
  }),
  WMAP9: Object.freeze({
    name: 'WMAP 9-Year',
    H0: 70.0,
    h: 0.700,
    Omega_m: 0.279,
    Omega_Lambda: 0.721,
    Omega_b: 0.0463,
    Omega_r: 8.5e-5,
    Omega_k: 0.0,
    w0: -1.0,
    wa: 0.0,
    sigma8: 0.821,
    ns: 0.972,
    growthIndexGamma: 0.55
  }),
  EDS: Object.freeze({
    name: 'Einstein-de Sitter (Fiducial Flat Matter)',
    H0: 70.0,
    h: 0.700,
    Omega_m: 1.0,
    Omega_Lambda: 0.0,
    Omega_b: 0.05,
    Omega_r: 0.0,
    Omega_k: 0.0,
    w0: -1.0,
    wa: 0.0,
    sigma8: 0.80,
    ns: 1.0,
    growthIndexGamma: 0.60
  })
});

/**
 * CosmologicalCalculator provides high-precision cosmological distance integrals,
 * background evolution, linear growth factors, and peculiar velocity continuity coefficients.
 */
export class CosmologicalCalculator {
  /**
   * @param {Object} [params] Cosmological parameters (defaults to CF4 standard).
   */
  constructor(params = COSMOLOGICAL_MODELS.CF4) {
    this.H0 = params.H0 ?? 74.6;
    this.h = params.h ?? (this.H0 / 100.0);
    this.Omega_m = params.Omega_m ?? 0.31;
    this.Omega_Lambda = params.Omega_Lambda ?? (1.0 - this.Omega_m);
    this.Omega_r = params.Omega_r ?? 0.0;
    this.Omega_k = params.Omega_k ?? (1.0 - this.Omega_m - this.Omega_Lambda - this.Omega_r);
    this.w0 = params.w0 ?? -1.0;
    this.wa = params.wa ?? 0.0;
    this.gamma = params.growthIndexGamma ?? 0.55;
    this.sigma8 = params.sigma8 ?? 0.81;
    this.ns = params.ns ?? 0.965;

    // Derived Hubble distance c / H0 in Mpc and Mpc/h
    this.HubbleDistanceMpc = PHYSICAL_CONSTANTS.SPEED_OF_LIGHT_KMS / this.H0;
    this.HubbleDistanceMpcOverH = PHYSICAL_CONSTANTS.SPEED_OF_LIGHT_KMS / 100.0; // exactly c/100 = 2997.92458 Mpc/h
    this.HubbleTimeGyr = (1.0 / (this.H0 * 1.0e3 / PHYSICAL_CONSTANTS.METERS_PER_MPC)) / PHYSICAL_CONSTANTS.SECONDS_PER_GYR;
    this.criticalDensity0_Msun_Mpc3 = PHYSICAL_CONSTANTS.CRITICAL_DENSITY_COEFF_MSUN_MPC3 * (this.h * this.h);
    this.meanMatterDensity0_Msun_Mpc3 = this.criticalDensity0_Msun_Mpc3 * this.Omega_m;
  }

  /**
   * Dark energy equation of state w(a) = w0 + wa * (1 - a) (Chevallier-Polarski-Linder param).
   * @param {number} a Scale factor.
   * @returns {number}
   */
  w(a) {
    return this.w0 + this.wa * (1.0 - a);
  }

  /**
   * Dark energy density scaling X_de(a) = exp(-3 int_1^a (1 + w(a'))/a' da').
   * @param {number} a Scale factor.
   * @returns {number}
   */
  darkEnergyScaling(a) {
    if (a <= 0.0) return 0.0;
    if (this.w0 === -1.0 && this.wa === 0.0) {
      return 1.0; // Pure cosmological constant Lambda
    }
    // For CPL: X_de(a) = a^(-3*(1 + w0 + wa)) * exp(-3 * wa * (1 - a))
    return Math.pow(a, -3.0 * (1.0 + this.w0 + this.wa)) * Math.exp(-3.0 * this.wa * (1.0 - a));
  }

  /**
   * Dimensionless Hubble expansion rate E(z) = H(z) / H0.
   * @param {number} z Redshift.
   * @returns {number}
   */
  E_z(z) {
    const a = 1.0 / (1.0 + z);
    return this.E_a(a);
  }

  /**
   * Dimensionless Hubble expansion rate E(a) = H(a) / H0.
   * @param {number} a Scale factor (a = 1 / (1 + z)).
   * @returns {number}
   */
  E_a(a) {
    if (a <= 0.0) return Infinity;
    const invA = 1.0 / a;
    const invA2 = invA * invA;
    const invA3 = invA2 * invA;
    const invA4 = invA3 * invA;

    const termM = this.Omega_m * invA3;
    const termR = this.Omega_r * invA4;
    const termK = this.Omega_k * invA2;
    const termDE = this.Omega_Lambda * this.darkEnergyScaling(a);

    const E2 = termM + termR + termK + termDE;
    return Math.sqrt(Math.max(1e-12, E2));
  }

  /**
   * Hubble parameter at redshift z in km/s/Mpc.
   * @param {number} z Redshift.
   * @returns {number}
   */
  H(z) {
    return this.H0 * this.E_z(z);
  }

  /**
   * Matter density parameter Omega_m(a) at scale factor a.
   * @param {number} a Scale factor.
   * @returns {number}
   */
  Omega_m_a(a) {
    const E = this.E_a(a);
    return (this.Omega_m * Math.pow(a, -3.0)) / (E * E);
  }

  /**
   * Matter density parameter Omega_m(z) at redshift z.
   * @param {number} z Redshift.
   * @returns {number}
   */
  Omega_m_z(z) {
    return this.Omega_m_a(1.0 / (1.0 + z));
  }

  /**
   * Linear growth rate f(z) = d ln D / d ln a ~= [Omega_m(z)]^gamma.
   * @param {number} z Redshift.
   * @returns {number}
   */
  growthRate_f(z = 0.0) {
    const om = this.Omega_m_z(z);
    return Math.pow(Math.max(0.0, om), this.gamma);
  }

  /**
   * Exact linear growth factor D(a) normalized to D(a=1) = 1.
   * Computed via Simpson numerical integration of Heath (1977) integral:
   *   D(a) propto E(a) * int_0^a da' / [a' * E(a')]^3
   * @param {number} a Scale factor in (0, 1].
   * @param {number} [numSteps=500] Quadrature subdivisions.
   * @returns {number} Normalized growth factor D(a) / D(1).
   */
  growthFactor_D(a = 1.0, numSteps = 500) {
    if (a <= 0.0) return 0.0;
    if (Math.abs(a - 1.0) < 1e-7) return 1.0;

    const unnormalizedD = (scaleA) => {
      if (scaleA <= 1e-6) return scaleA; // Matter-dominated asymptotic D(a) -> a
      const n = Math.max(100, numSteps);
      const h = scaleA / n;
      let sum = 0.0;

      for (let i = 0; i <= n; i++) {
        const aPrime = Math.max(1e-6, i * h);
        const EPrime = this.E_a(aPrime);
        const denom = Math.pow(aPrime * EPrime, 3.0);
        const integrand = denom > 1e-30 ? 1.0 / denom : 0.0;

        let weight = 2.0;
        if (i === 0 || i === n) weight = 1.0;
        else if (i % 2 === 1) weight = 4.0;

        sum += weight * integrand;
      }
      const integral = (h / 3.0) * sum;
      return (5.0 / 2.0) * this.Omega_m * this.E_a(scaleA) * integral;
    };

    const dA = unnormalizedD(a);
    const d1 = unnormalizedD(1.0);
    return d1 > 0.0 ? dA / d1 : a;
  }

  /**
   * Linear velocity-density continuity scale factor:
   *   H0 * f(z=0) in (km/s) / (Mpc/h)
   * 
   * Astrophysical note:
   * For H0 = 74.6 km/s/Mpc and Omega_m = 0.31:
   *   h = 0.746
   *   f = (0.31)^0.55 = 0.5257
   *   H0 * f = 74.6 * 0.5257 = 39.22 km/s/Mpc = (39.22 / 0.746) = 52.57 km/s / (Mpc/h) ~= 52.0 km/s/(Mpc/h)
   * 
   * This is the exact physical origin of the CF4 x52.0 velocity scale factor!
   * @param {number} [z=0.0] Redshift.
   * @returns {number} Velocity continuity factor in km/s / (Mpc/h).
   */
  linearContinuityFactor(z = 0.0) {
    const f = this.growthRate_f(z);
    const H_z = this.H(z); // km/s / Mpc
    // In units of (km/s) / (Mpc/h) = (H_z * f) / h
    return (H_z * f) / this.h;
  }

  /**
   * Comoving line-of-sight distance D_C(z) in Mpc and Mpc/h.
   * @param {number} z Redshift >= 0.
   * @param {number} [numSteps=500] Quadrature subdivisions.
   * @returns {{ comovingDistanceMpc: number, comovingDistanceMpcOverH: number }}
   */
  comovingDistance(z, numSteps = 500) {
    if (z <= 0.0) {
      return { comovingDistanceMpc: 0.0, comovingDistanceMpcOverH: 0.0 };
    }

    const n = Math.max(50, numSteps);
    const dz = z / n;
    let sum = 0.0;

    for (let i = 0; i <= n; i++) {
      const zPrime = i * dz;
      const Ez = this.E_z(zPrime);
      const integrand = 1.0 / Ez;

      let weight = 2.0;
      if (i === 0 || i === n) weight = 1.0;
      else if (i % 2 === 1) weight = 4.0;

      sum += weight * integrand;
    }

    const integral = (dz / 3.0) * sum;
    const dMpc = this.HubbleDistanceMpc * integral;
    const dMpcOverH = dMpc * this.h;

    return {
      comovingDistanceMpc: dMpc,
      comovingDistanceMpcOverH: dMpcOverH
    };
  }

  /**
   * Transverse comoving distance D_M(z) in Mpc.
   * @param {number} z Redshift.
   * @returns {number}
   */
  transverseComovingDistance(z) {
    const { comovingDistanceMpc: Dc } = this.comovingDistance(z);
    const Ok = this.Omega_k;
    const Dh = this.HubbleDistanceMpc;

    if (Math.abs(Ok) < 1e-6) {
      return Dc; // Flat universe
    } else if (Ok > 0.0) {
      // Open universe
      const sqrtOk = Math.sqrt(Ok);
      return Dh / sqrtOk * Math.sinh(sqrtOk * Dc / Dh);
    } else {
      // Closed universe
      const sqrtAbsOk = Math.sqrt(-Ok);
      return Dh / sqrtAbsOk * Math.sin(sqrtAbsOk * Dc / Dh);
    }
  }

  /**
   * Angular diameter distance D_A(z) = D_M(z) / (1 + z) in Mpc.
   * @param {number} z Redshift.
   * @returns {number}
   */
  angularDiameterDistance(z) {
    if (z <= 0.0) return 0.0;
    return this.transverseComovingDistance(z) / (1.0 + z);
  }

  /**
   * Luminosity distance D_L(z) = (1 + z) * D_M(z) in Mpc.
   * @param {number} z Redshift.
   * @returns {number}
   */
  luminosityDistance(z) {
    if (z <= 0.0) return 0.0;
    return (1.0 + z) * this.transverseComovingDistance(z);
  }

  /**
   * Distance modulus mu = 5 * log10(D_L / 10 pc) = 5 * log10(D_L_Mpc) + 25.
   * @param {number} z Redshift.
   * @returns {number}
   */
  distanceModulus(z) {
    const dL = this.luminosityDistance(z);
    if (dL <= 0.0) return 0.0;
    return 5.0 * Math.log10(dL) + 25.0;
  }

  /**
   * Inverts distance modulus mu to luminosity distance D_L in Mpc.
   * @param {number} mu Distance modulus.
   * @returns {number}
   */
  distanceFromModulus(mu) {
    return Math.pow(10.0, (mu - 25.0) / 5.0);
  }

  /**
   * Lookback time t_L(z) in Gyr.
   * t_L(z) = t_H * int_0^z dz' / ((1 + z') * E(z'))
   * @param {number} z Redshift.
   * @param {number} [numSteps=500]
   * @returns {number} Lookback time in Gyr.
   */
  lookbackTime(z, numSteps = 500) {
    if (z <= 0.0) return 0.0;

    const n = Math.max(50, numSteps);
    const dz = z / n;
    let sum = 0.0;

    for (let i = 0; i <= n; i++) {
      const zPrime = i * dz;
      const Ez = this.E_z(zPrime);
      const integrand = 1.0 / ((1.0 + zPrime) * Ez);

      let weight = 2.0;
      if (i === 0 || i === n) weight = 1.0;
      else if (i % 2 === 1) weight = 4.0;

      sum += weight * integrand;
    }

    const integral = (dz / 3.0) * sum;
    return this.HubbleTimeGyr * integral;
  }

  /**
   * Age of the universe at redshift z in Gyr.
   * @param {number} z Redshift.
   * @returns {number}
   */
  ageOfUniverse(z = 0.0) {
    const t0 = this.lookbackTime(1000.0); // Total age today
    const tLookback = this.lookbackTime(z);
    return Math.max(0.0, t0 - tLookback);
  }

  /**
   * Inverts pure Hubble recessional velocity cz (km/s) into comoving distance (Mpc/h).
   * In the local cosmological universe (cz < 30,000 km/s), d = cz / 100 in Mpc/h.
   * @param {number} cz Recessional velocity in km/s.
   * @returns {number} Comoving distance in Mpc/h.
   */
  velocityToComovingDistanceMpcOverH(cz) {
    const z = cz / PHYSICAL_CONSTANTS.SPEED_OF_LIGHT_KMS;
    if (z <= 0.0) return 0.0;
    if (z < 0.05) {
      // Local linear Hubble approximation: d (Mpc/h) = cz / 100
      return cz / 100.0;
    }
    return this.comovingDistance(z).comovingDistanceMpcOverH;
  }
}
