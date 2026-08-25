/**
 * @file cosmic_variance_calculator.js
 * @description Theoretical Lambda-CDM Cosmic Variance and Bulk Flow Expectation Calculator.
 *
 * Implements:
 * 1. Eisenstein & Hu (1998) analytical transfer function T(k) for Lambda-CDM matter power spectrum P(k).
 * 2. Spherical top-hat window function W(kR) = 3 j_1(kR) / (kR).
 * 3. Velocity power spectrum P_v(k) = (H_0 f)^2 P(k) / k^2 with growth rate f = Omega_m^0.55.
 * 4. Theoretical bulk flow variance sigma_bulk^2(R) = (H_0 f)^2 / (2 pi^2) \int P(k) W^2(kR) dk.
 * 5. Maxwell-Boltzmann 3D bulk flow distribution and tension p-value testing.
 *
 * @module bulk-flow/cosmic_variance_calculator
 */

/**
 * Standard Cosmological Parameters
 */
export const DEFAULT_COSMOLOGY = Object.freeze({
  Omega_m: 0.315,
  Omega_Lambda: 0.685,
  Omega_b: 0.049,
  h: 0.674, // H0 = 67.4 km/s / Mpc
  H0: 67.4,
  sigma8: 0.811,
  ns: 0.965
});

/**
 * Computes spherical top-hat window function W(x) = 3 (sin x - x cos x) / x^3.
 * @param {number} x - Dimensionless wave-number product k * R
 * @returns {number}
 */
export function sphericalTopHatWindow(x) {
  if (Math.abs(x) < 1e-4) {
    // Taylor expansion: 1 - x^2/10 + x^4/280
    const x2 = x * x;
    return 1.0 - 0.1 * x2 + (1.0 / 280.0) * x2 * x2;
  }
  return 3.0 * (Math.sin(x) - x * Math.cos(x)) / (x * x * x);
}

/**
 * Computes Eisenstein & Hu (1998) zero-baryon / fitting transfer function T(k).
 * @param {number} k - Wavenumber in h/Mpc
 * @param {Object} [cosmo=DEFAULT_COSMOLOGY]
 * @returns {number} Transfer function T(k)
 */
export function eisensteinHuTransfer(k, cosmo = DEFAULT_COSMOLOGY) {
  const { Omega_m, Omega_b, h } = cosmo;
  const theta2p7 = 2.7255 / 2.7; // T_CMB / 2.7 K

  // Shape parameter Gamma with baryon correction
  const s = 44.5 * Math.log(9.83 / (Omega_m * h * h)) / Math.sqrt(1.0 + 10.0 * Math.pow(Omega_b * h * h, 0.75));
  const alphaGamma = 1.0 - 0.328 * Math.log(431.0 * Omega_m * h * h) * (Omega_b / Omega_m) + 0.38 * Math.log(22.3 * Omega_m * h * h) * Math.pow(Omega_b / Omega_m, 2);
  const GammaEff = Omega_m * h * (alphaGamma + (1.0 - alphaGamma) / (1.0 + Math.pow(0.43 * k * s, 4)));

  const q = (k * theta2p7 * theta2p7) / GammaEff;
  const L0 = Math.log(2.0 * Math.E + 1.8 * q);
  const C0 = 14.2 + 731.0 / (1.0 + 62.5 * q);

  return L0 / (L0 + C0 * q * q);
}

/**
 * Theoretical Lambda-CDM Cosmic Variance Calculator.
 */
export class CosmicVarianceCalculator {
  /**
   * @param {Object} [cosmology=DEFAULT_COSMOLOGY]
   */
  constructor(cosmology = DEFAULT_COSMOLOGY) {
    this.cosmo = { ...DEFAULT_COSMOLOGY, ...cosmology };
    this.f = Math.pow(this.cosmo.Omega_m, 0.55); // Linear growth rate f ~ Omega_m^0.55
    this.h0f = this.cosmo.H0 * this.f;

    this.normAmplitude = this._computeSigma8Normalization();
  }

  /**
   * Normalizes matter power spectrum P(k) such that \sigma(8 Mpc/h) == \sigma_8.
   * @private
   */
  _computeSigma8Normalization() {
    const R8 = 8.0;
    const numPoints = 1000;
    const kMin = 1e-4;
    const kMax = 10.0;
    const dLogK = (Math.log(kMax) - Math.log(kMin)) / (numPoints - 1);

    let integral = 0.0;
    for (let i = 0; i < numPoints; i++) {
      const k = kMin * Math.exp(i * dLogK);
      const T = eisensteinHuTransfer(k, this.cosmo);
      const W = sphericalTopHatWindow(k * R8);
      const unnormPk = Math.pow(k, this.cosmo.ns) * T * T;

      // Integrand: k^2 / (2 pi^2) * P(k) * W^2(kR) * k dk
      const integrand = (k * k * k) / (2.0 * Math.PI * Math.PI) * unnormPk * W * W;
      integral += integrand * dLogK;
    }

    const sigma8Sq = this.cosmo.sigma8 * this.cosmo.sigma8;
    return sigma8Sq / Math.max(1e-12, integral);
  }

  /**
   * Evaluates normalized linear matter power spectrum P(k) in (Mpc/h)^3.
   * @param {number} k - Wavenumber in h/Mpc
   * @returns {number} P(k)
   */
  matterPowerSpectrum(k) {
    const T = eisensteinHuTransfer(k, this.cosmo);
    return this.normAmplitude * Math.pow(k, this.cosmo.ns) * T * T;
  }

  /**
   * Computes 1D bulk flow velocity dispersion sigma_bulk(R) in km/s inside sphere of radius R (Mpc/h).
   *
   * sigma_1D^2(R) = (H0 f)^2 / (2 pi^2) \int_0^\infty P(k) W^2(kR) dk.
   *
   * @param {number} R - Radius in Mpc/h
   * @returns {number} 1D velocity dispersion sigma_1D in km/s
   */
  computeBulkFlowDispersion1D(R) {
    const numPoints = 1500;
    const kMin = 1e-4;
    const kMax = 20.0;
    const dLogK = (Math.log(kMax) - Math.log(kMin)) / (numPoints - 1);

    let integral = 0.0;
    for (let i = 0; i < numPoints; i++) {
      const k = kMin * Math.exp(i * dLogK);
      const Pk = this.matterPowerSpectrum(k);
      const W = sphericalTopHatWindow(k * R);

      // Integrand: P(k) * W^2(kR) * k dk
      integral += Pk * W * W * k * dLogK;
    }

    const factor = (this.h0f * this.h0f) / (2.0 * Math.PI * Math.PI);
    const var1D = factor * integral;

    return Math.sqrt(Math.max(0.0, var1D));
  }

  /**
   * Computes expected 3D bulk flow magnitude and confidence intervals for sphere of radius R.
   * For 3 Gaussian components with 1D variance sigma^2, magnitude V follows Maxwell-Boltzmann:
   * E[V] = sigma * sqrt(8 / pi) ~ 1.5958 * sigma.
   *
   * @param {number} R - Sphere radius in Mpc/h
   * @returns {{
   *   radiusMpc: number,
   *   sigma1D: number,
   *   sigma3D: number,
   *   expectedMagnitude: number,
   *   ci68: [number, number],
   *   ci95: [number, number]
   * }}
   */
  computeBulkFlowExpectation(R) {
    const sigma1D = this.computeBulkFlowDispersion1D(R);
    const sigma3D = Math.sqrt(3.0) * sigma1D;
    const expectedMag = sigma1D * Math.sqrt(8.0 / Math.PI);

    // Chi(3) quantiles:
    // 68.3% CI: [0.88 * sigma_1D, 2.38 * sigma_1D]
    // 95.4% CI: [0.48 * sigma_1D, 3.03 * sigma_1D]
    const ci68 = [0.88 * sigma1D, 2.38 * sigma1D];
    const ci95 = [0.48 * sigma1D, 3.03 * sigma1D];

    return {
      radiusMpc: R,
      sigma1D,
      sigma3D,
      expectedMagnitude: expectedMag,
      ci68,
      ci95
    };
  }

  /**
   * Computes statistical tension (chi^2 and p-value) between observed bulk flow and Lambda-CDM.
   *
   * @param {number} R - Sphere radius in Mpc/h
   * @param {number} observedMagnitudeKms - Measured ||V_bulk|| in km/s
   * @param {number} [measurementErrorKms=0] - 1-sigma observational error
   * @returns {{
   *   chi2: number,
   *   pValue: number,
   *   tensionSigma: number,
   *   isConsistent: boolean
   * }}
   */
  evaluateTension(R, observedMagnitudeKms, measurementErrorKms = 0) {
    const sigma1D = this.computeBulkFlowDispersion1D(R);
    const totalSigma1D = Math.sqrt(sigma1D * sigma1D + (measurementErrorKms / Math.sqrt(3.0)) ** 2);

    // Chi^2 with 3 degrees of freedom: chi^2 = (V_obs / sigma_1D)^2
    const chi2 = (observedMagnitudeKms / totalSigma1D) ** 2;

    // Survival function for Chi^2 with 3 DOF:
    // P(X > chi2) = 2 * (1 - Phi(sqrt(chi2))) + sqrt(2/pi) * sqrt(chi2) * exp(-chi2/2)
    const sqrtChi2 = Math.sqrt(chi2);
    const t = 1.0 / (1.0 + 0.2316419 * sqrtChi2);
    const poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + 1.330274429 * t))));
    const normCdf = 1.0 - (1.0 / Math.sqrt(2.0 * Math.PI)) * Math.exp(-0.5 * chi2) * poly;
    const pValue = 2.0 * (1.0 - normCdf) + Math.sqrt(2.0 / Math.PI) * sqrtChi2 * Math.exp(-0.5 * chi2);

    // Tension in sigmas
    const tensionSigma = Math.sqrt(Math.max(0.0, chi2 - 3.0));

    return {
      chi2,
      pValue: Math.min(1.0, Math.max(0.0, pValue)),
      tensionSigma,
      isConsistent: pValue > 0.01 // Consistent at 99% confidence
    };
  }
}
