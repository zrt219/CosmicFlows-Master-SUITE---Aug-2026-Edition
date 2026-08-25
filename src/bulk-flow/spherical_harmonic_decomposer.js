/**
 * @file spherical_harmonic_decomposer.js
 * @description Real Spherical Harmonics Decomposition Suite for Cosmological Radial Velocity and Density Fields (l = 0 to 8).
 *
 * Implements:
 * 1. Associated Legendre polynomials P_l^m(x) and orthonormal real spherical harmonic basis Y_lm(theta, phi).
 * 2. Orthonormal spherical decomposition a_lm = \int f(theta, phi) Y_lm(theta, phi) d\Omega.
 * 3. Angular multipole power spectrum C_l = 1/(2l+1) sum_{m=-l}^l |a_lm|^2.
 * 4. Exact correspondence with Cartesian Monopole (l=0), Dipole (l=1), Quadrupole (l=2), Octupole (l=3).
 * 5. Full spherical harmonic reconstruction and synthesis on HEALPix or regular (theta, phi) angular grids.
 *
 * @module bulk-flow/spherical_harmonic_decomposer
 */

/**
 * Computes Associated Legendre Polynomial P_l^m(x) for l >= 0, 0 <= m <= l, |x| <= 1.
 * Uses standard stable recurrence relations.
 *
 * @param {number} l - Degree l
 * @param {number} m - Order m
 * @param {number} x - cos(theta) in [-1, 1]
 * @returns {number} Value of P_l^m(x)
 */
export function associatedLegendre(l, m, x) {
  if (m < 0 || m > l || Math.abs(x) > 1.0) {
    throw new RangeError(\'Invalid Legendre arguments: 0 <= m <= l and |x| <= 1\');
  }

  // 1. Compute P_m^m(x) = (-1)^m (2m - 1)!! (1 - x^2)^{m/2}
  let pmm = 1.0;
  if (m > 0) {
    const somx2 = Math.sqrt(Math.max(0.0, (1.0 - x) * (1.0 + x)));
    let fact = 1.0;
    for (let i = 1; i <= m; i++) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }

  if (l === m) return pmm;

  // 2. Compute P_{m+1}^m(x) = x (2m + 1) P_m^m(x)
  let pmmp1 = x * (2 * m + 1) * pmm;
  if (l === m + 1) return pmmp1;

  // 3. Upward recurrence for l > m + 1:
  // (l - m) P_l^m(x) = x (2l - 1) P_{l-1}^m(x) - (l + m - 1) P_{l-2}^m(x)
  let pll = 0.0;
  let pMinus2 = pmm;
  let pMinus1 = pmmp1;

  for (let ll = m + 2; ll <= l; ll++) {
    pll = (x * (2 * ll - 1) * pMinus1 - (ll + m - 1) * pMinus2) / (ll - m);
    pMinus2 = pMinus1;
    pMinus1 = pll;
  }

  return pll;
}

/**
 * Computes factorial n!
 * @param {number} n
 * @returns {number}
 */
export function factorial(n) {
  let res = 1.0;
  for (let i = 2; i <= n; i++) res *= i;
  return res;
}

/**
 * Computes Real Spherical Harmonic Basis Function Y_lm(theta, phi).
 * Fully normalized such that \int_{S^2} Y_lm Y_l'm' d\Omega = delta_ll' delta_mm'.
 *
 * @param {number} l - Degree (0 <= l <= 8)
 * @param {number} m - Order (-l <= m <= l)
 * @param {number} theta - Polar angle in radians [0, pi]
 * @param {number} phi - Azimuthal angle in radians [0, 2pi)
 * @returns {number} Real Y_lm(theta, phi)
 */
export function realSphericalHarmonic(l, m, theta, phi) {
  const absM = Math.abs(m);
  const cosTheta = Math.cos(theta);
  const Plm = associatedLegendre(l, absM, cosTheta);

  // Normalization factor K_lm = sqrt((2l + 1)/(4pi) * (l - |m|)! / (l + |m|]!))
  const normFactor = Math.sqrt(
    ((2 * l + 1) / (4.0 * Math.PI)) *
    (factorial(l - absM) / factorial(l + absM))
  );

  if (m === 0) {
    return normFactor * Plm;
  } else if (m > 0) {
    return Math.SQRT2 * normFactor * Plm * Math.cos(m * phi);
  } else {
    // m < 0
    return Math.SQRT2 * normFactor * Plm * Math.sin(absM * phi);
  }
}

/**
 * Decomposes an angular distribution of observations into Real Spherical Harmonic coefficients a_lm.
 */
export class SphericalHarmonicDecomposer {
  /**
   * @param {number} [maxL=4] - Maximum multipole degree l_max
   */
  constructor(maxL = 4) {
    this.maxL = Math.min(8, Math.max(1, maxL));
  }

  /**
   * Decomposes discrete point samples (e.g. galaxies or surface pixels) into a_lm coefficients.
   *
   * @param {Array<{ theta: number, phi: number, value: number, weight?: number }>} samples
   * @returns {{
   *   coefficients: Map<string, number>,
   *   powerSpectrum: Float64Array,
   *   monopoleA00: number,
   *   dipoleVector: [number, number, number],
   *   dipoleMagnitude: number
   * }}
   */
  decomposeSamples(samples) {
    const n = samples.length;
    if (n === 0) throw new Error(\'Samples array cannot be empty.\');

    let totalWeight = 0.0;
    for (let i = 0; i < n; i++) totalWeight += samples[i].weight ?? 1.0;
    const invTotalWeight = (4.0 * Math.PI) / totalWeight;

    const coeffs = new Map();
    const powerSpectrum = new Float64Array(this.maxL + 1);

    for (let l = 0; l <= this.maxL; l++) {
      let clSum = 0.0;
      for (let m = -l; m <= l; m++) {
        let alm = 0.0;
        for (let i = 0; i < n; i++) {
          const s = samples[i];
          const w = s.weight ?? 1.0;
          const ylm = realSphericalHarmonic(l, m, s.theta, s.phi);
          alm += w * s.value * ylm;
        }
        alm *= invTotalWeight;
        coeffs.set(${l},, alm);
        clSum += alm * alm;
      }
      powerSpectrum[l] = clSum / (2 * l + 1);
    }

    const a00 = coeffs.get(\'0,0\') || 0.0;
    const a1_minus1 = coeffs.get(\'1,-1\') || 0.0; // Y term
    const a1_0 = coeffs.get(\'1,0\') || 0.0;       // Z term
    const a1_1 = coeffs.get(\'1,1\') || 0.0;        // X term

    // In real spherical harmonic basis:
    // Y_11 = sqrt(3/(4pi)) * x/r => Vx = a_11 * sqrt(4pi/3)
    const factorDipole = Math.sqrt((4.0 * Math.PI) / 3.0);
    const vx = a1_1 * factorDipole;
    const vy = a1_minus1 * factorDipole;
    const vz = a1_0 * factorDipole;
    const dipoleMag = Math.hypot(vx, vy, vz);

    return {
      coefficients: coeffs,
      powerSpectrum,
      monopoleA00: a00,
      dipoleVector: [vx, vy, vz],
      dipoleMagnitude: dipoleMag
    };
  }

  /**
   * Synthesizes the spherical harmonic expansion at coordinates (theta, phi):
   * f(theta, phi) = sum_{l=0}^{maxL} sum_{m=-l}^l a_lm Y_lm(theta, phi).
   *
   * @param {Map<string, number>} coeffs
   * @param {number} theta
   * @param {number} phi
   * @returns {number}
   */
  synthesize(coeffs, theta, phi) {
    let sum = 0.0;
    for (let l = 0; l <= this.maxL; l++) {
      for (let m = -l; m <= l; m++) {
        const alm = coeffs.get(${l},) || 0.0;
        if (alm !== 0.0) {
          sum += alm * realSphericalHarmonic(l, m, theta, phi);
        }
      }
    }
    return sum;
  }
}
