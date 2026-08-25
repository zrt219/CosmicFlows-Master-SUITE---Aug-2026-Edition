/**
 * @file equatorial_galactic_transforms.js
 * @description Astrometric epoch precession (B1950 <-> J2000), sexagesimal angle parsing/formatting,
 * Galactic Zone of Avoidance (ZoA) spatial selection masks, and Tully-Fisher distance calculators
 * for the CosmicFlows-4 Research Workbench.
 * 
 * @module coordinates/equatorial_galactic_transforms
 */

import { ASTROMETRIC_CONSTANTS } from './canonical_frame.js';

const DEG2RAD = Math.PI / 180.0;
const RAD2DEG = 180.0 / Math.PI;

/**
 * Parses sexagesimal Right Ascension string ("12h 34m 56.7s" or "12:34:56.7") into decimal degrees [0, 360).
 * @param {string} raStr
 * @returns {number} RA in degrees
 */
export function parseRA(raStr) {
  if (typeof raStr !== 'string') {
    throw new TypeError(`parseRA: expected string, received ${typeof raStr}`);
  }
  const clean = raStr.trim().replace(/[hms]/gi, ':').replace(/\s+/g, '');
  const parts = clean.split(':').filter(p => p.length > 0).map(Number);
  if (parts.length === 0 || parts.some(isNaN)) {
    throw new Error(`parseRA: Failed to parse RA string '${raStr}'`);
  }
  const h = parts[0] || 0;
  const m = parts[1] || 0;
  const s = parts[2] || 0;
  return (h + m / 60.0 + s / 3600.0) * 15.0;
}

/**
 * Parses sexagesimal Declination string ("+27d 07m 41.7s", "+27:07:41.7", "-05 30 12") into decimal degrees [-90, +90].
 * @param {string} decStr
 * @returns {number} Dec in degrees
 */
export function parseDec(decStr) {
  if (typeof decStr !== 'string') {
    throw new TypeError(`parseDec: expected string, received ${typeof decStr}`);
  }
  const trimmed = decStr.trim();
  const sign = trimmed.startsWith('-') ? -1 : 1;
  const clean = trimmed.replace(/^[+-]/, '').replace(/[dms°'"]/gi, ':').replace(/\s+/g, ':');
  const parts = clean.split(':').filter(p => p.length > 0).map(Number);
  if (parts.length === 0 || parts.some(isNaN)) {
    throw new Error(`parseDec: Failed to parse Dec string '${decStr}'`);
  }
  const d = parts[0] || 0;
  const m = parts[1] || 0;
  const s = parts[2] || 0;
  return sign * (d + m / 60.0 + s / 3600.0);
}

/**
 * Formats decimal Right Ascension degrees into standard sexagesimal string "HH:MM:SS.ss".
 * @param {number} raDeg [0, 360)
 * @param {number} [decimalPlaces=2]
 * @returns {string}
 */
export function formatRA(raDeg, decimalPlaces = 2) {
  let deg = ((raDeg % 360.0) + 360.0) % 360.0;
  const totalHours = deg / 15.0;
  const h = Math.floor(totalHours);
  const totalMinutes = (totalHours - h) * 60.0;
  const m = Math.floor(totalMinutes);
  const s = (totalMinutes - m) * 60.0;

  const hStr = String(h).padStart(2, '0');
  const mStr = String(m).padStart(2, '0');
  const sStr = s.toFixed(decimalPlaces).padStart(3 + decimalPlaces, '0');
  return `${hStr}h ${mStr}m ${sStr}s`;
}

/**
 * Formats decimal Declination degrees into standard sexagesimal string "+-DD:MM:SS.s".
 * @param {number} decDeg [-90, +90]
 * @param {number} [decimalPlaces=1]
 * @returns {string}
 */
export function formatDec(decDeg, decimalPlaces = 1) {
  const sign = decDeg >= 0 ? '+' : '-';
  const absDec = Math.abs(decDeg);
  const d = Math.floor(absDec);
  const totalMinutes = (absDec - d) * 60.0;
  const m = Math.floor(totalMinutes);
  const s = (totalMinutes - m) * 60.0;

  const dStr = String(d).padStart(2, '0');
  const mStr = String(m).padStart(2, '0');
  const sStr = s.toFixed(decimalPlaces).padStart(2 + decimalPlaces, '0');
  return `${sign}${dStr}° ${mStr}' ${sStr}"`;
}

/**
 * Precesses Equatorial coordinates from B1950 to J2000 (standard Lieske / IAU 1976 matrix).
 * @param {number} raB1950 RA in degrees
 * @param {number} decB1950 Dec in degrees
 * @returns {{ raJ2000: number, decJ2000: number }}
 */
export function precessB1950ToJ2000(raB1950, decB1950) {
  // Transformation matrix for B1950 to J2000
  const alpha = raB1950 * DEG2RAD;
  const delta = decB1950 * DEG2RAD;

  const cosD = Math.cos(delta);
  const r0 = [
    cosD * Math.cos(alpha),
    cosD * Math.sin(alpha),
    Math.sin(delta)
  ];

  // Standard IAU B1950 -> J2000 precession matrix
  const P = [
    0.9999256794956877, -0.0111814828119656, -0.0048590038153592,
    0.0111814828119656,  0.9999374784934914, -0.0000271625947142,
    0.0048590038153592, -0.0000271702937440,  0.9999881946023742
  ];

  const r1 = [
    P[0] * r0[0] + P[1] * r0[1] + P[2] * r0[2],
    P[3] * r0[0] + P[4] * r0[1] + P[5] * r0[2],
    P[6] * r0[0] + P[7] * r0[1] + P[8] * r0[2]
  ];

  const norm = Math.sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
  const decRad = Math.asin(Math.max(-1.0, Math.min(1.0, r1[2] / norm)));
  let raRad = Math.atan2(r1[1], r1[0]);
  if (raRad < 0.0) raRad += 2.0 * Math.PI;

  return {
    raJ2000: raRad * RAD2DEG,
    decJ2000: decRad * RAD2DEG
  };
}

/**
 * Zone of Avoidance (ZoA) Galactic Latitude Filter & Extinction Mask.
 */
export class ZoneOfAvoidanceMask {
  /**
   * @param {number} [cutoffLatitudeDeg=10.0] Galactic latitude boundary |b| <= cutoff (deg)
   */
  constructor(cutoffLatitudeDeg = 10.0) {
    this.cutoffLatitudeDeg = cutoffLatitudeDeg;
  }

  /**
   * Tests whether a given Galactic position is masked by the ZoA.
   * @param {number} lDeg Galactic longitude
   * @param {number} bDeg Galactic latitude
   * @returns {boolean} True if inside the obscured Zone of Avoidance
   */
  isMasked(lDeg, bDeg) {
    return Math.abs(bDeg) <= this.cutoffLatitudeDeg;
  }

  /**
   * Computes an empirical dust extinction weight W_ext(b) in [0, 1] (0 = fully obscured, 1 = unextincted).
   * @param {number} bDeg Galactic latitude
   * @returns {number}
   */
  transmissionWeight(bDeg) {
    const absB = Math.abs(bDeg);
    if (absB >= this.cutoffLatitudeDeg) return 1.0;
    if (absB <= 2.0) return 0.0;
    // Smooth cosine transition
    const t = (absB - 2.0) / (this.cutoffLatitudeDeg - 2.0);
    return 0.5 * (1.0 - Math.cos(t * Math.PI));
  }
}

/**
 * Tully-Fisher (TF) Distance and Peculiar Velocity Solver.
 */
export class TullyFisherSolver {
  /**
   * @param {Object} [params]
   * @param {number} [params.zeroPoint=-21.5] Absolute magnitude zero point M_0
   * @param {number} [params.slope=-7.5] Tully-Fisher slope a in M = M_0 + a * (log10(W_mx) - 2.5)
   * @param {number} [params.H0=74.6] Hubble constant in km/s/Mpc
   */
  constructor(params = {}) {
    this.M0 = params.zeroPoint ?? -21.5;
    this.slope = params.slope ?? -7.5;
    this.H0 = params.H0 ?? 74.6;
  }

  /**
   * Computes absolute magnitude M from 21cm line width W_mx (km/s).
   * @param {number} wMx Line width in km/s (typically 50-600 km/s)
   * @returns {number} Absolute magnitude M
   */
  absoluteMagnitude(wMx) {
    if (wMx <= 0.0) throw new RangeError(`TullyFisherSolver: Invalid W_mx = ${wMx}`);
    const logW = Math.log10(wMx);
    return this.M0 + this.slope * (logW - 2.5);
  }

  /**
   * Computes distance modulus mu = m - M and luminosity distance D_L (Mpc) from apparent magnitude m and line width W_mx.
   * @param {number} apparentMag Apparent magnitude m
   * @param {number} wMx Line width in km/s
   * @returns {{ distanceModulus: number, distanceMpc: number, distanceMpcOverH: number }}
   */
  computeDistance(apparentMag, wMx) {
    const M = this.absoluteMagnitude(wMx);
    const mu = apparentMag - M;
    const distanceMpc = Math.pow(10.0, (mu - 25.0) / 5.0);
    const h = this.H0 / 100.0;
    const distanceMpcOverH = distanceMpc * h;

    return {
      distanceModulus: mu,
      distanceMpc,
      distanceMpcOverH
    };
  }

  /**
   * Computes peculiar velocity v_pec = cz - H0 * d (km/s).
   * @param {number} cz Observed recessional velocity in km/s (CMB frame)
   * @param {number} distanceMpc Luminosity/Hubble distance in Mpc
   * @returns {number} Peculiar velocity in km/s
   */
  peculiarVelocity(cz, distanceMpc) {
    const vHubble = this.H0 * distanceMpc;
    return cz - vHubble;
  }
}
