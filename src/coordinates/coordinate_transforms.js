/**
 * @file coordinate_transforms.js
 * @description Comprehensive astrometric coordinate frame transformations and 3D Jacobian kinematics
 * for the CosmicFlows-4 Research Workbench.
 * 
 * Supports exact transformations between:
 * 1. Supergalactic Cartesian (SGX, SGY, SGZ in Mpc/h) <-> Supergalactic Spherical (SGL, SGB, d)
 * 2. Galactic Cartesian (X_gal, Y_gal, Z_gal in Mpc/h) <-> Galactic Spherical (l, b, d)
 * 3. Equatorial Cartesian (X_eq, Y_eq, Z_eq in Mpc/h) <-> Equatorial Spherical (RA, Dec, cz / d)
 * 4. Inter-frame rotations:
 *    - Galactic <-> Equatorial (IAU J2000 definition: NGP alpha = 192.85948 deg, delta = 27.12825 deg, l0 = 32.93192 deg)
 *    - Supergalactic <-> Galactic (de Vaucouleurs 1976 / 1991 RC3: SGP l = 47.37 deg, b = 6.32 deg, SG0 l = 137.37 deg, b = 0 deg)
 *    - Supergalactic <-> Equatorial (Composite matrix R_SG_to_Eq = R_Gal_to_Eq * R_SG_to_Gal)
 * 5. Peculiar velocity rest frame adjustments:
 *    - Heliocentric recessional velocity v_helio
 *    - Local Group (LG) rest frame velocity v_LG (Karachentsev & Makarov 1996)
 *    - Cosmic Microwave Background (CMB) rest frame velocity v_CMB (Fixsen et al. 1996: v_apex = 369.0 km/s towards l = 263.99 deg, b = 48.26 deg)
 * 6. Vectorized batch transformations for large astronomical catalogs with Float64Array.
 * 
 * @module coordinates/coordinate_transforms
 */

import { ASTROMETRIC_CONSTANTS } from './canonical_frame.js';

const DEG2RAD = Math.PI / 180.0;
const RAD2DEG = 180.0 / Math.PI;

/**
 * Standard CMB Dipole parameters (Fixsen et al. 1996, Kogut et al. 1993, Planck 2018).
 * @readonly
 */
export const CMB_DIPOLE = Object.freeze({
  VELOCITY_KMS: 369.0,         // km/s
  GALACTIC_L_DEG: 263.99,      // deg
  GALACTIC_B_DEG: 48.26,       // deg
  RA_J2000_DEG: 168.01,        // deg (11h 12m)
  DEC_J2000_DEG: -7.05         // deg (-7° 03')
});

/**
 * Standard Local Group (LG) motion relative to Heliocentre (Karachentsev & Makarov 1996).
 * @readonly
 */
export const LOCAL_GROUP_MOTION = Object.freeze({
  VELOCITY_KMS: 316.0,         // km/s
  GALACTIC_L_DEG: 93.0,        // deg
  GALACTIC_B_DEG: -4.0         // deg
});

/**
 * Builds 3x3 orthonormal rotation matrix from Galactic to Equatorial coordinates (J2000).
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createGalacticToEquatorialMatrix() {
  const raNGP = ASTROMETRIC_CONSTANTS.RA_NGP_DEG * DEG2RAD;
  const decNGP = ASTROMETRIC_CONSTANTS.DEC_NGP_DEG * DEG2RAD;
  const l0 = ASTROMETRIC_CONSTANTS.L0_GAL_DEG * DEG2RAD;

  const sinDec = Math.sin(decNGP);
  const cosDec = Math.cos(decNGP);
  const sinRA = Math.sin(raNGP);
  const cosRA = Math.cos(raNGP);
  const sinL0 = Math.sin(l0);
  const cosL0 = Math.cos(l0);

  // Basis vectors of Galactic frame in Equatorial coordinates
  // ex_gal (towards l=0, b=0):
  const ex = [
    -sinDec * cosRA * sinL0 - sinRA * cosL0,
    -sinDec * sinRA * sinL0 + cosRA * cosL0,
    cosDec * sinL0
  ];

  // ey_gal (towards l=90, b=0):
  const ey = [
    -sinDec * cosRA * cosL0 + sinRA * sinL0,
    -sinDec * sinRA * cosL0 - cosRA * sinL0,
    cosDec * cosL0
  ];

  // ez_gal (towards NGP, b=90):
  const ez = [
    cosDec * cosRA,
    cosDec * sinRA,
    sinDec
  ];

  // Matrix where columns are ex, ey, ez
  return new Float64Array([
    ex[0], ey[0], ez[0],
    ex[1], ey[1], ez[1],
    ex[2], ey[2], ez[2]
  ]);
}

/**
 * Builds 3x3 orthonormal rotation matrix from Supergalactic to Galactic coordinates.
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createSupergalacticToGalacticMatrix() {
  const lSGP = ASTROMETRIC_CONSTANTS.L_SGP_DEG * DEG2RAD;
  const bSGP = ASTROMETRIC_CONSTANTS.B_SGP_DEG * DEG2RAD;
  const lSG0 = ASTROMETRIC_CONSTANTS.L_SG0_DEG * DEG2RAD;

  const sinB = Math.sin(bSGP);
  const cosB = Math.cos(bSGP);
  const sinL = Math.sin(lSGP);
  const cosL = Math.cos(lSGP);
  const sinL0 = Math.sin(lSG0);
  const cosL0 = Math.cos(lSG0);

  // Supergalactic North Pole in Galactic coordinates ez_sg
  const ez_sg = [cosB * cosL, cosB * sinL, sinB];

  // Supergalactic X-axis (SGL=0, SGB=0 in Galactic) ex_sg
  const ex_sg = [cosL0, sinL0, 0.0];

  // Supergalactic Y-axis (ey_sg = ez_sg x ex_sg)
  const ey_sg = [
    ez_sg[1] * ex_sg[2] - ez_sg[2] * ex_sg[1],
    ez_sg[2] * ex_sg[0] - ez_sg[0] * ex_sg[2],
    ez_sg[0] * ex_sg[1] - ez_sg[1] * ex_sg[0]
  ];

  return new Float64Array([
    ex_sg[0], ey_sg[0], ez_sg[0],
    ex_sg[1], ey_sg[1], ez_sg[1],
    ex_sg[2], ey_sg[2], ez_sg[2]
  ]);
}

/**
 * Transposes a 3x3 matrix (which is its exact inverse for orthogonal matrices).
 * @param {Float64Array} M
 * @returns {Float64Array}
 */
export function transpose3x3(M) {
  return new Float64Array([
    M[0], M[3], M[6],
    M[1], M[4], M[7],
    M[2], M[5], M[8]
  ]);
}

/**
 * Multiplies two 3x3 matrices: C = A * B.
 * @param {Float64Array} A
 * @param {Float64Array} B
 * @returns {Float64Array}
 */
export function multiply3x3(A, B) {
  const C = new Float64Array(9);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let sum = 0;
      for (let k = 0; k < 3; k++) {
        sum += A[i * 3 + k] * B[k * 3 + j];
      }
      C[i * 3 + j] = sum;
    }
  }
  return C;
}

/**
 * Multiplies a 3x3 matrix by a 3D vector: y = M * x.
 * @param {Float64Array} M
 * @param {Array<number>|Float64Array} v
 * @returns {Float64Array}
 */
export function applyMatrix3x3(M, v) {
  const x = v[0], y = v[1], z = v[2];
  return new Float64Array([
    M[0] * x + M[1] * y + M[2] * z,
    M[3] * x + M[4] * y + M[5] * z,
    M[6] * x + M[7] * y + M[8] * z
  ]);
}

// Pre-computed static orthonormal transformation matrices
export const R_GAL_TO_EQ = createGalacticToEquatorialMatrix();
export const R_EQ_TO_GAL = transpose3x3(R_GAL_TO_EQ);

export const R_SG_TO_GAL = createSupergalacticToGalacticMatrix();
export const R_GAL_TO_SG = transpose3x3(R_SG_TO_GAL);

export const R_SG_TO_EQ = multiply3x3(R_GAL_TO_EQ, R_SG_TO_GAL);
export const R_EQ_TO_SG = transpose3x3(R_SG_TO_EQ);

/**
 * Converts Spherical angles (longitude, latitude, distance) to Cartesian 3D (X, Y, Z).
 * @param {number} lonDeg Longitude in degrees [0, 360)
 * @param {number} latDeg Latitude in degrees [-90, +90]
 * @param {number} dist Radial distance
 * @returns {{ x: number, y: number, z: number }}
 */
export function sphericalToCartesian(lonDeg, latDeg, dist = 1.0) {
  const phi = lonDeg * DEG2RAD;
  const theta = latDeg * DEG2RAD;
  const cosTheta = Math.cos(theta);
  return {
    x: dist * cosTheta * Math.cos(phi),
    y: dist * cosTheta * Math.sin(phi),
    z: dist * Math.sin(theta)
  };
}

/**
 * Converts Cartesian 3D (X, Y, Z) to Spherical angles (longitude, latitude, distance).
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @returns {{ lonDeg: number, latDeg: number, dist: number }}
 */
export function cartesianToSpherical(x, y, z) {
  const r2 = x * x + y * y + z * z;
  const dist = Math.sqrt(r2);
  if (dist < 1e-12) {
    return { lonDeg: 0.0, latDeg: 0.0, dist: 0.0 };
  }

  const latRad = Math.asin(Math.max(-1.0, Math.min(1.0, z / dist)));
  let lonRad = Math.atan2(y, x);
  if (lonRad < 0.0) lonRad += 2.0 * Math.PI;

  return {
    lonDeg: lonRad * RAD2DEG,
    latDeg: latRad * RAD2DEG,
    dist
  };
}

/**
 * CoordinateTransforms provides high-level coordinate transformations.
 */
export class CoordinateTransforms {
  /**
   * Supergalactic Cartesian (SGX, SGY, SGZ) to Galactic Cartesian.
   * @param {number} sgx
   * @param {number} sgy
   * @param {number} sgz
   * @returns {{ x: number, y: number, z: number }}
   */
  static supergalacticToGalactic(sgx, sgy, sgz) {
    const v = applyMatrix3x3(R_SG_TO_GAL, [sgx, sgy, sgz]);
    return { x: v[0], y: v[1], z: v[2] };
  }

  /**
   * Galactic Cartesian to Supergalactic Cartesian.
   * @param {number} gx
   * @param {number} gy
   * @param {number} gz
   * @returns {{ sgx: number, sgy: number, sgz: number }}
   */
  static galacticToSupergalactic(gx, gy, gz) {
    const v = applyMatrix3x3(R_GAL_TO_SG, [gx, gy, gz]);
    return { sgx: v[0], sgy: v[1], sgz: v[2] };
  }

  /**
   * Supergalactic Cartesian (SGX, SGY, SGZ) to Equatorial Cartesian (J2000).
   * @param {number} sgx
   * @param {number} sgy
   * @param {number} sgz
   * @returns {{ x: number, y: number, z: number }}
   */
  static supergalacticToEquatorial(sgx, sgy, sgz) {
    const v = applyMatrix3x3(R_SG_TO_EQ, [sgx, sgy, sgz]);
    return { x: v[0], y: v[1], z: v[2] };
  }

  /**
   * Equatorial Cartesian (J2000) to Supergalactic Cartesian.
   * @param {number} eqx
   * @param {number} eqy
   * @param {number} eqz
   * @returns {{ sgx: number, sgy: number, sgz: number }}
   */
  static equatorialToSupergalactic(eqx, eqy, eqz) {
    const v = applyMatrix3x3(R_EQ_TO_SG, [eqx, eqy, eqz]);
    return { sgx: v[0], sgy: v[1], sgz: v[2] };
  }

  /**
   * Converts Equatorial celestial coordinates (RA, Dec, distance) to Supergalactic Cartesian (SGX, SGY, SGZ).
   * @param {number} raDeg Right Ascension in degrees [0, 360)
   * @param {number} decDeg Declination in degrees [-90, +90]
   * @param {number} dist Radial distance in Mpc/h
   * @returns {{ sgx: number, sgy: number, sgz: number }}
   */
  static equatorialSphericalToSupergalactic(raDeg, decDeg, dist) {
    const { x, y, z } = sphericalToCartesian(raDeg, decDeg, dist);
    return this.equatorialToSupergalactic(x, y, z);
  }

  /**
   * Converts Supergalactic Cartesian (SGX, SGY, SGZ) to Equatorial celestial coordinates (RA, Dec, distance).
   * @param {number} sgx
   * @param {number} sgy
   * @param {number} sgz
   * @returns {{ raDeg: number, decDeg: number, dist: number }}
   */
  static supergalacticToEquatorialSpherical(sgx, sgy, sgz) {
    const { x, y, z } = this.supergalacticToEquatorial(sgx, sgy, sgz);
    const sph = cartesianToSpherical(x, y, z);
    return { raDeg: sph.lonDeg, decDeg: sph.latDeg, dist: sph.dist };
  }

  /**
   * Converts Galactic celestial coordinates (l, b, distance) to Supergalactic Cartesian (SGX, SGY, SGZ).
   * @param {number} lDeg Galactic longitude in degrees
   * @param {number} bDeg Galactic latitude in degrees
   * @param {number} dist Radial distance in Mpc/h
   * @returns {{ sgx: number, sgy: number, sgz: number }}
   */
  static galacticSphericalToSupergalactic(lDeg, bDeg, dist) {
    const { x, y, z } = sphericalToCartesian(lDeg, bDeg, dist);
    return this.galacticToSupergalactic(x, y, z);
  }

  /**
   * Converts Supergalactic Cartesian (SGX, SGY, SGZ) to Galactic celestial coordinates (l, b, distance).
   * @param {number} sgx
   * @param {number} sgy
   * @param {number} sgz
   * @returns {{ lDeg: number, bDeg: number, dist: number }}
   */
  static supergalacticToGalacticSpherical(sgx, sgy, sgz) {
    const { x, y, z } = this.supergalacticToGalactic(sgx, sgy, sgz);
    const sph = cartesianToSpherical(x, y, z);
    return { lDeg: sph.lonDeg, bDeg: sph.latDeg, dist: sph.dist };
  }

  /**
   * Adjusts heliocentric recessional velocity to CMB rest frame.
   * v_CMB = v_helio + v_apex * [sin(b) sin(b_apex) + cos(b) cos(b_apex) cos(l - l_apex)]
   * @param {number} vHelio Recessional velocity in km/s
   * @param {number} lDeg Galactic longitude in degrees
   * @param {number} bDeg Galactic latitude in degrees
   * @returns {number} v_CMB in km/s
   */
  static heliocentricToCmbVelocity(vHelio, lDeg, bDeg) {
    const l = lDeg * DEG2RAD;
    const b = bDeg * DEG2RAD;
    const lApex = CMB_DIPOLE.GALACTIC_L_DEG * DEG2RAD;
    const bApex = CMB_DIPOLE.GALACTIC_B_DEG * DEG2RAD;

    const cosTheta = Math.sin(b) * Math.sin(bApex) + Math.cos(b) * Math.cos(bApex) * Math.cos(l - lApex);
    return vHelio + CMB_DIPOLE.VELOCITY_KMS * cosTheta;
  }

  /**
   * Adjusts heliocentric recessional velocity to Local Group rest frame.
   * @param {number} vHelio Recessional velocity in km/s
   * @param {number} lDeg Galactic longitude in degrees
   * @param {number} bDeg Galactic latitude in degrees
   * @returns {number} v_LG in km/s
   */
  static heliocentricToLgVelocity(vHelio, lDeg, bDeg) {
    const l = lDeg * DEG2RAD;
    const b = bDeg * DEG2RAD;
    const lApex = LOCAL_GROUP_MOTION.GALACTIC_L_DEG * DEG2RAD;
    const bApex = LOCAL_GROUP_MOTION.GALACTIC_B_DEG * DEG2RAD;

    const cosTheta = Math.sin(b) * Math.sin(bApex) + Math.cos(b) * Math.cos(bApex) * Math.cos(l - lApex);
    return vHelio + LOCAL_GROUP_MOTION.VELOCITY_KMS * cosTheta;
  }

  /**
   * Vectorized batch transformation of Float64Array catalog coordinates from Galactic (l, b, d) to Supergalactic (SGX, SGY, SGZ).
   * @param {Float64Array} lArray
   * @param {Float64Array} bArray
   * @param {Float64Array} dArray
   * @param {Float64Array} [outBuffer] Optional output buffer of size 3*N
   * @returns {Float64Array} Flat interleaved [sgx0, sgy0, sgz0, sgx1, sgy1, sgz1, ...]
   */
  static batchGalacticToSupergalactic(lArray, bArray, dArray, outBuffer) {
    const n = lArray.length;
    const out = outBuffer || new Float64Array(n * 3);
    const M = R_GAL_TO_SG;

    for (let i = 0; i < n; i++) {
      const phi = lArray[i] * DEG2RAD;
      const theta = bArray[i] * DEG2RAD;
      const d = dArray[i];

      const cosTheta = Math.cos(theta);
      const gx = d * cosTheta * Math.cos(phi);
      const gy = d * cosTheta * Math.sin(phi);
      const gz = d * Math.sin(theta);

      const idx = i * 3;
      out[idx]     = M[0] * gx + M[1] * gy + M[2] * gz;
      out[idx + 1] = M[3] * gx + M[4] * gy + M[5] * gz;
      out[idx + 2] = M[6] * gx + M[7] * gy + M[8] * gz;
    }
    return out;
  }
}
