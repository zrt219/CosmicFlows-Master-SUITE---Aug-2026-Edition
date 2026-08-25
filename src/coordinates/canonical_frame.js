/**
 * @fileoverview Canonical Astrometric Coordinate Frames and High-Performance TypedArray Transformations.
 * 
 * Provides exact IAU / de Vaucouleurs transformation matrices and conversions between:
 * - Supergalactic Cartesian (SGX, SGY, SGZ) & Spherical (SGL, SGB, d)
 * - Galactic Cartesian (X_gal, Y_gal, Z_gal) & Spherical (l, b, d)
 * - Equatorial Cartesian (X_eq, Y_eq, Z_eq) & Spherical (RA, Dec, cz / d)
 * 
 * Includes vectorized batch transformations with arbitrary array strides (e.g. interleaved Float32Array
 * position/velocity buffers for WebGL / GPU shaders).
 * 
 * @module coordinates/canonical_frame
 */

const DEG2RAD = Math.PI / 180.0;
const RAD2DEG = 180.0 / Math.PI;

/**
 * Standard IAU J2000 Astrometric Constants
 */
export const ASTROMETRIC_CONSTANTS = Object.freeze({
  // Galactic North Pole (J2000)
  RA_NGP_DEG: 192.85948,   // 12h 51m 26.275s
  DEC_NGP_DEG: 27.12825,   // +27° 07' 41.70"
  L0_GAL_DEG: 32.93192,    // Galactic center ascending node position angle

  // Galactic Center (J2000)
  RA_GC_DEG: 266.40510,    // 17h 45m 37.224s
  DEC_GC_DEG: -28.93617,   // -28° 56' 10.21"

  // Supergalactic North Pole in Galactic coordinates (de Vaucouleurs et al. 1976, 1991 RC3)
  L_SGP_DEG: 47.37,
  B_SGP_DEG: 6.32,

  // Supergalactic Plane / Galactic Plane ascending node intersection (SGL=0, SGB=0 in Galactic)
  L_SG0_DEG: 137.37,
  B_SG0_DEG: 0.0,

  // Default CF4 Hubble Constant (km/s/Mpc)
  DEFAULT_H0: 74.6
});

/**
 * Helper to multiply two 3x3 matrices (row-major flat arrays of length 9).
 * @param {Float64Array|number[]} A
 * @param {Float64Array|number[]} B
 * @returns {Float64Array} A * B
 */
function matMul3x3(A, B) {
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
 * Transpose a 3x3 matrix.
 * @param {Float64Array|number[]} M
 * @returns {Float64Array} M^T
 */
function matTranspose3x3(M) {
  const T = new Float64Array(9);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      T[j * 3 + i] = M[i * 3 + j];
    }
  }
  return T;
}

/**
 * Constructs Galactic to Equatorial (J2000) 3x3 Rotation Matrix.
 * Derived from standard IAU NGP angles:
 * alpha_GP = 192.85948°, delta_GP = 27.12825°, l_0 = 32.93192°
 * 
 * @returns {Float64Array} Row-major 3x3 matrix R_Gal2Eq such that X_eq = R_Gal2Eq * X_gal
 */
function buildGalacticToEquatorialMatrix() {
  const ra_gp = ASTROMETRIC_CONSTANTS.RA_NGP_DEG * DEG2RAD;
  const dec_gp = ASTROMETRIC_CONSTANTS.DEC_NGP_DEG * DEG2RAD;
  const ra_gc = ASTROMETRIC_CONSTANTS.RA_GC_DEG * DEG2RAD;
  const dec_gc = ASTROMETRIC_CONSTANTS.DEC_GC_DEG * DEG2RAD;

  // Z_gal axis in Equatorial coordinates (North Galactic Pole)
  const z_x = Math.cos(dec_gp) * Math.cos(ra_gp);
  const z_y = Math.cos(dec_gp) * Math.sin(ra_gp);
  const z_z = Math.sin(dec_gp);

  // X_gal target in Equatorial coordinates (Galactic Center)
  const x_raw_x = Math.cos(dec_gc) * Math.cos(ra_gc);
  const x_raw_y = Math.cos(dec_gc) * Math.sin(ra_gc);
  const x_raw_z = Math.sin(dec_gc);

  // Y_gal axis = Z_gal cross X_gal_raw
  const y_x = z_y * x_raw_z - z_z * x_raw_y;
  const y_y = z_z * x_raw_x - z_x * x_raw_z;
  const y_z = z_x * x_raw_y - z_y * x_raw_x;
  const y_norm = Math.hypot(y_x, y_y, y_z);

  const ny_x = y_x / y_norm;
  const ny_y = y_y / y_norm;
  const ny_z = y_z / y_norm;

  // X_gal axis = Y_gal cross Z_gal (strictly orthogonal and normalized)
  const nx_x = ny_y * z_z - ny_z * z_y;
  const nx_y = ny_z * z_x - ny_x * z_z;
  const nx_z = ny_x * z_y - ny_y * z_x;

  // R_Gal2Eq columns are nx, ny, z
  const R = new Float64Array([
    nx_x, ny_x, z_x,
    nx_y, ny_y, z_y,
    nx_z, ny_z, z_z
  ]);

  return R;
}

/**
 * Constructs Supergalactic to Galactic 3x3 Rotation Matrix.
 * Derived from standard de Vaucouleurs (1976 / 1991 RC3) definition:
 * SGP: (l=47.37°, b=6.32°), SG Origin: (l=137.37°, b=0.0°)
 * 
 * @returns {Float64Array} Row-major 3x3 matrix R_SG2Gal such that X_gal = R_SG2Gal * X_sg
 */
function buildSupergalacticToGalacticMatrix() {
  const l_sgp = ASTROMETRIC_CONSTANTS.L_SGP_DEG * DEG2RAD;
  const b_sgp = ASTROMETRIC_CONSTANTS.B_SGP_DEG * DEG2RAD;
  const l_sg0 = ASTROMETRIC_CONSTANTS.L_SG0_DEG * DEG2RAD;

  // Z_sg axis in Galactic Cartesian coordinates (Supergalactic North Pole)
  const z_x = Math.cos(b_sgp) * Math.cos(l_sgp);
  const z_y = Math.cos(b_sgp) * Math.sin(l_sgp);
  const z_z = Math.sin(b_sgp);

  // X_sg axis in Galactic Cartesian coordinates (SGL=0, SGB=0 origin)
  const x_x = Math.cos(l_sg0);
  const x_y = Math.sin(l_sg0);
  const x_z = 0.0;

  // Y_sg axis = Z_sg cross X_sg
  const y_x = z_y * x_z - z_z * x_y;
  const y_y = z_z * x_x - z_x * x_z;
  const y_z = z_x * x_y - z_y * x_x;

  // Normalize Y_sg to guarantee strict orthonormality
  const y_norm = Math.hypot(y_x, y_y, y_z);

  const R = new Float64Array([
    x_x, y_x / y_norm, z_x,
    x_y, y_y / y_norm, z_y,
    x_z, y_z / y_norm, z_z
  ]);

  return R;
}

// Precalculated Standard Rotation Matrices
export const ROT_GAL_TO_EQ = buildGalacticToEquatorialMatrix();
export const ROT_EQ_TO_GAL = matTranspose3x3(ROT_GAL_TO_EQ);

export const ROT_SG_TO_GAL = buildSupergalacticToGalacticMatrix();
export const ROT_GAL_TO_SG = matTranspose3x3(ROT_SG_TO_GAL);

// Direct Supergalactic <-> Equatorial transformations
export const ROT_SG_TO_EQ = matMul3x3(ROT_GAL_TO_EQ, ROT_SG_TO_GAL);
export const ROT_EQ_TO_SG = matTranspose3x3(ROT_SG_TO_EQ);

/**
 * Spherical coordinates to Cartesian 3D vector.
 * 
 * @param {number} lonDeg - Longitude / Right Ascension (degrees)
 * @param {number} latDeg - Latitude / Declination (degrees)
 * @param {number} [r=1.0] - Radial distance
 * @returns {[number, number, number]} Cartesian coordinates [x, y, z]
 */
export function sphericalToCartesian(lonDeg, latDeg, r = 1.0) {
  const phi = lonDeg * DEG2RAD;
  const theta = latDeg * DEG2RAD;
  const cosTheta = Math.cos(theta);

  const x = r * cosTheta * Math.cos(phi);
  const y = r * cosTheta * Math.sin(phi);
  const z = r * Math.sin(theta);
  return [x, y, z];
}

/**
 * Cartesian 3D vector to Spherical coordinates with robust pole singularity handling.
 * 
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @returns {[number, number, number]} [lonDeg (0 to 360), latDeg (-90 to +90), r]
 */
export function cartesianToSpherical(x, y, z) {
  const r = Math.hypot(x, y, z);
  if (r === 0.0) {
    return [0.0, 0.0, 0.0];
  }

  // Clamped latitude to [-1, 1] for asin stability
  const sinLat = Math.max(-1.0, Math.min(1.0, z / r));
  const lat = Math.asin(sinLat) * RAD2DEG;

  // Longitude in [0, 360)
  let lon = Math.atan2(y, x) * RAD2DEG;
  if (lon < 0.0) {
    lon += 360.0;
  }
  if (lon >= 360.0) {
    lon -= 360.0;
  }

  return [lon, lat, r];
}

/**
 * Multiplies a 3x3 matrix by a 3D vector.
 * @param {Float64Array|number[]} M
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @returns {[number, number, number]} Transformed vector
 */
export function transformPoint3D(M, x, y, z) {
  return [
    M[0] * x + M[1] * y + M[2] * z,
    M[3] * x + M[4] * y + M[5] * z,
    M[6] * x + M[7] * y + M[8] * z
  ];
}

/**
 * Converts Equatorial Spherical (RA, Dec, distance) to Galactic Spherical (l, b, distance).
 * 
 * @param {number} raDeg - Right Ascension (0 to 360 deg)
 * @param {number} decDeg - Declination (-90 to +90 deg)
 * @param {number} [dist=1.0] - Radial distance (Mpc/h or km/s)
 * @returns {[number, number, number]} [l (deg), b (deg), dist]
 */
export function equatorialToGalactic(raDeg, decDeg, dist = 1.0) {
  const [xEq, yEq, zEq] = sphericalToCartesian(raDeg, decDeg, dist);
  const [xGal, yGal, zGal] = transformPoint3D(ROT_EQ_TO_GAL, xEq, yEq, zEq);
  return cartesianToSpherical(xGal, yGal, zGal);
}

/**
 * Converts Galactic Spherical (l, b, distance) to Equatorial Spherical (RA, Dec, distance).
 * 
 * @param {number} lDeg - Galactic Longitude (0 to 360 deg)
 * @param {number} bDeg - Galactic Latitude (-90 to +90 deg)
 * @param {number} [dist=1.0] - Radial distance
 * @returns {[number, number, number]} [ra (deg), dec (deg), dist]
 */
export function galacticToEquatorial(lDeg, bDeg, dist = 1.0) {
  const [xGal, yGal, zGal] = sphericalToCartesian(lDeg, bDeg, dist);
  const [xEq, yEq, zEq] = transformPoint3D(ROT_GAL_TO_EQ, xGal, yGal, zGal);
  return cartesianToSpherical(xEq, yEq, zEq);
}

/**
 * Converts Galactic Spherical (l, b, distance) to Supergalactic Spherical (SGL, SGB, distance).
 * 
 * @param {number} lDeg - Galactic Longitude (0 to 360 deg)
 * @param {number} bDeg - Galactic Latitude (-90 to +90 deg)
 * @param {number} [dist=1.0] - Radial distance
 * @returns {[number, number, number]} [sgl (deg), sgb (deg), dist]
 */
export function galacticToSupergalactic(lDeg, bDeg, dist = 1.0) {
  const [xGal, yGal, zGal] = sphericalToCartesian(lDeg, bDeg, dist);
  const [sgx, sgy, sgz] = transformPoint3D(ROT_GAL_TO_SG, xGal, yGal, zGal);
  return cartesianToSpherical(sgx, sgy, sgz);
}

/**
 * Converts Supergalactic Spherical (SGL, SGB, distance) to Galactic Spherical (l, b, distance).
 * 
 * @param {number} sglDeg - Supergalactic Longitude (0 to 360 deg)
 * @param {number} sgbDeg - Supergalactic Latitude (-90 to +90 deg)
 * @param {number} [dist=1.0] - Radial distance
 * @returns {[number, number, number]} [l (deg), b (deg), dist]
 */
export function supergalacticToGalactic(sglDeg, sgbDeg, dist = 1.0) {
  const [sgx, sgy, sgz] = sphericalToCartesian(sglDeg, sgbDeg, dist);
  const [xGal, yGal, zGal] = transformPoint3D(ROT_SG_TO_GAL, sgx, sgy, sgz);
  return cartesianToSpherical(xGal, yGal, zGal);
}

/**
 * Converts Equatorial Spherical (RA, Dec, distance) directly to Supergalactic Cartesian (SGX, SGY, SGZ).
 * 
 * @param {number} raDeg - Right Ascension (deg)
 * @param {number} decDeg - Declination (deg)
 * @param {number} [dist=1.0] - Radial distance (Mpc/h)
 * @returns {[number, number, number]} Supergalactic Cartesian coordinates [SGX, SGY, SGZ]
 */
export function equatorialToSupergalacticCartesian(raDeg, decDeg, dist = 1.0) {
  const [xEq, yEq, zEq] = sphericalToCartesian(raDeg, decDeg, dist);
  return transformPoint3D(ROT_EQ_TO_SG, xEq, yEq, zEq);
}

/**
 * Converts Supergalactic Cartesian (SGX, SGY, SGZ) to Equatorial Spherical (RA, Dec, distance).
 * 
 * @param {number} sgx - Supergalactic X (Mpc/h)
 * @param {number} sgy - Supergalactic Y (Mpc/h)
 * @param {number} sgz - Supergalactic Z (Mpc/h)
 * @returns {[number, number, number]} Equatorial coordinates [RA (deg), Dec (deg), distance (Mpc/h)]
 */
export function supergalacticCartesianToEquatorial(sgx, sgy, sgz) {
  const [xEq, yEq, zEq] = transformPoint3D(ROT_SG_TO_EQ, sgx, sgy, sgz);
  return cartesianToSpherical(xEq, yEq, zEq);
}

/**
 * Converts Galactic Cartesian (X_gal, Y_gal, Z_gal) to Supergalactic Cartesian (SGX, SGY, SGZ).
 * 
 * @param {number} xGal
 * @param {number} yGal
 * @param {number} zGal
 * @returns {[number, number, number]} [SGX, SGY, SGZ]
 */
export function galacticCartesianToSupergalacticCartesian(xGal, yGal, zGal) {
  return transformPoint3D(ROT_GAL_TO_SG, xGal, yGal, zGal);
}

/**
 * Converts Supergalactic Cartesian (SGX, SGY, SGZ) to Galactic Cartesian (X_gal, Y_gal, Z_gal).
 * 
 * @param {number} sgx
 * @param {number} sgy
 * @param {number} sgz
 * @returns {[number, number, number]} [X_gal, Y_gal, Z_gal]
 */
export function supergalacticCartesianToGalacticCartesian(sgx, sgy, sgz) {
  return transformPoint3D(ROT_SG_TO_GAL, sgx, sgy, sgz);
}

/**
 * Rotates a 3D vector (e.g. peculiar velocity vx, vy, vz) between coordinate frames.
 * 
 * @param {[number, number, number]} vec - Input 3D vector [vx, vy, vz]
 * @param {'equatorial'|'galactic'|'supergalactic'} fromFrame - Source coordinate frame
 * @param {'equatorial'|'galactic'|'supergalactic'} toFrame - Target coordinate frame
 * @returns {[number, number, number]} Rotated 3D vector
 */
export function transformVelocityVector(vec, fromFrame, toFrame) {
  if (fromFrame === toFrame) return [...vec];

  let rotMat;
  if (fromFrame === 'equatorial' && toFrame === 'galactic') rotMat = ROT_EQ_TO_GAL;
  else if (fromFrame === 'galactic' && toFrame === 'equatorial') rotMat = ROT_GAL_TO_EQ;
  else if (fromFrame === 'galactic' && toFrame === 'supergalactic') rotMat = ROT_GAL_TO_SG;
  else if (fromFrame === 'supergalactic' && toFrame === 'galactic') rotMat = ROT_SG_TO_GAL;
  else if (fromFrame === 'equatorial' && toFrame === 'supergalactic') rotMat = ROT_EQ_TO_SG;
  else if (fromFrame === 'supergalactic' && toFrame === 'equatorial') rotMat = ROT_SG_TO_EQ;
  else throw new Error(`Unsupported frame transform from "${fromFrame}" to "${toFrame}"`);

  return transformPoint3D(rotMat, vec[0], vec[1], vec[2]);
}

/**
 * High-performance batch Array Stride transformation for WebGL / Float32Array buffers.
 * Converts coordinate tuples in-place or out-of-place across coordinate systems.
 * 
 * Supported modes:
 * - 'spherical_eq_to_cart_sg': (RA, Dec, Dist) -> (SGX, SGY, SGZ)
 * - 'spherical_gal_to_cart_sg': (l, b, Dist) -> (SGX, SGY, SGZ)
 * - 'cart_eq_to_cart_sg': (X_eq, Y_eq, Z_eq) -> (SGX, SGY, SGZ)
 * - 'cart_sg_to_cart_eq': (SGX, SGY, SGZ) -> (X_eq, Y_eq, Z_eq)
 * - 'cart_gal_to_cart_sg': (X_gal, Y_gal, Z_gal) -> (SGX, SGY, SGZ)
 * - 'cart_sg_to_cart_gal': (SGX, SGY, SGZ) -> (X_gal, Y_gal, Z_gal)
 * 
 * @param {Float32Array|Float64Array|number[]} inputBuffer - Source buffer
 * @param {Float32Array|Float64Array|number[]} [outputBuffer=inputBuffer] - Destination buffer (can be same)
 * @param {string} mode - Transformation mode
 * @param {object} [options]
 * @param {number} [options.stride=3] - Stride in floats between consecutive coordinate tuples
 * @param {number} [options.offset=0] - Offset in floats to the first coordinate component
 * @param {number} [options.count] - Number of points to transform (default: buffer length / stride)
 * @returns {Float32Array|Float64Array|number[]} The outputBuffer
 */
export function transformArrayStride(inputBuffer, outputBuffer = inputBuffer, mode, options = {}) {
  const stride = options.stride ?? 3;
  const offset = options.offset ?? 0;
  const totalPoints = options.count ?? Math.floor((inputBuffer.length - offset) / stride);

  if (totalPoints <= 0) return outputBuffer;

  let rotMat = null;
  let isSphericalEq = false;
  let isSphericalGal = false;

  switch (mode) {
    case 'spherical_eq_to_cart_sg':
      isSphericalEq = true;
      rotMat = ROT_EQ_TO_SG;
      break;
    case 'spherical_gal_to_cart_sg':
      isSphericalGal = true;
      rotMat = ROT_GAL_TO_SG;
      break;
    case 'cart_eq_to_cart_sg':
      rotMat = ROT_EQ_TO_SG;
      break;
    case 'cart_sg_to_cart_eq':
      rotMat = ROT_SG_TO_EQ;
      break;
    case 'cart_gal_to_cart_sg':
      rotMat = ROT_GAL_TO_SG;
      break;
    case 'cart_sg_to_cart_gal':
      rotMat = ROT_SG_TO_GAL;
      break;
    default:
      throw new Error(`Unknown transformation mode: "${mode}"`);
  }

  const m00 = rotMat[0], m01 = rotMat[1], m02 = rotMat[2];
  const m10 = rotMat[3], m11 = rotMat[4], m12 = rotMat[5];
  const m20 = rotMat[6], m21 = rotMat[7], m22 = rotMat[8];

  for (let i = 0; i < totalPoints; i++) {
    const idx = offset + i * stride;
    let x = inputBuffer[idx];
    let y = inputBuffer[idx + 1];
    let z = inputBuffer[idx + 2];

    if (isSphericalEq || isSphericalGal) {
      // Input is lonDeg, latDeg, dist
      const phi = x * DEG2RAD;
      const theta = y * DEG2RAD;
      const r = z;
      const cosTheta = Math.cos(theta);

      x = r * cosTheta * Math.cos(phi);
      y = r * cosTheta * Math.sin(phi);
      z = r * Math.sin(theta);
    }

    // Apply rotation matrix
    const outX = m00 * x + m01 * y + m02 * z;
    const outY = m10 * x + m11 * y + m12 * z;
    const outZ = m20 * x + m21 * y + m22 * z;

    outputBuffer[idx] = outX;
    outputBuffer[idx + 1] = outY;
    outputBuffer[idx + 2] = outZ;
  }

  return outputBuffer;
}

/**
 * Converts recessional velocity (cz in km/s) to metric distance (d in Mpc/h) via linear Hubble flow.
 * 
 * @param {number} cz - Recessional velocity in km/s
 * @param {number} [h0=74.6] - Hubble parameter in km/s/Mpc
 * @param {number} [h=0.746] - Dimensionless h (H0 / 100)
 * @returns {number} Distance in Mpc/h: d = (cz / H0) * h = cz / 100
 */
export function recessionalVelocityToDistance(cz, h0 = ASTROMETRIC_CONSTANTS.DEFAULT_H0, h = 0.746) {
  // Metric distance in Mpc: d_Mpc = cz / H0
  // Metric distance in Mpc/h: d_Mpch = d_Mpc * h = (cz / H0) * (H0 / 100) = cz / 100.0
  return (cz / h0) * h;
}

/**
 * Converts metric distance (d in Mpc/h) to recessional velocity (cz in km/s).
 * 
 * @param {number} distMpch - Metric distance in Mpc/h
 * @param {number} [h0=74.6] - Hubble parameter
 * @param {number} [h=0.746] - Dimensionless h
 * @returns {number} cz in km/s
 */
export function distanceToRecessionalVelocity(distMpch, h0 = ASTROMETRIC_CONSTANTS.DEFAULT_H0, h = 0.746) {
  return (distMpch / h) * h0;
}
