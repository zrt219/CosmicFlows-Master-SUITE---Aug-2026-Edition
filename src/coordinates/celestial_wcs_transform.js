/**
 * @file celestial_wcs_transform.js
 * @description High-Precision Astronomical Coordinate Frame Transformations, Redshift-Space Distortions (RSD / Kaiser Effect / FoG),
 * and Astrometric FITS World Coordinate System (WCS) Projection Parser for the CosmicFlows-4 Research Workbench.
 * 
 * Implements:
 * 1. High-Precision Astronomical Frame Conversions:
 *    - Supergalactic Cartesian (SGX, SGY, SGZ in Mpc/h or Mpc) <-> Galactic Spherical (l, b, d) & Cartesian (X, Y, Z)
 *    - Galactic (l, b, d) <-> ICRS Equatorial (RA, Dec, cz / d) [J2000 IAU Standards]
 *    - Supergalactic Longitude/Latitude (SGL, SGB) <-> Galactic & Equatorial
 *    - Heliocentric <-> CMB Rest Frame <-> Local Group Barycenter Frame conversions
 *      (Planck 2018 / Fixsen CMB Dipole: v_CMB = 369.8 km/s towards (l,b) = (264.021°, 48.253°))
 *    - Full 3D rotation matrices, Jacobians, and vectorized TypedArray stride transformations.
 * 
 * 2. Redshift-Space Distortions (RSD) & Kaiser / FoG Kinematics:
 *    - Linear Kaiser velocity mapping: s = x + (v_los / (a * H(z))) * \hat{r}
 *    - Non-linear dispersion damping (Lorentzian, Gaussian, Squared Lorentzian, Exponential)
 *    - Analytical & Quadrature Kaiser Power Spectrum Multipoles (P0, P2, P4)
 *    - Alcock-Paczynski (AP) geometric distortion factors (q_perp, q_para, F_AP)
 *    - Fingers-of-God (FoG) cluster identification, velocity dispersion estimation (Biweight/Gapper/Sigma-Clipped),
 *      and radial de-stretching / compression algorithms.
 *    - Iterative real-space coordinate and density field reconstruction.
 * 
 * 3. Astrometric FITS WCS CD/PC Matrix Projection Parser:
 *    - FITS Standard 4.0 compliant parser for CDi_j and PCi_j + CDELTi representations.
 *    - Full forward and inverse projections: TAN, SIN, ARC, STG, ZPN, CAR, CEA, MER, AIT, MOL, PAR.
 *    - SIP (Simple Imaging Polynomial) geometric distortion correction (HST/JWST/Spitzer).
 *    - Pixel-to-World (pix2world) and World-to-Pixel (world2pix) with LONPOLE/LATPOLE spherical rotations.
 *    - Footprint / sky polygon generator and pixel scale / FOV calculators.
 * 
 * @module coordinates/celestial_wcs_transform
 */

// ============================================================================
// 1. FUNDAMENTAL ASTROMETRIC, COSMOLOGICAL & PHYSICAL CONSTANTS
// ============================================================================

/** Degree to Radian conversion factor */
export const DEG2RAD = Math.PI / 180.0;
/** Radian to Degree conversion factor */
export const RAD2DEG = 180.0 / Math.PI;
/** Arcsecond to Degree conversion factor */
export const ARCSEC2DEG = 1.0 / 3600.0;
/** Degree to Arcsecond conversion factor */
export const DEG2ARCSEC = 3600.0;
/** Milliarcsecond to Degree conversion factor */
export const MAS2DEG = 1.0 / 3600000.0;
/** 2 * PI constant */
export const TWO_PI = 2.0 * Math.PI;
/** PI / 2 constant */
export const HALF_PI = Math.PI / 2.0;

/** Speed of light in vacuum (km/s) - CODATA 2018 */
export const SPEED_OF_LIGHT_KMS = 299792.458;

/**
 * Standard IAU J2000 Astrometric Orientation Constants
 */
export const ASTROMETRIC_CONSTANTS = Object.freeze({
  // North Galactic Pole (NGP) in J2000 ICRS Equatorial (IAU 1958 / 2000 definitions)
  RA_NGP_DEG: 192.85948,       // 12h 51m 26.2752s
  DEC_NGP_DEG: 27.12825,       // +27° 07' 41.70"
  L0_GAL_DEG: 32.93192,        // Galactic longitude of ascending node of galactic plane on equator

  // Galactic Center (Sgr A*) J2000
  RA_GC_DEG: 266.40510,        // 17h 45m 37.224s
  DEC_GC_DEG: -28.93617,       // -28° 56' 10.21"

  // Supergalactic North Pole (SGP) in Galactic coordinates (de Vaucouleurs et al. 1976, 1991 RC3)
  L_SGP_DEG: 47.37,
  B_SGP_DEG: 6.32,

  // Supergalactic Ascending Node (SG0) in Galactic coordinates (SGL=0°, SGB=0°)
  L_SG0_DEG: 137.37,
  B_SG0_DEG: 0.0,

  // Supergalactic North Pole in J2000 Equatorial coordinates (calculated)
  RA_SGP_DEG: 283.75429,
  DEC_SGP_DEG: 15.70034,

  // Default CosmicFlows-4 Hubble Constant
  DEFAULT_H0: 74.6,            // km/s/Mpc
  DEFAULT_H0_NORMALIZED: 100.0 // km/s / (Mpc/h)
});

/**
 * Velocity Reference Frame Dipole Parameters
 */
export const VELOCITY_REST_FRAMES = Object.freeze({
  // Cosmic Microwave Background (CMB) dipole (Planck 2018 / Fixsen et al. 1996)
  CMB_DIPOLE: Object.freeze({
    VELOCITY_KMS: 369.8,          // km/s
    GALACTIC_L_DEG: 264.021,      // deg
    GALACTIC_B_DEG: 48.253,       // deg
    RA_J2000_DEG: 168.01,         // deg (11h 12m 02s)
    DEC_J2000_DEG: -6.98          // deg (-6° 58' 48")
  }),

  // Local Group (LG) barycenter motion relative to Heliocentre (Courteau & van den Bergh 1999 / Karachentsev 1996)
  LOCAL_GROUP: Object.freeze({
    VELOCITY_KMS: 316.0,          // km/s
    GALACTIC_L_DEG: 93.0,         // deg
    GALACTIC_B_DEG: -4.0,         // deg
    RA_J2000_DEG: 326.69,         // deg
    DEC_J2000_DEG: 45.74          // deg
  }),

  // Galactic Standard of Rest (GSR) - Sun's orbital speed around Galactic Center
  GSR: Object.freeze({
    V_LSR_KMS: 220.0,             // km/s
    L_LSR_DEG: 90.0,              // deg
    B_LSR_DEG: 0.0,               // deg
    // Solar peculiar motion relative to LSR (Schönrich, Binney & Dehnen 2010): U=11.1, V=12.24, W=7.25 km/s
    U_PEC_KMS: 11.10,
    V_PEC_KMS: 12.24,
    W_PEC_KMS: 7.25
  })
});

/**
 * Standard Cosmological Parameters for RSD and Distance-Redshift Computations
 */
export const COSMOLOGY_DEFAULTS = Object.freeze({
  H0: 74.6,               // km/s/Mpc
  h: 0.746,               // dimensionless H0 / 100
  OMEGA_M: 0.315,         // Matter density parameter
  OMEGA_L: 0.685,         // Dark energy density parameter
  OMEGA_K: 0.0,           // Curvature density parameter (flat universe)
  OMEGA_R: 0.0,           // Radiation density parameter
  W0: -1.0,               // Dark energy equation of state parameter
  WA: 0.0,                // Dark energy CPL evolution parameter
  GROWTH_INDEX_GAMMA: 0.5454 // Linder growth index for LCDM
});


// ============================================================================
// 2. MATRIX & VECTOR LINEAR ALGEBRA PRIMITIVES (3D / 3x3)
// ============================================================================

/**
 * Multiplies two 3x3 row-major matrices: C = A * B.
 * @param {ArrayLike<number>} A - 9-element array
 * @param {ArrayLike<number>} B - 9-element array
 * @param {Float64Array} [out] - Optional destination array
 * @returns {Float64Array} 9-element row-major result
 */
export function matMul3x3(A, B, out = new Float64Array(9)) {
  for (let i = 0; i < 3; i++) {
    const rowOffset = i * 3;
    for (let j = 0; j < 3; j++) {
      out[rowOffset + j] =
        A[rowOffset] * B[j] +
        A[rowOffset + 1] * B[3 + j] +
        A[rowOffset + 2] * B[6 + j];
    }
  }
  return out;
}

/**
 * Multiplies a 3x3 matrix by a 3D column vector: u = M * v.
 * @param {ArrayLike<number>} M - 9-element row-major matrix
 * @param {ArrayLike<number>} v - 3-element vector [x, y, z]
 * @param {Float64Array} [out] - Optional destination array [x, y, z]
 * @returns {Float64Array} 3-element result vector
 */
export function matVecMul3x3(M, v, out = new Float64Array(3)) {
  const v0 = v[0], v1 = v[1], v2 = v[2];
  out[0] = M[0] * v0 + M[1] * v1 + M[2] * v2;
  out[1] = M[3] * v0 + M[4] * v1 + M[5] * v2;
  out[2] = M[6] * v0 + M[7] * v1 + M[8] * v2;
  return out;
}

/**
 * Computes transpose of a 3x3 matrix: T = M^T.
 * @param {ArrayLike<number>} M - 9-element row-major matrix
 * @param {Float64Array} [out] - Optional destination array
 * @returns {Float64Array} Transposed 3x3 matrix
 */
export function matTranspose3x3(M, out = new Float64Array(9)) {
  out[0] = M[0]; out[1] = M[3]; out[2] = M[6];
  out[3] = M[1]; out[4] = M[4]; out[5] = M[7];
  out[6] = M[2]; out[7] = M[5]; out[8] = M[8];
  return out;
}

/**
 * Computes the determinant of a 3x3 matrix.
 * @param {ArrayLike<number>} M - 9-element row-major matrix
 * @returns {number} Determinant det(M)
 */
export function matDet3x3(M) {
  return (
    M[0] * (M[4] * M[8] - M[5] * M[7]) -
    M[1] * (M[3] * M[8] - M[5] * M[6]) +
    M[2] * (M[3] * M[7] - M[4] * M[6])
  );
}

/**
 * Computes the trace of a 3x3 matrix.
 * @param {ArrayLike<number>} M - 9-element row-major matrix
 * @returns {number} Trace Tr(M)
 */
export function matTrace3x3(M) {
  return M[0] + M[4] + M[8];
}

/**
 * Computes exact inverse of a 3x3 matrix.
 * @param {ArrayLike<number>} M - 9-element row-major matrix
 * @param {Float64Array} [out] - Optional destination array
 * @returns {Float64Array} Inverse 3x3 matrix
 * @throws {Error} if matrix is singular (det === 0)
 */
export function matInverse3x3(M, out = new Float64Array(9)) {
  const det = matDet3x3(M);
  if (Math.abs(det) < 1e-15) {
    throw new Error(`matInverse3x3: Singular matrix encountered with det = ${det}`);
  }
  const invDet = 1.0 / det;

  out[0] = (M[4] * M[8] - M[5] * M[7]) * invDet;
  out[1] = (M[2] * M[7] - M[1] * M[8]) * invDet;
  out[2] = (M[1] * M[5] - M[2] * M[4]) * invDet;

  out[3] = (M[5] * M[6] - M[3] * M[8]) * invDet;
  out[4] = (M[0] * M[8] - M[2] * M[6]) * invDet;
  out[5] = (M[2] * M[3] - M[0] * M[5]) * invDet;

  out[6] = (M[3] * M[7] - M[4] * M[6]) * invDet;
  out[7] = (M[1] * M[6] - M[0] * M[7]) * invDet;
  out[8] = (M[0] * M[4] - M[1] * M[3]) * invDet;

  return out;
}

/**
 * Verifies if a 3x3 matrix is strictly orthogonal (M * M^T = I) and right-handed (det = +1).
 * @param {ArrayLike<number>} M - 9-element row-major matrix
 * @param {number} [tol=1e-7] - Maximum allowable tolerance
 * @returns {{isValid: boolean, maxOrthogonalityError: number, det: number}}
 */
export function verifyOrthonormality3x3(M, tol = 1e-7) {
  let maxErr = 0.0;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let sum = 0.0;
      for (let k = 0; k < 3; k++) {
        sum += M[i * 3 + k] * M[j * 3 + k];
      }
      const target = i === j ? 1.0 : 0.0;
      const err = Math.abs(sum - target);
      if (err > maxErr) maxErr = err;
    }
  }
  const det = matDet3x3(M);
  const detErr = Math.abs(det - 1.0);
  const isValid = maxErr <= tol && detErr <= tol;
  return { isValid, maxOrthogonalityError: maxErr, det };
}

/**
 * 3D Euclidean Vector Norm (magnitude).
 * @param {ArrayLike<number>} v
 * @returns {number}
 */
export function vec3Norm(v) {
  return Math.hypot(v[0], v[1], v[2]);
}

/**
 * 3D Dot Product.
 * @param {ArrayLike<number>} u
 * @param {ArrayLike<number>} v
 * @returns {number}
 */
export function vec3Dot(u, v) {
  return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

/**
 * 3D Cross Product: out = u x v.
 * @param {ArrayLike<number>} u
 * @param {ArrayLike<number>} v
 * @param {Float64Array} [out]
 * @returns {Float64Array}
 */
export function vec3Cross(u, v, out = new Float64Array(3)) {
  out[0] = u[1] * v[2] - u[2] * v[1];
  out[1] = u[2] * v[0] - u[0] * v[2];
  out[2] = u[0] * v[1] - u[1] * v[0];
  return out;
}

/**
 * Normalizes a 3D vector to unit length.
 * @param {ArrayLike<number>} v
 * @param {Float64Array} [out]
 * @returns {Float64Array} Unit vector
 */
export function vec3Normalize(v, out = new Float64Array(3)) {
  const norm = vec3Norm(v);
  if (norm === 0.0) {
    out[0] = 0.0; out[1] = 0.0; out[2] = 0.0;
    return out;
  }
  const invNorm = 1.0 / norm;
  out[0] = v[0] * invNorm;
  out[1] = v[1] * invNorm;
  out[2] = v[2] * invNorm;
  return out;
}

/**
 * Angle in radians between two 3D vectors.
 * @param {ArrayLike<number>} u
 * @param {ArrayLike<number>} v
 * @returns {number} Angle in [0, PI]
 */
export function vec3Angle(u, v) {
  const nu = vec3Norm(u);
  const nv = vec3Norm(v);
  if (nu === 0.0 || nv === 0.0) return 0.0;
  const cosTheta = Math.max(-1.0, Math.min(1.0, vec3Dot(u, v) / (nu * nv)));
  return Math.acos(cosTheta);
}

/**
 * Creates a 3D rotation matrix around an arbitrary unit axis by angle theta (Rodrigues formula).
 * @param {ArrayLike<number>} axis - 3-element axis vector (will be normalized)
 * @param {number} thetaRad - Rotation angle in radians
 * @param {Float64Array} [out]
 * @returns {Float64Array} 3x3 rotation matrix
 */
export function createRodriguesRotationMatrix(axis, thetaRad, out = new Float64Array(9)) {
  const u = vec3Normalize(axis);
  const c = Math.cos(thetaRad);
  const s = Math.sin(thetaRad);
  const C = 1.0 - c;
  const x = u[0], y = u[1], z = u[2];

  out[0] = x * x * C + c;     out[1] = x * y * C - z * s; out[2] = x * z * C + y * s;
  out[3] = y * x * C + z * s; out[4] = y * y * C + c;     out[5] = y * z * C - x * s;
  out[6] = z * x * C - y * s; out[7] = z * y * C + x * s; out[8] = z * z * C + c;

  return out;
}

/**
 * Computes 3x3 Jacobian Transformation Matrix for Spherical to Cartesian Coordinates.
 * J = [ dx/dr  dx/dphi  dx/dtheta ]
 *     [ dy/dr  dy/dphi  dy/dtheta ]
 *     [ dz/dr  dz/dphi  dz/dtheta ]
 * 
 * @param {number} lonDeg - Longitude in degrees
 * @param {number} latDeg - Latitude in degrees
 * @param {number} dist - Distance
 * @param {Float64Array} [out]
 * @returns {Float64Array} 9-element Jacobian matrix
 */
export function computeSphericalJacobian(lonDeg, latDeg, dist, out = new Float64Array(9)) {
  const phi = lonDeg * DEG2RAD;
  const theta = latDeg * DEG2RAD;
  const cosT = Math.cos(theta);
  const sinT = Math.sin(theta);
  const cosP = Math.cos(phi);
  const sinP = Math.sin(phi);

  out[0] = cosT * cosP;
  out[1] = -dist * cosT * sinP;
  out[2] = -dist * sinT * cosP;

  out[3] = cosT * sinP;
  out[4] = dist * cosT * cosP;
  out[5] = -dist * sinT * sinP;

  out[6] = sinT;
  out[7] = 0.0;
  out[8] = dist * cosT;

  return out;
}


// ============================================================================
// 3. ASTRONOMICAL ROTATION MATRICES (IAU J2000 & DE VAUCOULEURS 1976)
// ============================================================================

/**
 * Constructs Galactic to Equatorial (J2000) 3x3 Rotation Matrix R_Gal2Eq.
 * Defined such that X_eq = R_Gal2Eq * X_gal.
 * Columns of R_Gal2Eq are the Galactic unit vectors (x_gal, y_gal, z_gal) expressed in Equatorial coordinates.
 * 
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createGalacticToEquatorialMatrix() {
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

  return new Float64Array([
    nx_x, ny_x, z_x,
    nx_y, ny_y, z_y,
    nx_z, ny_z, z_z
  ]);
}

/**
 * Constructs Equatorial to Galactic (J2000) 3x3 Rotation Matrix R_Eq2Gal.
 * R_Eq2Gal = (R_Gal2Eq)^T.
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createEquatorialToGalacticMatrix() {
  const R_Gal2Eq = createGalacticToEquatorialMatrix();
  return matTranspose3x3(R_Gal2Eq);
}

/**
 * Constructs Supergalactic to Galactic 3x3 Rotation Matrix R_SG2Gal.
 * Defined by de Vaucouleurs et al. (1976 / 1991 RC3):
 * - SGP is at (l = 47.37°, b = 6.32°)
 * - Ascending node of SG plane on Galactic plane (SGL=0, SGB=0) is at (l = 137.37°, b = 0°)
 * 
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createSupergalacticToGalacticMatrix() {
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
  const y_norm = Math.hypot(y_x, y_y, y_z);

  return new Float64Array([
    x_x, y_x / y_norm, z_x,
    x_y, y_y / y_norm, z_y,
    x_z, y_z / y_norm, z_z
  ]);
}

/**
 * Constructs Galactic to Supergalactic 3x3 Rotation Matrix R_Gal2SG.
 * R_Gal2SG = (R_SG2Gal)^T.
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createGalacticToSupergalacticMatrix() {
  const R_SG2Gal = createSupergalacticToGalacticMatrix();
  return matTranspose3x3(R_SG2Gal);
}

/**
 * Constructs Supergalactic to Equatorial 3x3 Rotation Matrix R_SG2Eq.
 * R_SG2Eq = R_Gal2Eq * R_SG2Gal.
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createSupergalacticToEquatorialMatrix() {
  const R_Gal2Eq = createGalacticToEquatorialMatrix();
  const R_SG2Gal = createSupergalacticToGalacticMatrix();
  return matMul3x3(R_Gal2Eq, R_SG2Gal);
}

/**
 * Constructs Equatorial to Supergalactic 3x3 Rotation Matrix R_Eq2SG.
 * R_Eq2SG = (R_SG2Eq)^T = R_Gal2SG * R_Eq2Gal.
 * @returns {Float64Array} 9-element row-major matrix
 */
export function createEquatorialToSupergalacticMatrix() {
  const R_SG2Eq = createSupergalacticToEquatorialMatrix();
  return matTranspose3x3(R_SG2Eq);
}

// Pre-computed static rotation matrices
export const ROT_GAL_TO_EQ = createGalacticToEquatorialMatrix();
export const ROT_EQ_TO_GAL = createEquatorialToGalacticMatrix();
export const ROT_SG_TO_GAL = createSupergalacticToGalacticMatrix();
export const ROT_GAL_TO_SG = createGalacticToSupergalacticMatrix();
export const ROT_SG_TO_EQ = createSupergalacticToEquatorialMatrix();
export const ROT_EQ_TO_SG = createEquatorialToSupergalacticMatrix();


// ============================================================================
// 4. COORDINATE TRANSFORMATIONS (SPHERICAL & CARTESIAN)
// ============================================================================

/**
 * Wraps angle in degrees to [minDeg, maxDeg).
 * @param {number} deg
 * @param {number} [minDeg=0]
 * @param {number} [maxDeg=360]
 * @returns {number}
 */
export function wrapAngleDeg(deg, minDeg = 0.0, maxDeg = 360.0) {
  const range = maxDeg - minDeg;
  let wrapped = (deg - minDeg) % range;
  if (wrapped < 0.0) wrapped += range;
  return wrapped + minDeg;
}

/**
 * Clamps latitude / declination to [-90°, +90°].
 * @param {number} latDeg
 * @returns {number}
 */
export function clampLatitudeDeg(latDeg) {
  if (latDeg > 90.0) return 90.0;
  if (latDeg < -90.0) return -90.0;
  return latDeg;
}

/**
 * Converts Spherical Coordinates (Longitude, Latitude, Distance) to 3D Cartesian (X, Y, Z).
 * Convention: X = d * cos(lat) * cos(lon), Y = d * cos(lat) * sin(lon), Z = d * sin(lat).
 * 
 * @param {number} lonDeg - Longitude (RA, l, or SGL) in degrees
 * @param {number} latDeg - Latitude (Dec, b, or SGB) in degrees
 * @param {number} [dist=1.0] - Radial distance
 * @param {Float64Array} [out] - Optional destination array [x, y, z]
 * @returns {Float64Array} Cartesian [x, y, z]
 */
export function sphericalToCartesian(lonDeg, latDeg, dist = 1.0, out = new Float64Array(3)) {
  const phi = lonDeg * DEG2RAD;
  const theta = clampLatitudeDeg(latDeg) * DEG2RAD;
  const cosTheta = Math.cos(theta);

  out[0] = dist * cosTheta * Math.cos(phi);
  out[1] = dist * cosTheta * Math.sin(phi);
  out[2] = dist * Math.sin(theta);
  return out;
}

/**
 * Converts 3D Cartesian Coordinates (X, Y, Z) to Spherical (Longitude, Latitude, Distance).
 * Longitude in [0°, 360°), Latitude in [-90°, +90°].
 * 
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {Float64Array} [out] - Optional destination array [lonDeg, latDeg, dist]
 * @returns {Float64Array} Spherical [lonDeg, latDeg, dist]
 */
export function cartesianToSpherical(x, y, z, out = new Float64Array(3)) {
  const dist = Math.hypot(x, y, z);
  if (dist === 0.0) {
    out[0] = 0.0; out[1] = 0.0; out[2] = 0.0;
    return out;
  }
  const lonRad = Math.atan2(y, x);
  const latRad = Math.asin(Math.max(-1.0, Math.min(1.0, z / dist)));

  let lonDeg = lonRad * RAD2DEG;
  if (lonDeg < 0.0) lonDeg += 360.0;

  out[0] = lonDeg;
  out[1] = latRad * RAD2DEG;
  out[2] = dist;
  return out;
}

/**
 * Converts Equatorial Spherical (RA, Dec, d) to Galactic Spherical (l, b, d).
 * @param {number} raDeg - Right Ascension in degrees [0, 360)
 * @param {number} decDeg - Declination in degrees [-90, +90]
 * @param {number} [dist=1.0] - Radial distance
 * @returns {{l: number, b: number, dist: number}}
 */
export function equatorialToGalactic(raDeg, decDeg, dist = 1.0) {
  const vEq = sphericalToCartesian(raDeg, decDeg, dist);
  const vGal = matVecMul3x3(ROT_EQ_TO_GAL, vEq);
  const sph = cartesianToSpherical(vGal[0], vGal[1], vGal[2]);
  return { l: sph[0], b: sph[1], dist: sph[2] };
}

/**
 * Converts Galactic Spherical (l, b, d) to Equatorial Spherical (RA, Dec, d).
 * @param {number} lDeg - Galactic Longitude in degrees [0, 360)
 * @param {number} bDeg - Galactic Latitude in degrees [-90, +90]
 * @param {number} [dist=1.0] - Radial distance
 * @returns {{ra: number, dec: number, dist: number}}
 */
export function galacticToEquatorial(lDeg, bDeg, dist = 1.0) {
  const vGal = sphericalToCartesian(lDeg, bDeg, dist);
  const vEq = matVecMul3x3(ROT_GAL_TO_EQ, vGal);
  const sph = cartesianToSpherical(vEq[0], vEq[1], vEq[2]);
  return { ra: sph[0], dec: sph[1], dist: sph[2] };
}

/**
 * Converts Galactic Spherical (l, b, d) to Supergalactic Spherical (SGL, SGB, d).
 * @param {number} lDeg - Galactic Longitude in degrees
 * @param {number} bDeg - Galactic Latitude in degrees
 * @param {number} [dist=1.0] - Radial distance
 * @returns {{sgl: number, sgb: number, dist: number}}
 */
export function galacticToSupergalactic(lDeg, bDeg, dist = 1.0) {
  const vGal = sphericalToCartesian(lDeg, bDeg, dist);
  const vSG = matVecMul3x3(ROT_GAL_TO_SG, vGal);
  const sph = cartesianToSpherical(vSG[0], vSG[1], vSG[2]);
  return { sgl: sph[0], sgb: sph[1], dist: sph[2] };
}

/**
 * Converts Supergalactic Spherical (SGL, SGB, d) to Galactic Spherical (l, b, d).
 * @param {number} sglDeg - Supergalactic Longitude in degrees
 * @param {number} sgbDeg - Supergalactic Latitude in degrees
 * @param {number} [dist=1.0] - Radial distance
 * @returns {{l: number, b: number, dist: number}}
 */
export function supergalacticToGalactic(sglDeg, sgbDeg, dist = 1.0) {
  const vSG = sphericalToCartesian(sglDeg, sgbDeg, dist);
  const vGal = matVecMul3x3(ROT_SG_TO_GAL, vSG);
  const sph = cartesianToSpherical(vGal[0], vGal[1], vGal[2]);
  return { l: sph[0], b: sph[1], dist: sph[2] };
}

/**
 * Converts Equatorial Spherical (RA, Dec, d) to Supergalactic Spherical (SGL, SGB, d).
 * @param {number} raDeg
 * @param {number} decDeg
 * @param {number} [dist=1.0]
 * @returns {{sgl: number, sgb: number, dist: number}}
 */
export function equatorialToSupergalactic(raDeg, decDeg, dist = 1.0) {
  const vEq = sphericalToCartesian(raDeg, decDeg, dist);
  const vSG = matVecMul3x3(ROT_EQ_TO_SG, vEq);
  const sph = cartesianToSpherical(vSG[0], vSG[1], vSG[2]);
  return { sgl: sph[0], sgb: sph[1], dist: sph[2] };
}

/**
 * Converts Supergalactic Spherical (SGL, SGB, d) to Equatorial Spherical (RA, Dec, d).
 * @param {number} sglDeg
 * @param {number} sgbDeg
 * @param {number} [dist=1.0]
 * @returns {{ra: number, dec: number, dist: number}}
 */
export function supergalacticToEquatorial(sglDeg, sgbDeg, dist = 1.0) {
  const vSG = sphericalToCartesian(sglDeg, sgbDeg, dist);
  const vEq = matVecMul3x3(ROT_SG_TO_EQ, vSG);
  const sph = cartesianToSpherical(vEq[0], vEq[1], vEq[2]);
  return { ra: sph[0], dec: sph[1], dist: sph[2] };
}

/**
 * Direct Cartesian transformations between Galactic and Equatorial.
 */
export function cartesianGalToEq(xGal, yGal, zGal, out = new Float64Array(3)) {
  const v = [xGal, yGal, zGal];
  return matVecMul3x3(ROT_GAL_TO_EQ, v, out);
}

export function cartesianEqToGal(xEq, yEq, zEq, out = new Float64Array(3)) {
  const v = [xEq, yEq, zEq];
  return matVecMul3x3(ROT_EQ_TO_GAL, v, out);
}

/**
 * Direct Cartesian transformations between Galactic and Supergalactic.
 */
export function cartesianGalToSG(xGal, yGal, zGal, out = new Float64Array(3)) {
  const v = [xGal, yGal, zGal];
  return matVecMul3x3(ROT_GAL_TO_SG, v, out);
}

export function cartesianSGToGal(sgx, sgy, sgz, out = new Float64Array(3)) {
  const v = [sgx, sgy, sgz];
  return matVecMul3x3(ROT_SG_TO_GAL, v, out);
}

/**
 * Direct Cartesian transformations between Equatorial and Supergalactic.
 */
export function cartesianEqToSG(xEq, yEq, zEq, out = new Float64Array(3)) {
  const v = [xEq, yEq, zEq];
  return matVecMul3x3(ROT_EQ_TO_SG, v, out);
}

export function cartesianSGToEq(sgx, sgy, sgz, out = new Float64Array(3)) {
  const v = [sgx, sgy, sgz];
  return matVecMul3x3(ROT_SG_TO_EQ, v, out);
}


// ============================================================================
// 5. VELOCITY REFERENCE FRAME CONVERSIONS (HELIOCENTRIC, CMB, LOCAL GROUP, GSR)
// ============================================================================

/**
 * Computes the projection factor (cosine of angular separation) along line-of-sight to an apex.
 * @param {number} lDeg - Target Galactic longitude
 * @param {number} bDeg - Target Galactic latitude
 * @param {number} lApexDeg - Apex Galactic longitude
 * @param {number} bApexDeg - Apex Galactic latitude
 * @returns {number} cos(theta)
 */
export function computeApexProjection(lDeg, bDeg, lApexDeg, bApexDeg) {
  const l = lDeg * DEG2RAD;
  const b = bDeg * DEG2RAD;
  const lApex = lApexDeg * DEG2RAD;
  const bApex = bApexDeg * DEG2RAD;

  return Math.sin(b) * Math.sin(bApex) + Math.cos(b) * Math.cos(bApex) * Math.cos(l - lApex);
}

/**
 * Converts Heliocentric line-of-sight velocity to CMB Rest Frame velocity.
 * v_CMB = v_helio + v_apex * cos(theta_apex)
 * 
 * Default dipole (Planck 2018): v_apex = 369.8 km/s towards (l,b) = (264.021°, 48.253°)
 * 
 * @param {number} vHelioKms - Heliocentric velocity in km/s (cz)
 * @param {number} lDeg - Galactic longitude of object
 * @param {number} bDeg - Galactic latitude of object
 * @param {Object} [dipole=VELOCITY_REST_FRAMES.CMB_DIPOLE]
 * @returns {number} v_CMB in km/s
 */
export function heliocentricToCmb(vHelioKms, lDeg, bDeg, dipole = VELOCITY_REST_FRAMES.CMB_DIPOLE) {
  const proj = computeApexProjection(lDeg, bDeg, dipole.GALACTIC_L_DEG, dipole.GALACTIC_B_DEG);
  return vHelioKms + dipole.VELOCITY_KMS * proj;
}

/**
 * Converts CMB Rest Frame velocity to Heliocentric velocity.
 * v_helio = v_CMB - v_apex * cos(theta_apex)
 * 
 * @param {number} vCmbKms - CMB frame velocity in km/s
 * @param {number} lDeg
 * @param {number} bDeg
 * @param {Object} [dipole=VELOCITY_REST_FRAMES.CMB_DIPOLE]
 * @returns {number} v_helio in km/s
 */
export function cmbToHeliocentric(vCmbKms, lDeg, bDeg, dipole = VELOCITY_REST_FRAMES.CMB_DIPOLE) {
  const proj = computeApexProjection(lDeg, bDeg, dipole.GALACTIC_L_DEG, dipole.GALACTIC_B_DEG);
  return vCmbKms - dipole.VELOCITY_KMS * proj;
}

/**
 * Converts Heliocentric velocity to Local Group Barycenter frame velocity.
 * v_LG = v_helio + v_LG_apex * cos(theta_LG)
 * Default LG motion (Karachentsev & Makarov 1996): 316 km/s towards (93.0°, -4.0°)
 * 
 * @param {number} vHelioKms
 * @param {number} lDeg
 * @param {number} bDeg
 * @param {Object} [lgMotion=VELOCITY_REST_FRAMES.LOCAL_GROUP]
 * @returns {number} v_LG in km/s
 */
export function heliocentricToLocalGroup(vHelioKms, lDeg, bDeg, lgMotion = VELOCITY_REST_FRAMES.LOCAL_GROUP) {
  const proj = computeApexProjection(lDeg, bDeg, lgMotion.GALACTIC_L_DEG, lgMotion.GALACTIC_B_DEG);
  return vHelioKms + lgMotion.VELOCITY_KMS * proj;
}

/**
 * Converts Local Group Barycenter frame velocity to Heliocentric velocity.
 * @param {number} vLgKms
 * @param {number} lDeg
 * @param {number} bDeg
 * @param {Object} [lgMotion=VELOCITY_REST_FRAMES.LOCAL_GROUP]
 * @returns {number} v_helio in km/s
 */
export function localGroupToHeliocentric(vLgKms, lDeg, bDeg, lgMotion = VELOCITY_REST_FRAMES.LOCAL_GROUP) {
  const proj = computeApexProjection(lDeg, bDeg, lgMotion.GALACTIC_L_DEG, lgMotion.GALACTIC_B_DEG);
  return vLgKms - lgMotion.VELOCITY_KMS * proj;
}

/**
 * Converts CMB Rest Frame velocity directly to Local Group Barycenter velocity.
 * @param {number} vCmbKms
 * @param {number} lDeg
 * @param {number} bDeg
 * @param {Object} [dipole=VELOCITY_REST_FRAMES.CMB_DIPOLE]
 * @param {Object} [lgMotion=VELOCITY_REST_FRAMES.LOCAL_GROUP]
 * @returns {number} v_LG in km/s
 */
export function cmbToLocalGroup(vCmbKms, lDeg, bDeg, dipole = VELOCITY_REST_FRAMES.CMB_DIPOLE, lgMotion = VELOCITY_REST_FRAMES.LOCAL_GROUP) {
  const vHelio = cmbToHeliocentric(vCmbKms, lDeg, bDeg, dipole);
  return heliocentricToLocalGroup(vHelio, lDeg, bDeg, lgMotion);
}

/**
 * Converts Local Group Barycenter velocity to CMB Rest Frame velocity.
 * @param {number} vLgKms
 * @param {number} lDeg
 * @param {number} bDeg
 * @param {Object} [dipole=VELOCITY_REST_FRAMES.CMB_DIPOLE]
 * @param {Object} [lgMotion=VELOCITY_REST_FRAMES.LOCAL_GROUP]
 * @returns {number} v_CMB in km/s
 */
export function localGroupToCmb(vLgKms, lDeg, bDeg, dipole = VELOCITY_REST_FRAMES.CMB_DIPOLE, lgMotion = VELOCITY_REST_FRAMES.LOCAL_GROUP) {
  const vHelio = localGroupToHeliocentric(vLgKms, lDeg, bDeg, lgMotion);
  return heliocentricToCmb(vHelio, lDeg, bDeg, dipole);
}

/**
 * Converts Heliocentric line-of-sight velocity to Galactic Standard of Rest (GSR).
 * @param {number} vHelioKms
 * @param {number} lDeg
 * @param {number} bDeg
 * @param {Object} [gsrParams=VELOCITY_REST_FRAMES.GSR]
 * @returns {number} v_GSR in km/s
 */
export function heliocentricToGSR(vHelioKms, lDeg, bDeg, gsrParams = VELOCITY_REST_FRAMES.GSR) {
  const l = lDeg * DEG2RAD;
  const b = bDeg * DEG2RAD;
  const cosB = Math.cos(b);
  const sinB = Math.sin(b);
  const cosL = Math.cos(l);
  const sinL = Math.sin(l);

  const v_sun_los = (gsrParams.V_LSR_KMS + gsrParams.V_PEC_KMS) * sinL * cosB +
                    gsrParams.U_PEC_KMS * cosL * cosB +
                    gsrParams.W_PEC_KMS * sinB;

  return vHelioKms + v_sun_los;
}


// ============================================================================
// 6. VECTORIZED TYPEDARRAY BATCH TRANSFORMS (GPU/SHADER COMPATIBLE)
// ============================================================================

/**
 * Transforms an array/buffer of 3D Cartesian coordinates between frames.
 * @param {Float32Array|Float64Array} inBuffer
 * @param {Float32Array|Float64Array} outBuffer
 * @param {number} count
 * @param {'EQ'|'GAL'|'SG'} fromFrame
 * @param {'EQ'|'GAL'|'SG'} toFrame
 * @param {number} [stride=3]
 * @param {number} [offset=0]
 */
export function transformCartesianBatch(inBuffer, outBuffer, count, fromFrame, toFrame, stride = 3, offset = 0) {
  if (fromFrame === toFrame) {
    if (inBuffer !== outBuffer) {
      outBuffer.set(inBuffer);
    }
    return;
  }

  let R;
  if (fromFrame === 'GAL' && toFrame === 'EQ') R = ROT_GAL_TO_EQ;
  else if (fromFrame === 'EQ' && toFrame === 'GAL') R = ROT_EQ_TO_GAL;
  else if (fromFrame === 'GAL' && toFrame === 'SG') R = ROT_GAL_TO_SG;
  else if (fromFrame === 'SG' && toFrame === 'GAL') R = ROT_SG_TO_GAL;
  else if (fromFrame === 'EQ' && toFrame === 'SG') R = ROT_EQ_TO_SG;
  else if (fromFrame === 'SG' && toFrame === 'EQ') R = ROT_SG_TO_EQ;
  else throw new Error(`transformCartesianBatch: Unsupported frame pair ${fromFrame} -> ${toFrame}`);

  const r00 = R[0], r01 = R[1], r02 = R[2];
  const r10 = R[3], r11 = R[4], r12 = R[5];
  const r20 = R[6], r21 = R[7], r22 = R[8];

  for (let i = 0; i < count; i++) {
    const idx = i * stride + offset;
    const x = inBuffer[idx];
    const y = inBuffer[idx + 1];
    const z = inBuffer[idx + 2];

    outBuffer[idx] = r00 * x + r01 * y + r02 * z;
    outBuffer[idx + 1] = r10 * x + r11 * y + r12 * z;
    outBuffer[idx + 2] = r20 * x + r21 * y + r22 * z;
  }
}

/**
 * Transforms a full 6D phase space batch [X, Y, Z, Vx, Vy, Vz] between astronomical frames.
 * @param {Float32Array|Float64Array} inBuffer
 * @param {Float32Array|Float64Array} outBuffer
 * @param {number} count
 * @param {'EQ'|'GAL'|'SG'} fromFrame
 * @param {'EQ'|'GAL'|'SG'} toFrame
 * @param {number} [stride=6]
 */
export function transformPhaseSpace6DBatch(inBuffer, outBuffer, count, fromFrame, toFrame, stride = 6) {
  transformCartesianBatch(inBuffer, outBuffer, count, fromFrame, toFrame, stride, 0);
  transformCartesianBatch(inBuffer, outBuffer, count, fromFrame, toFrame, stride, 3);
}


// ============================================================================
// 7. REDSHIFT-SPACE DISTORTIONS (RSD) & KAISER / FOG ENGINE
// ============================================================================

/**
 * Cosmological Background & Expansion Rate Engine
 */
export class CosmologyEngine {
  /**
   * @param {Object} [params]
   * @param {number} [params.H0=74.6] - Hubble constant in km/s/Mpc
   * @param {number} [params.OmegaM=0.315] - Matter density parameter
   * @param {number} [params.OmegaL=0.685] - Dark energy density parameter
   * @param {number} [params.OmegaK=0.0] - Curvature density
   * @param {number} [params.OmegaR=0.0] - Radiation density
   * @param {number} [params.w0=-1.0] - Dark energy equation of state
   * @param {number} [params.wa=0.0] - Dark energy time evolution
   */
  constructor(params = {}) {
    this.H0 = params.H0 !== undefined ? params.H0 : COSMOLOGY_DEFAULTS.H0;
    this.h = this.H0 / 100.0;
    this.OmegaM = params.OmegaM !== undefined ? params.OmegaM : COSMOLOGY_DEFAULTS.OMEGA_M;
    this.OmegaL = params.OmegaL !== undefined ? params.OmegaL : COSMOLOGY_DEFAULTS.OMEGA_L;
    this.OmegaR = params.OmegaR !== undefined ? params.OmegaR : 0.0;
    this.OmegaK = params.OmegaK !== undefined ? params.OmegaK : (1.0 - this.OmegaM - this.OmegaL - this.OmegaR);
    this.w0 = params.w0 !== undefined ? params.w0 : COSMOLOGY_DEFAULTS.W0;
    this.wa = params.wa !== undefined ? params.wa : COSMOLOGY_DEFAULTS.WA;
  }

  /**
   * Dimensionless Hubble expansion parameter E(z) = H(z) / H0.
   * @param {number} z - Redshift
   * @returns {number} E(z)
   */
  E(z) {
    if (z < -0.999) throw new Error(`CosmologyEngine: Invalid negative redshift z = ${z}`);
    const zp1 = 1.0 + z;
    const zp1_2 = zp1 * zp1;
    const zp1_3 = zp1_2 * zp1;
    const zp1_4 = zp1_2 * zp1_2;

    let de_evol = 1.0;
    if (this.w0 !== -1.0 || this.wa !== 0.0) {
      de_evol = Math.pow(zp1, 3.0 * (1.0 + this.w0 + this.wa)) * Math.exp(-3.0 * this.wa * z / zp1);
    }

    const E2 = this.OmegaR * zp1_4 +
               this.OmegaM * zp1_3 +
               this.OmegaK * zp1_2 +
               this.OmegaL * de_evol;

    return Math.sqrt(Math.max(1e-12, E2));
  }

  /**
   * Hubble parameter H(z) in km/s/Mpc.
   * @param {number} z
   * @returns {number}
   */
  H(z) {
    return this.H0 * this.E(z);
  }

  /**
   * Scale factor a(z) = 1 / (1 + z).
   * @param {number} z
   * @returns {number}
   */
  a(z) {
    return 1.0 / (1.0 + z);
  }

  /**
   * Matter density parameter at redshift z.
   * @param {number} z
   * @returns {number}
   */
  OmegaM_z(z) {
    const zp1 = 1.0 + z;
    const Ez = this.E(z);
    return (this.OmegaM * zp1 * zp1 * zp1) / (Ez * Ez);
  }

  /**
   * Linear growth rate f(z) = d ln D / d ln a ≈ [Omega_m(z)]^gamma.
   * @param {number} z
   * @returns {number}
   */
  growthRateF(z) {
    const omz = this.OmegaM_z(z);
    const gamma = 0.5454 + 0.02 * (1.0 + this.w0);
    return Math.pow(Math.max(1e-6, omz), gamma);
  }

  /**
   * Linear growth factor D(z) normalized to D(z=0) = 1.
   * @param {number} z
   * @returns {number}
   */
  growthFactorD(z) {
    const zMax = 200.0;
    const nSteps = 128;
    const h = (zMax - z) / nSteps;

    let sum = 0.0;
    for (let i = 0; i < nSteps; i++) {
      const zA = z + i * h;
      const zB = z + (i + 1) * h;
      const zMid = 0.5 * (zA + zB);

      const fA = (1.0 + zA) / Math.pow(this.E(zA), 3);
      const fB = (1.0 + zB) / Math.pow(this.E(zB), 3);
      const fMid = (1.0 + zMid) / Math.pow(this.E(zMid), 3);

      sum += (h / 6.0) * (fA + 4.0 * fMid + fB);
    }

    const D_unnorm = 2.5 * this.OmegaM * this.E(z) * sum;

    if (this._D0 === undefined) {
      let sum0 = 0.0;
      const h0 = zMax / nSteps;
      for (let i = 0; i < nSteps; i++) {
        const zA = i * h0;
        const zB = (i + 1) * h0;
        const zMid = 0.5 * (zA + zB);
        const fA = (1.0 + zA) / Math.pow(this.E(zA), 3);
        const fB = (1.0 + zB) / Math.pow(this.E(zB), 3);
        const fMid = (1.0 + zMid) / Math.pow(this.E(zMid), 3);
        sum0 += (h0 / 6.0) * (fA + 4.0 * fMid + fB);
      }
      this._D0 = 2.5 * this.OmegaM * this.E(0.0) * sum0;
    }

    return D_unnorm / this._D0;
  }

  /**
   * Co-moving line-of-sight distance D_C(z) in Mpc.
   * @param {number} z
   * @returns {number}
   */
  comovingDistanceMpc(z) {
    if (z <= 0.0) return 0.0;
    const nSteps = 64;
    const h = z / nSteps;
    let sum = 0.0;
    for (let i = 0; i < nSteps; i++) {
      const zA = i * h;
      const zB = (i + 1) * h;
      const zMid = 0.5 * (zA + zB);
      const fA = 1.0 / this.E(zA);
      const fB = 1.0 / this.E(zB);
      const fMid = 1.0 / this.E(zMid);
      sum += (h / 6.0) * (fA + 4.0 * fMid + fB);
    }
    return (SPEED_OF_LIGHT_KMS / this.H0) * sum;
  }

  /**
   * Angular diameter distance D_A(z) in Mpc.
   * @param {number} z
   * @returns {number}
   */
  angularDiameterDistanceMpc(z) {
    if (z <= 0.0) return 0.0;
    const dC = this.comovingDistanceMpc(z);
    return dC / (1.0 + z);
  }

  /**
   * Luminosity distance D_L(z) in Mpc.
   * @param {number} z
   * @returns {number}
   */
  luminosityDistanceMpc(z) {
    if (z <= 0.0) return 0.0;
    const dC = this.comovingDistanceMpc(z);
    return (1.0 + z) * dC;
  }

  /**
   * Computes Alcock-Paczynski distortion parameters relative to a fiducial cosmology.
   * @param {number} z
   * @param {CosmologyEngine} fiducialCosmology
   * @returns {{qPerp: number, qPara: number, F_AP: number}}
   */
  alcockPaczynskiParameters(z, fiducialCosmology) {
    const da = this.angularDiameterDistanceMpc(z);
    const daFid = fiducialCosmology.angularDiameterDistanceMpc(z);
    const hz = this.H(z);
    const hzFid = fiducialCosmology.H(z);

    const qPerp = daFid > 0 ? da / daFid : 1.0;
    const qPara = hz > 0 ? hzFid / hz : 1.0;
    const F_AP = qPerp > 0 ? qPara / qPerp : 1.0;

    return { qPerp, qPara, F_AP };
  }

  /**
   * Inverts comoving distance to redshift z via Newton-Raphson.
   * @param {number} distMpc
   * @param {number} [tol=1e-8]
   * @returns {number}
   */
  redshiftFromComovingDistance(distMpc, tol = 1e-8) {
    if (distMpc <= 0.0) return 0.0;
    let z = (this.H0 * distMpc) / SPEED_OF_LIGHT_KMS;

    for (let iter = 0; iter < 50; iter++) {
      const d = this.comovingDistanceMpc(z);
      const diff = d - distMpc;
      if (Math.abs(diff) < tol * distMpc) return z;

      const dDdz = SPEED_OF_LIGHT_KMS / this.H(z);
      z -= diff / dDdz;
      if (z < 0.0) z = 1e-6;
    }
    return z;
  }
}

/**
 * Forward Real-Space to Redshift-Space mapping:
 * s = x + (v_los / (a * H(z))) * \hat{r}
 * 
 * @param {ArrayLike<number>} xReal - Real-space 3D position [x, y, z] (in Mpc or Mpc/h)
 * @param {ArrayLike<number>} vPec - Peculiar velocity 3D vector [vx, vy, vz] (in km/s)
 * @param {Object} [options]
 * @param {CosmologyEngine} [options.cosmology]
 * @param {boolean} [options.isMpcOverH=true]
 * @param {Float64Array} [out]
 * @returns {Float64Array} Redshift-space position s [sx, sy, sz]
 */
export function realToRedshiftSpace(xReal, vPec, options = {}, out = new Float64Array(3)) {
  const cosmo = options.cosmology || new CosmologyEngine();
  const isMpcOverH = options.isMpcOverH !== undefined ? options.isMpcOverH : true;

  const r = vec3Norm(xReal);
  if (r < 1e-10) {
    out[0] = xReal[0]; out[1] = xReal[1]; out[2] = xReal[2];
    return out;
  }

  const rHat_x = xReal[0] / r;
  const rHat_y = xReal[1] / r;
  const rHat_z = xReal[2] / r;

  const v_los = vPec[0] * rHat_x + vPec[1] * rHat_y + vPec[2] * rHat_z;

  const distMpc = isMpcOverH ? r / cosmo.h : r;
  const z = cosmo.redshiftFromComovingDistance(distMpc);
  const a = cosmo.a(z);
  const H_z = cosmo.H(z);

  let delta_r = v_los / (a * H_z);
  if (isMpcOverH) {
    delta_r *= cosmo.h;
  }

  out[0] = xReal[0] + delta_r * rHat_x;
  out[1] = xReal[1] + delta_r * rHat_y;
  out[2] = xReal[2] + delta_r * rHat_z;

  return out;
}

/**
 * Inverse Redshift-Space to Real-Space mapping (Iterative Kaiser Inversion):
 * Solves s = x + (v_los / (a(z_x) H(z_x))) \hat{x}.
 * 
 * @param {ArrayLike<number>} sRedshift - Redshift-space 3D position [sx, sy, sz]
 * @param {ArrayLike<number>} vPec - Peculiar velocity 3D vector [vx, vy, vz] (in km/s)
 * @param {Object} [options]
 * @param {Float64Array} [out]
 * @returns {Float64Array} Real-space position x [x, y, z]
 */
export function redshiftToRealSpace(sRedshift, vPec, options = {}, out = new Float64Array(3)) {
  const cosmo = options.cosmology || new CosmologyEngine();
  const isMpcOverH = options.isMpcOverH !== undefined ? options.isMpcOverH : true;

  const s = vec3Norm(sRedshift);
  if (s < 1e-10) {
    out[0] = sRedshift[0]; out[1] = sRedshift[1]; out[2] = sRedshift[2];
    return out;
  }

  const sHat_x = sRedshift[0] / s;
  const sHat_y = sRedshift[1] / s;
  const sHat_z = sRedshift[2] / s;

  const v_los = vPec[0] * sHat_x + vPec[1] * sHat_y + vPec[2] * sHat_z;

  // Fixed-point iteration to evaluate a(z) and H(z) at real-space distance r
  let rEst = s;
  for (let iter = 0; iter < 4; iter++) {
    const distMpc = isMpcOverH ? rEst / cosmo.h : rEst;
    const z = cosmo.redshiftFromComovingDistance(distMpc);
    const a = cosmo.a(z);
    const H_z = cosmo.H(z);
    let delta_s = v_los / (a * H_z);
    if (isMpcOverH) delta_s *= cosmo.h;
    rEst = s - delta_s;
  }

  const distMpc = isMpcOverH ? rEst / cosmo.h : rEst;
  const z = cosmo.redshiftFromComovingDistance(distMpc);
  const a = cosmo.a(z);
  const H_z = cosmo.H(z);

  let delta_s = v_los / (a * H_z);
  if (isMpcOverH) {
    delta_s *= cosmo.h;
  }

  out[0] = sRedshift[0] - delta_s * sHat_x;
  out[1] = sRedshift[1] - delta_s * sHat_y;
  out[2] = sRedshift[2] - delta_s * sHat_z;

  return out;
}

/**
 * Anisotropic Redshift-Space Power Spectrum Multipoles (Kaiser 1987 + FoG damping).
 */
export class KaiserRSDModel {
  /**
   * @param {Object} [params]
   * @param {number} [params.growthRateF=0.55]
   * @param {number} [params.biasB=1.0]
   * @param {number} [params.sigmaV=300.0]
   * @param {'lorentzian'|'gaussian'|'squared_lorentzian'|'exponential'|'none'} [params.fogModel='lorentzian']
   * @param {number} [params.H0=100.0]
   */
  constructor(params = {}) {
    this.growthRateF = params.growthRateF !== undefined ? params.growthRateF : 0.55;
    this.biasB = params.biasB !== undefined ? params.biasB : 1.0;
    this.beta = this.growthRateF / Math.max(1e-4, this.biasB);
    this.sigmaV = params.sigmaV !== undefined ? params.sigmaV : 300.0;
    this.fogModel = params.fogModel || 'lorentzian';
    this.H0 = params.H0 || 100.0;
  }

  /**
   * Evaluates the non-linear Fingers-of-God damping factor D_FoG(k_los).
   * @param {number} k
   * @param {number} mu
   * @returns {number}
   */
  fogDamping(k, mu) {
    if (this.fogModel === 'none' || this.sigmaV === 0.0) return 1.0;
    const k_los = k * mu;
    const sigmaDist = this.sigmaV / this.H0;
    const kSigma = Math.abs(k_los * sigmaDist);
    const kSigma2 = kSigma * kSigma;

    switch (this.fogModel) {
      case 'lorentzian':
        return 1.0 / (1.0 + 0.5 * kSigma2);
      case 'squared_lorentzian': {
        const l = 1.0 / (1.0 + 0.5 * kSigma2);
        return l * l;
      }
      case 'gaussian':
        return Math.exp(-kSigma2);
      case 'exponential':
        return Math.exp(-Math.SQRT2 * kSigma);
      default:
        return 1.0;
    }
  }

  /**
   * Evaluates anisotropic 2D power spectrum P(k, mu).
   * @param {number} k
   * @param {number} mu
   * @param {number} PkReal
   * @returns {number}
   */
  anisotropicPower(k, mu, PkReal) {
    const mu2 = mu * mu;
    const kaiserFactor = 1.0 + this.beta * mu2;
    const kaiserFactor2 = kaiserFactor * kaiserFactor;
    const damping = this.fogDamping(k, mu);
    return kaiserFactor2 * PkReal * damping;
  }

  /**
   * Exact analytical linear Kaiser multipoles.
   * @param {number} PkReal
   * @returns {{P0: number, P2: number, P4: number}}
   */
  analyticalLinearMultipoles(PkReal) {
    const b = this.beta;
    const b2 = b * b;

    const P0 = (1.0 + (2.0 / 3.0) * b + (1.0 / 5.0) * b2) * PkReal;
    const P2 = ((4.0 / 3.0) * b + (4.0 / 7.0) * b2) * PkReal;
    const P4 = (8.0 / 35.0) * b2 * PkReal;

    return { P0, P2, P4 };
  }

  /**
   * Numerically computes Legendre multipoles P_l(k) via Simpson's quadrature.
   * @param {number} k
   * @param {number} PkReal
   * @returns {{P0: number, P2: number, P4: number}}
   */
  computeMultipolesNumerical(k, PkReal) {
    const n = 120;
    const dMu = 1.0 / n;

    let sum0 = 0.0, sum2 = 0.0, sum4 = 0.0;

    for (let i = 0; i <= n; i++) {
      const mu = i * dMu;
      const weight = (i === 0 || i === n) ? 1.0 : (i % 2 === 1 ? 4.0 : 2.0);

      const mu2 = mu * mu;
      const mu4 = mu2 * mu2;

      const L0 = 1.0;
      const L2 = 0.5 * (3.0 * mu2 - 1.0);
      const L4 = 0.125 * (35.0 * mu4 - 30.0 * mu2 + 3.0);

      const P_s = this.anisotropicPower(k, mu, PkReal);

      sum0 += weight * P_s * L0;
      sum2 += weight * P_s * L2;
      sum4 += weight * P_s * L4;
    }

    const norm = dMu / 3.0;
    const P0 = (2 * 0 + 1) * sum0 * norm;
    const P2 = (2 * 2 + 1) * sum2 * norm;
    const P4 = (2 * 4 + 1) * sum4 * norm;

    return { P0, P2, P4 };
  }
}

/**
 * Fingers-of-God (FoG) Virialized Velocity Dispersion & Cluster Radial Compressor.
 */
export class ClusterFoGCompressor {
  /**
   * Estimates robust 1D velocity dispersion using Beers, Flynn & Gebhardt (1990) Biweight Scale.
   * @param {ArrayLike<number>} velocities
   * @param {number} [tuningC=6.0]
   * @returns {{median: number, sigmaBiweight: number, count: number}}
   */
  static estimateBiweightDispersion(velocities, tuningC = 6.0) {
    const n = velocities.length;
    if (n < 2) {
      return { median: velocities[0] || 0.0, sigmaBiweight: 0.0, count: n };
    }

    const sorted = Array.from(velocities).sort((a, b) => a - b);
    const mid = Math.floor(n / 2);
    const median = n % 2 === 1 ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);

    const absDevs = sorted.map(v => Math.abs(v - median)).sort((a, b) => a - b);
    const mad = n % 2 === 1 ? absDevs[mid] : 0.5 * (absDevs[mid - 1] + absDevs[mid]);

    if (mad < 1e-12) {
      return { median, sigmaBiweight: 0.0, count: n };
    }

    const cMAD = tuningC * mad;
    let numSum = 0.0;
    let denSum = 0.0;

    for (let i = 0; i < n; i++) {
      const u = (velocities[i] - median) / cMAD;
      if (Math.abs(u) < 1.0) {
        const u2 = u * u;
        const oneMinusU2 = 1.0 - u2;
        const oneMinus5U2 = 1.0 - 5.0 * u2;
        numSum += Math.pow(velocities[i] - median, 2) * Math.pow(oneMinusU2, 4);
        denSum += oneMinusU2 * oneMinus5U2;
      }
    }

    const sigmaBiweight = Math.sqrt(n) * Math.sqrt(numSum) / Math.abs(denSum);
    return { median, sigmaBiweight, count: n };
  }

  /**
   * Gapper Estimator (Wainer & Thissen 1976 / Beers et al. 1990) for small cluster sample sizes (n < 15).
   * @param {ArrayLike<number>} velocities
   * @returns {{sigmaGapper: number, count: number}}
   */
  static estimateGapperDispersion(velocities) {
    const n = velocities.length;
    if (n < 2) return { sigmaGapper: 0.0, count: n };

    const sorted = Array.from(velocities).sort((a, b) => a - b);
    let sum = 0.0;

    for (let i = 1; i < n; i++) {
      const weight = i * (n - i);
      const gap = sorted[i] - sorted[i - 1];
      sum += weight * gap;
    }

    const factor = Math.sqrt(Math.PI) / (n * (n - 1));
    return { sigmaGapper: factor * sum, count: n };
  }

  /**
   * Standard iterative sigma-clipped velocity dispersion.
   * @param {ArrayLike<number>} velocities
   * @param {number} [nSigma=3.0]
   * @param {number} [maxIter=5]
   * @returns {{mean: number, std: number, count: number}}
   */
  static estimateSigmaClippedDispersion(velocities, nSigma = 3.0, maxIter = 5) {
    let current = Array.from(velocities);
    for (let iter = 0; iter < maxIter; iter++) {
      const n = current.length;
      if (n < 2) break;
      const mean = current.reduce((a, b) => a + b, 0.0) / n;
      const variance = current.reduce((a, b) => a + Math.pow(b - mean, 2), 0.0) / (n - 1);
      const std = Math.sqrt(variance);

      const filtered = current.filter(v => Math.abs(v - mean) <= nSigma * std);
      if (filtered.length === current.length) break;
      current = filtered;
    }

    const n = current.length;
    const mean = current.reduce((a, b) => a + b, 0.0) / Math.max(1, n);
    const variance = n > 1 ? current.reduce((a, b) => a + Math.pow(b - mean, 2), 0.0) / (n - 1) : 0.0;
    return { mean, std: Math.sqrt(variance), count: n };
  }

  /**
   * Compresses radial elongation of galaxies in a cluster FoG towards the cluster centroid.
   * @param {Array<{x: number, y: number, z: number, cz?: number}>} memberPoints
   * @param {{x: number, y: number, z: number}} clusterCenter
   * @param {number} sigmaV
   * @param {number} rVirMpc
   * @param {Object} [options]
   * @returns {Array<{x: number, y: number, z: number}>}
   */
  static compressClusterFoG(memberPoints, clusterCenter, sigmaV, rVirMpc, options = {}) {
    const H0 = options.H0 || 100.0;
    const rC = Math.hypot(clusterCenter.x, clusterCenter.y, clusterCenter.z);
    if (rC < 1e-10) return memberPoints.map(p => ({ ...p }));

    const rHat_x = clusterCenter.x / rC;
    const rHat_y = clusterCenter.y / rC;
    const rHat_z = clusterCenter.z / rC;

    const elongationFactor = Math.max(1.0, 1.0 + (sigmaV / Math.max(1e-4, H0 * rVirMpc)));
    const compressionRatio = 1.0 / elongationFactor;

    return memberPoints.map(p => {
      const dx = p.x - clusterCenter.x;
      const dy = p.y - clusterCenter.y;
      const dz = p.z - clusterCenter.z;

      const dr_los = dx * rHat_x + dy * rHat_y + dz * rHat_z;

      const dx_perp = dx - dr_los * rHat_x;
      const dy_perp = dy - dr_los * rHat_y;
      const dz_perp = dz - dr_los * rHat_z;

      const dr_los_compressed = dr_los * compressionRatio;

      return {
        x: clusterCenter.x + dx_perp + dr_los_compressed * rHat_x,
        y: clusterCenter.y + dy_perp + dr_los_compressed * rHat_y,
        z: clusterCenter.z + dz_perp + dr_los_compressed * rHat_z
      };
    });
  }
}


// ============================================================================
// 8. ASTROMETRIC FITS WCS CD/PC MATRIX PROJECTION PARSER
// ============================================================================

/**
 * Astrometric FITS WCS Header Data Representation.
 */
export class WCSHeader {
  constructor(params = {}) {
    this.crpix1 = params.CRPIX1 !== undefined ? Number(params.CRPIX1) : 1.0;
    this.crpix2 = params.CRPIX2 !== undefined ? Number(params.CRPIX2) : 1.0;

    this.crval1 = params.CRVAL1 !== undefined ? Number(params.CRVAL1) : 0.0;
    this.crval2 = params.CRVAL2 !== undefined ? Number(params.CRVAL2) : 0.0;

    this.ctype1 = (params.CTYPE1 || 'RA---TAN').trim().toUpperCase();
    this.ctype2 = (params.CTYPE2 || 'DEC--TAN').trim().toUpperCase();

    this.cdelt1 = params.CDELT1 !== undefined ? Number(params.CDELT1) : 1.0;
    this.cdelt2 = params.CDELT2 !== undefined ? Number(params.CDELT2) : 1.0;

    this.lonpole = params.LONPOLE !== undefined ? Number(params.LONPOLE) : 180.0;
    this.latpole = params.LATPOLE !== undefined ? Number(params.LATPOLE) : 90.0;

    this.naxis1 = params.NAXIS1 !== undefined ? Number(params.NAXIS1) : 1024;
    this.naxis2 = params.NAXIS2 !== undefined ? Number(params.NAXIS2) : 1024;

    this.cd = new Float64Array(4);
    this._initializeCDMatrix(params);

    this.sip = this._initializeSIP(params);
    this.projectionCode = this._extractProjectionCode();
  }

  _initializeCDMatrix(params) {
    if (params.CD1_1 !== undefined || params.CD1_2 !== undefined || params.CD2_1 !== undefined || params.CD2_2 !== undefined) {
      this.cd[0] = params.CD1_1 !== undefined ? Number(params.CD1_1) : this.cdelt1;
      this.cd[1] = params.CD1_2 !== undefined ? Number(params.CD1_2) : 0.0;
      this.cd[2] = params.CD2_1 !== undefined ? Number(params.CD2_1) : 0.0;
      this.cd[3] = params.CD2_2 !== undefined ? Number(params.CD2_2) : this.cdelt2;
    } else if (params.PC1_1 !== undefined || params.PC1_2 !== undefined || params.PC2_1 !== undefined || params.PC2_2 !== undefined) {
      const pc11 = params.PC1_1 !== undefined ? Number(params.PC1_1) : 1.0;
      const pc12 = params.PC1_2 !== undefined ? Number(params.PC1_2) : 0.0;
      const pc21 = params.PC2_1 !== undefined ? Number(params.PC2_1) : 0.0;
      const pc22 = params.PC2_2 !== undefined ? Number(params.PC2_2) : 1.0;

      this.cd[0] = this.cdelt1 * pc11;
      this.cd[1] = this.cdelt1 * pc12;
      this.cd[2] = this.cdelt2 * pc21;
      this.cd[3] = this.cdelt2 * pc22;
    } else if (params.CROTA2 !== undefined || params.CROTA1 !== undefined) {
      const rotDeg = params.CROTA2 !== undefined ? Number(params.CROTA2) : Number(params.CROTA1);
      const rotRad = rotDeg * DEG2RAD;
      const cosR = Math.cos(rotRad);
      const sinR = Math.sin(rotRad);

      this.cd[0] = this.cdelt1 * cosR;
      this.cd[1] = -this.cdelt2 * sinR;
      this.cd[2] = this.cdelt1 * sinR;
      this.cd[3] = this.cdelt2 * cosR;
    } else {
      this.cd[0] = this.cdelt1;
      this.cd[1] = 0.0;
      this.cd[2] = 0.0;
      this.cd[3] = this.cdelt2;
    }
  }

  _initializeSIP(params) {
    const hasSip = params.A_ORDER !== undefined || params.B_ORDER !== undefined ||
                   this.ctype1.includes('-SIP') || this.ctype2.includes('-SIP');
    if (!hasSip) return null;

    const aOrder = params.A_ORDER !== undefined ? Number(params.A_ORDER) : 0;
    const bOrder = params.B_ORDER !== undefined ? Number(params.B_ORDER) : 0;
    const apOrder = params.AP_ORDER !== undefined ? Number(params.AP_ORDER) : 0;
    const bpOrder = params.BP_ORDER !== undefined ? Number(params.BP_ORDER) : 0;

    const aCoeffs = {};
    const bCoeffs = {};
    const apCoeffs = {};
    const bpCoeffs = {};

    for (const [key, val] of Object.entries(params)) {
      const uKey = key.toUpperCase();
      let match;
      if ((match = uKey.match(/^A_(\d+)_(\d+)$/))) {
        aCoeffs[`${match[1]}_${match[2]}`] = Number(val);
      } else if ((match = uKey.match(/^B_(\d+)_(\d+)$/))) {
        bCoeffs[`${match[1]}_${match[2]}`] = Number(val);
      } else if ((match = uKey.match(/^AP_(\d+)_(\d+)$/))) {
        apCoeffs[`${match[1]}_${match[2]}`] = Number(val);
      } else if ((match = uKey.match(/^BP_(\d+)_(\d+)$/))) {
        bpCoeffs[`${match[1]}_${match[2]}`] = Number(val);
      }
    }

    return { aOrder, bOrder, apOrder, bpOrder, aCoeffs, bCoeffs, apCoeffs, bpCoeffs };
  }

  _extractProjectionCode() {
    if (this.ctype1.length >= 8 && this.ctype1.charAt(4) === '-') {
      return this.ctype1.substring(5, 8).trim();
    }
    if (this.ctype2.length >= 8 && this.ctype2.charAt(4) === '-') {
      return this.ctype2.substring(5, 8).trim();
    }
    return 'TAN';
  }

  get cdDeterminant() {
    return this.cd[0] * this.cd[3] - this.cd[1] * this.cd[2];
  }

  get cdInverse() {
    const det = this.cdDeterminant;
    if (Math.abs(det) < 1e-20) {
      throw new Error(`WCSHeader: Singular CD matrix with det = ${det}`);
    }
    const invDet = 1.0 / det;
    return new Float64Array([
      this.cd[3] * invDet,
      -this.cd[1] * invDet,
      -this.cd[2] * invDet,
      this.cd[0] * invDet
    ]);
  }

  get pixelScaleArcsec() {
    const det = Math.abs(this.cdDeterminant);
    return Math.sqrt(det) * 3600.0;
  }
}

/**
 * Astrometric FITS WCS Projections and Coordinate Transformations Engine.
 */
export class FITS_WCS_Parser {
  constructor(headerInput) {
    if (headerInput instanceof WCSHeader) {
      this.wcs = headerInput;
    } else if (typeof headerInput === 'string') {
      this.wcs = FITS_WCS_Parser.parseFITSHeaderCards(headerInput);
    } else if (typeof headerInput === 'object') {
      this.wcs = new WCSHeader(headerInput);
    } else {
      this.wcs = new WCSHeader({});
    }
  }

  static parseFITSHeaderCards(rawHeaderString) {
    const params = {};
    const lines = rawHeaderString.includes('\n')
      ? rawHeaderString.split('\n')
      : rawHeaderString.match(/.{1,80}/g) || [];

    for (const rawLine of lines) {
      const line = rawLine.trimEnd();
      if (!line || line.startsWith('COMMENT') || line.startsWith('HISTORY')) continue;
      if (line.startsWith('END')) break;

      const eqIdx = line.indexOf('=');
      if (eqIdx > 0) {
        const key = line.substring(0, eqIdx).trim().toUpperCase();
        let valStr = line.substring(eqIdx + 1);

        const slashIdx = valStr.indexOf('/');
        if (slashIdx >= 0) {
          valStr = valStr.substring(0, slashIdx);
        }
        valStr = valStr.trim();

        if (valStr.startsWith("'")) {
          const endQuote = valStr.indexOf("'", 1);
          if (endQuote > 0) {
            params[key] = valStr.substring(1, endQuote).trim();
          } else {
            params[key] = valStr.replace(/'/g, '').trim();
          }
        } else if (valStr === 'T' || valStr === 'F') {
          params[key] = valStr === 'T';
        } else {
          const num = Number(valStr);
          params[key] = isNaN(num) ? valStr : num;
        }
      }
    }

    return new WCSHeader(params);
  }

  pixelToIntermediate(xPix, yPix, origin = 1) {
    const u = xPix - this.wcs.crpix1 + (origin === 0 ? 1 : 0);
    const v = yPix - this.wcs.crpix2 + (origin === 0 ? 1 : 0);

    let uDistorted = u;
    let vDistorted = v;

    if (this.wcs.sip) {
      const sip = this.wcs.sip;
      let du = 0.0;
      let dv = 0.0;

      for (const [key, coeff] of Object.entries(sip.aCoeffs)) {
        const [p, q] = key.split('_').map(Number);
        du += coeff * Math.pow(u, p) * Math.pow(v, q);
      }
      for (const [key, coeff] of Object.entries(sip.bCoeffs)) {
        const [p, q] = key.split('_').map(Number);
        dv += coeff * Math.pow(u, p) * Math.pow(v, q);
      }

      uDistorted += du;
      vDistorted += dv;
    }

    const cd = this.wcs.cd;
    const x_int = cd[0] * uDistorted + cd[1] * vDistorted;
    const y_int = cd[2] * uDistorted + cd[3] * vDistorted;

    return { x_int, y_int };
  }

  intermediateToPixel(xInt, yInt, origin = 1) {
    const invCD = this.wcs.cdInverse;
    const uPrime = invCD[0] * xInt + invCD[1] * yInt;
    const vPrime = invCD[2] * xInt + invCD[3] * yInt;

    let u = uPrime;
    let v = vPrime;

    if (this.wcs.sip && (this.wcs.sip.apOrder > 0 || this.wcs.sip.bpOrder > 0)) {
      const sip = this.wcs.sip;
      let du = 0.0;
      let dv = 0.0;

      for (const [key, coeff] of Object.entries(sip.apCoeffs)) {
        const [p, q] = key.split('_').map(Number);
        du += coeff * Math.pow(uPrime, p) * Math.pow(vPrime, q);
      }
      for (const [key, coeff] of Object.entries(sip.bpCoeffs)) {
        const [p, q] = key.split('_').map(Number);
        dv += coeff * Math.pow(uPrime, p) * Math.pow(vPrime, q);
      }

      u += du;
      v += dv;
    }

    const xPix = u + this.wcs.crpix1 - (origin === 0 ? 1 : 0);
    const yPix = v + this.wcs.crpix2 - (origin === 0 ? 1 : 0);

    return { xPix, yPix };
  }

  intermediateToNativeSpherical(xDeg, yDeg, projCode) {
    const x = xDeg * DEG2RAD;
    const y = yDeg * DEG2RAD;
    const r = Math.hypot(x, y);

    let phiRad = 0.0;
    let thetaRad = 0.0;

    switch (projCode) {
      case 'TAN': {
        phiRad = Math.atan2(x, -y);
        thetaRad = r === 0.0 ? HALF_PI : Math.atan2(1.0, r);
        break;
      }
      case 'SIN': {
        if (r > 1.0) throw new Error(`WCS Projection SIN: Point out of sky projection radius (r = ${r})`);
        phiRad = Math.atan2(x, -y);
        thetaRad = Math.acos(r);
        break;
      }
      case 'ARC': {
        if (r > Math.PI) throw new Error(`WCS Projection ARC: Point beyond antipodal limit (r = ${r})`);
        phiRad = Math.atan2(x, -y);
        thetaRad = HALF_PI - r;
        break;
      }
      case 'STG': {
        phiRad = Math.atan2(x, -y);
        thetaRad = HALF_PI - 2.0 * Math.atan(0.5 * r);
        break;
      }
      case 'CAR': {
        phiRad = x;
        thetaRad = y;
        break;
      }
      case 'CEA': {
        if (Math.abs(y) > 1.0) throw new Error(`WCS Projection CEA: Invalid y = ${y}`);
        phiRad = x;
        thetaRad = Math.asin(Math.max(-1.0, Math.min(1.0, y)));
        break;
      }
      case 'MER': {
        phiRad = x;
        thetaRad = 2.0 * Math.atan(Math.exp(y)) - HALF_PI;
        break;
      }
      case 'AIT': {
        const z2 = 1.0 - Math.pow(x / 4.0, 2) - Math.pow(y / 2.0, 2);
        const z = Math.sqrt(Math.max(0, z2));
        phiRad = 2.0 * Math.atan2(z * x / 2.0, 2.0 * z2 - 1.0);
        thetaRad = Math.asin(Math.max(-1.0, Math.min(1.0, z * y)));
        break;
      }
      case 'MOL': {
        const sinGamma = Math.max(-1.0, Math.min(1.0, y / Math.SQRT2));
        const gamma = Math.asin(sinGamma);
        const sinTheta = (2.0 * gamma + Math.sin(2.0 * gamma)) / Math.PI;
        thetaRad = Math.asin(Math.max(-1.0, Math.min(1.0, sinTheta)));
        const cosGamma = Math.cos(gamma);
        phiRad = Math.abs(cosGamma) > 1e-10 ? (Math.PI * x) / (2.0 * Math.SQRT2 * cosGamma) : 0.0;
        break;
      }
      case 'PAR': {
        const sinThird = Math.max(-1.0, Math.min(1.0, y / 3.0));
        thetaRad = 3.0 * Math.asin(sinThird);
        const denom = 2.0 * Math.cos(2.0 * thetaRad / 3.0) - 1.0;
        phiRad = Math.abs(denom) > 1e-10 ? (Math.PI * x / 2.0) / denom : 0.0;
        break;
      }
      default:
        phiRad = Math.atan2(x, -y);
        thetaRad = r === 0.0 ? HALF_PI : Math.atan2(1.0, r);
        break;
    }

    return { phiDeg: phiRad * RAD2DEG, thetaDeg: thetaRad * RAD2DEG };
  }

  nativeSphericalToIntermediate(phiDeg, thetaDeg, projCode) {
    const phi = phiDeg * DEG2RAD;
    const theta = clampLatitudeDeg(thetaDeg) * DEG2RAD;

    let xRad = 0.0;
    let yRad = 0.0;

    switch (projCode) {
      case 'TAN': {
        if (theta <= 0.0) throw new Error(`WCS Projection TAN: Coordinate theta <= 0° cannot be projected.`);
        const r = 1.0 / Math.tan(theta);
        xRad = r * Math.sin(phi);
        yRad = -r * Math.cos(phi);
        break;
      }
      case 'SIN': {
        const r = Math.cos(theta);
        xRad = r * Math.sin(phi);
        yRad = -r * Math.cos(phi);
        break;
      }
      case 'ARC': {
        const r = HALF_PI - theta;
        xRad = r * Math.sin(phi);
        yRad = -r * Math.cos(phi);
        break;
      }
      case 'STG': {
        const r = 2.0 * Math.tan(0.5 * (HALF_PI - theta));
        xRad = r * Math.sin(phi);
        yRad = -r * Math.cos(phi);
        break;
      }
      case 'CAR': {
        xRad = phi;
        yRad = theta;
        break;
      }
      case 'CEA': {
        xRad = phi;
        yRad = Math.sin(theta);
        break;
      }
      case 'MER': {
        xRad = phi;
        yRad = Math.log(Math.tan(0.25 * Math.PI + 0.5 * theta));
        break;
      }
      case 'AIT': {
        const cosTheta = Math.cos(theta);
        const cosHalfPhi = Math.cos(0.5 * phi);
        const gamma = Math.sqrt(2.0 / Math.max(1e-12, 1.0 + cosTheta * cosHalfPhi));
        xRad = 2.0 * gamma * cosTheta * Math.sin(0.5 * phi);
        yRad = gamma * Math.sin(theta);
        break;
      }
      case 'MOL': {
        const target = Math.PI * Math.sin(theta);
        let gamma = theta;
        for (let iter = 0; iter < 10; iter++) {
          const f = 2.0 * gamma + Math.sin(2.0 * gamma) - target;
          const fPrime = 2.0 + 2.0 * Math.cos(2.0 * gamma);
          if (Math.abs(f) < 1e-12) break;
          gamma -= f / Math.max(1e-8, fPrime);
        }
        xRad = (2.0 * Math.SQRT2 / Math.PI) * phi * Math.cos(gamma);
        yRad = Math.SQRT2 * Math.sin(gamma);
        break;
      }
      case 'PAR': {
        xRad = (2.0 / Math.PI) * phi * (2.0 * Math.cos(2.0 * theta / 3.0) - 1.0);
        yRad = 3.0 * Math.sin(theta / 3.0);
        break;
      }
      default: {
        const r = 1.0 / Math.tan(theta);
        xRad = r * Math.sin(phi);
        yRad = -r * Math.cos(phi);
        break;
      }
    }

    return { xDeg: xRad * RAD2DEG, yDeg: yRad * RAD2DEG };
  }

  nativeToCelestial(phiDeg, thetaDeg) {
    const alpha0 = this.wcs.crval1 * DEG2RAD;
    const delta0 = this.wcs.crval2 * DEG2RAD;
    const phi = phiDeg * DEG2RAD;
    const theta = thetaDeg * DEG2RAD;

    const isZenithal = ['TAN', 'SIN', 'ARC', 'STG', 'ZPN'].includes(this.wcs.projectionCode);

    if (isZenithal) {
      const sinTheta = Math.sin(theta);
      const cosTheta = Math.cos(theta);
      const sinDelta0 = Math.sin(delta0);
      const cosDelta0 = Math.cos(delta0);
      const cosPhi = Math.cos(phi);
      const sinPhi = Math.sin(phi);

      const sinDelta = sinTheta * sinDelta0 + cosTheta * cosDelta0 * cosPhi;
      const delta = Math.asin(Math.max(-1.0, Math.min(1.0, sinDelta)));

      const y = -cosTheta * sinPhi;
      const x = sinTheta * cosDelta0 - cosTheta * sinDelta0 * cosPhi;
      const alpha = alpha0 + Math.atan2(y, x);

      return {
        worldLonDeg: wrapAngleDeg(alpha * RAD2DEG, 0, 360),
        worldLatDeg: clampLatitudeDeg(delta * RAD2DEG)
      };
    } else {
      const sinTheta = Math.sin(theta);
      const cosTheta = Math.cos(theta);
      const sinDelta0 = Math.sin(delta0);
      const cosDelta0 = Math.cos(delta0);
      const cosPhi = Math.cos(phi);
      const sinPhi = Math.sin(phi);

      const sinDelta = sinTheta * cosDelta0 + cosTheta * sinDelta0 * cosPhi;
      const delta = Math.asin(Math.max(-1.0, Math.min(1.0, sinDelta)));

      const y = cosTheta * sinPhi;
      const x = cosTheta * cosDelta0 * cosPhi - sinTheta * sinDelta0;
      const alpha = alpha0 + Math.atan2(y, x);

      return {
        worldLonDeg: wrapAngleDeg(alpha * RAD2DEG, 0, 360),
        worldLatDeg: clampLatitudeDeg(delta * RAD2DEG)
      };
    }
  }

  celestialToNative(worldLonDeg, worldLatDeg) {
    const alpha0 = this.wcs.crval1 * DEG2RAD;
    const delta0 = this.wcs.crval2 * DEG2RAD;
    const alpha = worldLonDeg * DEG2RAD;
    const delta = worldLatDeg * DEG2RAD;

    const isZenithal = ['TAN', 'SIN', 'ARC', 'STG', 'ZPN'].includes(this.wcs.projectionCode);

    const dAlpha = alpha - alpha0;
    const sinDelta = Math.sin(delta);
    const cosDelta = Math.cos(delta);
    const sinDelta0 = Math.sin(delta0);
    const cosDelta0 = Math.cos(delta0);
    const cosDAlpha = Math.cos(dAlpha);
    const sinDAlpha = Math.sin(dAlpha);

    if (isZenithal) {
      const sinTheta = sinDelta * sinDelta0 + cosDelta * cosDelta0 * cosDAlpha;
      const theta = Math.asin(Math.max(-1.0, Math.min(1.0, sinTheta)));

      const y = -cosDelta * sinDAlpha;
      const x = sinDelta * cosDelta0 - cosDelta * sinDelta0 * cosDAlpha;
      const phi = Math.atan2(y, x);

      return {
        phiDeg: wrapAngleDeg(phi * RAD2DEG, -180, 180),
        thetaDeg: clampLatitudeDeg(theta * RAD2DEG)
      };
    } else {
      const sinTheta = sinDelta * cosDelta0 - cosDelta * sinDelta0 * cosDAlpha;
      const theta = Math.asin(Math.max(-1.0, Math.min(1.0, sinTheta)));

      const y = cosDelta * sinDAlpha;
      const x = cosDelta * cosDelta0 * cosDAlpha + sinDelta * sinDelta0;
      const phi = Math.atan2(y, x);

      return {
        phiDeg: wrapAngleDeg(phi * RAD2DEG, -180, 180),
        thetaDeg: clampLatitudeDeg(theta * RAD2DEG)
      };
    }
  }

  pixelToWorld(xPix, yPix, origin = 1) {
    const { x_int, y_int } = this.pixelToIntermediate(xPix, yPix, origin);
    const { phiDeg, thetaDeg } = this.intermediateToNativeSpherical(x_int, y_int, this.wcs.projectionCode);
    const { worldLonDeg, worldLatDeg } = this.nativeToCelestial(phiDeg, thetaDeg);
    return { worldLon: worldLonDeg, worldLat: worldLatDeg };
  }

  worldToPixel(worldLonDeg, worldLatDeg, origin = 1) {
    const { phiDeg, thetaDeg } = this.celestialToNative(worldLonDeg, worldLatDeg);
    const { xDeg, yDeg } = this.nativeSphericalToIntermediate(phiDeg, thetaDeg, this.wcs.projectionCode);
    return this.intermediateToPixel(xDeg, yDeg, origin);
  }

  computeFootprint(origin = 1) {
    const w = this.wcs.naxis1;
    const h = this.wcs.naxis2;
    const minP = origin === 0 ? 0 : 1;
    const maxPx = origin === 0 ? w : w;
    const maxPy = origin === 0 ? h : h;

    const corners = [
      { xPix: minP, yPix: minP },
      { xPix: maxPx, yPix: minP },
      { xPix: maxPx, yPix: maxPy },
      { xPix: minP, yPix: maxPy }
    ];

    return corners.map(c => {
      const world = this.pixelToWorld(c.xPix, c.yPix, origin);
      return { ...c, worldLon: world.worldLon, worldLat: world.worldLat };
    });
  }

  exportToFITSHeaderCards() {
    const formatCard = (key, value, comment = '') => {
      const paddedKey = (key + '        ').substring(0, 8);
      let valStr = '';
      if (typeof value === 'string') {
        valStr = `'${value}'`;
        valStr = (valStr + '                    ').substring(0, 20);
      } else if (typeof value === 'boolean') {
        valStr = value ? 'T' : 'F';
        valStr = ('                   ' + valStr).slice(-20);
      } else if (typeof value === 'number') {
        valStr = Number.isInteger(value) ? value.toString() : value.toExponential(10);
        valStr = ('                   ' + valStr).slice(-20);
      }
      const commentStr = comment ? ` / ${comment}` : '';
      const card = `${paddedKey}= ${valStr}${commentStr}`;
      return (card + ' '.repeat(80)).substring(0, 80);
    };

    const cards = [
      formatCard('WCSAXES', 2, 'Number of WCS axes'),
      formatCard('CRPIX1', this.wcs.crpix1, 'X reference pixel'),
      formatCard('CRPIX2', this.wcs.crpix2, 'Y reference pixel'),
      formatCard('CRVAL1', this.wcs.crval1, 'Reference longitude (deg)'),
      formatCard('CRVAL2', this.wcs.crval2, 'Reference latitude (deg)'),
      formatCard('CTYPE1', this.wcs.ctype1, 'Coordinate type for axis 1'),
      formatCard('CTYPE2', this.wcs.ctype2, 'Coordinate type for axis 2'),
      formatCard('CD1_1', this.wcs.cd[0], 'Linear transformation matrix'),
      formatCard('CD1_2', this.wcs.cd[1], 'Linear transformation matrix'),
      formatCard('CD2_1', this.wcs.cd[2], 'Linear transformation matrix'),
      formatCard('CD2_2', this.wcs.cd[3], 'Linear transformation matrix'),
      formatCard('LONPOLE', this.wcs.lonpole, 'Native longitude of celestial pole'),
      formatCard('LATPOLE', this.wcs.latpole, 'Native latitude of celestial pole')
    ];

    return cards;
  }
}


// ============================================================================
// 9. MAJOR COSMOGRAPHY HUBS BENCHMARK VALIDATION CATALOG
// ============================================================================

/**
 * Standard Cosmography Benchmark Hubs (Tully et al. 2014, 2023 CF4 Catalog).
 */
export const COSMOGRAPHY_HUBS = Object.freeze([
  {
    name: 'Virgo Cluster Core (M87)',
    eq: { ra: 187.70593, dec: 12.39112, distMpc: 16.5, cz: 1307.0 },
    gal: { l: 283.778, b: 74.491 },
    sg: { sgx: -3.675, sgy: 16.071, sgz: -0.676 }
  },
  {
    name: 'Coma Cluster (Abell 1656)',
    eq: { ra: 194.898, dec: 27.959, distMpc: 100.0, cz: 6925.0 },
    gal: { l: 58.085, b: 88.011 },
    sg: { sgx: 0.645, sgy: 98.957, sgz: 14.391 }
  },
  {
    name: 'Centaurus Cluster (Abell 3526)',
    eq: { ra: 192.20, dec: -41.31, distMpc: 49.5, cz: 3041.0 },
    gal: { l: 302.399, b: 21.559 },
    sg: { sgx: -44.474, sgy: 19.388, sgz: -9.818 }
  },
  {
    name: 'Norma Cluster (Abell 3627 / Great Attractor Core)',
    eq: { ra: 243.83, dec: -60.83, distMpc: 68.0, cz: 4871.0 },
    gal: { l: 325.366, b: -7.180 },
    sg: { sgx: -66.811, sgy: -9.481, sgz: 8.392 }
  },
  {
    name: 'Perseus Cluster (Abell 426)',
    eq: { ra: 49.95, dec: 41.51, distMpc: 75.0, cz: 5366.0 },
    gal: { l: 150.576, b: -13.263 },
    sg: { sgx: 71.069, sgy: -15.266, sgz: -18.470 }
  },
  {
    name: 'Shapley Supercluster Core (Abell 3558)',
    eq: { ra: 201.99, dec: -31.50, distMpc: 200.0, cz: 14450.0 },
    gal: { l: 311.988, b: 30.729 },
    sg: { sgx: -171.161, sgy: 103.349, sgz: -4.778 }
  }
]);
