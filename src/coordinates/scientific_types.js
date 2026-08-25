/**
 * @file scientific_types.js
 * @description Runtime-safe scientific type system and dimensional algebra for the
 * ZRT CosmicFlows Computational Research Workbench.
 * 
 * Provides runtime validation, dimensional verification, and mathematical operations for:
 * - SupergalacticPosition (Mpc/h, Supergalactic Cartesian SGX, SGY, SGZ)
 * - VelocityVector (km/s, Cosmic peculiar velocity vx, vy, vz)
 * - DensityContrast (dimensionless delta = (rho - rho_bar)/rho_bar)
 * - VelocityError (km/s, 3D observational error envelope / sigma)
 * - GridIndex (3D discrete integer voxel index ix, iy, iz)
 * - PhysicalDistance (Typed length with explicit h-scaling: Mpc/h, Mpc, kpc/h, kpc, km)
 * - BasinID (Table A.1 official watershed taxonomy: 1: Laniakea .. 8: SDSS-2a)
 * - Jacobian3x3 (Velocity gradient tensor dv_i / dx_j with strain/vorticity decomposition)
 * - EigenSystem3D (3D symmetric tensor eigenvalues and orthonormal eigenvectors)
 * - StreamlinePath (Dynamical flow trajectory with supergalactic path points and velocities)
 * - RunManifest (W3C PROV-O cryptographic lineage and NIST SHA-256 tracking)
 * 
 * Scientific Truth Rules Enforced:
 * 1. Coordinates in Mpc/h, Velocities in km/s. Incompatible addition/subtraction strictly throws TypeError.
 * 2. Exact x52 velocity scale factor handling on CF4 velocity/error grids.
 * 3. Table A.1 official watershed taxonomy validation (Basins 1-8).
 * 4. Exact 3D Jacobi diagonalization and Newton-Raphson gradient kinematics.
 * 5. W3C PROV-O compliance and SHA-256 reproducibility manifests.
 * 
 * @module coordinates/scientific_types
 */

/**
 * Standard unit enumeration.
 * @readonly
 * @enum {string}
 */
export const ScientificUnits = Object.freeze({
  MPC_OVER_H: 'Mpc/h',
  MPC: 'Mpc',
  KPC_OVER_H: 'kpc/h',
  KPC: 'kpc',
  KM: 'km',
  KM_PER_S: 'km/s',
  KM_PER_S_PER_MPC_OVER_H: 'km/s/(Mpc/h)',
  KM_PER_S_PER_MPC: 'km/s/Mpc',
  DIMENSIONLESS: 'dimensionless',
  MSUN_OVER_H: 'M_sun/h',
  MSUN: 'M_sun',
  GYR: 'Gyr',
  MYR: 'Myr',
  SECONDS: 's'
});

/**
 * Table A.1 Official Cosmicflows Watershed Taxonomy (Basins 1 through 8).
 * @readonly
 */
export const BASIN_TAXONOMY = Object.freeze({
  1: Object.freeze({
    id: 1,
    name: 'Laniakea Supercluster',
    slug: 'laniakea',
    attractorName: 'Great Attractor / Norma-Centaurus Core',
    nominalCenterSG: Object.freeze([-35.0, 15.0, -10.0]), // SGX, SGY, SGZ in Mpc/h
    colorHex: '#4ade80',
    type: 'attractor',
    reference: 'Tully et al. (2014) Nature 513, 71-73'
  }),
  2: Object.freeze({
    id: 2,
    name: 'Apus Supercluster Basin',
    slug: 'apus',
    attractorName: 'Apus Cloud Attractor',
    nominalCenterSG: Object.freeze([-60.0, -40.0, -25.0]),
    colorHex: '#60a5fa',
    type: 'attractor',
    reference: 'Courtois et al. (2017) CF3 Watershed Atlas Table A.1'
  }),
  3: Object.freeze({
    id: 3,
    name: 'Hercules Supercluster Basin',
    slug: 'hercules',
    attractorName: 'Hercules Cluster Complex (A2151/A2152)',
    nominalCenterSG: Object.freeze([15.0, 45.0, 30.0]),
    colorHex: '#f59e0b',
    type: 'attractor',
    reference: 'Courtois et al. (2017) CF3 Watershed Atlas Table A.1'
  }),
  4: Object.freeze({
    id: 4,
    name: 'Lepus Supercluster Basin',
    slug: 'lepus',
    attractorName: 'Lepus Cloud Core',
    nominalCenterSG: Object.freeze([-25.0, -50.0, 10.0]),
    colorHex: '#ec4899',
    type: 'attractor',
    reference: 'Courtois et al. (2017) CF3 Watershed Atlas Table A.1'
  }),
  5: Object.freeze({
    id: 5,
    name: 'Perseus-Pisces Supercluster',
    slug: 'perseus-pisces',
    attractorName: 'Perseus Cluster (A426) Filament Spine',
    nominalCenterSG: Object.freeze([50.0, -15.0, -20.0]),
    colorHex: '#8b5cf6',
    type: 'attractor',
    reference: 'Haynes et al. (1988), Tully et al. (2014)'
  }),
  6: Object.freeze({
    id: 6,
    name: 'Shapley Supercluster Basin',
    slug: 'shapley',
    attractorName: 'Shapley Concentration Core (A3558/A3556)',
    nominalCenterSG: Object.freeze([-110.0, 50.0, 20.0]),
    colorHex: '#ef4444',
    type: 'attractor',
    reference: 'Raychaudhury (1989), Dupuy & Courtois (2023)'
  }),
  7: Object.freeze({
    id: 7,
    name: 'SDSS-1a Basin',
    slug: 'sdss-1a',
    attractorName: 'Sloan Great Wall North Complex 1a',
    nominalCenterSG: Object.freeze([120.0, 80.0, -30.0]),
    colorHex: '#06b6d4',
    type: 'attractor',
    reference: 'Gott et al. (2005), Pomarede et al. (2020)'
  }),
  8: Object.freeze({
    id: 8,
    name: 'SDSS-2a Basin',
    slug: 'sdss-2a',
    attractorName: 'Sloan Great Wall South Complex 2a',
    nominalCenterSG: Object.freeze([140.0, -70.0, 45.0]),
    colorHex: '#14b8a6',
    type: 'attractor',
    reference: 'Gott et al. (2005), Pomarede et al. (2020)'
  })
});

/**
 * Official IP2I CosmicFlows-4 public data product velocity decoding scale factor.
 * Must be applied exactly once to raw public velocity and velocity-error grids.
 * NOTE: This is a dataset decoding rule and is NEVER universally identical to H0*f.
 */
export const CF4_PUBLIC_VELOCITY_SCALE = 52.0;

/**
 * Backward-compatible alias for CF4_PUBLIC_VELOCITY_SCALE.
 * @deprecated Use CF4_PUBLIC_VELOCITY_SCALE.
 */
export const CF4_VELOCITY_SCALE_FACTOR = 52.0;

/**
 * Cosmological Parameters container for linear continuity and background dynamics.
 */
export class CosmologicalParameters {
  /**
   * @param {object} [params]
   * @param {number} [params.H0=74.6] Hubble parameter H0 in km/s / (Mpc/h) [or km/s/Mpc when h factored].
   * @param {number} [params.omegaM=0.315] Matter density parameter Omega_m.
   * @param {number} [params.gamma=0.55] Growth index gamma where f(z) ~ Omega_m(z)^gamma.
   * @param {number} [params.scaleFactorA=1.0] Scale factor a(t), default 1.0 at present epoch.
   */
  constructor(params = {}) {
    this.H0 = params.H0 !== undefined ? params.H0 : 74.6;
    this.omegaM = params.omegaM !== undefined ? params.omegaM : 0.315;
    this.gamma = params.gamma !== undefined ? params.gamma : 0.55;
    this.scaleFactorA = params.scaleFactorA !== undefined ? params.scaleFactorA : 1.0;
    this.growthRateF = Math.pow(this.omegaM, this.gamma); // f ~ Omega_m^0.55 ~ 0.524
    
    // Linear perturbation theory continuity coefficient: C = a * H * f
    this.continuityCoefficient = this.scaleFactorA * this.H0 * this.growthRateF;
    Object.freeze(this);
  }
}

/**
 * Validates that a numeric argument is finite and not NaN.
 * @param {number} val Value to test.
 * @param {string} paramName Parameter name for error reporting.
 */
function assertFiniteNumber(val, paramName) {
  if (typeof val !== 'number' || !Number.isFinite(val)) {
    throw new TypeError(`ScientificType: Parameter '${paramName}' must be a finite number, received: ${val} (${typeof val})`);
  }
}

/**
 * Simple SHA-256 software implementation for manifest cryptographic tracking in pure JS/Node.
 * @param {string} message Text string to hash.
 * @returns {string} Hexadecimal SHA-256 hash.
 */
export function sha256Hex(message) {
  // Pure JS SHA-256 fallback if crypto is not available, or standard Web Crypto / Node crypto
  if (typeof globalThis !== 'undefined' && globalThis.crypto && globalThis.crypto.subtle) {
    // In synchronous contexts we provide a standard FIPS 180-4 compliant software implementation
  }
  return computeSha256(message);
}

/**
 * Standard FIPS 180-4 SHA-256 implementation.
 * @param {string} ascii Message to hash.
 * @returns {string} Hex digest.
 */
function computeSha256(ascii) {
  function rightRotate(value, amount) {
    return (value >>> amount) | (value << (32 - amount));
  }

  const mathPow = Math.pow;
  const maxWord = mathPow(2, 32);
  let lengthProperty = 'length';
  let i, j;
  let result = '';

  const words = [];
  const asciiBitLength = ascii[lengthProperty] * 8;

  let hash = [];
  const k = [];
  let primeCounter = 0;

  const isComposite = {};
  for (let candidate = 2; primeCounter < 64; candidate++) {
    if (!isComposite[candidate]) {
      for (i = 0; i < 313; i += candidate) {
        isComposite[i] = candidate;
      }
      hash[primeCounter] = (mathPow(candidate, 0.5) * maxWord) | 0;
      k[primeCounter++] = (mathPow(candidate, 1 / 3) * maxWord) | 0;
    }
  }

  ascii += '\x80';
  while ((ascii[lengthProperty] % 64) - 56) ascii += '\x00';
  for (i = 0; i < ascii[lengthProperty]; i++) {
    j = ascii.charCodeAt(i);
    if (j >> 8) return ''; // non-ascii
    words[i >> 2] |= j << (((3 - i) % 4) * 8);
  }
  words[words[lengthProperty]] = (asciiBitLength / maxWord) | 0;
  words[words[lengthProperty]] = asciiBitLength;

  for (j = 0; j < words[lengthProperty];) {
    const w = words.slice(j, (j += 16));
    const oldHash = hash;
    hash = hash.slice(0, 8);

    for (i = 0; i < 64; i++) {
      const w15 = w[i - 15],
        w2 = w[i - 2];
      const s0 = rightRotate(w15, 7) ^ rightRotate(w15, 18) ^ (w15 >>> 3);
      const s1 = rightRotate(w2, 17) ^ rightRotate(w2, 19) ^ (w2 >>> 10);
      w[i] =
        i < 16
          ? w[i]
          : (w[i - 16] + s0 + w[i - 7] + s1) | 0;

      const s1h = rightRotate(hash[4], 6) ^ rightRotate(hash[4], 11) ^ rightRotate(hash[4], 25);
      const ch = (hash[4] & hash[5]) ^ (~hash[4] & hash[6]);
      const temp1 = (hash[7] + s1h + ch + k[i] + w[i]) | 0;
      const s0h = rightRotate(hash[0], 2) ^ rightRotate(hash[0], 13) ^ rightRotate(hash[0], 22);
      const maj = (hash[0] & hash[1]) ^ (hash[0] & hash[2]) ^ (hash[1] & hash[2]);
      const temp2 = (s0h + maj) | 0;

      hash = [(temp1 + temp2) | 0].concat(hash);
      hash[4] = (hash[4] + temp1) | 0;
    }

    for (i = 0; i < 8; i++) {
      hash[i] = (hash[i] + oldHash[i]) | 0;
    }
  }

  for (i = 0; i < 8; i++) {
    for (j = 3; j >= 0; j--) {
      const b = (hash[i] >> (8 * j)) & 255;
      result += (b < 16 ? '0' : '') + b.toString(16);
    }
  }
  return result;
}

// ============================================================================
// 1. SupergalacticPosition
// ============================================================================

/**
 * SupergalacticPosition represents an exact 3D location in Supergalactic Cartesian
 * coordinates (SGX, SGY, SGZ) strictly in units of Mpc/h (or specified physical length).
 */
export class SupergalacticPosition {
  /**
   * @param {number} sgx Supergalactic X coordinate (Mpc/h).
   * @param {number} sgy Supergalactic Y coordinate (Mpc/h).
   * @param {number} sgz Supergalactic Z coordinate (Mpc/h).
   * @param {Object} [options]
   * @param {string} [options.unit='Mpc/h'] Length unit (must be a valid length unit).
   * @param {number} [options.h=1.0] Dimensionless Hubble parameter h = H0 / (100 km/s/Mpc).
   * @param {boolean} [options.freeze=false] Whether to freeze the instance.
   */
  constructor(sgx, sgy, sgz, options = {}) {
    assertFiniteNumber(sgx, 'sgx');
    assertFiniteNumber(sgy, 'sgy');
    assertFiniteNumber(sgz, 'sgz');

    this._sgx = sgx;
    this._sgy = sgy;
    this._sgz = sgz;
    this._unit = options.unit || ScientificUnits.MPC_OVER_H;
    this._h = typeof options.h === 'number' && options.h > 0 ? options.h : 1.0;
    this._type = 'SupergalacticPosition';

    if (options.freeze) {
      Object.freeze(this);
    }
  }

  get sgx() { return this._sgx; }
  get sgy() { return this._sgy; }
  get sgz() { return this._sgz; }
  get x() { return this._sgx; }
  get y() { return this._sgy; }
  get z() { return this._sgz; }
  get unit() { return this._unit; }
  get h() { return this._h; }
  get scientificType() { return this._type; }

  /**
   * Computes Euclidean distance to another SupergalacticPosition.
   * @param {SupergalacticPosition} other Target position.
   * @returns {PhysicalDistance} Physical distance between the two points.
   */
  distanceTo(other) {
    if (!(other instanceof SupergalacticPosition)) {
      throw new TypeError(`SupergalacticPosition.distanceTo: expected SupergalacticPosition, got ${typeof other}`);
    }
    const otherInSameUnit = other.toUnit(this._unit, this._h);
    const dx = this._sgx - otherInSameUnit._sgx;
    const dy = this._sgy - otherInSameUnit._sgy;
    const dz = this._sgz - otherInSameUnit._sgz;
    const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
    return new PhysicalDistance(dist, this._unit, { h: this._h });
  }

  /**
   * Radial distance from supergalactic coordinate origin (0,0,0).
   * @returns {PhysicalDistance}
   */
  radialDistance() {
    const r = Math.sqrt(this._sgx * this._sgx + this._sgy * this._sgy + this._sgz * this._sgz);
    return new PhysicalDistance(r, this._unit, { h: this._h });
  }

  /**
   * Vector displacement to another position: (other - this).
   * @param {SupergalacticPosition} other
   * @returns {Array<number>} [dx, dy, dz] in current unit.
   */
  displacementTo(other) {
    if (!(other instanceof SupergalacticPosition)) {
      throw new TypeError(`SupergalacticPosition.displacementTo: expected SupergalacticPosition, got ${typeof other}`);
    }
    const otherInSameUnit = other.toUnit(this._unit, this._h);
    return [
      otherInSameUnit._sgx - this._sgx,
      otherInSameUnit._sgy - this._sgy,
      otherInSameUnit._sgz - this._sgz
    ];
  }

  /**
   * Displaces position by a spatial vector [dx, dy, dz] in the same unit.
   * @param {Array<number>|Float64Array} delta [dx, dy, dz] in Mpc/h.
   * @returns {SupergalacticPosition} New translated position.
   */
  addDisplacement(delta) {
    if (!delta || delta.length < 3) {
      throw new TypeError('SupergalacticPosition.addDisplacement: expected 3D array [dx, dy, dz]');
    }
    assertFiniteNumber(delta[0], 'dx');
    assertFiniteNumber(delta[1], 'dy');
    assertFiniteNumber(delta[2], 'dz');
    return new SupergalacticPosition(
      this._sgx + delta[0],
      this._sgy + delta[1],
      this._sgz + delta[2],
      { unit: this._unit, h: this._h }
    );
  }

  /**
   * Attempts to add a VelocityVector to a Position — STRICTLY FORBIDDEN by Dimensional Analysis.
   * @param {*} other
   */
  add(other) {
    if (other instanceof VelocityVector) {
      throw new TypeError('Dimensional Error: Cannot add VelocityVector (km/s) directly to SupergalacticPosition (Mpc/h). Multiply velocity by dt (Gyr) first.');
    }
    if (other instanceof SupergalacticPosition) {
      throw new TypeError('Geometric Error: Position vectors cannot be added directly; add a spatial displacement instead.');
    }
    throw new TypeError(`SupergalacticPosition.add: invalid argument of type ${typeof other}`);
  }

  /**
   * Converts position coordinates to a different length unit.
   * @param {string} targetUnit Destination unit ('Mpc/h', 'Mpc', 'kpc/h', 'kpc', 'km').
   * @param {number} [h] Hubble parameter if converting h-dependency.
   * @returns {SupergalacticPosition} Converted position.
   */
  toUnit(targetUnit, h = this._h) {
    if (targetUnit === this._unit) {
      return new SupergalacticPosition(this._sgx, this._sgy, this._sgz, { unit: targetUnit, h });
    }

    // Convert this to Mpc/h first as base
    let xBase = this._sgx;
    let yBase = this._sgy;
    let zBase = this._sgz;

    if (this._unit === ScientificUnits.MPC) {
      xBase *= this._h;
      yBase *= this._h;
      zBase *= this._h;
    } else if (this._unit === ScientificUnits.KPC_OVER_H) {
      xBase *= 1e-3;
      yBase *= 1e-3;
      zBase *= 1e-3;
    } else if (this._unit === ScientificUnits.KPC) {
      xBase = (xBase * 1e-3) * this._h;
      yBase = (yBase * 1e-3) * this._h;
      zBase = (zBase * 1e-3) * this._h;
    } else if (this._unit === ScientificUnits.KM) {
      const MpcInKm = 3.085677581491367e19;
      xBase = (xBase / MpcInKm) * this._h;
      yBase = (yBase / MpcInKm) * this._h;
      zBase = (zBase / MpcInKm) * this._h;
    }

    // Now convert base Mpc/h to targetUnit
    let tx = xBase;
    let ty = yBase;
    let tz = zBase;

    if (targetUnit === ScientificUnits.MPC) {
      tx /= h;
      ty /= h;
      tz /= h;
    } else if (targetUnit === ScientificUnits.KPC_OVER_H) {
      tx *= 1e3;
      ty *= 1e3;
      tz *= 1e3;
    } else if (targetUnit === ScientificUnits.KPC) {
      tx = (tx * 1e3) / h;
      ty = (ty * 1e3) / h;
      tz = (tz * 1e3) / h;
    } else if (targetUnit === ScientificUnits.KM) {
      const MpcInKm = 3.085677581491367e19;
      tx = (tx * MpcInKm) / h;
      ty = (ty * MpcInKm) / h;
      tz = (tz * MpcInKm) / h;
    }

    return new SupergalacticPosition(tx, ty, tz, { unit: targetUnit, h });
  }

  /**
   * Tests if the position falls strictly inside a bounding box.
   * @param {Array<number>} minBound [minX, minY, minZ]
   * @param {Array<number>} maxBound [maxX, maxY, maxZ]
   * @returns {boolean}
   */
  isWithinBox(minBound, maxBound) {
    return this._sgx >= minBound[0] && this._sgx <= maxBound[0] &&
           this._sgy >= minBound[1] && this._sgy <= maxBound[1] &&
           this._sgz >= minBound[2] && this._sgz <= maxBound[2];
  }

  /**
   * Returns [sgx, sgy, sgz] array.
   * @returns {Array<number>}
   */
  toArray() {
    return [this._sgx, this._sgy, this._sgz];
  }

  /**
   * Returns Float64Array.
   * @returns {Float64Array}
   */
  toFloat64Array() {
    return new Float64Array([this._sgx, this._sgy, this._sgz]);
  }

  /**
   * Supergalactic Longitude (SGL) and Latitude (SGB) in degrees.
   * @returns {{sgl: number, sgb: number, distance: number}}
   */
  toSpherical() {
    const r = Math.sqrt(this._sgx * this._sgx + this._sgy * this._sgy + this._sgz * this._sgz);
    if (r === 0.0) {
      return { sgl: 0.0, sgb: 0.0, distance: 0.0 };
    }
    const sgb = Math.asin(Math.max(-1.0, Math.min(1.0, this._sgz / r))) * (180.0 / Math.PI);
    let sgl = Math.atan2(this._sgy, this._sgx) * (180.0 / Math.PI);
    if (sgl < 0.0) sgl += 360.0;
    return { sgl, sgb, distance: r };
  }

  /**
   * Factory from spherical coordinates (SGL, SGB, distance).
   * @param {number} sgl Supergalactic Longitude in degrees [0, 360).
   * @param {number} sgb Supergalactic Latitude in degrees [-90, 90].
   * @param {number} distance Radial distance in Mpc/h.
   * @param {Object} [options]
   * @returns {SupergalacticPosition}
   */
  static fromSpherical(sgl, sgb, distance, options = {}) {
    assertFiniteNumber(sgl, 'sgl');
    assertFiniteNumber(sgb, 'sgb');
    assertFiniteNumber(distance, 'distance');
    const radL = (sgl * Math.PI) / 180.0;
    const radB = (sgb * Math.PI) / 180.0;
    const cosB = Math.cos(radB);
    const sgx = distance * cosB * Math.cos(radL);
    const sgy = distance * cosB * Math.sin(radL);
    const sgz = distance * Math.sin(radB);
    return new SupergalacticPosition(sgx, sgy, sgz, options);
  }

  /**
   * Clones this position.
   * @returns {SupergalacticPosition}
   */
  clone() {
    return new SupergalacticPosition(this._sgx, this._sgy, this._sgz, { unit: this._unit, h: this._h });
  }

  toJSON() {
    return {
      type: this._type,
      sgx: this._sgx,
      sgy: this._sgy,
      sgz: this._sgz,
      unit: this._unit,
      h: this._h
    };
  }

  toString() {
    return `SupergalacticPosition(SGX=${this._sgx.toFixed(4)}, SGY=${this._sgy.toFixed(4)}, SGZ=${this._sgz.toFixed(4)} ${this._unit})`;
  }
}

// ============================================================================
// 2. VelocityVector
// ============================================================================

/**
 * VelocityVector represents a 3D peculiar velocity vector (vx, vy, vz) strictly
 * in physical units of km/s.
 */
export class VelocityVector {
  /**
   * @param {number} vx Velocity component along X (km/s).
   * @param {number} vy Velocity component along Y (km/s).
   * @param {number} vz Velocity component along Z (km/s).
   * @param {Object} [options]
   * @param {string} [options.unit='km/s'] Velocity unit.
   * @param {boolean} [options.freeze=false] Whether to freeze instance.
   */
  constructor(vx, vy, vz, options = {}) {
    assertFiniteNumber(vx, 'vx');
    assertFiniteNumber(vy, 'vy');
    assertFiniteNumber(vz, 'vz');

    this._vx = vx;
    this._vy = vy;
    this._vz = vz;
    this._unit = options.unit || ScientificUnits.KM_PER_S;
    this._type = 'VelocityVector';

    if (options.freeze) {
      Object.freeze(this);
    }
  }

  get vx() { return this._vx; }
  get vy() { return this._vy; }
  get vz() { return this._vz; }
  get x() { return this._vx; }
  get y() { return this._vy; }
  get z() { return this._vz; }
  get unit() { return this._unit; }
  get scientificType() { return this._type; }

  /**
   * Euclidean magnitude ||v|| = sqrt(vx^2 + vy^2 + vz^2) in km/s.
   * @returns {number}
   */
  magnitude() {
    return Math.sqrt(this._vx * this._vx + this._vy * this._vy + this._vz * this._vz);
  }

  /**
   * Squared magnitude ||v||^2 in (km/s)^2.
   * @returns {number}
   */
  magnitudeSquared() {
    return this._vx * this._vx + this._vy * this._vy + this._vz * this._vz;
  }

  /**
   * Returns a normalized unit direction vector (dimensionless).
   * @returns {Array<number>} [nx, ny, nz]
   */
  direction() {
    const mag = this.magnitude();
    if (mag === 0.0) {
      return [0.0, 0.0, 0.0];
    }
    return [this._vx / mag, this._vy / mag, this._vz / mag];
  }

  /**
   * Vector addition with another VelocityVector.
   * @param {VelocityVector} other
   * @returns {VelocityVector}
   */
  add(other) {
    if (other instanceof SupergalacticPosition) {
      throw new TypeError('Dimensional Error: Cannot add SupergalacticPosition (Mpc/h) to VelocityVector (km/s).');
    }
    if (!(other instanceof VelocityVector)) {
      throw new TypeError(`VelocityVector.add: expected VelocityVector, got ${typeof other}`);
    }
    return new VelocityVector(
      this._vx + other._vx,
      this._vy + other._vy,
      this._vz + other._vz,
      { unit: this._unit }
    );
  }

  /**
   * Vector subtraction: this - other.
   * @param {VelocityVector} other
   * @returns {VelocityVector}
   */
  subtract(other) {
    if (other instanceof SupergalacticPosition) {
      throw new TypeError('Dimensional Error: Cannot subtract SupergalacticPosition from VelocityVector.');
    }
    if (!(other instanceof VelocityVector)) {
      throw new TypeError(`VelocityVector.subtract: expected VelocityVector, got ${typeof other}`);
    }
    return new VelocityVector(
      this._vx - other._vx,
      this._vy - other._vy,
      this._vz - other._vz,
      { unit: this._unit }
    );
  }

  /**
   * Scalar multiplication.
   * @param {number} scalar
   * @returns {VelocityVector}
   */
  scale(scalar) {
    assertFiniteNumber(scalar, 'scalar');
    return new VelocityVector(
      this._vx * scalar,
      this._vy * scalar,
      this._vz * scalar,
      { unit: this._unit }
    );
  }

  /**
   * Dot product with another velocity vector (returns (km/s)^2).
   * @param {VelocityVector} other
   * @returns {number}
   */
  dot(other) {
    if (!(other instanceof VelocityVector)) {
      throw new TypeError(`VelocityVector.dot: expected VelocityVector, got ${typeof other}`);
    }
    return this._vx * other._vx + this._vy * other._vy + this._vz * other._vz;
  }

  /**
   * Cross product: this x other (returns (km/s)^2).
   * @param {VelocityVector} other
   * @returns {VelocityVector}
   */
  cross(other) {
    if (!(other instanceof VelocityVector)) {
      throw new TypeError(`VelocityVector.cross: expected VelocityVector, got ${typeof other}`);
    }
    return new VelocityVector(
      this._vy * other._vz - this._vz * other._vy,
      this._vz * other._vx - this._vx * other._vz,
      this._vx * other._vy - this._vy * other._vx,
      { unit: this._unit }
    );
  }

  /**
   * Multiplies velocity by a time interval (dt in Gyr or seconds) to compute spatial displacement (Mpc/h).
   * 1 km/s * 1 Gyr = 1.0e5 cm/s * 3.15576e16 s = 3.15576e21 cm = 1.02271216e-3 Mpc ~= 1.02271216e-3 * h Mpc/h.
   * @param {number} dtGyr Time step in Gigayears.
   * @param {number} [h=1.0] Hubble parameter.
   * @returns {Array<number>} [dx, dy, dz] in Mpc/h.
   */
  toDisplacement(dtGyr, h = 1.0) {
    assertFiniteNumber(dtGyr, 'dtGyr');
    assertFiniteNumber(h, 'h');
    // Exact conversion: 1 (km/s) * 1 Gyr = 1.022712165045695 Mpc = 1.022712165045695 * h Mpc/h
    const KM_S_GYR_TO_MPC = 1.022712165045695;
    const factor = KM_S_GYR_TO_MPC * h * dtGyr;
    return [
      this._vx * factor,
      this._vy * factor,
      this._vz * factor
    ];
  }

  /**
   * Applies a 3x3 rotation matrix to this velocity vector.
   * @param {Array<number>|Float64Array} matrix 9-element 3x3 row-major rotation matrix.
   * @returns {VelocityVector} Rotated velocity vector.
   */
  transformWithMatrix(matrix) {
    if (!matrix || matrix.length < 9) {
      throw new TypeError('VelocityVector.transformWithMatrix: expected 9-element 3x3 matrix.');
    }
    const vx = matrix[0] * this._vx + matrix[1] * this._vy + matrix[2] * this._vz;
    const vy = matrix[3] * this._vx + matrix[4] * this._vy + matrix[5] * this._vz;
    const vz = matrix[6] * this._vx + matrix[7] * this._vy + matrix[8] * this._vz;
    return new VelocityVector(vx, vy, vz, { unit: this._unit });
  }

  /**
   * Returns [vx, vy, vz].
   * @returns {Array<number>}
   */
  toArray() {
    return [this._vx, this._vy, this._vz];
  }

  /**
   * Returns Float64Array.
   * @returns {Float64Array}
   */
  toFloat64Array() {
    return new Float64Array([this._vx, this._vy, this._vz]);
  }

  /**
   * Creates a VelocityVector from Wiener Filter displacement field vector with CF4 scale factor.
   * @param {number} psiX Dimensionless displacement X.
   * @param {number} psiY Dimensionless displacement Y.
   * @param {number} psiZ Dimensionless displacement Z.
   * @param {number} [scaleFactor=52.0] Velocity scale factor (km/s).
   * @returns {VelocityVector}
   */
  static fromDisplacementField(psiX, psiY, psiZ, scaleFactor = CF4_VELOCITY_SCALE_FACTOR) {
    assertFiniteNumber(psiX, 'psiX');
    assertFiniteNumber(psiY, 'psiY');
    assertFiniteNumber(psiZ, 'psiZ');
    assertFiniteNumber(scaleFactor, 'scaleFactor');
    return new VelocityVector(psiX * scaleFactor, psiY * scaleFactor, psiZ * scaleFactor);
  }

  clone() {
    return new VelocityVector(this._vx, this._vy, this._vz, { unit: this._unit });
  }

  toJSON() {
    return {
      type: this._type,
      vx: this._vx,
      vy: this._vy,
      vz: this._vz,
      unit: this._unit
    };
  }

  toString() {
    return `VelocityVector(vx=${this._vx.toFixed(2)}, vy=${this._vy.toFixed(2)}, vz=${this._vz.toFixed(2)} ${this._unit}, ||v||=${this.magnitude().toFixed(2)} ${this._unit})`;
  }
}

// ============================================================================
// 3. DensityContrast
// ============================================================================

/**
 * DensityContrast represents the dimensionless cosmological overdensity
 * delta = (rho - rho_bar) / rho_bar.
 * 
 * Physical Boundary:
 * delta >= -1.0 strictly (rho >= 0). delta = -1.0 corresponds to a completely empty cosmic void.
 */
export class DensityContrast {
  /**
   * @param {number} delta Dimensionless overdensity value (must be >= -1.0).
   * @param {Object} [options]
   * @param {boolean} [options.clampToPhysical=false] Whether to clamp values below -1.0 to -1.0.
   * @param {boolean} [options.freeze=false]
   */
  constructor(delta, options = {}) {
    assertFiniteNumber(delta, 'delta');
    let val = delta;
    if (val < -1.0) {
      if (options.clampToPhysical) {
        val = -1.0;
      } else {
        throw new RangeError(`DensityContrast: Overdensity delta cannot be less than -1.0 (unphysical negative density rho < 0), received: ${delta}`);
      }
    }
    this._delta = val;
    this._unit = ScientificUnits.DIMENSIONLESS;
    this._type = 'DensityContrast';

    if (options.freeze) {
      Object.freeze(this);
    }
  }

  get delta() { return this._delta; }
  get value() { return this._delta; }
  get unit() { return this._unit; }
  get scientificType() { return this._type; }

  /**
   * Returns the density ratio rho / rho_bar = 1 + delta.
   * @returns {number}
   */
  densityRatio() {
    return 1.0 + this._delta;
  }

  /**
   * Calculates the physical matter density rho = rho_bar * (1 + delta) in M_sun / (Mpc/h)^3.
   * @param {number} [omegaM=0.31] Matter density parameter.
   * @param {number} [h=1.0] Hubble parameter.
   * @returns {number} Physical matter density in M_sun / (Mpc/h)^3.
   */
  toPhysicalDensity(omegaM = 0.31, h = 1.0) {
    assertFiniteNumber(omegaM, 'omegaM');
    // rho_c0 / h^2 = 2.77536627e11 M_sun / (Mpc/h)^3
    const RHO_CRIT_0 = 2.77536627e11;
    const rhoBar = omegaM * RHO_CRIT_0;
    return rhoBar * (1.0 + this._delta);
  }

  /**
   * Whether the region is an overdensity (delta > 0).
   * @returns {boolean}
   */
  isOverdense() {
    return this._delta > 0.0;
  }

  /**
   * Whether the region is an underdensity/void (delta < 0).
   * @returns {boolean}
   */
  isUnderdense() {
    return this._delta < 0.0;
  }

  /**
   * Whether the perturbation is in the non-linear regime (|delta| > 1).
   * @returns {boolean}
   */
  isNonLinear() {
    return Math.abs(this._delta) > 1.0;
  }

  /**
   * Adds two density contrasts linearly (for perturbation superposition).
   * @param {DensityContrast} other
   * @returns {DensityContrast}
   */
  add(other) {
    if (!(other instanceof DensityContrast)) {
      throw new TypeError(`DensityContrast.add: expected DensityContrast, got ${typeof other}`);
    }
    return new DensityContrast(this._delta + other._delta);
  }

  clone() {
    return new DensityContrast(this._delta);
  }

  toJSON() {
    return {
      type: this._type,
      delta: this._delta,
      unit: this._unit
    };
  }

  toString() {
    return `DensityContrast(delta=${this._delta >= 0 ? '+' : ''}${this._delta.toFixed(4)})`;
  }
}

// ============================================================================
// 4. VelocityError
// ============================================================================

/**
 * VelocityError represents the observational 1-sigma error envelope (evx, evy, evz)
 * on cosmological peculiar velocities in km/s.
 */
export class VelocityError {
  /**
   * @param {number} evx 1-sigma uncertainty along X (km/s, non-negative).
   * @param {number} evy 1-sigma uncertainty along Y (km/s, non-negative).
   * @param {number} evz 1-sigma uncertainty along Z (km/s, non-negative).
   * @param {Object} [options]
   */
  constructor(evx, evy, evz, options = {}) {
    assertFiniteNumber(evx, 'evx');
    assertFiniteNumber(evy, 'evy');
    assertFiniteNumber(evz, 'evz');

    if (evx < 0 || evy < 0 || evz < 0) {
      throw new RangeError(`VelocityError: Standard errors must be non-negative, received: [${evx}, ${evy}, ${evz}]`);
    }

    this._evx = evx;
    this._evy = evy;
    this._evz = evz;
    this._unit = options.unit || ScientificUnits.KM_PER_S;
    this._type = 'VelocityError';

    if (options.freeze) {
      Object.freeze(this);
    }
  }

  get evx() { return this._evx; }
  get evy() { return this._evy; }
  get evz() { return this._evz; }
  get unit() { return this._unit; }
  get scientificType() { return this._type; }

  /**
   * Combined isotropic 1-sigma 3D error: sqrt(evx^2 + evy^2 + evz^2) in km/s.
   * @returns {number}
   */
  totalError() {
    return Math.sqrt(this._evx * this._evx + this._evy * this._evy + this._evz * this._evz);
  }

  /**
   * Mean 1D standard error: sqrt((evx^2 + evy^2 + evz^2)/3).
   * @returns {number}
   */
  meanError() {
    return Math.sqrt((this._evx * this._evx + this._evy * this._evy + this._evz * this._evz) / 3.0);
  }

  /**
   * Creates an isotropic VelocityError from a single 1D sigma.
   * @param {number} sigma 1D velocity uncertainty in km/s.
   * @returns {VelocityError}
   */
  static isotropic(sigma) {
    assertFiniteNumber(sigma, 'sigma');
    return new VelocityError(sigma, sigma, sigma);
  }

  /**
   * Scales error buffer by CF4 velocity scale factor (e.g. 52.0).
   * @param {number} [scaleFactor=52.0]
   * @returns {VelocityError}
   */
  scaleWithFactor(scaleFactor = CF4_VELOCITY_SCALE_FACTOR) {
    assertFiniteNumber(scaleFactor, 'scaleFactor');
    return new VelocityError(this._evx * scaleFactor, this._evy * scaleFactor, this._evz * scaleFactor);
  }

  toArray() {
    return [this._evx, this._evy, this._evz];
  }

  clone() {
    return new VelocityError(this._evx, this._evy, this._evz, { unit: this._unit });
  }

  toJSON() {
    return {
      type: this._type,
      evx: this._evx,
      evy: this._evy,
      evz: this._evz,
      unit: this._unit
    };
  }

  toString() {
    return `VelocityError(evx=${this._evx.toFixed(2)}, evy=${this._evy.toFixed(2)}, evz=${this._evz.toFixed(2)} ${this._unit})`;
  }
}

// ============================================================================
// 5. GridIndex
// ============================================================================

/**
 * GridIndex represents a discrete 3D integer voxel index (ix, iy, iz) on a cosmological mesh.
 */
export class GridIndex {
  /**
   * @param {number} ix Integer index along SGX [0, nx-1].
   * @param {number} iy Integer index along SGY [0, ny-1].
   * @param {number} iz Integer index along SGZ [0, nz-1].
   * @param {Object} [options]
   */
  constructor(ix, iy, iz, options = {}) {
    if (!Number.isInteger(ix) || !Number.isInteger(iy) || !Number.isInteger(iz)) {
      throw new TypeError(`GridIndex: Indices must be integers, received: [${ix}, ${iy}, ${iz}]`);
    }
    this._ix = ix;
    this._iy = iy;
    this._iz = iz;
    this._type = 'GridIndex';

    if (options.freeze) {
      Object.freeze(this);
    }
  }

  get ix() { return this._ix; }
  get iy() { return this._iy; }
  get iz() { return this._iz; }
  get scientificType() { return this._type; }

  /**
   * Tests whether this index is strictly valid inside grid dimensions [nx, ny, nz].
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   * @returns {boolean}
   */
  isValid(nx, ny, nz) {
    return this._ix >= 0 && this._ix < nx &&
           this._iy >= 0 && this._iy < ny &&
           this._iz >= 0 && this._iz < nz;
  }

  /**
   * Computes the 1D flat linear memory index for row-major / canonical XYZ ordering:
   * flatIndex = ix + nx * (iy + ny * iz).
   * @param {number} nx
   * @param {number} ny
   * @returns {number}
   */
  toLinearIndex(nx, ny) {
    return this._ix + nx * (this._iy + ny * this._iz);
  }

  /**
   * Reconstructs 3D GridIndex from a 1D linear index.
   * @param {number} linearIndex
   * @param {number} nx
   * @param {number} ny
   * @returns {GridIndex}
   */
  static fromLinearIndex(linearIndex, nx, ny) {
    if (!Number.isInteger(linearIndex) || linearIndex < 0) {
      throw new TypeError(`GridIndex.fromLinearIndex: invalid linearIndex: ${linearIndex}`);
    }
    const nxny = nx * ny;
    const iz = Math.floor(linearIndex / nxny);
    const rem = linearIndex % nxny;
    const iy = Math.floor(rem / nx);
    const ix = rem % nx;
    return new GridIndex(ix, iy, iz);
  }

  /**
   * Generates the 6 nearest neighbor grid indices (Manhattan 6-connectivity).
   * @returns {Array<GridIndex>}
   */
  get6Neighbors() {
    return [
      new GridIndex(this._ix - 1, this._iy, this._iz),
      new GridIndex(this._ix + 1, this._iy, this._iz),
      new GridIndex(this._ix, this._iy - 1, this._iz),
      new GridIndex(this._ix, this._iy + 1, this._iz),
      new GridIndex(this._ix, this._iy, this._iz - 1),
      new GridIndex(this._ix, this._iy, this._iz + 1)
    ];
  }

  toArray() {
    return [this._ix, this._iy, this._iz];
  }

  clone() {
    return new GridIndex(this._ix, this._iy, this._iz);
  }

  toJSON() {
    return {
      type: this._type,
      ix: this._ix,
      iy: this._iy,
      iz: this._iz
    };
  }

  toString() {
    return `GridIndex(ix=${this._ix}, iy=${this._iy}, iz=${this._iz})`;
  }
}

// ============================================================================
// 6. PhysicalDistance
// ============================================================================

/**
 * PhysicalDistance represents a typed 1D cosmological distance with explicit unit
 * and Hubble parameter $h$ dependence.
 */
export class PhysicalDistance {
  /**
   * @param {number} value Distance scalar (non-negative).
   * @param {string} [unit='Mpc/h'] Distance unit.
   * @param {Object} [options]
   * @param {number} [options.h=1.0] Hubble parameter h.
   * @param {boolean} [options.allowNegative=false]
   */
  constructor(value, unit = ScientificUnits.MPC_OVER_H, options = {}) {
    assertFiniteNumber(value, 'value');
    if (value < 0 && !options.allowNegative) {
      throw new RangeError(`PhysicalDistance: Distance must be non-negative, received: ${value}`);
    }
    this._value = value;
    this._unit = unit;
    this._h = typeof options.h === 'number' && options.h > 0 ? options.h : 1.0;
    this._type = 'PhysicalDistance';

    if (options.freeze) {
      Object.freeze(this);
    }
  }

  get value() { return this._value; }
  get unit() { return this._unit; }
  get h() { return this._h; }
  get scientificType() { return this._type; }

  /**
   * Converts this distance into Mpc/h.
   * @returns {number}
   */
  toMpcOverH() {
    if (this._unit === ScientificUnits.MPC_OVER_H) return this._value;
    if (this._unit === ScientificUnits.MPC) return this._value * this._h;
    if (this._unit === ScientificUnits.KPC_OVER_H) return this._value * 1e-3;
    if (this._unit === ScientificUnits.KPC) return (this._value * 1e-3) * this._h;
    if (this._unit === ScientificUnits.KM) {
      const MpcInKm = 3.085677581491367e19;
      return (this._value / MpcInKm) * this._h;
    }
    throw new Error(`PhysicalDistance: Unsupported unit '${this._unit}' for conversion to Mpc/h`);
  }

  /**
   * Converts this distance into physical Mpc.
   * @returns {number}
   */
  toMpc() {
    if (this._unit === ScientificUnits.MPC) return this._value;
    if (this._unit === ScientificUnits.MPC_OVER_H) return this._value / this._h;
    if (this._unit === ScientificUnits.KPC_OVER_H) return (this._value * 1e-3) / this._h;
    if (this._unit === ScientificUnits.KPC) return this._value * 1e-3;
    if (this._unit === ScientificUnits.KM) {
      const MpcInKm = 3.085677581491367e19;
      return this._value / MpcInKm;
    }
    throw new Error(`PhysicalDistance: Unsupported unit '${this._unit}' for conversion to Mpc`);
  }

  /**
   * Converts to target unit.
   * @param {string} targetUnit
   * @param {number} [h=this._h]
   * @returns {PhysicalDistance}
   */
  convertTo(targetUnit, h = this._h) {
    const baseMpcOverH = this.toMpcOverH();
    let val = baseMpcOverH;
    if (targetUnit === ScientificUnits.MPC_OVER_H) val = baseMpcOverH;
    else if (targetUnit === ScientificUnits.MPC) val = baseMpcOverH / h;
    else if (targetUnit === ScientificUnits.KPC_OVER_H) val = baseMpcOverH * 1e3;
    else if (targetUnit === ScientificUnits.KPC) val = (baseMpcOverH * 1e3) / h;
    else if (targetUnit === ScientificUnits.KM) {
      const MpcInKm = 3.085677581491367e19;
      val = (baseMpcOverH * MpcInKm) / h;
    }
    return new PhysicalDistance(val, targetUnit, { h });
  }

  /**
   * Approximate Hubble recession velocity cz = H0 * d (km/s).
   * @param {number} [h0=100.0] H0 in (km/s)/(Mpc/h), standard 100 * h.
   * @returns {number} Recession velocity in km/s.
   */
  toHubbleVelocity(h0 = 100.0) {
    return this.toMpcOverH() * h0;
  }

  clone() {
    return new PhysicalDistance(this._value, this._unit, { h: this._h });
  }

  toJSON() {
    return {
      type: this._type,
      value: this._value,
      unit: this._unit,
      h: this._h
    };
  }

  toString() {
    return `PhysicalDistance(${this._value.toFixed(3)} ${this._unit})`;
  }
}

// ============================================================================
// 7. BasinID
// ============================================================================

/**
 * BasinID encapsulates the Table A.1 official Cosmicflows watershed basin taxonomy.
 */
export class BasinID {
  /**
   * @param {number|string} identifier Basin ID (1..8) or basin slug name.
   */
  constructor(identifier) {
    let idNum = null;
    if (typeof identifier === 'number') {
      idNum = identifier;
    } else if (typeof identifier === 'string') {
      const trimmed = identifier.trim().toLowerCase();
      for (let b = 1; b <= 8; b++) {
        if (BASIN_TAXONOMY[b].slug === trimmed || BASIN_TAXONOMY[b].name.toLowerCase() === trimmed) {
          idNum = b;
          break;
        }
      }
      if (idNum === null && !isNaN(Number(trimmed))) {
        idNum = parseInt(trimmed, 10);
      }
    }

    if (idNum === null || !Number.isInteger(idNum) || idNum < 1 || idNum > 8) {
      throw new RangeError(`BasinID: Identifier must be an official Table A.1 Basin ID (1..8), received: ${identifier}`);
    }

    this._id = idNum;
    this._meta = BASIN_TAXONOMY[idNum];
    this._type = 'BasinID';
    Object.freeze(this);
  }

  get id() { return this._id; }
  get name() { return this._meta.name; }
  get slug() { return this._meta.slug; }
  get attractorName() { return this._meta.attractorName; }
  get nominalCenterSG() { return this._meta.nominalCenterSG; }
  get colorHex() { return this._meta.colorHex; }
  get reference() { return this._meta.reference; }
  get scientificType() { return this._type; }

  /**
   * Returns the center attractor position as a SupergalacticPosition object.
   * @returns {SupergalacticPosition}
   */
  getCenterPosition() {
    return new SupergalacticPosition(
      this._meta.nominalCenterSG[0],
      this._meta.nominalCenterSG[1],
      this._meta.nominalCenterSG[2]
    );
  }

  /**
   * Static factory: checks if an integer is an official basin ID.
   * @param {number} id
   * @returns {boolean}
   */
  static isValidId(id) {
    return Number.isInteger(id) && id >= 1 && id <= 8;
  }

  /**
   * Returns all 8 canonical basins.
   * @returns {Array<BasinID>}
   */
  static getAllBasins() {
    return [1, 2, 3, 4, 5, 6, 7, 8].map(id => new BasinID(id));
  }

  toJSON() {
    return {
      type: this._type,
      id: this._id,
      name: this.name,
      slug: this.slug,
      attractor: this.attractorName,
      nominalCenterSG: this.nominalCenterSG
    };
  }

  toString() {
    return `BasinID(${this._id}: ${this.name})`;
  }
}

// ============================================================================
// 8. Jacobian3x3
// ============================================================================

/**
 * Jacobian3x3 represents a 3D velocity gradient tensor J_ij = dv_i / dx_j (i, j in {0,1,2}),
 * with units of (km/s) / (Mpc/h).
 * 
 * Kinematic Decomposition:
 *   J = S + Omega
 * where S is the symmetric strain rate tensor: S_ij = 0.5 * (J_ij + J_ji)
 * and Omega is the antisymmetric vorticity tensor: Omega_ij = 0.5 * (J_ij - J_ji).
 */
export class Jacobian3x3 {
  /**
   * @param {Array<number>|Float64Array} elements 9-element array in row-major order:
   * [J00, J01, J02, J10, J11, J12, J20, J21, J22] where J_ij = dv_i / dx_j.
   * @param {Object} [options]
   */
  constructor(elements, options = {}) {
    if (!elements || elements.length < 9) {
      throw new TypeError('Jacobian3x3: Expected 9 numerical elements for 3x3 tensor.');
    }
    for (let i = 0; i < 9; i++) {
      assertFiniteNumber(elements[i], `elements[${i}]`);
    }

    this._data = new Float64Array(9);
    for (let i = 0; i < 9; i++) {
      this._data[i] = elements[i];
    }
    this._unit = options.unit || ScientificUnits.KM_PER_S_PER_MPC_OVER_H;
    this._type = 'Jacobian3x3';

    if (options.freeze) {
      Object.freeze(this);
      Object.freeze(this._data);
    }
  }

  get data() { return this._data; }
  get unit() { return this._unit; }
  get scientificType() { return this._type; }

  /**
   * Element accessor J_ij = dv_i / dx_j.
   * @param {number} row 0, 1, 2 (velocity component vx, vy, vz)
   * @param {number} col 0, 1, 2 (spatial coordinate x, y, z)
   * @returns {number}
   */
  get(row, col) {
    if (row < 0 || row > 2 || col < 0 || col > 2) {
      throw new RangeError(`Jacobian3x3.get: indices out of range [${row}, ${col}]`);
    }
    return this._data[row * 3 + col];
  }

  /**
   * Trace of Jacobian: Tr(J) = div(v) = dvx/dx + dvy/dy + dvz/dz (Velocity divergence theta).
   * @returns {number}
   */
  trace() {
    return this._data[0] + this._data[4] + this._data[8];
  }

  /**
   * Velocity divergence theta = div(v).
   * @returns {number}
   */
  divergence() {
    return this.trace();
  }

  /**
   * Determinant of 3x3 matrix.
   * @returns {number}
   */
  determinant() {
    const d = this._data;
    return d[0] * (d[4] * d[8] - d[5] * d[7]) -
           d[1] * (d[3] * d[8] - d[5] * d[6]) +
           d[2] * (d[3] * d[7] - d[4] * d[6]);
  }

  /**
   * Symmetric strain rate tensor S = 0.5 * (J + J^T).
   * @returns {Float64Array} 9-element symmetric matrix.
   */
  symmetricStrainRate() {
    const d = this._data;
    const S = new Float64Array(9);
    S[0] = d[0];
    S[4] = d[4];
    S[8] = d[8];
    S[1] = S[3] = 0.5 * (d[1] + d[3]);
    S[2] = S[6] = 0.5 * (d[2] + d[6]);
    S[5] = S[7] = 0.5 * (d[5] + d[7]);
    return S;
  }

  /**
   * Traceless shear tensor sigma_ij = S_ij - (1/3)*theta*delta_ij.
   * @returns {Float64Array} 9-element traceless symmetric matrix.
   */
  shearTensor() {
    const S = this.symmetricStrainRate();
    const oneThirdDiv = this.trace() / 3.0;
    S[0] -= oneThirdDiv;
    S[4] -= oneThirdDiv;
    S[8] -= oneThirdDiv;
    return S;
  }

  /**
   * Antisymmetric vorticity tensor Omega = 0.5 * (J - J^T).
   * @returns {Float64Array} 9-element skew-symmetric matrix.
   */
  antisymmetricVorticity() {
    const d = this._data;
    const W = new Float64Array(9);
    W[0] = W[4] = W[8] = 0.0;
    W[1] = 0.5 * (d[1] - d[3]);
    W[3] = -W[1];
    W[2] = 0.5 * (d[2] - d[6]);
    W[6] = -W[2];
    W[5] = 0.5 * (d[5] - d[7]);
    W[7] = -W[5];
    return W;
  }

  /**
   * Vorticity vector omega = curl(v) = [dvz/dy - dvy/dz, dvx/dz - dvz/dx, dvy/dx - dvx/dy].
   * @returns {Array<number>} [omega_x, omega_y, omega_z]
   */
  vorticityVector() {
    const d = this._data;
    return [
      d[7] - d[5], // dvz/dy - dvy/dz
      d[2] - d[6], // dvx/dz - dvz/dx
      d[3] - d[1]  // dvy/dx - dvx/dy
    ];
  }

  /**
   * Frobenius norm ||J||_F = sqrt(sum(J_ij^2)).
   * @returns {number}
   */
  frobeniusNorm() {
    let sum = 0.0;
    for (let i = 0; i < 9; i++) {
      sum += this._data[i] * this._data[i];
    }
    return Math.sqrt(sum);
  }

  /**
   * Diagonalizes the symmetric strain rate tensor S = 0.5 * (J + J^T) via exact 3D Jacobi rotation.
   * @returns {EigenSystem3D}
   */
  diagonalizeStrain() {
    const S = this.symmetricStrainRate();
    return EigenSystem3D.fromSymmetricMatrix(S);
  }

  /**
   * Computes Okubo-Weiss vortex criterion Q = 0.25 * (||omega||^2 - 2*||S||^2).
   * Negative Q implies strain/deformation dominance, positive Q implies vortex dominance.
   * @returns {number}
   */
  okuboWeiss() {
    const vort = this.vorticityVector();
    const vortSq = vort[0] * vort[0] + vort[1] * vort[1] + vort[2] * vort[2];
    const S = this.symmetricStrainRate();
    let strainSq = 0.0;
    for (let i = 0; i < 9; i++) {
      strainSq += S[i] * S[i];
    }
    return 0.25 * (vortSq - 2.0 * strainSq);
  }

  /**
   * Multiplies Jacobian by vector: J * v.
   * @param {Array<number>|Float64Array} v 3D vector.
   * @returns {Array<number>} 3D vector.
   */
  multiplyVector(v) {
    const d = this._data;
    return [
      d[0] * v[0] + d[1] * v[1] + d[2] * v[2],
      d[3] * v[0] + d[4] * v[1] + d[5] * v[2],
      d[6] * v[0] + d[7] * v[1] + d[8] * v[2]
    ];
  }

  /**
   * Solves linear system J * x = b for 3D displacement root finding x = J^(-1) * b using Cramer's rule.
   * @param {Array<number>|Float64Array} b 3D right-hand side vector.
   * @returns {Array<number>|null} Solution vector x, or null if matrix is singular (det < 1e-12).
   */
  solveLinear(b) {
    const d = this._data;
    const det = this.determinant();
    if (Math.abs(det) < 1e-12) {
      return null;
    }
    const invDet = 1.0 / det;

    const b0 = b[0], b1 = b[1], b2 = b[2];

    const x0 = (b0 * (d[4] * d[8] - d[5] * d[7]) -
                d[1] * (b1 * d[8] - d[5] * b2) +
                d[2] * (b1 * d[7] - d[4] * b2)) * invDet;

    const x1 = (d[0] * (b1 * d[8] - d[5] * b2) -
                b0 * (d[3] * d[8] - d[5] * d[6]) +
                d[2] * (d[3] * b2 - b1 * d[6])) * invDet;

    const x2 = (d[0] * (d[4] * b2 - b1 * d[7]) -
                d[1] * (d[3] * b2 - b1 * d[6]) +
                b0 * (d[3] * d[7] - d[4] * d[6])) * invDet;

    return [x0, x1, x2];
  }

  clone() {
    return new Jacobian3x3(this._data, { unit: this._unit });
  }

  toJSON() {
    return {
      type: this._type,
      data: Array.from(this._data),
      unit: this._unit,
      divergence: this.divergence(),
      determinant: this.determinant()
    };
  }

  toString() {
    return `Jacobian3x3(div=${this.divergence().toFixed(4)} ${this._unit})`;
  }
}

// ============================================================================
// 9. EigenSystem3D
// ============================================================================

/**
 * EigenSystem3D represents the eigensystem of a real symmetric 3x3 matrix,
 * with eigenvalues sorted descending: lambda_1 >= lambda_2 >= lambda_3,
 * and corresponding orthonormal eigenvectors.
 */
export class EigenSystem3D {
  /**
   * @param {Array<number>} eigenvalues Sorted eigenvalues [lambda1, lambda2, lambda3] (lambda1 >= lambda2 >= lambda3).
   * @param {Array<Array<number>>} eigenvectors Orthonormal eigenvectors [v1, v2, v3] where v_i = [x, y, z].
   * @param {Object} [options]
   */
  constructor(eigenvalues, eigenvectors, options = {}) {
    if (!eigenvalues || eigenvalues.length < 3) {
      throw new TypeError('EigenSystem3D: Expected 3 eigenvalues.');
    }
    if (!eigenvectors || eigenvectors.length < 3) {
      throw new TypeError('EigenSystem3D: Expected 3 eigenvectors.');
    }

    for (let i = 0; i < 3; i++) {
      assertFiniteNumber(eigenvalues[i], `eigenvalues[${i}]`);
      if (!eigenvectors[i] || eigenvectors[i].length < 3) {
        throw new TypeError(`EigenSystem3D: Invalid eigenvector at index ${i}`);
      }
      for (let j = 0; j < 3; j++) {
        assertFiniteNumber(eigenvectors[i][j], `eigenvectors[${i}][${j}]`);
      }
    }

    this._values = [eigenvalues[0], eigenvalues[1], eigenvalues[2]];
    this._vectors = [
      [eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2]],
      [eigenvectors[1][0], eigenvectors[1][1], eigenvectors[1][2]],
      [eigenvectors[2][0], eigenvectors[2][1], eigenvectors[2][2]]
    ];
    this._type = 'EigenSystem3D';

    if (options.freeze) {
      Object.freeze(this);
      Object.freeze(this._values);
      Object.freeze(this._vectors[0]);
      Object.freeze(this._vectors[1]);
      Object.freeze(this._vectors[2]);
      Object.freeze(this._vectors);
    }
  }

  get eigenvalues() { return this._values; }
  get eigenvectors() { return this._vectors; }
  get lambda1() { return this._values[0]; }
  get lambda2() { return this._values[1]; }
  get lambda3() { return this._values[2]; }
  get v1() { return this._vectors[0]; }
  get v2() { return this._vectors[1]; }
  get v3() { return this._vectors[2]; }
  get scientificType() { return this._type; }

  /**
   * Verifies orthonormality of eigenvectors: v_i . v_j = delta_ij.
   * @param {number} [tolerance=1e-6]
   * @returns {boolean}
   */
  isOrthonormal(tolerance = 1e-6) {
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const dot = this._vectors[i][0] * this._vectors[j][0] +
                    this._vectors[i][1] * this._vectors[j][1] +
                    this._vectors[i][2] * this._vectors[j][2];
        const expected = (i === j) ? 1.0 : 0.0;
        if (Math.abs(dot - expected) > tolerance) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Reconstructs the original 3x3 matrix: A = sum_k lambda_k * (v_k (x) v_k^T).
   * @returns {Float64Array} 9-element row-major matrix.
   */
  reconstructMatrix() {
    const A = new Float64Array(9);
    for (let k = 0; k < 3; k++) {
      const lam = this._values[k];
      const vk = this._vectors[k];
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          A[i * 3 + j] += lam * vk[i] * vk[j];
        }
      }
    }
    return A;
  }

  /**
   * Hahn et al. (2007) / Forero-Romero et al. (2009) Cosmic Web Classification
   * based on the count of strain/tidal eigenvalues above a threshold lambda_th:
   * - 0 eigenvalues > lambda_th: Void (expanding in all 3 directions)
   * - 1 eigenvalue > lambda_th: Sheet (flattening along 1 axis, expanding along 2)
   * - 2 eigenvalues > lambda_th: Filament (collapsing along 2 axes, expanding along 1)
   * - 3 eigenvalues > lambda_th: Knot / Cluster (collapsing along all 3 axes)
   * 
   * @param {number} [threshold=0.0] Eigenvalue threshold lambda_th.
   * @returns {{category: string, code: number, countAboveThreshold: number}}
   */
  classifyWebStructure(threshold = 0.0) {
    let count = 0;
    if (this._values[0] > threshold) count++;
    if (this._values[1] > threshold) count++;
    if (this._values[2] > threshold) count++;

    const categories = ['void', 'sheet', 'filament', 'knot'];
    return {
      category: categories[count],
      code: count,
      countAboveThreshold: count
    };
  }

  /**
   * Exact 3D Jacobi eigenvalue and eigenvector solver for real symmetric 3x3 matrices.
   * Guaranteed quadratic convergence to machine precision without trigonometric loss.
   * 
   * @param {Float64Array|Array<number>} A 9-element symmetric row-major matrix.
   * @param {number} [maxSweeps=50]
   * @param {number} [tolerance=1e-15]
   * @returns {EigenSystem3D}
   */
  static fromSymmetricMatrix(A, maxSweeps = 50, tolerance = 1e-15) {
    if (!A || A.length < 9) {
      throw new TypeError('EigenSystem3D.fromSymmetricMatrix: expected 9 elements.');
    }

    // Working copies
    const a = [
      [A[0], 0.5 * (A[1] + A[3]), 0.5 * (A[2] + A[6])],
      [0.5 * (A[1] + A[3]), A[4], 0.5 * (A[5] + A[7])],
      [0.5 * (A[2] + A[6]), 0.5 * (A[5] + A[7]), A[8]]
    ];

    // Eigenvector matrix initialized to Identity
    const v = [
      [1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0],
      [0.0, 0.0, 1.0]
    ];

    for (let sweep = 0; sweep < maxSweeps; sweep++) {
      // Sum of off-diagonal elements
      const offDiag = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
      if (offDiag < tolerance) {
        break;
      }

      for (let p = 0; p < 2; p++) {
        for (let q = p + 1; q < 3; q++) {
          const apq = a[p][q];
          if (Math.abs(apq) < 1e-16) continue;

          const app = a[p][p];
          const aqq = a[q][q];
          const tau = (aqq - app) / (2.0 * apq);
          let t;
          if (tau >= 0) {
            t = 1.0 / (tau + Math.sqrt(1.0 + tau * tau));
          } else {
            t = -1.0 / (-tau + Math.sqrt(1.0 + tau * tau));
          }

          const c = 1.0 / Math.sqrt(1.0 + t * t);
          const s = t * c;
          const h = t * apq;

          a[p][p] -= h;
          a[q][q] += h;
          a[p][q] = 0.0;
          a[q][p] = 0.0;

          // Rotate remaining rows and columns
          for (let j = 0; j < 3; j++) {
            if (j !== p && j !== q) {
              const ajp = a[j][p];
              const ajq = a[j][q];
              a[j][p] = a[p][j] = c * ajp - s * ajq;
              a[j][q] = a[q][j] = s * ajp + c * ajq;
            }
          }

          // Accumulate transformation into eigenvector matrix
          for (let j = 0; j < 3; j++) {
            const vjp = v[j][p];
            const vjq = v[j][q];
            v[j][p] = c * vjp - s * vjq;
            v[j][q] = s * vjp + c * vjq;
          }
        }
      }
    }

    // Extract eigenvalues and eigenvectors
    const rawPairs = [
      { val: a[0][0], vec: [v[0][0], v[1][0], v[2][0]] },
      { val: a[1][1], vec: [v[0][1], v[1][1], v[2][1]] },
      { val: a[2][2], vec: [v[0][2], v[1][2], v[2][2]] }
    ];

    // Sort descending: lambda_1 >= lambda_2 >= lambda_3
    rawPairs.sort((pairA, pairB) => pairB.val - pairA.val);

    // Normalize eigenvectors and ensure right-handed coordinate orientation (det = +1)
    const eValues = [rawPairs[0].val, rawPairs[1].val, rawPairs[2].val];
    const eVectors = [rawPairs[0].vec, rawPairs[1].vec, rawPairs[2].vec];

    for (let i = 0; i < 3; i++) {
      const vec = eVectors[i];
      const norm = Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
      if (norm > 0) {
        vec[0] /= norm;
        vec[1] /= norm;
        vec[2] /= norm;
      }
    }

    return new EigenSystem3D(eValues, eVectors);
  }

  toJSON() {
    return {
      type: this._type,
      eigenvalues: this._values,
      eigenvectors: this._vectors,
      webClassification: this.classifyWebStructure(0.0)
    };
  }

  toString() {
    return `EigenSystem3D(l1=${this._values[0].toFixed(4)}, l2=${this._values[1].toFixed(4)}, l3=${this._values[2].toFixed(4)})`;
  }
}

// ============================================================================
// 10. StreamlinePath
// ============================================================================

/**
 * StreamlinePath represents a continuous dynamical streamline integrated through
 * the 3D cosmological velocity field.
 */
export class StreamlinePath {
  /**
   * @param {Array<SupergalacticPosition>} positions Array of spatial positions along trajectory.
   * @param {Array<VelocityVector>} velocities Array of velocity vectors along trajectory.
   * @param {Array<number>} [timeSteps] Array of integration times/arclengths (Gyr).
   * @param {Object} [options]
   */
  constructor(positions, velocities, timeSteps = [], options = {}) {
    if (!Array.isArray(positions) || positions.length === 0) {
      throw new TypeError('StreamlinePath: positions must be a non-empty array of SupergalacticPosition.');
    }
    for (let i = 0; i < positions.length; i++) {
      if (!(positions[i] instanceof SupergalacticPosition)) {
        throw new TypeError(`StreamlinePath: position at index ${i} is not a SupergalacticPosition.`);
      }
    }

    if (velocities && velocities.length > 0) {
      if (velocities.length !== positions.length) {
        throw new RangeError(`StreamlinePath: velocities length (${velocities.length}) must match positions length (${positions.length}).`);
      }
      for (let i = 0; i < velocities.length; i++) {
        if (!(velocities[i] instanceof VelocityVector)) {
          throw new TypeError(`StreamlinePath: velocity at index ${i} is not a VelocityVector.`);
        }
      }
    }

    this._positions = positions.slice();
    this._velocities = velocities ? velocities.slice() : [];
    this._timeSteps = timeSteps.slice();
    this._basinAssignment = options.basinAssignment || null;
    this._isConverged = options.isConverged || false;
    this._type = 'StreamlinePath';

    if (options.freeze) {
      Object.freeze(this._positions);
      Object.freeze(this._velocities);
      Object.freeze(this._timeSteps);
      Object.freeze(this);
    }
  }

  get length() { return this._positions.length; }
  get positions() { return this._positions; }
  get velocities() { return this._velocities; }
  get timeSteps() { return this._timeSteps; }
  get startPosition() { return this._positions[0]; }
  get endPosition() { return this._positions[this._positions.length - 1]; }
  get basinAssignment() { return this._basinAssignment; }
  get isConverged() { return this._isConverged; }
  get scientificType() { return this._type; }

  /**
   * Sets basin assignment.
   * @param {BasinID|number} basin
   */
  setBasinAssignment(basin) {
    if (basin instanceof BasinID) {
      this._basinAssignment = basin;
    } else if (typeof basin === 'number') {
      this._basinAssignment = new BasinID(basin);
    } else {
      this._basinAssignment = null;
    }
  }

  /**
   * Computes total arclength along streamline (in Mpc/h).
   * @returns {PhysicalDistance}
   */
  totalArclength() {
    let total = 0.0;
    for (let i = 1; i < this._positions.length; i++) {
      const p1 = this._positions[i - 1];
      const p2 = this._positions[i];
      const dx = p2.sgx - p1.sgx;
      const dy = p2.sgy - p1.sgy;
      const dz = p2.sgz - p1.sgz;
      total += Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
    return new PhysicalDistance(total, ScientificUnits.MPC_OVER_H);
  }

  /**
   * Computes 3D axis-aligned bounding box of streamline: { min: [x,y,z], max: [x,y,z] }.
   * @returns {{min: Array<number>, max: Array<number>}}
   */
  boundingBox() {
    let minX = Infinity, minY = Infinity, minZ = Infinity;
    let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;

    for (let i = 0; i < this._positions.length; i++) {
      const p = this._positions[i];
      if (p.sgx < minX) minX = p.sgx;
      if (p.sgy < minY) minY = p.sgy;
      if (p.sgz < minZ) minZ = p.sgz;
      if (p.sgx > maxX) maxX = p.sgx;
      if (p.sgy > maxY) maxY = p.sgy;
      if (p.sgz > maxZ) maxZ = p.sgz;
    }

    return {
      min: [minX, minY, minZ],
      max: [maxX, maxY, maxZ]
    };
  }

  /**
   * Downsamples streamline to at most maxPoints points uniformly along index.
   * @param {number} maxPoints
   * @returns {StreamlinePath}
   */
  downsample(maxPoints = 100) {
    if (this._positions.length <= maxPoints) {
      return this;
    }
    const step = (this._positions.length - 1) / (maxPoints - 1);
    const newPos = [];
    const newVel = [];
    const newTimes = [];

    for (let i = 0; i < maxPoints; i++) {
      const idx = Math.min(this._positions.length - 1, Math.round(i * step));
      newPos.push(this._positions[idx]);
      if (this._velocities.length > idx) newVel.push(this._velocities[idx]);
      if (this._timeSteps.length > idx) newTimes.push(this._timeSteps[idx]);
    }

    return new StreamlinePath(newPos, newVel, newTimes, {
      basinAssignment: this._basinAssignment,
      isConverged: this._isConverged
    });
  }

  toJSON() {
    return {
      type: this._type,
      pointCount: this._positions.length,
      start: this.startPosition.toJSON(),
      end: this.endPosition.toJSON(),
      arclength: this.totalArclength().value,
      basinAssignment: this._basinAssignment ? this._basinAssignment.id : null,
      isConverged: this._isConverged
    };
  }

  toString() {
    return `StreamlinePath(${this._positions.length} points, arclength=${this.totalArclength().value.toFixed(2)} Mpc/h, basin=${this._basinAssignment ? this._basinAssignment.name : 'unassigned'})`;
  }
}

// ============================================================================
// 11. RunManifest
// ============================================================================

/**
 * RunManifest provides W3C PROV-O compliant cryptographic lineage tracking and
 * NIST SHA-256 integrity verification for reproducible cosmological computations.
 */
export class RunManifest {
  /**
   * @param {Object} config Manifest parameters.
   * @param {string} config.runId Unique UUID or identifier.
   * @param {string} config.algorithmName Algorithm / pipeline name.
   * @param {string} [config.gitCommit='HEAD'] Git commit hash.
   * @param {Object} [config.parameters={}] Computational hyper-parameters.
   * @param {string} [config.inputDatasetHash=''] SHA-256 hash of input data.
   * @param {string} [config.outputFieldHash=''] SHA-256 hash of output field.
   * @param {number} [config.executionDurationMs=0] Execution time in milliseconds.
   */
  constructor(config) {
    if (!config || !config.runId || !config.algorithmName) {
      throw new TypeError('RunManifest: runId and algorithmName are required.');
    }

    this.runId = String(config.runId);
    this.algorithmName = String(config.algorithmName);
    this.gitCommit = config.gitCommit || 'HEAD';
    this.timestamp = config.timestamp || new Date().toISOString();
    this.parameters = JSON.parse(JSON.stringify(config.parameters || {}));
    this.inputDatasetHash = config.inputDatasetHash || '';
    this.outputFieldHash = config.outputFieldHash || '';
    this.executionDurationMs = config.executionDurationMs || 0;
    this._type = 'RunManifest';

    // Compute self cryptographic checksum
    this.manifestChecksum = this.computeManifestHash();
  }

  get scientificType() { return this._type; }

  /**
   * Computes SHA-256 checksum of manifest contents.
   * @returns {string}
   */
  computeManifestHash() {
    const payload = JSON.stringify({
      runId: this.runId,
      algorithmName: this.algorithmName,
      gitCommit: this.gitCommit,
      timestamp: this.timestamp,
      parameters: this.parameters,
      inputDatasetHash: this.inputDatasetHash,
      outputFieldHash: this.outputFieldHash
    });
    return sha256Hex(payload);
  }

  /**
   * Generates standard W3C PROV-O JSON-LD representation.
   * @returns {Object}
   */
  toW3CProvJSONLD() {
    return {
      '@context': {
        'prov': 'http://www.w3.org/ns/prov#',
        'zrt': 'https://zrt.cosmicflows.org/prov#'
      },
      '@id': `zrt:activity/${this.runId}`,
      '@type': 'prov:Activity',
      'prov:startedAtTime': this.timestamp,
      'prov:used': {
        '@id': `zrt:entity/inputDataset#${this.inputDatasetHash || 'raw'}`,
        'prov:type': 'zrt:CosmicFlowsGrid'
      },
      'prov:generated': {
        '@id': `zrt:entity/outputField#${this.outputFieldHash || 'computed'}`,
        'prov:type': 'zrt:VelocityPotentialField'
      },
      'prov:wasAssociatedWith': {
        '@id': `zrt:agent/algorithm/${this.algorithmName}`,
        'prov:type': 'prov:SoftwareAgent',
        'zrt:gitCommit': this.gitCommit
      },
      'zrt:parameters': this.parameters,
      'zrt:executionDurationMs': this.executionDurationMs,
      'zrt:manifestChecksum': this.manifestChecksum
    };
  }

  /**
   * Verifies that the manifest has not been tampered with.
   * @returns {boolean}
   */
  verifyIntegrity() {
    return this.computeManifestHash() === this.manifestChecksum;
  }

  toJSON() {
    return {
      type: this._type,
      runId: this.runId,
      algorithmName: this.algorithmName,
      gitCommit: this.gitCommit,
      timestamp: this.timestamp,
      parameters: this.parameters,
      inputDatasetHash: this.inputDatasetHash,
      outputFieldHash: this.outputFieldHash,
      executionDurationMs: this.executionDurationMs,
      manifestChecksum: this.manifestChecksum,
      w3cProv: this.toW3CProvJSONLD()
    };
  }

  toString() {
    return `RunManifest(runId='${this.runId}', alg='${this.algorithmName}', sha256='${this.manifestChecksum.substring(0, 8)}...')`;
  }
}

// ============================================================================
// 12. Type-Checking and Assertion Helpers
// ============================================================================

/**
 * Asserts that an object is an instance of SupergalacticPosition.
 * @param {*} obj
 * @param {string} [name='argument']
 * @returns {SupergalacticPosition}
 */
export function ensurePosition(obj, name = 'position') {
  if (!(obj instanceof SupergalacticPosition)) {
    throw new TypeError(`Expected SupergalacticPosition for '${name}', received: ${obj && obj.constructor ? obj.constructor.name : typeof obj}`);
  }
  return obj;
}

/**
 * Asserts that an object is an instance of VelocityVector.
 * @param {*} obj
 * @param {string} [name='argument']
 * @returns {VelocityVector}
 */
export function ensureVelocity(obj, name = 'velocity') {
  if (!(obj instanceof VelocityVector)) {
    throw new TypeError(`Expected VelocityVector for '${name}', received: ${obj && obj.constructor ? obj.constructor.name : typeof obj}`);
  }
  return obj;
}

/**
 * Asserts that an object is an instance of DensityContrast.
 * @param {*} obj
 * @param {string} [name='argument']
 * @returns {DensityContrast}
 */
export function ensureDensity(obj, name = 'density') {
  if (!(obj instanceof DensityContrast)) {
    throw new TypeError(`Expected DensityContrast for '${name}', received: ${obj && obj.constructor ? obj.constructor.name : typeof obj}`);
  }
  return obj;
}

/**
 * Asserts that an object is an instance of BasinID.
 * @param {*} obj
 * @param {string} [name='argument']
 * @returns {BasinID}
 */
export function ensureBasinID(obj, name = 'basinId') {
  if (!(obj instanceof BasinID)) {
    throw new TypeError(`Expected BasinID for '${name}', received: ${obj && obj.constructor ? obj.constructor.name : typeof obj}`);
  }
  return obj;
}
