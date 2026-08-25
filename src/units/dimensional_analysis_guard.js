/**
 * @file dimensional_analysis_guard.js
 * @description Strict 7-Base Physical Dimensional Analysis Guard and Dynamic Unit Converter
 * for Cosmological Field Processing, CF4 Kinematics, and Large-Scale Structure Analysis.
 * 
 * Enforces Scientific Invariants (AGENTS.md Pillar 1 & Rule 1-4):
 * 1. Strict 7-Base Physical Dimension Tensor: [L, M, T, I, Theta, N, J] + h-exponent + a-exponent.
 * 2. Static & Runtime Dimensional Consistency Checker:
 *    - Spatial Position [L=1, h^-1] (Mpc/h)
 *    - Peculiar Velocity [L=1, T=-1] (km/s)
 *    - Density Contrast delta [dimensionless]
 *    - Hubble Parameter H(z) [T=-1] (km/s/Mpc)
 * 3. Explicit DimensionalityMismatchError (extends TypeError) thrown upon unphysical operations
 *    (e.g., position + velocity, density + potential).
 * 4. Exact high-precision unit conversions across cosmic distances (Mpc/h, kpc, pc, km, m, AU, ly),
 *    velocities (km/s, m/s, c, redshift z, relativistic rapidity), masses (M_sun/h, kg),
 *    and cosmological density/potential/shear tensors.
 * 5. Symbolic AST dimensional expression parser and compile-time contract guards.
 * 
 * @module units/dimensional_analysis_guard
 */

// ============================================================================
// 1. FUNDAMENTAL PHYSICAL AND COSMOLOGICAL CONSTANTS (CODATA 2022 / IAU 2015)
// ============================================================================

export const CONSTANTS = Object.freeze({
  // Speed of light in vacuum (exact SI)
  SPEED_OF_LIGHT_SI: 299792458.0, // m / s
  SPEED_OF_LIGHT_KMS: 299792.458, // km / s
  SPEED_OF_LIGHT_CGS: 2.99792458e10, // cm / s

  // Newtonian Gravitational Constant (CODATA 2022)
  GRAVITATIONAL_CONSTANT_SI: 6.67430e-11, // m^3 / (kg * s^2)
  GRAVITATIONAL_CONSTANT_CGS: 6.67430e-8, // cm^3 / (g * s^2)
  // G in (km/s)^2 * Mpc / M_sun
  GRAVITATIONAL_CONSTANT_ASTRO: 4.30091e-9, // (km/s)^2 * Mpc / M_sun
  // G in (km/s)^2 * (Mpc/h) / (10^10 M_sun/h)
  GRAVITATIONAL_CONSTANT_GADGET: 4.30091e-4,

  // Solar Mass (IAU 2015 Nominal Solar Mass)
  SOLAR_MASS_KG: 1.98847e30, // kg
  SOLAR_MASS_G: 1.98847e33, // g
  EARTH_MASS_KG: 5.9722e24, // kg
  JUPITER_MASS_KG: 1.89813e27, // kg

  // Solar Luminosity
  SOLAR_LUMINOSITY_WATTS: 3.828e26, // W
  SOLAR_LUMINOSITY_ERG_S: 3.828e33, // erg / s

  // Astrophysical Length Standards (IAU 2012 / IAU 2015)
  METERS_PER_AU: 149597870700.0, // m (exact IAU 2012)
  METERS_PER_PC: 3.085677581491367e16, // m (exact IAU 2015: 648000 / pi * AU)
  METERS_PER_KPC: 3.085677581491367e19, // m
  METERS_PER_MPC: 3.085677581491367e22, // m
  METERS_PER_GPC: 3.085677581491367e25, // m
  METERS_PER_LY: 9.4607304725808e15, // m (Julian light year: 365.25 * 86400 * c)

  KM_PER_MPC: 3.085677581491367e19, // km
  KM_PER_KPC: 3.085677581491367e16, // km
  KM_PER_PC: 3.085677581491367e13, // km

  // Time Standards (Julian Calendar)
  SECONDS_PER_MINUTE: 60.0,
  SECONDS_PER_HOUR: 3600.0,
  SECONDS_PER_DAY: 86400.0,
  SECONDS_PER_JULIAN_YEAR: 31557600.0, // s (365.25 days)
  SECONDS_PER_KYR: 3.15576e10, // s
  SECONDS_PER_MYR: 3.15576e13, // s
  SECONDS_PER_GYR: 3.15576e16, // s

  // Thermodynamics & Microphysics
  BOLTZMANN_CONSTANT_SI: 1.380649e-23, // J / K (exact)
  PLANCK_CONSTANT_SI: 6.62607015e-34, // J * s (exact)
  ELEMENTARY_CHARGE_SI: 1.602176634e-19, // C (exact)
  ELECTRON_VOLT_JOULES: 1.602176634e-19, // J

  // Standard Hubble Benchmark (H100 = 100 km/s/Mpc)
  H100_KMS_MPC: 100.0,
  // 1 km/s / Mpc in SI units (s^-1) = 1.0e3 / METERS_PER_MPC
  KMS_MPC_TO_SEC_INV: 3.240779289444365e-18, // s^-1

  // Cosmological Critical Density at z=0: rho_c0 / h^2 = 3 * (100 km/s/Mpc)^2 / (8 * pi * G)
  // in M_sun / Mpc^3
  CRITICAL_DENSITY_COEFF_MSUN_MPC3: 2.77536627245708e11, // (M_sun / Mpc^3) / h^2
  // in kg / m^3
  CRITICAL_DENSITY_COEFF_KG_M3: 1.8784711e-26, // (kg / m^3) / h^2

  // CosmicFlows-4 Official Scaling Factor (Rule 3)
  CF4_VELOCITY_SCALE_FACTOR: 52.0 // km/s per raw unit
});

/**
 * Standard cosmological parameter benchmarks.
 */
export const COSMOLOGY_BENCHMARKS = Object.freeze({
  CF4: Object.freeze({
    name: 'CosmicFlows-4 Standard',
    H0: 74.6, // km/s/Mpc
    h: 0.746,
    Omega_m: 0.31,
    Omega_Lambda: 0.69,
    Omega_r: 8.48e-5,
    Omega_k: 0.0,
    growthRate_f0: 0.52573, // (0.31)^0.55
    continuityScale_H0f: 52.57 // km/s / (Mpc/h)
  }),
  PLANCK18: Object.freeze({
    name: 'Planck 2018 Final',
    H0: 67.4,
    h: 0.674,
    Omega_m: 0.315,
    Omega_Lambda: 0.685,
    Omega_r: 9.16e-5,
    Omega_k: 0.0,
    growthRate_f0: 0.5307,
    continuityScale_H0f: 53.07
  }),
  WMAP9: Object.freeze({
    name: 'WMAP 9-Year',
    H0: 70.0,
    h: 0.700,
    Omega_m: 0.279,
    Omega_Lambda: 0.721,
    Omega_r: 8.5e-5,
    Omega_k: 0.0,
    growthRate_f0: 0.4950,
    continuityScale_H0f: 49.50
  })
});

// ============================================================================
// 2. EXPLICIT ERROR TYPES
// ============================================================================

/**
 * Thrown when an illegal dimensional arithmetic or conversion operation is attempted.
 * Subclasses TypeError to satisfy both TypeError catch blocks and named inspection.
 */
export class DimensionalityMismatchError extends TypeError {
  /**
   * @param {string} message - Descriptive error explanation.
   * @param {object} [details] - Diagnostic payload.
   */
  constructor(message, details = {}) {
    const fullMessage = message.startsWith('DimensionalityMismatchError:')
      ? message
      : `DimensionalityMismatchError: ${message}`;
    super(fullMessage);
    this.name = 'DimensionalityMismatchError';
    this.sourceDimension = details.sourceDimension ?? null;
    this.targetDimension = details.targetDimension ?? null;
    this.sourceUnit = details.sourceUnit ?? null;
    this.targetUnit = details.targetUnit ?? null;
    this.operation = details.operation ?? null;
    this.context = details.context ?? null;
  }
}

/**
 * Thrown when a fundamental cosmological or dataset invariant is violated.
 */
export class InvariantViolationError extends Error {
  constructor(message, details = {}) {
    super(`[SCIENTIFIC INVARIANT VIOLATION] ${message}`);
    this.name = 'InvariantViolationError';
    this.details = details;
  }
}

/**
 * Thrown when a unit expression cannot be parsed or resolved.
 */
export class UnitParseError extends Error {
  constructor(message, expression = '') {
    super(`UnitParseError: ${message}${expression ? ` for expression "${expression}"` : ''}`);
    this.name = 'UnitParseError';
    this.expression = expression;
  }
}

// ============================================================================
// 3. STRICT 7-BASE PHYSICAL DIMENSION TENSOR CLASS
// ============================================================================

/**
 * Base Physical Dimension Indices in the 7-SI Base Tensor:
 * Index 0: L  - Length [meter]
 * Index 1: M  - Mass [kilogram]
 * Index 2: T  - Time [second]
 * Index 3: I  - Electric Current [ampere]
 * Index 4: Θ  - Thermodynamic Temperature [kelvin]
 * Index 5: N  - Amount of Substance [mole]
 * Index 6: J  - Luminous Intensity [candela]
 * 
 * Plus cosmological scale extensions:
 * - hExp: Exponent of dimensionless Hubble parameter h (where H0 = 100 h km/s/Mpc)
 * - aExp: Exponent of cosmological scale factor a = 1/(1+z) (comoving vs physical)
 */
export class DimensionTensor {
  /**
   * Constructs an immutable 7-base physical dimension tensor.
   * 
   * @param {number} [L=0] - Length exponent
   * @param {number} [M=0] - Mass exponent
   * @param {number} [T=0] - Time exponent
   * @param {number} [I=0] - Electric Current exponent
   * @param {number} [Theta=0] - Thermodynamic Temperature exponent
   * @param {number} [N=0] - Amount of Substance exponent
   * @param {number} [J=0] - Luminous Intensity exponent
   * @param {number} [hExp=0] - Cosmological h exponent
   * @param {number} [aExp=0] - Cosmological scale factor a exponent
   */
  constructor(L = 0, M = 0, T = 0, I = 0, Theta = 0, N = 0, J = 0, hExp = 0, aExp = 0) {
    this.L = Number(L) || 0;
    this.M = Number(M) || 0;
    this.T = Number(T) || 0;
    this.I = Number(I) || 0;
    this.Theta = Number(Theta) || 0;
    this.N = Number(N) || 0;
    this.J = Number(J) || 0;
    this.hExp = Number(hExp) || 0;
    this.aExp = Number(aExp) || 0;
    Object.freeze(this);
  }

  /**
   * Returns array of 7 SI base dimension exponents: [L, M, T, I, Theta, N, J].
   * @returns {number[]}
   */
  toBaseArray() {
    return [this.L, this.M, this.T, this.I, this.Theta, this.N, this.J];
  }

  /**
   * Returns array of all 9 exponents: [L, M, T, I, Theta, N, J, hExp, aExp].
   * @returns {number[]}
   */
  toArray() {
    return [this.L, this.M, this.T, this.I, this.Theta, this.N, this.J, this.hExp, this.aExp];
  }

  /**
   * Dimensionless check (all 7 SI base dimensions must be zero).
   * @returns {boolean}
   */
  isDimensionless() {
    return this.L === 0 && this.M === 0 && this.T === 0 &&
           this.I === 0 && this.Theta === 0 && this.N === 0 && this.J === 0;
  }

  /**
   * Strictly dimensionless check (including h and a exponents).
   * @returns {boolean}
   */
  isStrictlyDimensionless() {
    return this.isDimensionless() && this.hExp === 0 && this.aExp === 0;
  }

  /**
   * Checks exact physical equivalence with another dimension tensor.
   * Compares 7 SI base dimensions and optionally h/a exponents.
   * 
   * @param {DimensionTensor} other
   * @param {boolean} [strictCosmo=true] - If true, hExp and aExp must also match.
   * @returns {boolean}
   */
  equals(other, strictCosmo = true) {
    if (!(other instanceof DimensionTensor)) return false;
    const baseEqual = (
      this.L === other.L &&
      this.M === other.M &&
      this.T === other.T &&
      this.I === other.I &&
      this.Theta === other.Theta &&
      this.N === other.N &&
      this.J === other.J
    );
    if (!baseEqual) return false;
    if (strictCosmo) {
      return this.hExp === other.hExp && this.aExp === other.aExp;
    }
    return true;
  }

  /**
   * Multiplies two dimension tensors (adds exponents).
   * @param {DimensionTensor} other
   * @returns {DimensionTensor}
   */
  multiply(other) {
    if (!(other instanceof DimensionTensor)) {
      throw new DimensionalityMismatchError(`Expected DimensionTensor in multiply, received ${typeof other}`);
    }
    return new DimensionTensor(
      this.L + other.L,
      this.M + other.M,
      this.T + other.T,
      this.I + other.I,
      this.Theta + other.Theta,
      this.N + other.N,
      this.J + other.J,
      this.hExp + other.hExp,
      this.aExp + other.aExp
    );
  }

  /**
   * Divides this dimension tensor by another (subtracts exponents).
   * @param {DimensionTensor} other
   * @returns {DimensionTensor}
   */
  divide(other) {
    if (!(other instanceof DimensionTensor)) {
      throw new DimensionalityMismatchError(`Expected DimensionTensor in divide, received ${typeof other}`);
    }
    return new DimensionTensor(
      this.L - other.L,
      this.M - other.M,
      this.T - other.T,
      this.I - other.I,
      this.Theta - other.Theta,
      this.N - other.N,
      this.J - other.J,
      this.hExp - other.hExp,
      this.aExp - other.aExp
    );
  }

  /**
   * Raises this dimension tensor to a real power.
   * @param {number} power
   * @returns {DimensionTensor}
   */
  pow(power) {
    const p = Number(power);
    if (isNaN(p)) throw new TypeError(`Power must be a valid number, received ${power}`);
    return new DimensionTensor(
      this.L * p,
      this.M * p,
      this.T * p,
      this.I * p,
      this.Theta * p,
      this.N * p,
      this.J * p,
      this.hExp * p,
      this.aExp * p
    );
  }

  /**
   * Computes the nth root of this dimension tensor.
   * @param {number} n
   * @returns {DimensionTensor}
   */
  root(n) {
    if (n === 0) throw new RangeError('Cannot take 0th root of DimensionTensor');
    return this.pow(1.0 / n);
  }

  /**
   * Inverts this dimension tensor (negates all exponents).
   * @returns {DimensionTensor}
   */
  invert() {
    return this.pow(-1.0);
  }

  /**
   * Formats tensor as human-readable scientific dimension string (e.g. "[L^1 T^-1 h^-1]").
   * @returns {string}
   */
  toString() {
    if (this.isStrictlyDimensionless()) return '[Dimensionless]';
    const parts = [];
    if (this.L !== 0) parts.push(this.L === 1 ? 'L' : `L^${this.L}`);
    if (this.M !== 0) parts.push(this.M === 1 ? 'M' : `M^${this.M}`);
    if (this.T !== 0) parts.push(this.T === 1 ? 'T' : `T^${this.T}`);
    if (this.I !== 0) parts.push(this.I === 1 ? 'I' : `I^${this.I}`);
    if (this.Theta !== 0) parts.push(this.Theta === 1 ? 'Θ' : `Θ^${this.Theta}`);
    if (this.N !== 0) parts.push(this.N === 1 ? 'N' : `N^${this.N}`);
    if (this.J !== 0) parts.push(this.J === 1 ? 'J' : `J^${this.J}`);
    if (this.hExp !== 0) parts.push(this.hExp === 1 ? 'h' : `h^${this.hExp}`);
    if (this.aExp !== 0) parts.push(this.aExp === 1 ? 'a' : `a^${this.aExp}`);
    return `[${parts.join(' ')}]`;
  }

  /**
   * Formats tensor as LaTeX expression (e.g. "\\mathrm{L}\\mathrm{T}^{-1}h^{-1}").
   * @returns {string}
   */
  toLaTeX() {
    if (this.isStrictlyDimensionless()) return '1';
    const parts = [];
    const syms = [
      ['\\mathrm{L}', this.L],
      ['\\mathrm{M}', this.M],
      ['\\mathrm{T}', this.T],
      ['\\mathrm{I}', this.I],
      ['\\Theta', this.Theta],
      ['\\mathrm{N}', this.N],
      ['\\mathrm{J}', this.J],
      ['h', this.hExp],
      ['a', this.aExp]
    ];
    for (const [sym, exp] of syms) {
      if (exp === 1) parts.push(sym);
      else if (exp !== 0) parts.push(`${sym}^{${exp}}`);
    }
    return parts.join('');
  }

  /**
   * Creates a DimensionTensor from a raw object or vector.
   * @param {Object|Array<number>} input
   * @returns {DimensionTensor}
   */
  static from(input) {
    if (input instanceof DimensionTensor) return input;
    if (Array.isArray(input)) {
      return new DimensionTensor(
        input[0] || 0,
        input[1] || 0,
        input[2] || 0,
        input[3] || 0,
        input[4] || 0,
        input[5] || 0,
        input[6] || 0,
        input[7] || 0,
        input[8] || 0
      );
    }
    if (input && typeof input === 'object') {
      return new DimensionTensor(
        input.L ?? input.l ?? 0,
        input.M ?? input.m ?? 0,
        input.T ?? input.t ?? 0,
        input.I ?? input.i ?? 0,
        input.Theta ?? input.theta ?? 0,
        input.N ?? input.n ?? 0,
        input.J ?? input.j ?? 0,
        input.hExp ?? input.h ?? 0,
        input.aExp ?? input.a ?? 0
      );
    }
    throw new TypeError(`Cannot create DimensionTensor from ${typeof input}`);
  }
}

// ============================================================================
// 4. STANDARD CANONICAL PHYSICAL & COSMOLOGICAL DIMENSIONS
// ============================================================================

export const PHYSICAL_DIMENSIONS = Object.freeze({
  // Dimensionless
  DIMENSIONLESS: Object.freeze(new DimensionTensor(0, 0, 0, 0, 0, 0, 0, 0, 0)),

  // Spatial & Distance
  LENGTH: Object.freeze(new DimensionTensor(1, 0, 0, 0, 0, 0, 0, 0, 0)),
  LENGTH_H_INV: Object.freeze(new DimensionTensor(1, 0, 0, 0, 0, 0, 0, -1, 0)), // Mpc/h, kpc/h, pc/h
  AREA: Object.freeze(new DimensionTensor(2, 0, 0, 0, 0, 0, 0, 0, 0)),
  AREA_H_INV2: Object.freeze(new DimensionTensor(2, 0, 0, 0, 0, 0, 0, -2, 0)),
  VOLUME: Object.freeze(new DimensionTensor(3, 0, 0, 0, 0, 0, 0, 0, 0)),
  VOLUME_H_INV3: Object.freeze(new DimensionTensor(3, 0, 0, 0, 0, 0, 0, -3, 0)), // (Mpc/h)^3
  WAVENUMBER: Object.freeze(new DimensionTensor(-1, 0, 0, 0, 0, 0, 0, 0, 0)), // 1/Mpc, 1/m
  WAVENUMBER_H: Object.freeze(new DimensionTensor(-1, 0, 0, 0, 0, 0, 0, 1, 0)), // h/Mpc

  // Kinematic & Dynamic
  TIME: Object.freeze(new DimensionTensor(0, 0, 1, 0, 0, 0, 0, 0, 0)),
  FREQUENCY: Object.freeze(new DimensionTensor(0, 0, -1, 0, 0, 0, 0, 0, 0)), // s^-1
  EXPANSION_RATE: Object.freeze(new DimensionTensor(0, 0, -1, 0, 0, 0, 0, 0, 0)), // km/s/Mpc, s^-1
  EXPANSION_RATE_H: Object.freeze(new DimensionTensor(0, 0, -1, 0, 0, 0, 0, 1, 0)), // km/s/(Mpc/h)
  VELOCITY: Object.freeze(new DimensionTensor(1, 0, -1, 0, 0, 0, 0, 0, 0)), // km/s, m/s, c
  ACCELERATION: Object.freeze(new DimensionTensor(1, 0, -2, 0, 0, 0, 0, 0, 0)), // m/s^2, cm/s^2
  ACCELERATION_ASTRO: Object.freeze(new DimensionTensor(1, 0, -2, 0, 0, 0, 0, 1, 0)), // (km/s)^2 / (Mpc/h)
  SPECIFIC_ENERGY: Object.freeze(new DimensionTensor(2, 0, -2, 0, 0, 0, 0, 0, 0)), // (km/s)^2, m^2/s^2, J/kg
  GRAVITATIONAL_POTENTIAL: Object.freeze(new DimensionTensor(2, 0, -2, 0, 0, 0, 0, 0, 0)),

  // Mass & Matter
  MASS: Object.freeze(new DimensionTensor(0, 1, 0, 0, 0, 0, 0, 0, 0)), // kg, M_sun
  MASS_H_INV: Object.freeze(new DimensionTensor(0, 1, 0, 0, 0, 0, 0, -1, 0)), // M_sun/h, 1e10 M_sun/h
  MASS_DENSITY: Object.freeze(new DimensionTensor(-3, 1, 0, 0, 0, 0, 0, 0, 0)), // kg/m^3, g/cm^3, M_sun/Mpc^3
  MASS_DENSITY_H2: Object.freeze(new DimensionTensor(-3, 1, 0, 0, 0, 0, 0, 2, 0)), // (M_sun/h)/(Mpc/h)^3, rho_crit,0
  MASS_DENSITY_H3: Object.freeze(new DimensionTensor(-3, 1, 0, 0, 0, 0, 0, 3, 0)), // M_sun/(Mpc/h)^3
  SURFACE_DENSITY: Object.freeze(new DimensionTensor(-2, 1, 0, 0, 0, 0, 0, 0, 0)), // kg/m^2, M_sun/pc^2
  SURFACE_DENSITY_H: Object.freeze(new DimensionTensor(-2, 1, 0, 0, 0, 0, 0, 1, 0)), // M_sun/(Mpc/h)^2

  // Force, Energy & Power
  FORCE: Object.freeze(new DimensionTensor(1, 1, -2, 0, 0, 0, 0, 0, 0)), // N, dyne
  ENERGY: Object.freeze(new DimensionTensor(2, 1, -2, 0, 0, 0, 0, 0, 0)), // J, erg, eV
  POWER: Object.freeze(new DimensionTensor(2, 1, -3, 0, 0, 0, 0, 0, 0)), // W, erg/s, L_sun
  PRESSURE: Object.freeze(new DimensionTensor(-1, 1, -2, 0, 0, 0, 0, 0, 0)), // Pa, bar, J/m^3
  ENERGY_DENSITY: Object.freeze(new DimensionTensor(-1, 1, -2, 0, 0, 0, 0, 0, 0)),

  // Electromagnetism & Thermodynamics
  CURRENT: Object.freeze(new DimensionTensor(0, 0, 0, 1, 0, 0, 0, 0, 0)), // A
  ELECTRIC_CHARGE: Object.freeze(new DimensionTensor(0, 0, 1, 1, 0, 0, 0, 0, 0)), // C = A*s
  VOLTAGE: Object.freeze(new DimensionTensor(2, 1, -3, -1, 0, 0, 0, 0, 0)), // V
  MAGNETIC_FIELD: Object.freeze(new DimensionTensor(0, 1, -2, -1, 0, 0, 0, 0, 0)), // Tesla, Gauss
  TEMPERATURE: Object.freeze(new DimensionTensor(0, 0, 0, 0, 1, 0, 0, 0, 0)), // K
  AMOUNT_OF_SUBSTANCE: Object.freeze(new DimensionTensor(0, 0, 0, 0, 0, 1, 0, 0, 0)), // mol
  LUMINOUS_INTENSITY: Object.freeze(new DimensionTensor(0, 0, 0, 0, 0, 0, 1, 0, 0)), // cd

  // Universal Constants Dimensions
  GRAVITATIONAL_CONSTANT: Object.freeze(new DimensionTensor(3, -1, -2, 0, 0, 0, 0, 0, 0)), // m^3 kg^-1 s^-2
  POWER_SPECTRUM: Object.freeze(new DimensionTensor(3, 0, 0, 0, 0, 0, 0, -3, 0)) // P(k) in (Mpc/h)^3
});

// ============================================================================
// 5. COMPREHENSIVE UNIT DEFINITIONS & CANONICAL REGISTRY
// ============================================================================

/**
 * UnitDefinition encapsulates unit metadata, conversion to SI base, and cosmological scaling.
 */
export class UnitDefinition {
  /**
   * @param {Object} config
   * @param {string} config.symbol - Primary symbol (e.g. 'Mpc/h', 'km/s')
   * @param {DimensionTensor} config.dimension - 7-base dimension tensor
   * @param {number} config.siFactor - Multiplier to canonical SI base units (assuming h=1, a=1)
   * @param {number} [config.hExponent=0] - Cosmological h exponent
   * @param {number} [config.aExponent=0] - Cosmological scale factor a exponent
   * @param {string} [config.name] - Human-readable name
   * @param {string} [config.category] - Category (length, velocity, mass, etc.)
   * @param {string[]} [config.aliases=[]] - Alternate symbols
   */
  constructor({
    symbol,
    dimension,
    siFactor,
    hExponent = 0,
    aExponent = 0,
    name = symbol,
    category = 'custom',
    aliases = []
  }) {
    if (!symbol || typeof symbol !== 'string') {
      throw new Error('UnitDefinition requires a valid string symbol');
    }
    this.symbol = symbol;
    this.dimension = DimensionTensor.from(dimension);
    this.siFactor = Number(siFactor);
    this.hExponent = Number(hExponent) || 0;
    this.aExponent = Number(aExponent) || 0;
    this.name = name;
    this.category = category;
    this.aliases = Object.freeze([...aliases]);
    Object.freeze(this);
  }

  /**
   * Computes the effective conversion factor to SI base units for given cosmological parameters.
   * @param {object} [context={}]
   * @param {number} [context.h=0.746]
   * @param {number} [context.a=1.0]
   * @returns {number} Multiplier to SI base
   */
  getEffectiveMultiplierToSI(context = {}) {
    const h = typeof context.h === 'number' ? context.h : 0.746;
    const a = typeof context.a === 'number' ? context.a : 1.0;
    let factor = this.siFactor;
    if (this.hExponent !== 0) {
      factor *= Math.pow(h, this.hExponent);
    }
    if (this.aExponent !== 0) {
      factor *= Math.pow(a, this.aExponent);
    }
    return factor;
  }
}

/**
 * UnitRegistry manages standard physical, cosmological, and user-defined units.
 */
export class UnitRegistry {
  constructor() {
    /** @type {Map<string, UnitDefinition>} */
    this._units = new Map();
    /** @type {Map<string, string>} */
    this._aliases = new Map();
    this._initStandardCatalog();
  }

  /**
   * Registers a UnitDefinition into the registry.
   * @param {UnitDefinition} unitDef
   * @returns {UnitRegistry} this
   */
  register(unitDef) {
    if (!(unitDef instanceof UnitDefinition)) {
      throw new TypeError('UnitRegistry.register requires a UnitDefinition instance');
    }
    this._units.set(unitDef.symbol, unitDef);
    this._aliases.set(unitDef.symbol.toLowerCase(), unitDef.symbol);

    // Register aliases
    for (const alias of unitDef.aliases) {
      this._aliases.set(alias, unitDef.symbol);
      this._aliases.set(alias.toLowerCase(), unitDef.symbol);
    }

    // Register whitespace-stripped aliases
    const stripped = unitDef.symbol.replace(/\s+/g, '');
    if (stripped !== unitDef.symbol) {
      this._aliases.set(stripped, unitDef.symbol);
      this._aliases.set(stripped.toLowerCase(), unitDef.symbol);
    }

    return this;
  }

  /**
   * Resolves and retrieves a UnitDefinition.
   * @param {string} symbol
   * @returns {UnitDefinition}
   */
  get(symbol) {
    if (!symbol || typeof symbol !== 'string') {
      throw new UnitParseError(`Invalid unit symbol: ${symbol}`, String(symbol));
    }
    const trimmed = symbol.trim();
    if (this._units.has(trimmed)) {
      return this._units.get(trimmed);
    }

    // Check alias maps
    const aliasMatch = this._aliases.get(trimmed) || this._aliases.get(trimmed.toLowerCase());
    if (aliasMatch && this._units.has(aliasMatch)) {
      return this._units.get(aliasMatch);
    }

    // Try without internal whitespace
    const noSpace = trimmed.replace(/\s+/g, '');
    if (this._units.has(noSpace)) {
      return this._units.get(noSpace);
    }
    const noSpaceAlias = this._aliases.get(noSpace) || this._aliases.get(noSpace.toLowerCase());
    if (noSpaceAlias && this._units.has(noSpaceAlias)) {
      return this._units.get(noSpaceAlias);
    }

    // Fallback: parse compound unit dynamically
    return this._parseCompoundUnit(trimmed);
  }

  /**
   * Checks if unit symbol is registered or resolvable.
   * @param {string} symbol
   * @returns {boolean}
   */
  has(symbol) {
    try {
      this.get(symbol);
      return true;
    } catch {
      return false;
    }
  }

  /**
   * Converts a numeric value from sourceUnit to targetUnit.
   * 
   * @param {number} value - Magnitude in sourceUnit
   * @param {string|UnitDefinition} fromUnit - Source unit
   * @param {string|UnitDefinition} toUnit - Target unit
   * @param {object} [context] - Cosmological context { h, a, H0, z }
   * @returns {number} Converted magnitude in toUnit
   */
  convert(value, fromUnit, toUnit, context = {}) {
    if (typeof value !== 'number' || isNaN(value)) {
      throw new TypeError(`UnitRegistry.convert value must be a valid number, received ${value}`);
    }

    const src = fromUnit instanceof UnitDefinition ? fromUnit : this.get(fromUnit);
    const dst = toUnit instanceof UnitDefinition ? toUnit : this.get(toUnit);

    if (src.symbol === dst.symbol) return value;

    // Check 7-base dimensional consistency (ignoring h/a for conversion equation check,
    // since h and a are explicitly scaled in getEffectiveMultiplierToSI)
    if (!src.dimension.equals(dst.dimension, false)) {
      throw new DimensionalityMismatchError(
        `Cannot convert incompatible dimensions: "${src.symbol}" (${src.dimension}) to "${dst.symbol}" (${dst.dimension})`,
        {
          sourceDimension: src.dimension,
          targetDimension: dst.dimension,
          sourceUnit: src.symbol,
          targetUnit: dst.symbol,
          operation: 'convert',
          context
        }
      );
    }

    const srcMultiplier = src.getEffectiveMultiplierToSI(context);
    const dstMultiplier = dst.getEffectiveMultiplierToSI(context);

    if (dstMultiplier === 0) {
      throw new RangeError(`Target unit "${dst.symbol}" has zero SI multiplier`);
    }

    const valSI = value * srcMultiplier;
    return valSI / dstMultiplier;
  }

  /**
   * Initializes the exhaustive scientific and cosmological unit catalog.
   * @private
   */
  _initStandardCatalog() {
    // ------------------------------------------------------------------------
    // DIMENSIONLESS
    // ------------------------------------------------------------------------
    const dimLess = PHYSICAL_DIMENSIONS.DIMENSIONLESS;
    this.register(new UnitDefinition({ symbol: '1', dimension: dimLess, siFactor: 1.0, name: 'Dimensionless Scalar', category: 'dimensionless', aliases: ['', 'scalar', 'none', 'dimless'] }));
    this.register(new UnitDefinition({ symbol: 'dimensionless', dimension: dimLess, siFactor: 1.0, name: 'Dimensionless', category: 'dimensionless' }));
    this.register(new UnitDefinition({ symbol: 'delta', dimension: dimLess, siFactor: 1.0, name: 'Density Contrast (delta)', category: 'dimensionless', aliases: ['density_contrast'] }));
    this.register(new UnitDefinition({ symbol: 'percent', dimension: dimLess, siFactor: 0.01, name: 'Percentage', category: 'dimensionless', aliases: ['%'] }));
    this.register(new UnitDefinition({ symbol: 'ppm', dimension: dimLess, siFactor: 1e-6, name: 'Parts Per Million', category: 'dimensionless' }));
    this.register(new UnitDefinition({ symbol: 'mag', dimension: dimLess, siFactor: 1.0, name: 'Astronomical Magnitude', category: 'dimensionless', aliases: ['magnitude'] }));

    // ------------------------------------------------------------------------
    // LENGTH & DISTANCE
    // ------------------------------------------------------------------------
    const dLen = PHYSICAL_DIMENSIONS.LENGTH;
    const dLenH = PHYSICAL_DIMENSIONS.LENGTH_H_INV;

    // Standard metric lengths
    this.register(new UnitDefinition({ symbol: 'm', dimension: dLen, siFactor: 1.0, name: 'Meter', category: 'length', aliases: ['meter', 'meters'] }));
    this.register(new UnitDefinition({ symbol: 'km', dimension: dLen, siFactor: 1.0e3, name: 'Kilometer', category: 'length', aliases: ['kilometer', 'kilometers'] }));
    this.register(new UnitDefinition({ symbol: 'cm', dimension: dLen, siFactor: 1.0e-2, name: 'Centimeter', category: 'length', aliases: ['centimeter'] }));
    this.register(new UnitDefinition({ symbol: 'mm', dimension: dLen, siFactor: 1.0e-3, name: 'Millimeter', category: 'length', aliases: ['millimeter'] }));
    this.register(new UnitDefinition({ symbol: 'um', dimension: dLen, siFactor: 1.0e-6, name: 'Micrometer / Micron', category: 'length', aliases: ['micron', 'micrometer'] }));
    this.register(new UnitDefinition({ symbol: 'nm', dimension: dLen, siFactor: 1.0e-9, name: 'Nanometer', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'pm', dimension: dLen, siFactor: 1.0e-12, name: 'Picometer', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'fm', dimension: dLen, siFactor: 1.0e-15, name: 'Femtometer / Fermi', category: 'length', aliases: ['fermi'] }));
    this.register(new UnitDefinition({ symbol: 'angstrom', dimension: dLen, siFactor: 1.0e-10, name: 'Angstrom', category: 'length', aliases: ['Å', 'Angstrom', 'A_ring'] }));

    // Astronomical lengths
    this.register(new UnitDefinition({ symbol: 'au', dimension: dLen, siFactor: CONSTANTS.METERS_PER_AU, name: 'Astronomical Unit', category: 'length', aliases: ['AU', 'Au'] }));
    this.register(new UnitDefinition({ symbol: 'ly', dimension: dLen, siFactor: CONSTANTS.METERS_PER_LY, name: 'Light Year', category: 'length', aliases: ['light-year', 'lightyear', 'lyr'] }));
    this.register(new UnitDefinition({ symbol: 'pc', dimension: dLen, siFactor: CONSTANTS.METERS_PER_PC, name: 'Parsec', category: 'length', aliases: ['parsec'] }));
    this.register(new UnitDefinition({ symbol: 'kpc', dimension: dLen, siFactor: CONSTANTS.METERS_PER_KPC, name: 'Kiloparsec', category: 'length', aliases: ['kiloparsec'] }));
    this.register(new UnitDefinition({ symbol: 'Mpc', dimension: dLen, siFactor: CONSTANTS.METERS_PER_MPC, name: 'Megaparsec', category: 'length', aliases: ['megaparsec', 'mpc'] }));
    this.register(new UnitDefinition({ symbol: 'Gpc', dimension: dLen, siFactor: CONSTANTS.METERS_PER_GPC, name: 'Gigaparsec', category: 'length', aliases: ['gigaparsec', 'gpc'] }));

    // Light-travel times as lengths
    this.register(new UnitDefinition({ symbol: 'ls', dimension: dLen, siFactor: CONSTANTS.SPEED_OF_LIGHT_SI * 1.0, name: 'Light Second', category: 'length', aliases: ['light-second'] }));
    this.register(new UnitDefinition({ symbol: 'lm', dimension: dLen, siFactor: CONSTANTS.SPEED_OF_LIGHT_SI * 60.0, name: 'Light Minute', category: 'length', aliases: ['light-minute'] }));
    this.register(new UnitDefinition({ symbol: 'lh', dimension: dLen, siFactor: CONSTANTS.SPEED_OF_LIGHT_SI * 3600.0, name: 'Light Hour', category: 'length', aliases: ['light-hour'] }));
    this.register(new UnitDefinition({ symbol: 'ld', dimension: dLen, siFactor: CONSTANTS.SPEED_OF_LIGHT_SI * 86400.0, name: 'Light Day', category: 'length', aliases: ['light-day'] }));

    // Cosmological h-dependent Lengths: 1 Mpc/h = (1/h) * Mpc -> siFactor = METERS_PER_MPC, hExponent = -1
    this.register(new UnitDefinition({ symbol: 'Mpc/h', dimension: dLenH, siFactor: CONSTANTS.METERS_PER_MPC, hExponent: -1, name: 'Megaparsec per h (CF4 Grid Distance)', category: 'length', aliases: ['mpc/h', 'h^-1 Mpc', 'Mpc*h^-1', 'Mpc h^-1'] }));
    this.register(new UnitDefinition({ symbol: 'kpc/h', dimension: dLenH, siFactor: CONSTANTS.METERS_PER_KPC, hExponent: -1, name: 'Kiloparsec per h', category: 'length', aliases: ['kpc h^-1', 'h^-1 kpc'] }));
    this.register(new UnitDefinition({ symbol: 'pc/h', dimension: dLenH, siFactor: CONSTANTS.METERS_PER_PC, hExponent: -1, name: 'Parsec per h', category: 'length', aliases: ['pc h^-1', 'h^-1 pc'] }));
    this.register(new UnitDefinition({ symbol: 'Gpc/h', dimension: dLenH, siFactor: CONSTANTS.METERS_PER_GPC, hExponent: -1, name: 'Gigaparsec per h', category: 'length', aliases: ['Gpc h^-1', 'h^-1 Gpc'] }));

    // ------------------------------------------------------------------------
    // VELOCITY
    // ------------------------------------------------------------------------
    const dVel = PHYSICAL_DIMENSIONS.VELOCITY;
    this.register(new UnitDefinition({ symbol: 'm/s', dimension: dVel, siFactor: 1.0, name: 'Meter per second (SI)', category: 'velocity', aliases: ['m s^-1', 'mps'] }));
    this.register(new UnitDefinition({ symbol: 'km/s', dimension: dVel, siFactor: 1.0e3, name: 'Kilometer per second (CF4 Peculiar Velocity)', category: 'velocity', aliases: ['kms', 'km s^-1', 'km/sec', 'kms^-1'] }));
    this.register(new UnitDefinition({ symbol: 'cm/s', dimension: dVel, siFactor: 1.0e-2, name: 'Centimeter per second (CGS)', category: 'velocity', aliases: ['cm/sec', 'cm s^-1'] }));
    this.register(new UnitDefinition({ symbol: 'mm/s', dimension: dVel, siFactor: 1.0e-3, name: 'Millimeter per second', category: 'velocity' }));
    this.register(new UnitDefinition({ symbol: 'km/h', dimension: dVel, siFactor: 1.0e3 / 3600.0, name: 'Kilometer per hour', category: 'velocity', aliases: ['kph'] }));
    this.register(new UnitDefinition({ symbol: 'mph', dimension: dVel, siFactor: 1609.344 / 3600.0, name: 'Mile per hour', category: 'velocity' }));
    this.register(new UnitDefinition({ symbol: 'c', dimension: dVel, siFactor: CONSTANTS.SPEED_OF_LIGHT_SI, name: 'Speed of Light in Vacuum', category: 'velocity', aliases: ['speed_of_light'] }));
    this.register(new UnitDefinition({ symbol: 'pc/Myr', dimension: dVel, siFactor: CONSTANTS.METERS_PER_PC / CONSTANTS.SECONDS_PER_MYR, name: 'Parsec per Megayear (~0.9778 km/s)', category: 'velocity', aliases: ['pc/myr'] }));
    this.register(new UnitDefinition({ symbol: 'AU/yr', dimension: dVel, siFactor: CONSTANTS.METERS_PER_AU / CONSTANTS.SECONDS_PER_JULIAN_YEAR, name: 'Astronomical Unit per year (~4.7404 km/s)', category: 'velocity', aliases: ['au/yr'] }));

    // ------------------------------------------------------------------------
    // TIME & FREQUENCY & HUBBLE EXPANSION RATE
    // ------------------------------------------------------------------------
    const dTime = PHYSICAL_DIMENSIONS.TIME;
    const dFreq = PHYSICAL_DIMENSIONS.FREQUENCY;
    const dExpH = PHYSICAL_DIMENSIONS.EXPANSION_RATE_H;

    this.register(new UnitDefinition({ symbol: 's', dimension: dTime, siFactor: 1.0, name: 'Second (SI)', category: 'time', aliases: ['sec', 'second', 'seconds'] }));
    this.register(new UnitDefinition({ symbol: 'ms', dimension: dTime, siFactor: 1.0e-3, name: 'Millisecond', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'us', dimension: dTime, siFactor: 1.0e-6, name: 'Microsecond', category: 'time', aliases: ['μs'] }));
    this.register(new UnitDefinition({ symbol: 'ns', dimension: dTime, siFactor: 1.0e-9, name: 'Nanosecond', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'min', dimension: dTime, siFactor: CONSTANTS.SECONDS_PER_MINUTE, name: 'Minute', category: 'time', aliases: ['minute', 'minutes'] }));
    this.register(new UnitDefinition({ symbol: 'hr', dimension: dTime, siFactor: CONSTANTS.SECONDS_PER_HOUR, name: 'Hour', category: 'time', aliases: ['hour', 'hours', 'h_time'] }));
    this.register(new UnitDefinition({ symbol: 'day', dimension: dTime, siFactor: CONSTANTS.SECONDS_PER_DAY, name: 'Day', category: 'time', aliases: ['days', 'd'] }));
    this.register(new UnitDefinition({ symbol: 'yr', dimension: dTime, siFactor: CONSTANTS.SECONDS_PER_JULIAN_YEAR, name: 'Julian Year', category: 'time', aliases: ['year', 'years', 'a_year'] }));
    this.register(new UnitDefinition({ symbol: 'kyr', dimension: dTime, siFactor: CONSTANTS.SECONDS_PER_KYR, name: 'Kiloyear', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'Myr', dimension: dTime, siFactor: CONSTANTS.SECONDS_PER_MYR, name: 'Megayear', category: 'time', aliases: ['myr'] }));
    this.register(new UnitDefinition({ symbol: 'Gyr', dimension: dTime, siFactor: CONSTANTS.SECONDS_PER_GYR, name: 'Gigayear', category: 'time', aliases: ['gyr'] }));

    // Frequencies & Expansion Rates
    this.register(new UnitDefinition({ symbol: 'Hz', dimension: dFreq, siFactor: 1.0, name: 'Hertz (1/s)', category: 'frequency', aliases: ['1/s', 's^-1'] }));
    this.register(new UnitDefinition({ symbol: 'kHz', dimension: dFreq, siFactor: 1.0e3, name: 'Kilohertz', category: 'frequency' }));
    this.register(new UnitDefinition({ symbol: 'MHz', dimension: dFreq, siFactor: 1.0e6, name: 'Megahertz', category: 'frequency' }));
    this.register(new UnitDefinition({ symbol: 'GHz', dimension: dFreq, siFactor: 1.0e9, name: 'Gigahertz', category: 'frequency' }));
    this.register(new UnitDefinition({ symbol: '1/yr', dimension: dFreq, siFactor: 1.0 / CONSTANTS.SECONDS_PER_JULIAN_YEAR, name: 'Per Year', category: 'frequency', aliases: ['yr^-1'] }));
    this.register(new UnitDefinition({ symbol: '1/Gyr', dimension: dFreq, siFactor: 1.0 / CONSTANTS.SECONDS_PER_GYR, name: 'Per Gigayear', category: 'frequency', aliases: ['Gyr^-1'] }));

    // Hubble Parameters:
    // km/s/Mpc in SI: (1000 m/s) / (3.085677581491367e22 m) = 3.240779289e-18 s^-1
    this.register(new UnitDefinition({
      symbol: 'km/s/Mpc',
      dimension: dFreq,
      siFactor: CONSTANTS.KMS_MPC_TO_SEC_INV,
      hExponent: 0,
      name: 'Kilometers per second per Megaparsec',
      category: 'expansion_rate',
      aliases: ['km/s/mpc', 'km s^-1 Mpc^-1', 'km/s Mpc^-1']
    }));

    // km/s/(Mpc/h) = h * (km/s/Mpc) -> siFactor = KMS_MPC_TO_SEC_INV, hExponent = +1
    this.register(new UnitDefinition({
      symbol: 'km/s/(Mpc/h)',
      dimension: dExpH,
      siFactor: CONSTANTS.KMS_MPC_TO_SEC_INV,
      hExponent: 1,
      name: 'Kilometers per second per (Mpc/h) [100 h km/s/Mpc]',
      category: 'expansion_rate',
      aliases: ['km/s/(mpc/h)', 'km/s/Mpc/h', '100 h km/s/Mpc']
    }));

    // ------------------------------------------------------------------------
    // MASS
    // ------------------------------------------------------------------------
    const dMass = PHYSICAL_DIMENSIONS.MASS;
    const dMassH = PHYSICAL_DIMENSIONS.MASS_H_INV;

    this.register(new UnitDefinition({ symbol: 'kg', dimension: dMass, siFactor: 1.0, name: 'Kilogram (SI)', category: 'mass', aliases: ['kilogram', 'kilograms'] }));
    this.register(new UnitDefinition({ symbol: 'g', dimension: dMass, siFactor: 1.0e-3, name: 'Gram (CGS)', category: 'mass', aliases: ['gram', 'grams'] }));
    this.register(new UnitDefinition({ symbol: 'mg', dimension: dMass, siFactor: 1.0e-6, name: 'Milligram', category: 'mass' }));
    this.register(new UnitDefinition({ symbol: 'tonne', dimension: dMass, siFactor: 1.0e3, name: 'Metric Tonne', category: 'mass', aliases: ['t', 'ton'] }));

    // Astronomical masses
    this.register(new UnitDefinition({ symbol: 'M_sun', dimension: dMass, siFactor: CONSTANTS.SOLAR_MASS_KG, name: 'Solar Mass', category: 'mass', aliases: ['Msun', 'MSun', 'M_Sun', 'solMass'] }));
    this.register(new UnitDefinition({ symbol: 'M_earth', dimension: dMass, siFactor: CONSTANTS.EARTH_MASS_KG, name: 'Earth Mass', category: 'mass', aliases: ['Mearth'] }));
    this.register(new UnitDefinition({ symbol: 'M_jupiter', dimension: dMass, siFactor: CONSTANTS.JUPITER_MASS_KG, name: 'Jupiter Mass', category: 'mass', aliases: ['Mjup'] }));
    this.register(new UnitDefinition({ symbol: '1e10 M_sun', dimension: dMass, siFactor: 1.0e10 * CONSTANTS.SOLAR_MASS_KG, name: '10^10 Solar Masses', category: 'mass', aliases: ['10^10 M_sun', '10^10 Msun'] }));
    this.register(new UnitDefinition({ symbol: '1e12 M_sun', dimension: dMass, siFactor: 1.0e12 * CONSTANTS.SOLAR_MASS_KG, name: '10^12 Solar Masses (Milky Way scale)', category: 'mass', aliases: ['10^12 M_sun', '10^12 Msun'] }));

    // Cosmological h-dependent masses: M_sun/h = (1/h) * M_sun -> siFactor = SOLAR_MASS_KG, hExponent = -1
    this.register(new UnitDefinition({ symbol: 'M_sun/h', dimension: dMassH, siFactor: CONSTANTS.SOLAR_MASS_KG, hExponent: -1, name: 'Solar Mass per h', category: 'mass', aliases: ['Msun/h', 'h^-1 M_sun', 'M_sun h^-1'] }));
    this.register(new UnitDefinition({ symbol: '1e10 M_sun/h', dimension: dMassH, siFactor: 1.0e10 * CONSTANTS.SOLAR_MASS_KG, hExponent: -1, name: '10^10 Solar Masses per h (Standard GADGET Mass)', category: 'mass', aliases: ['10^10 M_sun/h', '10^10 Msun/h', '1e10 Msun/h'] }));
    this.register(new UnitDefinition({ symbol: '1e12 M_sun/h', dimension: dMassH, siFactor: 1.0e12 * CONSTANTS.SOLAR_MASS_KG, hExponent: -1, name: '10^12 Solar Masses per h', category: 'mass', aliases: ['10^12 M_sun/h'] }));

    // ------------------------------------------------------------------------
    // DENSITY (MASS DENSITY)
    // ------------------------------------------------------------------------
    const dDens = PHYSICAL_DIMENSIONS.MASS_DENSITY;
    const dDensH2 = PHYSICAL_DIMENSIONS.MASS_DENSITY_H2;
    const dDensH3 = PHYSICAL_DIMENSIONS.MASS_DENSITY_H3;
    const msunMpc3SI = CONSTANTS.SOLAR_MASS_KG / Math.pow(CONSTANTS.METERS_PER_MPC, 3);

    this.register(new UnitDefinition({ symbol: 'kg/m^3', dimension: dDens, siFactor: 1.0, name: 'Kilograms per cubic meter (SI)', category: 'density', aliases: ['kg m^-3', 'kg/m3'] }));
    this.register(new UnitDefinition({ symbol: 'g/cm^3', dimension: dDens, siFactor: 1.0e3, name: 'Grams per cubic centimeter (CGS)', category: 'density', aliases: ['g cm^-3', 'g/cm3', 'g/cc'] }));
    this.register(new UnitDefinition({ symbol: 'M_sun/Mpc^3', dimension: dDens, siFactor: msunMpc3SI, hExponent: 0, name: 'Solar Masses per cubic Megaparsec', category: 'density', aliases: ['Msun/Mpc^3', 'M_sun Mpc^-3'] }));
    this.register(new UnitDefinition({ symbol: 'M_sun/kpc^3', dimension: dDens, siFactor: CONSTANTS.SOLAR_MASS_KG / Math.pow(CONSTANTS.METERS_PER_KPC, 3), hExponent: 0, name: 'Solar Masses per cubic Kiloparsec', category: 'density' }));

    // h-dependent densities
    // (M_sun/h) / (Mpc/h)^3 = (M_sun * h^-1) / (Mpc^3 * h^-3) = h^2 * (M_sun / Mpc^3) -> hExponent = +2
    this.register(new UnitDefinition({
      symbol: '(M_sun/h)/(Mpc/h)^3',
      dimension: dDensH2,
      siFactor: msunMpc3SI,
      hExponent: 2,
      name: 'Solar Mass over h per (Mpc/h)^3',
      category: 'density',
      aliases: ['(Msun/h)/(Mpc/h)^3', 'M_sun/h/(Mpc/h)^3', 'h^2 M_sun/Mpc^3']
    }));

    // M_sun / (Mpc/h)^3 = M_sun / (Mpc^3 * h^-3) = h^3 * (M_sun / Mpc^3) -> hExponent = +3
    this.register(new UnitDefinition({
      symbol: 'M_sun/(Mpc/h)^3',
      dimension: dDensH3,
      siFactor: msunMpc3SI,
      hExponent: 3,
      name: 'Solar Mass per cubic (Mpc/h)',
      category: 'density',
      aliases: ['Msun/(Mpc/h)^3']
    }));

    // Critical density unit at z=0: rho_crit,0 = 2.77536627e11 * h^2 M_sun / Mpc^3
    this.register(new UnitDefinition({
      symbol: 'rho_crit',
      dimension: dDensH2,
      siFactor: CONSTANTS.CRITICAL_DENSITY_COEFF_MSUN_MPC3 * msunMpc3SI,
      hExponent: 2,
      name: 'Cosmological Critical Density rho_c,0',
      category: 'density',
      aliases: ['rho_crit_0', 'rho_c0', 'rho_critical']
    }));

    // ------------------------------------------------------------------------
    // GRAVITATIONAL POTENTIAL / SPECIFIC ENERGY
    // ------------------------------------------------------------------------
    const dPot = PHYSICAL_DIMENSIONS.SPECIFIC_ENERGY;
    this.register(new UnitDefinition({ symbol: 'm^2/s^2', dimension: dPot, siFactor: 1.0, name: 'Square meter per square second (SI)', category: 'potential', aliases: ['J/kg'] }));
    this.register(new UnitDefinition({ symbol: '(km/s)^2', dimension: dPot, siFactor: 1.0e6, name: 'Square kilometer per square second (Astrophysical Potential)', category: 'potential', aliases: ['km^2/s^2', '(kms)^2', 'km2/s2'] }));
    this.register(new UnitDefinition({ symbol: 'erg/g', dimension: dPot, siFactor: 1.0e-4, name: 'Ergs per gram (CGS)', category: 'potential' }));
    this.register(new UnitDefinition({ symbol: 'c^2', dimension: dPot, siFactor: CONSTANTS.SPEED_OF_LIGHT_SI * CONSTANTS.SPEED_OF_LIGHT_SI, name: 'Speed of Light Squared', category: 'potential' }));

    // ------------------------------------------------------------------------
    // ACCELERATION & GRAVITATIONAL FIELD
    // ------------------------------------------------------------------------
    const dAcc = PHYSICAL_DIMENSIONS.ACCELERATION;
    const dAccH = PHYSICAL_DIMENSIONS.ACCELERATION_ASTRO;
    this.register(new UnitDefinition({ symbol: 'm/s^2', dimension: dAcc, siFactor: 1.0, name: 'Meter per square second (SI)', category: 'acceleration', aliases: ['m s^-2', 'm/s2'] }));
    this.register(new UnitDefinition({ symbol: 'cm/s^2', dimension: dAcc, siFactor: 1.0e-2, name: 'Gal / Centimeter per square second (CGS)', category: 'acceleration', aliases: ['Gal', 'gal'] }));
    this.register(new UnitDefinition({ symbol: 'g_earth', dimension: dAcc, siFactor: 9.80665, name: 'Standard Earth Gravity (g0)', category: 'acceleration', aliases: ['g0', 'gee'] }));
    // (km/s)^2 / (Mpc/h) = 1.0e6 / (METERS_PER_MPC / h) = (1.0e6 / METERS_PER_MPC) * h m/s^2 -> hExponent = +1
    this.register(new UnitDefinition({
      symbol: '(km/s)^2/(Mpc/h)',
      dimension: dAccH,
      siFactor: 1.0e6 / CONSTANTS.METERS_PER_MPC,
      hExponent: 1,
      name: '(km/s)^2 per (Mpc/h) [Astrophysical Acceleration]',
      category: 'acceleration',
      aliases: ['(km/s)^2/(mpc/h)', '(km/s)2/(Mpc/h)']
    }));
    this.register(new UnitDefinition({
      symbol: '(km/s)^2/Mpc',
      dimension: dAcc,
      siFactor: 1.0e6 / CONSTANTS.METERS_PER_MPC,
      hExponent: 0,
      name: '(km/s)^2 per Mpc',
      category: 'acceleration'
    }));

    // ------------------------------------------------------------------------
    // FORCE, ENERGY, POWER & PRESSURE
    // ------------------------------------------------------------------------
    const dForce = PHYSICAL_DIMENSIONS.FORCE;
    const dEnergy = PHYSICAL_DIMENSIONS.ENERGY;
    const dPower = PHYSICAL_DIMENSIONS.POWER;
    const dPress = PHYSICAL_DIMENSIONS.PRESSURE;

    // Force
    this.register(new UnitDefinition({ symbol: 'N', dimension: dForce, siFactor: 1.0, name: 'Newton (SI)', category: 'force', aliases: ['newton'] }));
    this.register(new UnitDefinition({ symbol: 'dyn', dimension: dForce, siFactor: 1.0e-5, name: 'Dyne (CGS)', category: 'force', aliases: ['dyne'] }));

    // Energy
    this.register(new UnitDefinition({ symbol: 'J', dimension: dEnergy, siFactor: 1.0, name: 'Joule (SI)', category: 'energy', aliases: ['joule', 'joules'] }));
    this.register(new UnitDefinition({ symbol: 'erg', dimension: dEnergy, siFactor: 1.0e-7, name: 'Erg (CGS)', category: 'energy', aliases: ['ergs'] }));
    this.register(new UnitDefinition({ symbol: 'eV', dimension: dEnergy, siFactor: CONSTANTS.ELECTRON_VOLT_JOULES, name: 'Electronvolt', category: 'energy' }));
    this.register(new UnitDefinition({ symbol: 'keV', dimension: dEnergy, siFactor: 1.0e3 * CONSTANTS.ELECTRON_VOLT_JOULES, name: 'Kiloelectronvolt', category: 'energy' }));
    this.register(new UnitDefinition({ symbol: 'MeV', dimension: dEnergy, siFactor: 1.0e6 * CONSTANTS.ELECTRON_VOLT_JOULES, name: 'Megaelectronvolt', category: 'energy' }));
    this.register(new UnitDefinition({ symbol: 'GeV', dimension: dEnergy, siFactor: 1.0e9 * CONSTANTS.ELECTRON_VOLT_JOULES, name: 'Gigaelectronvolt', category: 'energy' }));

    // Power
    this.register(new UnitDefinition({ symbol: 'W', dimension: dPower, siFactor: 1.0, name: 'Watt (SI)', category: 'power', aliases: ['watt', 'watts', 'J/s'] }));
    this.register(new UnitDefinition({ symbol: 'kW', dimension: dPower, siFactor: 1.0e3, name: 'Kilowatt', category: 'power' }));
    this.register(new UnitDefinition({ symbol: 'MW', dimension: dPower, siFactor: 1.0e6, name: 'Megawatt', category: 'power' }));
    this.register(new UnitDefinition({ symbol: 'erg/s', dimension: dPower, siFactor: 1.0e-7, name: 'Ergs per second (CGS)', category: 'power' }));
    this.register(new UnitDefinition({ symbol: 'L_sun', dimension: dPower, siFactor: CONSTANTS.SOLAR_LUMINOSITY_WATTS, name: 'Solar Luminosity', category: 'power', aliases: ['Lsun', 'L_Sun'] }));

    // Pressure & Energy Density
    this.register(new UnitDefinition({ symbol: 'Pa', dimension: dPress, siFactor: 1.0, name: 'Pascal (SI)', category: 'pressure', aliases: ['pascal', 'N/m^2'] }));
    this.register(new UnitDefinition({ symbol: 'bar', dimension: dPress, siFactor: 1.0e5, name: 'Bar', category: 'pressure' }));
    this.register(new UnitDefinition({ symbol: 'atm', dimension: dPress, siFactor: 101325.0, name: 'Standard Atmosphere', category: 'pressure' }));
    this.register(new UnitDefinition({ symbol: 'dyn/cm^2', dimension: dPress, siFactor: 0.1, name: 'Dynes per square centimeter (CGS)', category: 'pressure', aliases: ['erg/cm^3', 'erg/cm3'] }));
    this.register(new UnitDefinition({ symbol: 'J/m^3', dimension: dPress, siFactor: 1.0, name: 'Joules per cubic meter', category: 'pressure', aliases: ['J/m3'] }));

    // ------------------------------------------------------------------------
    // THERMODYNAMICS & ELECTROMAGNETISM
    // ------------------------------------------------------------------------
    const dTemp = PHYSICAL_DIMENSIONS.TEMPERATURE;
    const dCurr = PHYSICAL_DIMENSIONS.CURRENT;
    const dCharge = PHYSICAL_DIMENSIONS.ELECTRIC_CHARGE;
    const dVolt = PHYSICAL_DIMENSIONS.VOLTAGE;
    const dMag = PHYSICAL_DIMENSIONS.MAGNETIC_FIELD;

    this.register(new UnitDefinition({ symbol: 'K', dimension: dTemp, siFactor: 1.0, name: 'Kelvin (SI)', category: 'temperature', aliases: ['kelvin'] }));
    this.register(new UnitDefinition({ symbol: 'A', dimension: dCurr, siFactor: 1.0, name: 'Ampere (SI)', category: 'current', aliases: ['amp', 'ampere'] }));
    this.register(new UnitDefinition({ symbol: 'mA', dimension: dCurr, siFactor: 1.0e-3, name: 'Milliampere', category: 'current' }));
    this.register(new UnitDefinition({ symbol: 'uA', dimension: dCurr, siFactor: 1.0e-6, name: 'Microampere', category: 'current', aliases: ['μA'] }));
    this.register(new UnitDefinition({ symbol: 'C', dimension: dCharge, siFactor: 1.0, name: 'Coulomb (SI)', category: 'charge', aliases: ['coulomb'] }));
    this.register(new UnitDefinition({ symbol: 'V', dimension: dVolt, siFactor: 1.0, name: 'Volt (SI)', category: 'voltage', aliases: ['volt'] }));
    this.register(new UnitDefinition({ symbol: 'T', dimension: dMag, siFactor: 1.0, name: 'Tesla (SI)', category: 'magnetic_field', aliases: ['tesla'] }));
    this.register(new UnitDefinition({ symbol: 'G', dimension: dMag, siFactor: 1.0e-4, name: 'Gauss (CGS)', category: 'magnetic_field', aliases: ['gauss'] }));
    this.register(new UnitDefinition({ symbol: 'uG', dimension: dMag, siFactor: 1.0e-10, name: 'Microgauss', category: 'magnetic_field', aliases: ['μG'] }));
    this.register(new UnitDefinition({ symbol: 'nG', dimension: dMag, siFactor: 1.0e-13, name: 'Nanogauss', category: 'magnetic_field' }));

    // ------------------------------------------------------------------------
    // ANGULAR MEASURES (Dimensionless in SI base, but dimensionally tracked)
    // ------------------------------------------------------------------------
    this.register(new UnitDefinition({ symbol: 'rad', dimension: dimLess, siFactor: 1.0, name: 'Radian', category: 'angle', aliases: ['radian', 'radians'] }));
    this.register(new UnitDefinition({ symbol: 'deg', dimension: dimLess, siFactor: Math.PI / 180.0, name: 'Degree', category: 'angle', aliases: ['degree', 'degrees', '°'] }));
    this.register(new UnitDefinition({ symbol: 'arcmin', dimension: dimLess, siFactor: Math.PI / (180.0 * 60.0), name: 'Arcminute', category: 'angle', aliases: ['arcminute', 'arcminutes', 'amin', "'"] }));
    this.register(new UnitDefinition({ symbol: 'arcsec', dimension: dimLess, siFactor: Math.PI / (180.0 * 3600.0), name: 'Arcsecond', category: 'angle', aliases: ['arcsecond', 'arcseconds', 'asec', '"'] }));
    this.register(new UnitDefinition({ symbol: 'mas', dimension: dimLess, siFactor: Math.PI / (180.0 * 3600.0 * 1.0e3), name: 'Milliarcsecond', category: 'angle' }));
    this.register(new UnitDefinition({ symbol: 'uas', dimension: dimLess, siFactor: Math.PI / (180.0 * 3600.0 * 1.0e6), name: 'Microarcsecond', category: 'angle', aliases: ['μas'] }));
    this.register(new UnitDefinition({ symbol: 'sr', dimension: dimLess, siFactor: 1.0, name: 'Steradian', category: 'solid_angle', aliases: ['steradian'] }));
    this.register(new UnitDefinition({ symbol: 'sq_deg', dimension: dimLess, siFactor: Math.pow(Math.PI / 180.0, 2), name: 'Square Degree', category: 'solid_angle', aliases: ['deg^2', 'deg2'] }));
  }

  /**
   * Dynamically parses compound unit expressions (e.g. "km/s / (Mpc/h)", "M_sun / kpc^3", "kg * m / s^2").
   * @private
   * @param {string} expr
   * @returns {UnitDefinition}
   */
  _parseCompoundUnit(expr) {
    const trimmed = expr.trim();
    if (!trimmed) {
      return this.get('1');
    }

    // Tokenize expression
    const tokens = this._tokenizeCompoundUnit(trimmed);
    const parsed = this._evaluateUnitTokens(tokens, expr);

    // Create and cache ad-hoc UnitDefinition
    const customDef = new UnitDefinition({
      symbol: trimmed,
      dimension: parsed.dimension,
      siFactor: parsed.siFactor,
      hExponent: parsed.hExponent,
      aExponent: parsed.aExponent,
      name: `Compound (${trimmed})`,
      category: 'compound'
    });

    this.register(customDef);
    return customDef;
  }

  /**
   * Tokenizes a compound unit string into symbols, operators (*, /, ^), parentheses, and numbers.
   * @private
   */
  _tokenizeCompoundUnit(expr) {
    const tokens = [];
    let i = 0;
    const len = expr.length;

    while (i < len) {
      const c = expr[i];

      if (/\s/.test(c)) {
        i++;
        continue;
      }

      if (c === '*' || c === '·' || c === '•') {
        tokens.push({ type: 'OP', value: '*' });
        i++;
      } else if (c === '/') {
        tokens.push({ type: 'OP', value: '/' });
        i++;
      } else if (c === '^') {
        tokens.push({ type: 'OP', value: '^' });
        i++;
      } else if (c === '(' || c === '[' || c === '{') {
        tokens.push({ type: 'LPAREN', value: '(' });
        i++;
      } else if (c === ')' || c === ']' || c === '}') {
        tokens.push({ type: 'RPAREN', value: ')' });
        i++;
      } else if (c === '-' || c === '+' || /\d/.test(c)) {
        // Number / signed integer (exponent or multiplier)
        let numStr = c;
        i++;
        while (i < len && /[\d.]/.test(expr[i])) {
          numStr += expr[i];
          i++;
        }
        const val = parseFloat(numStr);
        if (!isNaN(val)) {
          tokens.push({ type: 'NUMBER', value: val });
        } else {
          tokens.push({ type: 'IDENT', value: numStr });
        }
      } else {
        // Identifier / Unit symbol token
        let ident = '';
        while (i < len && !/[\s*·•/\^()[\]{}]/.test(expr[i])) {
          ident += expr[i];
          i++;
        }
        tokens.push({ type: 'IDENT', value: ident });
      }
    }

    return tokens;
  }

  /**
   * Evaluates tokenized unit expression using shunting-yard / recursive descent.
   * @private
   */
  _evaluateUnitTokens(tokens, originalExpr) {
    let pos = 0;

    const parsePrimary = () => {
      if (pos >= tokens.length) {
        throw new UnitParseError('Unexpected end of unit expression', originalExpr);
      }
      const tok = tokens[pos];

      if (tok.type === 'LPAREN') {
        pos++;
        const inner = parseExpr();
        if (pos >= tokens.length || tokens[pos].type !== 'RPAREN') {
          throw new UnitParseError('Mismatched parentheses in unit expression', originalExpr);
        }
        pos++;
        return inner;
      }

      if (tok.type === 'NUMBER') {
        pos++;
        return {
          dimension: PHYSICAL_DIMENSIONS.DIMENSIONLESS,
          siFactor: tok.value,
          hExponent: 0,
          aExponent: 0
        };
      }

      if (tok.type === 'IDENT') {
        pos++;
        const symbol = tok.value;

        // Check registered
        if (this._units.has(symbol) || this._aliases.has(symbol) || this._aliases.has(symbol.toLowerCase())) {
          const resolved = this.get(symbol);
          return {
            dimension: resolved.dimension,
            siFactor: resolved.siFactor,
            hExponent: resolved.hExponent,
            aExponent: resolved.aExponent
          };
        }

        // Check for 'h' or 'h^-1'
        if (symbol === 'h') {
          return {
            dimension: PHYSICAL_DIMENSIONS.DIMENSIONLESS,
            siFactor: 1.0,
            hExponent: 1,
            aExponent: 0
          };
        }

        // Check for 'a'
        if (symbol === 'a') {
          return {
            dimension: PHYSICAL_DIMENSIONS.DIMENSIONLESS,
            siFactor: 1.0,
            hExponent: 0,
            aExponent: 1
          };
        }

        // Standard SI prefix matching
        const prefixParsed = this._tryMatchPrefix(symbol);
        if (prefixParsed) {
          return prefixParsed;
        }

        throw new UnitParseError(`Unknown unit symbol "${symbol}"`, originalExpr);
      }

      throw new UnitParseError(`Unexpected token "${tok.value}" in unit expression`, originalExpr);
    };

    const parsePower = () => {
      let left = parsePrimary();
      while (pos < tokens.length && tokens[pos].type === 'OP' && tokens[pos].value === '^') {
        pos++;
        if (pos >= tokens.length) {
          throw new UnitParseError('Expected exponent number after "^"', originalExpr);
        }
        const expTok = tokens[pos];
        if (expTok.type !== 'NUMBER') {
          throw new UnitParseError(`Expected numeric exponent after "^", received "${expTok.value}"`, originalExpr);
        }
        pos++;
        const p = expTok.value;
        left = {
          dimension: left.dimension.pow(p),
          siFactor: Math.pow(left.siFactor, p),
          hExponent: left.hExponent * p,
          aExponent: left.aExponent * p
        };
      }
      return left;
    };

    const parseExpr = () => {
      let left = parsePower();

      while (pos < tokens.length && tokens[pos].type !== 'RPAREN') {
        const tok = tokens[pos];
        if (tok.type === 'OP' && (tok.value === '*' || tok.value === '/')) {
          pos++;
          const right = parsePower();
          if (tok.value === '*') {
            left = {
              dimension: left.dimension.multiply(right.dimension),
              siFactor: left.siFactor * right.siFactor,
              hExponent: left.hExponent + right.hExponent,
              aExponent: left.aExponent + right.aExponent
            };
          } else {
            left = {
              dimension: left.dimension.divide(right.dimension),
              siFactor: left.siFactor / right.siFactor,
              hExponent: left.hExponent - right.hExponent,
              aExponent: left.aExponent - right.aExponent
            };
          }
        } else if (tok.type === 'IDENT' || tok.type === 'LPAREN') {
          // Implicit multiplication
          const right = parsePower();
          left = {
            dimension: left.dimension.multiply(right.dimension),
            siFactor: left.siFactor * right.siFactor,
            hExponent: left.hExponent + right.hExponent,
            aExponent: left.aExponent + right.aExponent
          };
        } else {
          break;
        }
      }

      return left;
    };

    const result = parseExpr();
    if (pos < tokens.length) {
      throw new UnitParseError(`Unexpected trailing token "${tokens[pos].value}"`, originalExpr);
    }
    return result;
  }

  /**
   * Helper to decompose unit with standard SI prefixes.
   * @private
   */
  _tryMatchPrefix(symbol) {
    const prefixes = [
      { name: 'Y', factor: 1e24 },
      { name: 'Z', factor: 1e21 },
      { name: 'E', factor: 1e18 },
      { name: 'P', factor: 1e15 },
      { name: 'T', factor: 1e12 },
      { name: 'G', factor: 1e9 },
      { name: 'M', factor: 1e6 },
      { name: 'k', factor: 1e3 },
      { name: 'h', factor: 1e2 },
      { name: 'da', factor: 10.0 },
      { name: 'd', factor: 1e-1 },
      { name: 'c', factor: 1e-2 },
      { name: 'm', factor: 1e-3 },
      { name: 'u', factor: 1e-6 },
      { name: 'μ', factor: 1e-6 },
      { name: 'n', factor: 1e-9 },
      { name: 'p', factor: 1e-12 },
      { name: 'f', factor: 1e-15 },
      { name: 'a', factor: 1e-18 },
      { name: 'z', factor: 1e-21 },
      { name: 'y', factor: 1e-24 }
    ];

    for (const p of prefixes) {
      if (symbol.startsWith(p.name) && symbol.length > p.name.length) {
        const baseSymbol = symbol.slice(p.name.length);
        if (this._units.has(baseSymbol) || this._aliases.has(baseSymbol)) {
          const baseUnit = this.get(baseSymbol);
          return {
            dimension: baseUnit.dimension,
            siFactor: baseUnit.siFactor * p.factor,
            hExponent: baseUnit.hExponent,
            aExponent: baseUnit.aExponent
          };
        }
      }
    }
    return null;
  }
}

/** Global default UnitRegistry singleton */
export const defaultUnitRegistry = new UnitRegistry();

// ============================================================================
// 6. STRONGLY-TYPED SCALAR PHYSICAL QUANTITY (GUARDED QUANTITY)
// ============================================================================

/**
 * DimensionalQuantity wraps a numerical magnitude with immutable dimensional typing,
 * strict arithmetic assertion guards, and cosmological h/a transformations.
 */
export class DimensionalQuantity {
  /**
   * @param {number} value - Numerical scalar magnitude
   * @param {string|UnitDefinition} [unit='1'] - Unit symbol or definition
   * @param {UnitRegistry} [registry=defaultUnitRegistry] - Associated registry
   */
  constructor(value, unit = '1', registry = defaultUnitRegistry) {
    if (typeof value !== 'number' || isNaN(value)) {
      throw new TypeError(`DimensionalQuantity magnitude must be a valid number, received ${value}`);
    }
    this.value = value;
    this.registry = registry;

    if (unit instanceof UnitDefinition) {
      this.unitDef = unit;
      this.unit = unit.symbol;
    } else {
      this.unit = unit;
      this.unitDef = registry.get(unit);
    }

    this.dimension = this.unitDef.dimension;
    Object.freeze(this);
  }

  /**
   * Checks if quantity is dimensionless.
   * @returns {boolean}
   */
  isDimensionless() {
    return this.dimension.isDimensionless();
  }

  /**
   * Checks if quantity is finite.
   * @returns {boolean}
   */
  isFinite() {
    return Number.isFinite(this.value);
  }

  /**
   * Returns whether the magnitude is exactly zero.
   * @returns {boolean}
   */
  isZero() {
    return this.value === 0.0;
  }

  /**
   * Returns true if magnitude > 0.
   * @returns {boolean}
   */
  isPositive() {
    return this.value > 0.0;
  }

  /**
   * Returns true if magnitude < 0.
   * @returns {boolean}
   */
  isNegative() {
    return this.value < 0.0;
  }

  /**
   * Converts this quantity to target unit.
   * Throws DimensionalityMismatchError if dimensions are physically incompatible.
   * 
   * @param {string|UnitDefinition} targetUnit
   * @param {object} [context={}] Cosmological parameters (e.g. { h: 0.746, a: 1.0 })
   * @returns {DimensionalQuantity} Converted quantity
   */
  to(targetUnit, context = {}) {
    const convertedVal = this.registry.convert(this.value, this.unitDef, targetUnit, context);
    const dstDef = targetUnit instanceof UnitDefinition ? targetUnit : this.registry.get(targetUnit);
    return new DimensionalQuantity(convertedVal, dstDef, this.registry);
  }

  /**
   * Converts this quantity to canonical SI base units.
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  toSI(context = {}) {
    const h = typeof context.h === 'number' ? context.h : 0.746;
    const a = typeof context.a === 'number' ? context.a : 1.0;
    const siVal = this.value * this.unitDef.getEffectiveMultiplierToSI({ h, a });

    // Look for matching registered base SI unit
    for (const u of this.registry._units.values()) {
      if (u.dimension.equals(this.dimension, true) && u.hExponent === 0 && u.aExponent === 0 && u.siFactor === 1.0) {
        return new DimensionalQuantity(siVal, u, this.registry);
      }
    }

    // Default to composite SI definition
    const compositeDef = new UnitDefinition({
      symbol: `[SI: ${this.dimension.toString()}]`,
      dimension: this.dimension,
      siFactor: 1.0,
      hExponent: 0,
      aExponent: 0
    });
    return new DimensionalQuantity(siVal, compositeDef, this.registry);
  }

  /**
   * Adds another quantity (strictly verified for dimensional consistency).
   * Throws DimensionalityMismatchError upon dimension mismatch (e.g. Position + Velocity).
   * 
   * @param {DimensionalQuantity} other
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  add(other, context = {}) {
    if (!(other instanceof DimensionalQuantity)) {
      throw new TypeError(`Cannot add non-DimensionalQuantity (${typeof other}) to DimensionalQuantity`);
    }

    if (!this.dimension.equals(other.dimension, false)) {
      throw new DimensionalityMismatchError(
        `Illegal addition: Cannot add incompatible dimensions: "${this.unit}" (${this.dimension}) + "${other.unit}" (${other.dimension})`,
        {
          sourceDimension: this.dimension,
          targetDimension: other.dimension,
          sourceUnit: this.unit,
          targetUnit: other.unit,
          operation: 'add',
          context
        }
      );
    }

    const otherInThis = other.to(this.unitDef, context);
    return new DimensionalQuantity(this.value + otherInThis.value, this.unitDef, this.registry);
  }

  /**
   * Subtracts another quantity (strictly verified for dimensional consistency).
   * Throws DimensionalityMismatchError upon dimension mismatch.
   * 
   * @param {DimensionalQuantity} other
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  subtract(other, context = {}) {
    if (!(other instanceof DimensionalQuantity)) {
      throw new TypeError(`Cannot subtract non-DimensionalQuantity (${typeof other}) from DimensionalQuantity`);
    }

    if (!this.dimension.equals(other.dimension, false)) {
      throw new DimensionalityMismatchError(
        `Illegal subtraction: Cannot subtract incompatible dimensions: "${this.unit}" (${this.dimension}) - "${other.unit}" (${other.dimension})`,
        {
          sourceDimension: this.dimension,
          targetDimension: other.dimension,
          sourceUnit: this.unit,
          targetUnit: other.unit,
          operation: 'subtract',
          context
        }
      );
    }

    const otherInThis = other.to(this.unitDef, context);
    return new DimensionalQuantity(this.value - otherInThis.value, this.unitDef, this.registry);
  }

  /**
   * Multiplies by a scalar number or another DimensionalQuantity.
   * Computes product dimension tensor and combined SI multiplier.
   * 
   * @param {number|DimensionalQuantity} factor
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  multiply(factor, context = {}) {
    if (typeof factor === 'number') {
      return new DimensionalQuantity(this.value * factor, this.unitDef, this.registry);
    }

    if (factor instanceof DimensionalQuantity) {
      const newDim = this.dimension.multiply(factor.dimension);
      const newHExp = this.unitDef.hExponent + factor.unitDef.hExponent;
      const newAExp = this.unitDef.aExponent + factor.unitDef.aExponent;

      const compositeSymbol = `(${this.unit}*${factor.unit})`;
      const compositeDef = new UnitDefinition({
        symbol: compositeSymbol,
        dimension: newDim,
        siFactor: this.unitDef.siFactor * factor.unitDef.siFactor,
        hExponent: newHExp,
        aExponent: newAExp
      });

      return new DimensionalQuantity(this.value * factor.value, compositeDef, this.registry);
    }

    throw new TypeError(`DimensionalQuantity.multiply factor must be number or DimensionalQuantity, received ${typeof factor}`);
  }

  /**
   * Divides by a scalar number or another DimensionalQuantity.
   * 
   * @param {number|DimensionalQuantity} divisor
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  divide(divisor, context = {}) {
    if (typeof divisor === 'number') {
      if (divisor === 0.0) throw new RangeError('Division by zero in DimensionalQuantity.divide');
      return new DimensionalQuantity(this.value / divisor, this.unitDef, this.registry);
    }

    if (divisor instanceof DimensionalQuantity) {
      if (divisor.value === 0.0) throw new RangeError('Division by zero in DimensionalQuantity.divide');

      const newDim = this.dimension.divide(divisor.dimension);
      const newHExp = this.unitDef.hExponent - divisor.unitDef.hExponent;
      const newAExp = this.unitDef.aExponent - divisor.unitDef.aExponent;

      const compositeSymbol = `(${this.unit}/${divisor.unit})`;
      const compositeDef = new UnitDefinition({
        symbol: compositeSymbol,
        dimension: newDim,
        siFactor: this.unitDef.siFactor / divisor.unitDef.siFactor,
        hExponent: newHExp,
        aExponent: newAExp
      });

      return new DimensionalQuantity(this.value / divisor.value, compositeDef, this.registry);
    }

    throw new TypeError(`DimensionalQuantity.divide divisor must be number or DimensionalQuantity, received ${typeof divisor}`);
  }

  /**
   * Raises quantity to a real power.
   * @param {number} p
   * @returns {DimensionalQuantity}
   */
  pow(p) {
    const power = Number(p);
    if (isNaN(power)) throw new TypeError(`DimensionalQuantity.pow power must be a valid number, received ${p}`);

    const newDim = this.dimension.pow(power);
    const newHExp = this.unitDef.hExponent * power;
    const newAExp = this.unitDef.aExponent * power;

    const powDef = new UnitDefinition({
      symbol: `(${this.unit})^${power}`,
      dimension: newDim,
      siFactor: Math.pow(this.unitDef.siFactor, power),
      hExponent: newHExp,
      aExponent: newAExp
    });

    return new DimensionalQuantity(Math.pow(this.value, power), powDef, this.registry);
  }

  /**
   * Computes square root of quantity.
   * @returns {DimensionalQuantity}
   */
  sqrt() {
    if (this.value < 0.0) {
      throw new RangeError(`Cannot compute square root of negative quantity (${this.value} ${this.unit})`);
    }
    return this.pow(0.5);
  }

  /**
   * Computes absolute value.
   * @returns {DimensionalQuantity}
   */
  abs() {
    return new DimensionalQuantity(Math.abs(this.value), this.unitDef, this.registry);
  }

  /**
   * Negates the quantity.
   * @returns {DimensionalQuantity}
   */
  negate() {
    return new DimensionalQuantity(-this.value, this.unitDef, this.registry);
  }

  /**
   * Compares with another quantity.
   * Returns -1 if this < other, 0 if equal, +1 if this > other.
   * Throws DimensionalityMismatchError if dimensions mismatch.
   * 
   * @param {DimensionalQuantity} other
   * @param {object} [context={}]
   * @returns {number}
   */
  compareTo(other, context = {}) {
    if (!(other instanceof DimensionalQuantity)) {
      throw new TypeError(`Cannot compare DimensionalQuantity with ${typeof other}`);
    }
    if (!this.dimension.equals(other.dimension, false)) {
      throw new DimensionalityMismatchError(
        `Cannot compare incompatible dimensions: "${this.unit}" (${this.dimension}) vs "${other.unit}" (${other.dimension})`,
        { sourceDimension: this.dimension, targetDimension: other.dimension, operation: 'compareTo' }
      );
    }
    const otherVal = other.to(this.unitDef, context).value;
    if (this.value < otherVal) return -1;
    if (this.value > otherVal) return 1;
    return 0;
  }

  /**
   * Tests numerical equality within absolute / relative tolerance.
   * 
   * @param {DimensionalQuantity} other
   * @param {number} [relTol=1e-7]
   * @param {number} [absTol=1e-12]
   * @param {object} [context={}]
   * @returns {boolean}
   */
  equals(other, relTol = 1e-7, absTol = 1e-12, context = {}) {
    if (!(other instanceof DimensionalQuantity)) return false;
    if (!this.dimension.equals(other.dimension, false)) return false;
    const otherInThis = other.to(this.unitDef, context).value;
    const diff = Math.abs(this.value - otherInThis);
    const maxVal = Math.max(Math.abs(this.value), Math.abs(otherInThis), 1e-15);
    return diff <= absTol || (diff / maxVal) <= relTol;
  }

  /**
   * Formats for scientific publication display.
   * @param {number} [precision=6]
   * @returns {string}
   */
  format(precision = 6) {
    const absV = Math.abs(this.value);
    const valStr = (absV >= 1e6 || (absV < 1e-3 && absV !== 0.0))
      ? this.value.toExponential(precision - 1)
      : parseFloat(this.value.toPrecision(precision)).toString();
    return `${valStr} ${this.unit}`;
  }

  toString() {
    return `${this.value} ${this.unit}`;
  }

  valueOf() {
    return this.value;
  }

  // --------------------------------------------------------------------------
  // STATISTICAL REDUCTIONS OVER HOMOGENEOUS ARRAYS
  // --------------------------------------------------------------------------

  /**
   * Sums an array of homogeneous DimensionalQuantity elements.
   * @param {DimensionalQuantity[]} items
   * @param {string|UnitDefinition} [targetUnit]
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  static sum(items, targetUnit, context = {}) {
    if (!Array.isArray(items) || items.length === 0) {
      throw new RangeError('DimensionalQuantity.sum requires a non-empty array');
    }
    const baseUnit = targetUnit ?? items[0].unitDef;
    let acc = new DimensionalQuantity(0.0, baseUnit, items[0].registry);
    for (const item of items) {
      acc = acc.add(item, context);
    }
    return acc;
  }

  /**
   * Computes the arithmetic mean of an array of homogeneous quantities.
   * @param {DimensionalQuantity[]} items
   * @param {string|UnitDefinition} [targetUnit]
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  static mean(items, targetUnit, context = {}) {
    const total = DimensionalQuantity.sum(items, targetUnit, context);
    return total.divide(items.length);
  }

  /**
   * Computes root-mean-square (RMS) of homogeneous quantities.
   * @param {DimensionalQuantity[]} items
   * @param {string|UnitDefinition} [targetUnit]
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  static rms(items, targetUnit, context = {}) {
    if (!Array.isArray(items) || items.length === 0) {
      throw new RangeError('DimensionalQuantity.rms requires a non-empty array');
    }
    const baseUnit = targetUnit ?? items[0].unitDef;
    let sumSq = 0.0;
    for (const item of items) {
      const v = item.to(baseUnit, context).value;
      sumSq += v * v;
    }
    const rmsVal = Math.sqrt(sumSq / items.length);
    return new DimensionalQuantity(rmsVal, baseUnit, items[0].registry);
  }

  /**
   * Finds minimum quantity.
   * @param {DimensionalQuantity[]} items
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  static min(items, context = {}) {
    if (!Array.isArray(items) || items.length === 0) {
      throw new RangeError('DimensionalQuantity.min requires a non-empty array');
    }
    let minItem = items[0];
    for (let i = 1; i < items.length; i++) {
      if (items[i].compareTo(minItem, context) < 0) {
        minItem = items[i];
      }
    }
    return minItem;
  }

  /**
   * Finds maximum quantity.
   * @param {DimensionalQuantity[]} items
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  static max(items, context = {}) {
    if (!Array.isArray(items) || items.length === 0) {
      throw new RangeError('DimensionalQuantity.max requires a non-empty array');
    }
    let maxItem = items[0];
    for (let i = 1; i < items.length; i++) {
      if (items[i].compareTo(maxItem, context) > 0) {
        maxItem = items[i];
      }
    }
    return maxItem;
  }
}

// ============================================================================
// 7. 3D VECTOR & 3X3 TENSOR DIMENSIONAL QUANTITIES
// ============================================================================

/**
 * 3D Physical Vector Quantity with strict dimensional integrity on components,
 * dot products, cross products, coordinate frame transformations, and norm operations.
 */
export class VectorDimensionalQuantity {
  /**
   * @param {number|number[]|{x:number, y:number, z:number}} x - X component or array/object
   * @param {number} [y=0] - Y component
   * @param {number} [z=0] - Z component
   * @param {string|UnitDefinition} [unit='1'] - Vector physical unit
   * @param {UnitRegistry} [registry=defaultUnitRegistry]
   */
  constructor(x, y = 0, z = 0, unit = '1', registry = defaultUnitRegistry) {
    let vx = 0;
    let vy = 0;
    let vz = 0;
    let vUnit = unit;

    if (Array.isArray(x)) {
      vx = Number(x[0]) || 0;
      vy = Number(x[1]) || 0;
      vz = Number(x[2]) || 0;
      if (typeof y === 'string' || y instanceof UnitDefinition) {
        vUnit = y;
      }
    } else if (x && typeof x === 'object') {
      vx = Number(x.x ?? x.sgx ?? x[0]) || 0;
      vy = Number(x.y ?? x.sgy ?? x[1]) || 0;
      vz = Number(x.z ?? x.sgz ?? x[2]) || 0;
      if (typeof y === 'string' || y instanceof UnitDefinition) {
        vUnit = y;
      }
    } else {
      vx = Number(x) || 0;
      vy = Number(y) || 0;
      vz = Number(z) || 0;
    }

    if (isNaN(vx) || isNaN(vy) || isNaN(vz)) {
      throw new TypeError(`Vector components must be valid numbers: [${vx}, ${vy}, ${vz}]`);
    }

    this.x = vx;
    this.y = vy;
    this.z = vz;
    this.registry = registry;

    if (vUnit instanceof UnitDefinition) {
      this.unitDef = vUnit;
      this.unit = vUnit.symbol;
    } else {
      this.unit = vUnit;
      this.unitDef = registry.get(vUnit);
    }

    this.dimension = this.unitDef.dimension;
    Object.freeze(this);
  }

  toArray() {
    return [this.x, this.y, this.z];
  }

  toObject() {
    return { x: this.x, y: this.y, z: this.z, unit: this.unit };
  }

  /**
   * Computes Euclidean norm / magnitude as a scalar DimensionalQuantity.
   * @returns {DimensionalQuantity}
   */
  magnitude() {
    const mag = Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    return new DimensionalQuantity(mag, this.unitDef, this.registry);
  }

  /**
   * Converts vector components to target unit.
   * @param {string|UnitDefinition} targetUnit
   * @param {object} [context={}]
   * @returns {VectorDimensionalQuantity}
   */
  to(targetUnit, context = {}) {
    const convX = this.registry.convert(this.x, this.unitDef, targetUnit, context);
    const convY = this.registry.convert(this.y, this.unitDef, targetUnit, context);
    const convZ = this.registry.convert(this.z, this.unitDef, targetUnit, context);
    const dstDef = targetUnit instanceof UnitDefinition ? targetUnit : this.registry.get(targetUnit);
    return new VectorDimensionalQuantity(convX, convY, convZ, dstDef, this.registry);
  }

  /**
   * Adds another 3D vector quantity (strictly dimensionally checked).
   * Throws DimensionalityMismatchError if dimensions do not match.
   * 
   * @param {VectorDimensionalQuantity} other
   * @param {object} [context={}]
   * @returns {VectorDimensionalQuantity}
   */
  add(other, context = {}) {
    if (!(other instanceof VectorDimensionalQuantity)) {
      throw new TypeError(`Cannot add non-VectorDimensionalQuantity (${typeof other}) to VectorDimensionalQuantity`);
    }

    if (!this.dimension.equals(other.dimension, false)) {
      throw new DimensionalityMismatchError(
        `Illegal vector addition: Cannot add incompatible dimensions: "${this.unit}" (${this.dimension}) + "${other.unit}" (${other.dimension})`,
        {
          sourceDimension: this.dimension,
          targetDimension: other.dimension,
          sourceUnit: this.unit,
          targetUnit: other.unit,
          operation: 'vector_add',
          context
        }
      );
    }

    const otherInThis = other.to(this.unitDef, context);
    return new VectorDimensionalQuantity(
      this.x + otherInThis.x,
      this.y + otherInThis.y,
      this.z + otherInThis.z,
      this.unitDef,
      this.registry
    );
  }

  /**
   * Subtracts another 3D vector quantity (strictly dimensionally checked).
   * @param {VectorDimensionalQuantity} other
   * @param {object} [context={}]
   * @returns {VectorDimensionalQuantity}
   */
  subtract(other, context = {}) {
    if (!(other instanceof VectorDimensionalQuantity)) {
      throw new TypeError(`Cannot subtract non-VectorDimensionalQuantity (${typeof other}) from VectorDimensionalQuantity`);
    }

    if (!this.dimension.equals(other.dimension, false)) {
      throw new DimensionalityMismatchError(
        `Illegal vector subtraction: Cannot subtract incompatible dimensions: "${this.unit}" (${this.dimension}) - "${other.unit}" (${other.dimension})`,
        {
          sourceDimension: this.dimension,
          targetDimension: other.dimension,
          sourceUnit: this.unit,
          targetUnit: other.unit,
          operation: 'vector_subtract',
          context
        }
      );
    }

    const otherInThis = other.to(this.unitDef, context);
    return new VectorDimensionalQuantity(
      this.x - otherInThis.x,
      this.y - otherInThis.y,
      this.z - otherInThis.z,
      this.unitDef,
      this.registry
    );
  }

  /**
   * Multiplies vector by a scalar or scalar DimensionalQuantity.
   * @param {number|DimensionalQuantity} scalar
   * @returns {VectorDimensionalQuantity}
   */
  multiply(scalar) {
    if (typeof scalar === 'number') {
      return new VectorDimensionalQuantity(this.x * scalar, this.y * scalar, this.z * scalar, this.unitDef, this.registry);
    }
    if (scalar instanceof DimensionalQuantity) {
      const newDim = this.dimension.multiply(scalar.dimension);
      const compositeDef = new UnitDefinition({
        symbol: `(${this.unit}*${scalar.unit})`,
        dimension: newDim,
        siFactor: this.unitDef.siFactor * scalar.unitDef.siFactor,
        hExponent: this.unitDef.hExponent + scalar.unitDef.hExponent,
        aExponent: this.unitDef.aExponent + scalar.unitDef.aExponent
      });
      return new VectorDimensionalQuantity(
        this.x * scalar.value,
        this.y * scalar.value,
        this.z * scalar.value,
        compositeDef,
        this.registry
      );
    }
    throw new TypeError(`Vector multiply scalar must be number or DimensionalQuantity, received ${typeof scalar}`);
  }

  /**
   * Computes the vector dot product (returns scalar DimensionalQuantity with dim1 * dim2).
   * @param {VectorDimensionalQuantity} other
   * @param {object} [context={}]
   * @returns {DimensionalQuantity}
   */
  dot(other, context = {}) {
    if (!(other instanceof VectorDimensionalQuantity)) {
      throw new TypeError(`Vector dot product requires VectorDimensionalQuantity, received ${typeof other}`);
    }

    const newDim = this.dimension.multiply(other.dimension);
    const compositeDef = new UnitDefinition({
      symbol: `(${this.unit}*${other.unit})`,
      dimension: newDim,
      siFactor: this.unitDef.siFactor * other.unitDef.siFactor,
      hExponent: this.unitDef.hExponent + other.unitDef.hExponent,
      aExponent: this.unitDef.aExponent + other.unitDef.aExponent
    });

    const val = (this.x * other.x) + (this.y * other.y) + (this.z * other.z);
    return new DimensionalQuantity(val, compositeDef, this.registry);
  }

  /**
   * Computes the vector cross product (returns 3D vector with dim1 * dim2).
   * @param {VectorDimensionalQuantity} other
   * @param {object} [context={}]
   * @returns {VectorDimensionalQuantity}
   */
  cross(other, context = {}) {
    if (!(other instanceof VectorDimensionalQuantity)) {
      throw new TypeError(`Vector cross product requires VectorDimensionalQuantity, received ${typeof other}`);
    }

    const newDim = this.dimension.multiply(other.dimension);
    const compositeDef = new UnitDefinition({
      symbol: `(${this.unit}*${other.unit})`,
      dimension: newDim,
      siFactor: this.unitDef.siFactor * other.unitDef.siFactor,
      hExponent: this.unitDef.hExponent + other.unitDef.hExponent,
      aExponent: this.unitDef.aExponent + other.unitDef.aExponent
    });

    const cx = this.y * other.z - this.z * other.y;
    const cy = this.z * other.x - this.x * other.z;
    const cz = this.x * other.y - this.y * other.x;

    return new VectorDimensionalQuantity(cx, cy, cz, compositeDef, this.registry);
  }

  /**
   * Returns dimensionless unit direction vector.
   * @returns {VectorDimensionalQuantity}
   */
  normalize() {
    const mag = Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    if (mag === 0.0) throw new RangeError('Cannot normalize zero vector');
    return new VectorDimensionalQuantity(
      this.x / mag,
      this.y / mag,
      this.z / mag,
      '1',
      this.registry
    );
  }

  toString() {
    return `[${this.x}, ${this.y}, ${this.z}] ${this.unit}`;
  }
}

/**
 * 3x3 Second-Rank Physical Tensor Quantity (for strain rates, velocity shear, tidal field, quadrupole).
 */
export class TensorDimensionalQuantity {
  /**
   * @param {number[][]} matrix - 3x3 numerical array
   * @param {string|UnitDefinition} [unit='1']
   * @param {UnitRegistry} [registry=defaultUnitRegistry]
   */
  constructor(matrix, unit = '1', registry = defaultUnitRegistry) {
    if (!Array.isArray(matrix) || matrix.length !== 3 || !matrix.every(r => Array.isArray(r) && r.length === 3)) {
      throw new TypeError('TensorDimensionalQuantity requires a 3x3 nested number array');
    }

    this.m = [
      [Number(matrix[0][0]) || 0, Number(matrix[0][1]) || 0, Number(matrix[0][2]) || 0],
      [Number(matrix[1][0]) || 0, Number(matrix[1][1]) || 0, Number(matrix[1][2]) || 0],
      [Number(matrix[2][0]) || 0, Number(matrix[2][1]) || 0, Number(matrix[2][2]) || 0]
    ];

    this.registry = registry;

    if (unit instanceof UnitDefinition) {
      this.unitDef = unit;
      this.unit = unit.symbol;
    } else {
      this.unit = unit;
      this.unitDef = registry.get(unit);
    }

    this.dimension = this.unitDef.dimension;
    Object.freeze(this);
  }

  /**
   * Computes trace of the tensor as scalar DimensionalQuantity.
   * @returns {DimensionalQuantity}
   */
  trace() {
    const tr = this.m[0][0] + this.m[1][1] + this.m[2][2];
    return new DimensionalQuantity(tr, this.unitDef, this.registry);
  }

  /**
   * Computes Frobenius norm of the tensor.
   * @returns {DimensionalQuantity}
   */
  frobeniusNorm() {
    let sumSq = 0.0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        sumSq += this.m[i][j] * this.m[i][j];
      }
    }
    return new DimensionalQuantity(Math.sqrt(sumSq), this.unitDef, this.registry);
  }

  /**
   * Computes exact Jacobi eigenvalue diagonalization for symmetric 3x3 tensors.
   * Returns eigenvalues sorted lambda_1 >= lambda_2 >= lambda_3 with physical dimensions preserved.
   * 
   * @param {number} [maxSweeps=50]
   * @param {number} [tolerance=1e-15]
   * @returns {{ eigenvalues: DimensionalQuantity[], eigenvectors: number[][] }}
   */
  diagonalizeSymmetric(maxSweeps = 50, tolerance = 1e-15) {
    const a = [
      [this.m[0][0], this.m[0][1], this.m[0][2]],
      [this.m[1][0], this.m[1][1], this.m[1][2]],
      [this.m[2][0], this.m[2][1], this.m[2][2]]
    ];

    const v = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ];

    for (let sweep = 0; sweep < maxSweeps; sweep++) {
      let maxOff = 0.0;
      for (let i = 0; i < 3; i++) {
        for (let j = i + 1; j < 3; j++) {
          maxOff = Math.max(maxOff, Math.abs(a[i][j]));
        }
      }
      if (maxOff < tolerance) break;

      for (let p = 0; p < 2; p++) {
        for (let q = p + 1; q < 3; q++) {
          const apq = a[p][q];
          if (Math.abs(apq) < tolerance) continue;

          const app = a[p][p];
          const aqq = a[q][q];
          const theta = (aqq - app) / (2.0 * apq);
          const t = theta >= 0
            ? 1.0 / (theta + Math.sqrt(theta * theta + 1.0))
            : -1.0 / (-theta + Math.sqrt(theta * theta + 1.0));
          const c = 1.0 / Math.sqrt(t * t + 1.0);
          const s = t * c;
          const tau = s / (1.0 + c);

          a[p][p] = app - t * apq;
          a[q][q] = aqq + t * apq;
          a[p][q] = 0.0;
          a[q][p] = 0.0;

          for (let r = 0; r < 3; r++) {
            if (r !== p && r !== q) {
              const arp = a[r][p];
              const arq = a[r][q];
              a[r][p] = arp - s * (arq + tau * arp);
              a[p][r] = a[r][p];
              a[r][q] = arq + s * (arp - tau * arq);
              a[q][r] = a[r][q];
            }
          }

          for (let r = 0; r < 3; r++) {
            const vrp = v[r][p];
            const vrq = v[r][q];
            v[r][p] = vrp - s * (vrq + tau * vrp);
            v[r][q] = vrq + s * (vrp - tau * vrq);
          }
        }
      }
    }

    const indices = [0, 1, 2];
    indices.sort((i, j) => a[j][j] - a[i][i]);

    const sortedVals = indices.map(i => new DimensionalQuantity(a[i][i], this.unitDef, this.registry));
    const sortedVecs = indices.map(col => [v[0][col], v[1][col], v[2][col]]);

    return {
      eigenvalues: sortedVals,
      eigenvectors: sortedVecs
    };
  }
}

// ============================================================================
// 8. COMPILE-TIME & RUNTIME SYMBOLIC DIMENSIONAL EXPRESSION CHECKER
// ============================================================================

/**
 * DimensionalExpressionChecker analyzes mathematical expression strings and ASTs
 * to verify dimensional homogeneity at compile-time or static initialization.
 */
export class DimensionalExpressionChecker {
  /**
   * @param {Map<string, DimensionTensor|string>|Object} [symbolTable={}]
   * @param {UnitRegistry} [registry=defaultUnitRegistry]
   */
  constructor(symbolTable = {}, registry = defaultUnitRegistry) {
    this.registry = registry;
    this.symbolTable = new Map();

    if (symbolTable instanceof Map) {
      for (const [k, v] of symbolTable.entries()) {
        this.setSymbol(k, v);
      }
    } else if (typeof symbolTable === 'object') {
      for (const [k, v] of Object.entries(symbolTable)) {
        this.setSymbol(k, v);
      }
    }
  }

  /**
   * Sets or updates a symbol's dimension type.
   * @param {string} symbol
   * @param {DimensionTensor|string} dimensionOrUnit
   */
  setSymbol(symbol, dimensionOrUnit) {
    if (dimensionOrUnit instanceof DimensionTensor) {
      this.symbolTable.set(symbol, dimensionOrUnit);
    } else if (typeof dimensionOrUnit === 'string') {
      const def = this.registry.get(dimensionOrUnit);
      this.symbolTable.set(symbol, def.dimension);
    } else {
      throw new TypeError(`Invalid dimension specifier for symbol "${symbol}"`);
    }
  }

  /**
   * Validates dimensional consistency of a mathematical string expression
   * (e.g. "0.5 * rho * v^2 + P", "G * M / r^2", "v_bulk + H0 * r").
   * 
   * Throws DimensionalityMismatchError if addition/subtraction operands have incompatible dimensions.
   * 
   * @param {string} expression
   * @returns {DimensionTensor} Resulting physical dimension of the expression
   */
  check(expression) {
    if (!expression || typeof expression !== 'string') {
      throw new TypeError('Expression must be a non-empty string');
    }

    const tokens = this._tokenize(expression);
    let pos = 0;

    const parsePrimary = () => {
      if (pos >= tokens.length) {
        throw new UnitParseError('Unexpected end of expression', expression);
      }
      const tok = tokens[pos];

      if (tok.type === 'LPAREN') {
        pos++;
        const inner = parseAddSub();
        if (pos >= tokens.length || tokens[pos].type !== 'RPAREN') {
          throw new UnitParseError('Mismatched parentheses in expression', expression);
        }
        pos++;
        return inner;
      }

      if (tok.type === 'NUMBER') {
        pos++;
        return PHYSICAL_DIMENSIONS.DIMENSIONLESS;
      }

      if (tok.type === 'IDENT') {
        pos++;
        const sym = tok.value;
        if (this.symbolTable.has(sym)) {
          return this.symbolTable.get(sym);
        }
        // Check known constants
        if (sym === 'c' || sym === 'C') return PHYSICAL_DIMENSIONS.VELOCITY;
        if (sym === 'G') return PHYSICAL_DIMENSIONS.GRAVITATIONAL_CONSTANT;
        if (sym === 'H0' || sym === 'H') return PHYSICAL_DIMENSIONS.EXPANSION_RATE;
        if (sym === 'delta' || sym === 'f' || sym === 'pi' || sym === 'PI') return PHYSICAL_DIMENSIONS.DIMENSIONLESS;

        throw new UnitParseError(`Unbound symbol "${sym}" in expression dimensional analysis`, expression);
      }

      throw new UnitParseError(`Unexpected token "${tok.value}"`, expression);
    };

    const parsePower = () => {
      let left = parsePrimary();
      while (pos < tokens.length && tokens[pos].type === 'OP' && tokens[pos].value === '^') {
        pos++;
        if (pos >= tokens.length || tokens[pos].type !== 'NUMBER') {
          throw new UnitParseError('Expected numeric exponent after "^"', expression);
        }
        const p = tokens[pos].value;
        pos++;
        left = left.pow(p);
      }
      return left;
    };

    const parseMulDiv = () => {
      let left = parsePower();
      while (pos < tokens.length && tokens[pos].type === 'OP' && (tokens[pos].value === '*' || tokens[pos].value === '/')) {
        const op = tokens[pos].value;
        pos++;
        const right = parsePower();
        if (op === '*') {
          left = left.multiply(right);
        } else {
          left = left.divide(right);
        }
      }
      return left;
    };

    const parseAddSub = () => {
      let left = parseMulDiv();
      while (pos < tokens.length && tokens[pos].type === 'OP' && (tokens[pos].value === '+' || tokens[pos].value === '-')) {
        const op = tokens[pos].value;
        pos++;
        const right = parseMulDiv();

        // STRICT CHECK: Addition and Subtraction require exact dimensional match!
        if (!left.equals(right, false)) {
          throw new DimensionalityMismatchError(
            `Illegal expression: Incompatible dimensions in "${op}" operation: ${left} ${op} ${right}`,
            { sourceDimension: left, targetDimension: right, operation: `expression_${op}` }
          );
        }
      }
      return left;
    };

    const result = parseAddSub();
    if (pos < tokens.length) {
      throw new UnitParseError(`Trailing characters after expression: "${tokens[pos].value}"`, expression);
    }
    return result;
  }

  _tokenize(expr) {
    const tokens = [];
    let i = 0;
    const len = expr.length;

    while (i < len) {
      const c = expr[i];
      if (/\s/.test(c)) {
        i++;
        continue;
      }
      if (c === '+' || c === '-' || c === '*' || c === '/' || c === '^') {
        tokens.push({ type: 'OP', value: c });
        i++;
      } else if (c === '(') {
        tokens.push({ type: 'LPAREN', value: '(' });
        i++;
      } else if (c === ')') {
        tokens.push({ type: 'RPAREN', value: ')' });
        i++;
      } else if (/\d/.test(c) || (c === '.' && i + 1 < len && /\d/.test(expr[i + 1]))) {
        let numStr = c;
        i++;
        while (i < len && /[\d.eE+-]/.test(expr[i])) {
          if ((expr[i] === '+' || expr[i] === '-') && !/[eE]/.test(expr[i - 1])) {
            break;
          }
          numStr += expr[i];
          i++;
        }
        tokens.push({ type: 'NUMBER', value: parseFloat(numStr) });
      } else if (/[a-zA-Z_]/.test(c)) {
        let ident = c;
        i++;
        while (i < len && /[a-zA-Z0-9_]/.test(expr[i])) {
          ident += expr[i];
          i++;
        }
        tokens.push({ type: 'IDENT', value: ident });
      } else {
        throw new UnitParseError(`Illegal character "${c}" in expression`, expr);
      }
    }
    return tokens;
  }
}

// ============================================================================
// 9. HIGHER-ORDER FUNCTION DIMENSIONAL CONTRACT GUARDS
// ============================================================================

/**
 * Creates a dimensionally guarded wrapper around a function.
 * Enforces input dimension signatures and return dimension assertions at runtime.
 * 
 * @param {Object} contract
 * @param {Array<DimensionTensor|string>} contract.inputs - Expected input dimensions
 * @param {DimensionTensor|string} contract.output - Expected output dimension
 * @param {string} [contract.functionName='guardedFunction'] - Descriptive name for errors
 * @param {Function} targetFn - The mathematical function to guard
 * @returns {Function}
 */
export function createGuardedFunction(contract, targetFn) {
  if (typeof targetFn !== 'function') {
    throw new TypeError('createGuardedFunction targetFn must be a callable function');
  }

  const inDims = (contract.inputs || []).map(d => d instanceof DimensionTensor ? d : defaultUnitRegistry.get(d).dimension);
  const outDim = contract.output instanceof DimensionTensor ? contract.output : defaultUnitRegistry.get(contract.output).dimension;
  const fnName = contract.functionName || targetFn.name || 'guardedFunction';

  return function (...args) {
    if (args.length < inDims.length) {
      throw new DimensionalityMismatchError(
        `${fnName}: Expected at least ${inDims.length} arguments, received ${args.length}`
      );
    }

    for (let i = 0; i < inDims.length; i++) {
      const arg = args[i];
      const expectedDim = inDims[i];

      if (arg instanceof DimensionalQuantity || arg instanceof VectorDimensionalQuantity) {
        if (!arg.dimension.equals(expectedDim, false)) {
          throw new DimensionalityMismatchError(
            `${fnName}: Argument ${i} dimension mismatch: expected ${expectedDim}, received ${arg.dimension} (${arg.unit})`,
            { sourceDimension: arg.dimension, targetDimension: expectedDim, operation: `${fnName}_arg_${i}` }
          );
        }
      } else if (typeof arg === 'number') {
        if (!expectedDim.isDimensionless()) {
          throw new DimensionalityMismatchError(
            `${fnName}: Argument ${i} is raw number but expected dimension ${expectedDim}`,
            { targetDimension: expectedDim, operation: `${fnName}_arg_${i}` }
          );
        }
      }
    }

    const result = targetFn.apply(this, args);

    if (result instanceof DimensionalQuantity || result instanceof VectorDimensionalQuantity) {
      if (!result.dimension.equals(outDim, false)) {
        throw new DimensionalityMismatchError(
          `${fnName}: Return value dimension mismatch: expected ${outDim}, received ${result.dimension} (${result.unit})`,
          { sourceDimension: result.dimension, targetDimension: outDim, operation: `${fnName}_return` }
        );
      }
    }

    return result;
  };
}

// ============================================================================
// 10. SCIENTIFIC INVARIANT & COSMOLOGICAL KINEMATIC GUARDS
// ============================================================================

/**
 * Asserts strict dimensional separation between spatial positions (Mpc/h) and velocities (km/s).
 * Enforces Rule 1: Never combine position and velocity directly without explicit dt interval.
 * 
 * @param {*} pos - Spatial position
 * @param {*} vel - Velocity vector
 * @param {DimensionalQuantity|number} [dt] - Time interval
 * @param {object} [context={}]
 * @returns {DimensionalQuantity|VectorDimensionalQuantity} Propagated position
 */
export function propagatePositionWithVelocity(pos, vel, dt, context = {}) {
  if (dt === undefined || dt === null) {
    throw new DimensionalityMismatchError(
      'Illegal kinematic combination: Cannot add or combine spatial position directly with velocity without an explicit time interval dt. ' +
      '(Rule 1 Invariant Violation: position_Mpc_h + velocity_km_s is forbidden).'
    );
  }

  const h = typeof context.h === 'number' ? context.h : 0.746;
  const timeQty = dt instanceof DimensionalQuantity
    ? dt.to('s', context)
    : new DimensionalQuantity(Number(dt), 's');

  if (vel instanceof VectorDimensionalQuantity) {
    const dx_m = vel.x * 1.0e3 * timeQty.value;
    const dy_m = vel.y * 1.0e3 * timeQty.value;
    const dz_m = vel.z * 1.0e3 * timeQty.value;

    const mpcPerMeter = 1.0 / CONSTANTS.METERS_PER_MPC;
    const dx_Mpc_h = dx_m * mpcPerMeter * h;
    const dy_Mpc_h = dy_m * mpcPerMeter * h;
    const dz_Mpc_h = dz_m * mpcPerMeter * h;

    const dispVec = new VectorDimensionalQuantity(dx_Mpc_h, dy_Mpc_h, dz_Mpc_h, 'Mpc/h');

    if (pos instanceof VectorDimensionalQuantity) {
      const posMpcH = pos.to('Mpc/h', context);
      return posMpcH.add(dispVec, context);
    }
    return dispVec;
  }

  if (vel instanceof DimensionalQuantity) {
    const vMps = vel.to('m/s', context).value;
    const distMeters = vMps * timeQty.value;
    const distMpcH = (distMeters / CONSTANTS.METERS_PER_MPC) * h;
    const dispQty = new DimensionalQuantity(distMpcH, 'Mpc/h');

    if (pos instanceof DimensionalQuantity) {
      const posMpcH = pos.to('Mpc/h', context);
      return posMpcH.add(dispQty, context);
    }
    return dispQty;
  }

  throw new TypeError('propagatePositionWithVelocity requires valid VectorDimensionalQuantity or DimensionalQuantity for vel');
}

/**
 * Verifies that the official CosmicFlows-4 x52.0 velocity scaling factor has been applied exactly once.
 * (Rule 3 Invariant Guard).
 * 
 * @param {boolean} currentAppliedState - Current state flag
 * @returns {boolean} New applied state (true)
 */
export function assertSingleCF4VelocityScaling(currentAppliedState) {
  if (currentAppliedState === true) {
    throw new InvariantViolationError(
      'Attempted double application of the official CosmicFlows-4 x52.0 velocity scale factor! ' +
      '(AGENTS.md Rule 3: Velocity scale factor must be applied exactly once).'
    );
  }
  return true;
}

/**
 * Verifies linear continuity diagnostic: div(v) approx -H0 * f * delta.
 * (Rule 4 Invariant Guard).
 * 
 * @param {number|DimensionalQuantity} divV - Velocity divergence in km/s / (Mpc/h)
 * @param {number|DimensionalQuantity} delta - Density contrast (dimensionless)
 * @param {object} [context={}] - Cosmological parameters
 * @param {number} [tolerance=0.35] - Fractional residual tolerance
 * @returns {{ satisfiesLinear: boolean, expectedDivV: number, actualDivV: number, relativeResidual: number }}
 */
export function evaluateLinearContinuityDiagnostic(divV, delta, context = {}, tolerance = 0.35) {
  const actualVal = divV instanceof DimensionalQuantity
    ? divV.to('km/s/(Mpc/h)', context).value
    : Number(divV);

  const deltaVal = delta instanceof DimensionalQuantity
    ? delta.to('1', context).value
    : Number(delta);

  const h = typeof context.h === 'number' ? context.h : 0.746;
  const H0 = typeof context.H0 === 'number' ? context.H0 : 74.6;
  const Omega_m = typeof context.Omega_m === 'number' ? context.Omega_m : 0.31;
  const f = Math.pow(Omega_m, 0.55);

  const continuityScale = (H0 * f) / h;
  const expectedVal = -continuityScale * deltaVal;

  const diff = Math.abs(actualVal - expectedVal);
  const norm = Math.max(1.0, Math.abs(expectedVal), Math.abs(actualVal));
  const relativeResidual = diff / norm;

  return {
    satisfiesLinear: relativeResidual <= tolerance,
    expectedDivV: expectedVal,
    actualDivV: actualVal,
    relativeResidual,
    continuityScale
  };
}

// ============================================================================
// 11. EXACT COSMOLOGICAL DISTANCE & VELOCITY / REDSHIFT CONVERSIONS
// ============================================================================

/**
 * Converts cosmological redshift z to pure Hubble recessional velocity v = c * z (km/s).
 * @param {number} z - Redshift (z >= 0)
 * @param {boolean} [relativistic=false] - If true, uses relativistic Doppler formula
 * @returns {DimensionalQuantity}
 */
export function redshiftToVelocity(z, relativistic = false) {
  const numZ = Number(z);
  if (isNaN(numZ) || numZ < 0.0) {
    throw new RangeError(`Redshift must be non-negative number, received ${z}`);
  }
  if (!relativistic) {
    const vKms = CONSTANTS.SPEED_OF_LIGHT_KMS * numZ;
    return new DimensionalQuantity(vKms, 'km/s');
  }
  const zp1Sq = Math.pow(1.0 + numZ, 2);
  const beta = (zp1Sq - 1.0) / (zp1Sq + 1.0);
  const vKms = CONSTANTS.SPEED_OF_LIGHT_KMS * beta;
  return new DimensionalQuantity(vKms, 'km/s');
}

/**
 * Converts recessional velocity (km/s) to cosmological redshift z.
 * @param {number|DimensionalQuantity} vel - Velocity
 * @param {boolean} [relativistic=false] - If true, inverts relativistic Doppler formula
 * @returns {DimensionalQuantity} Dimensionless redshift z
 */
export function velocityToRedshift(vel, relativistic = false) {
  const vKms = vel instanceof DimensionalQuantity ? vel.to('km/s').value : Number(vel);
  if (isNaN(vKms) || vKms < 0.0) {
    throw new RangeError(`Velocity must be non-negative number, received ${vel}`);
  }
  if (!relativistic) {
    const z = vKms / CONSTANTS.SPEED_OF_LIGHT_KMS;
    return new DimensionalQuantity(z, '1');
  }
  const beta = vKms / CONSTANTS.SPEED_OF_LIGHT_KMS;
  if (beta >= 1.0) {
    throw new RangeError(`Relativistic velocity cannot equal or exceed c (${vKms} km/s)`);
  }
  const z = Math.sqrt((1.0 + beta) / (1.0 - beta)) - 1.0;
  return new DimensionalQuantity(z, '1');
}

/**
 * Converts pure Hubble velocity (km/s) to local comoving distance in Mpc/h.
 * In linear Hubble regime (v < 30,000 km/s), d = v / 100 in Mpc/h.
 * 
 * @param {number|DimensionalQuantity} vel
 * @param {object} [context={}]
 * @returns {DimensionalQuantity} Distance in Mpc/h
 */
export function velocityToComovingDistance(vel, context = {}) {
  const vKms = vel instanceof DimensionalQuantity ? vel.to('km/s', context).value : Number(vel);
  if (vKms < 0.0) throw new RangeError(`Velocity must be non-negative, received ${vKms}`);

  const dMpcH = vKms / 100.0;
  return new DimensionalQuantity(dMpcH, 'Mpc/h');
}

/**
 * Converts comoving distance in Mpc/h to pure Hubble velocity in km/s.
 * @param {number|DimensionalQuantity} dist
 * @param {object} [context={}]
 * @returns {DimensionalQuantity} Velocity in km/s
 */
export function comovingDistanceToVelocity(dist, context = {}) {
  const dMpcH = dist instanceof DimensionalQuantity ? dist.to('Mpc/h', context).value : Number(dist);
  if (dMpcH < 0.0) throw new RangeError(`Distance must be non-negative, received ${dMpcH}`);

  const vKms = dMpcH * 100.0;
  return new DimensionalQuantity(vKms, 'km/s');
}
