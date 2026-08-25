/**
 * @fileoverview Dimensionally Consistent Unit Registry and Typed Physical Quantities for Cosmography.
 * 
 * Provides rigorous dimensional analysis, unit conversions, physical constants,
 * and $h$-parameter scaling for astrophysical quantities (Mpc/h, km/s, M_sun/h, km/s/(Mpc/h), rho_crit).
 * 
 * @module units/unit_registry
 */

/**
 * Fundamental Physical and Cosmological Constants (SI and Astrophysical Units)
 * CODATA 2022 / IAU 2015 / Planck 2018 / CF4
 */
export const CONSTANTS = Object.freeze({
  // Speed of light in vacuum (exact)
  C: 299792458.0, // m/s
  C_KMS: 299792.458, // km/s

  // Newtonian Gravitational Constant (m^3 kg^-1 s^-2)
  G_SI: 6.67430e-11,
  // G in (km/s)^2 * Mpc / M_sun
  G_ASTRO: 4.30091e-9,

  // Solar Mass (kg)
  M_SUN: 1.98847e30, // kg
  M_EARTH: 5.9722e24, // kg

  // Length conversions (meters)
  AU: 1.495978707e11, // m
  PC: 3.085677581491367e16, // m
  KPC: 3.085677581491367e19, // m
  MPC: 3.085677581491367e22, // m
  GPC: 3.085677581491367e25, // m
  LY: 9.4607304725808e15, // m

  // Time conversions (seconds)
  MINUTE: 60.0,
  HOUR: 3600.0,
  DAY: 86400.0,
  YEAR: 31557600.0, // 365.25 days (Julian year)
  MYR: 31557600.0 * 1.0e6,
  GYR: 31557600.0 * 1.0e9,

  // Boltzmann Constant (J/K)
  K_B: 1.380649e-23,

  // Electronvolt (J)
  EV: 1.602176634e-19,
  KEV: 1.602176634e-16,

  // Solar Luminosity (W and erg/s)
  L_SUN: 3.828e26, // Watts
  L_SUN_ERG: 3.828e33, // erg/s

  // Standard Hubble constant reference H0 / h (km/s / Mpc)
  H100_KMS_MPC: 100.0,
  // 1 km/s / Mpc in s^-1
  KMS_MPC_TO_SEC_INV: 1.0e3 / 3.085677581491367e22, // ~3.240779289e-18 s^-1

  // Critical density coefficient rho_crit,0 / h^2 = 3 * (100 km/s/Mpc)^2 / (8 * pi * G) in M_sun / Mpc^3
  // rho_c0 = 2.77536627e11 * h^2 M_sun / Mpc^3 = 2.77536627e11 M_sun / (Mpc/h)^3
  RHO_CRIT_COEFF_MSUN_MPC3: 2.77536627e11
});

/**
 * 5D Dimension Vector: [Length (L), Mass (M), Time (T), Temperature (Theta), Angle (A)]
 */
export class Dimension {
  /**
   * @param {number} [l=0] - Length exponent
   * @param {number} [m=0] - Mass exponent
   * @param {number} [t=0] - Time exponent
   * @param {number} [theta=0] - Temperature exponent
   * @param {number} [angle=0] - Angle exponent
   */
  constructor(l = 0, m = 0, t = 0, theta = 0, angle = 0) {
    this.l = l;
    this.m = m;
    this.t = t;
    this.theta = theta;
    this.angle = angle;
  }

  /**
   * Checks if two dimensions are identical.
   * @param {Dimension} other
   * @returns {boolean}
   */
  equals(other) {
    return this.l === other.l &&
           this.m === other.m &&
           this.t === other.t &&
           this.theta === other.theta &&
           this.angle === other.angle;
  }

  /**
   * Dimensionless check.
   * @returns {boolean}
   */
  isDimensionless() {
    return this.l === 0 && this.m === 0 && this.t === 0 && this.theta === 0 && this.angle === 0;
  }

  /**
   * Multiplies dimensions (adds exponents).
   * @param {Dimension} other
   * @returns {Dimension}
   */
  multiply(other) {
    return new Dimension(
      this.l + other.l,
      this.m + other.m,
      this.t + other.t,
      this.theta + other.theta,
      this.angle + other.angle
    );
  }

  /**
   * Divides dimensions (subtracts exponents).
   * @param {Dimension} other
   * @returns {Dimension}
   */
  divide(other) {
    return new Dimension(
      this.l - other.l,
      this.m - other.m,
      this.t - other.t,
      this.theta - other.theta,
      this.angle - other.angle
    );
  }

  /**
   * Raises dimension to a power (multiplies exponents).
   * @param {number} power
   * @returns {Dimension}
   */
  pow(power) {
    return new Dimension(
      this.l * power,
      this.m * power,
      this.t * power,
      this.theta * power,
      this.angle * power
    );
  }

  /**
   * Formats dimension vector as a human-readable string (e.g. "L^1 T^-1").
   * @returns {string}
   */
  toString() {
    if (this.isDimensionless()) return 'Dimensionless';
    const parts = [];
    if (this.l !== 0) parts.push(`L^${this.l}`);
    if (this.m !== 0) parts.push(`M^${this.m}`);
    if (this.t !== 0) parts.push(`T^${this.t}`);
    if (this.theta !== 0) parts.push(`Θ^${this.theta}`);
    if (this.angle !== 0) parts.push(`A^${this.angle}`);
    return parts.join(' ');
  }
}

/** Standard Dimension Definitions */
export const DIMENSIONS = Object.freeze({
  DIMENSIONLESS: new Dimension(0, 0, 0, 0, 0),
  LENGTH: new Dimension(1, 0, 0, 0, 0),
  AREA: new Dimension(2, 0, 0, 0, 0),
  VOLUME: new Dimension(3, 0, 0, 0, 0),
  MASS: new Dimension(0, 1, 0, 0, 0),
  TIME: new Dimension(0, 0, 1, 0, 0),
  VELOCITY: new Dimension(1, 0, -1, 0, 0),
  ACCELERATION: new Dimension(1, 0, -2, 0, 0),
  FORCE: new Dimension(1, 1, -2, 0, 0),
  ENERGY: new Dimension(2, 1, -2, 0, 0),
  POWER: new Dimension(2, 1, -3, 0, 0),
  DENSITY: new Dimension(-3, 1, 0, 0, 0),
  SURFACE_DENSITY: new Dimension(-2, 1, 0, 0, 0),
  EXPANSION_RATE: new Dimension(0, 0, -1, 0, 0), // Hubble parameter 1/Time
  GRAVITATIONAL_POTENTIAL: new Dimension(2, 0, -2, 0, 0), // (km/s)^2 = m^2/s^2
  TEMPERATURE: new Dimension(0, 0, 0, 1, 0),
  ANGLE: new Dimension(0, 0, 0, 0, 1)
});

/**
 * Custom Error for dimension mismatches.
 */
export class DimensionMismatchError extends Error {
  /**
   * @param {string} fromUnit
   * @param {string} toUnit
   * @param {Dimension} fromDim
   * @param {Dimension} toDim
   */
  constructor(fromUnit, toUnit, fromDim, toDim) {
    super(`Cannot convert incompatible dimensions: "${fromUnit}" (${fromDim.toString()}) to "${toUnit}" (${toDim.toString()})`);
    this.name = 'DimensionMismatchError';
    this.fromUnit = fromUnit;
    this.toUnit = toUnit;
    this.fromDim = fromDim;
    this.toDim = toDim;
  }
}

/**
 * Definition of a single unit in the registry.
 */
export class UnitDefinition {
  /**
   * @param {object} params
   * @param {string} params.symbol - Unit symbol (e.g. 'Mpc/h', 'km/s')
   * @param {Dimension} params.dimension - 5D dimension vector
   * @param {number} params.siFactor - Conversion multiplier to canonical SI base
   * @param {number} [params.hExponent=0] - Exponent of cosmological h parameter (H0 = 100 h km/s/Mpc)
   * @param {string} [params.name] - Descriptive name
   * @param {string} [params.category] - Category (length, velocity, etc.)
   */
  constructor({ symbol, dimension, siFactor, hExponent = 0, name = symbol, category = 'custom' }) {
    this.symbol = symbol;
    this.dimension = dimension;
    this.siFactor = siFactor;
    this.hExponent = hExponent;
    this.name = name;
    this.category = category;
  }

  /**
   * Computes the effective scale to SI base given a cosmological h value.
   * @param {number} [h=1.0] - Cosmological h parameter
   * @returns {number} Effective SI multiplier
   */
  getEffectiveSiFactor(h = 1.0) {
    if (this.hExponent === 0 || h === 1.0) {
      return this.siFactor;
    }
    return this.siFactor * Math.pow(h, this.hExponent);
  }
}

/**
 * Registry of known physical and cosmological units.
 */
export class UnitRegistry {
  constructor() {
    /** @type {Map<string, UnitDefinition>} */
    this.units = new Map();
    this._initBaseUnits();
  }

  /**
   * Initializes canonical astrophysical and SI units.
   * @private
   */
  _initBaseUnits() {
    // ----------------------------------------------------
    // Dimensionless
    // ----------------------------------------------------
    this.register(new UnitDefinition({
      symbol: '1',
      dimension: DIMENSIONS.DIMENSIONLESS,
      siFactor: 1.0,
      name: 'Dimensionless Scalar',
      category: 'dimensionless'
    }));
    this.register(new UnitDefinition({
      symbol: 'dimensionless',
      dimension: DIMENSIONS.DIMENSIONLESS,
      siFactor: 1.0,
      name: 'Dimensionless',
      category: 'dimensionless'
    }));
    this.register(new UnitDefinition({
      symbol: 'percent',
      dimension: DIMENSIONS.DIMENSIONLESS,
      siFactor: 0.01,
      name: 'Percentage',
      category: 'dimensionless'
    }));
    this.register(new UnitDefinition({
      symbol: 'mag',
      dimension: DIMENSIONS.DIMENSIONLESS,
      siFactor: 1.0,
      name: 'Magnitude',
      category: 'dimensionless'
    }));

    // ----------------------------------------------------
    // Length
    // ----------------------------------------------------
    this.register(new UnitDefinition({ symbol: 'm', dimension: DIMENSIONS.LENGTH, siFactor: 1.0, name: 'Meter', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'cm', dimension: DIMENSIONS.LENGTH, siFactor: 1e-2, name: 'Centimeter', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'mm', dimension: DIMENSIONS.LENGTH, siFactor: 1e-3, name: 'Millimeter', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'km', dimension: DIMENSIONS.LENGTH, siFactor: 1e3, name: 'Kilometer', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'au', dimension: DIMENSIONS.LENGTH, siFactor: CONSTANTS.AU, name: 'Astronomical Unit', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'ly', dimension: DIMENSIONS.LENGTH, siFactor: CONSTANTS.LY, name: 'Light Year', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'pc', dimension: DIMENSIONS.LENGTH, siFactor: CONSTANTS.PC, name: 'Parsec', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'kpc', dimension: DIMENSIONS.LENGTH, siFactor: CONSTANTS.KPC, name: 'Kiloparsec', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'Mpc', dimension: DIMENSIONS.LENGTH, siFactor: CONSTANTS.MPC, name: 'Megaparsec', category: 'length' }));
    this.register(new UnitDefinition({ symbol: 'Gpc', dimension: DIMENSIONS.LENGTH, siFactor: CONSTANTS.GPC, name: 'Gigaparsec', category: 'length' }));

    // h-dependent Lengths: 1 Mpc/h = (1/h) Mpc -> siFactor = MPC, hExponent = -1
    this.register(new UnitDefinition({
      symbol: 'Mpc/h',
      dimension: DIMENSIONS.LENGTH,
      siFactor: CONSTANTS.MPC,
      hExponent: -1,
      name: 'Megaparsec per h',
      category: 'length'
    }));
    this.register(new UnitDefinition({
      symbol: 'kpc/h',
      dimension: DIMENSIONS.LENGTH,
      siFactor: CONSTANTS.KPC,
      hExponent: -1,
      name: 'Kiloparsec per h',
      category: 'length'
    }));
    this.register(new UnitDefinition({
      symbol: 'pc/h',
      dimension: DIMENSIONS.LENGTH,
      siFactor: CONSTANTS.PC,
      hExponent: -1,
      name: 'Parsec per h',
      category: 'length'
    }));

    // ----------------------------------------------------
    // Mass
    // ----------------------------------------------------
    this.register(new UnitDefinition({ symbol: 'kg', dimension: DIMENSIONS.MASS, siFactor: 1.0, name: 'Kilogram', category: 'mass' }));
    this.register(new UnitDefinition({ symbol: 'g', dimension: DIMENSIONS.MASS, siFactor: 1e-3, name: 'Gram', category: 'mass' }));
    this.register(new UnitDefinition({ symbol: 'M_sun', dimension: DIMENSIONS.MASS, siFactor: CONSTANTS.M_SUN, name: 'Solar Mass', category: 'mass' }));
    this.register(new UnitDefinition({ symbol: 'Msun', dimension: DIMENSIONS.MASS, siFactor: CONSTANTS.M_SUN, name: 'Solar Mass', category: 'mass' }));
    this.register(new UnitDefinition({ symbol: 'M_earth', dimension: DIMENSIONS.MASS, siFactor: CONSTANTS.M_EARTH, name: 'Earth Mass', category: 'mass' }));
    this.register(new UnitDefinition({ symbol: '1e10 M_sun', dimension: DIMENSIONS.MASS, siFactor: 1e10 * CONSTANTS.M_SUN, name: '10^10 Solar Masses', category: 'mass' }));
    this.register(new UnitDefinition({ symbol: '1e12 M_sun', dimension: DIMENSIONS.MASS, siFactor: 1e12 * CONSTANTS.M_SUN, name: '10^12 Solar Masses', category: 'mass' }));
    this.register(new UnitDefinition({
      symbol: 'M_sun/h',
      dimension: DIMENSIONS.MASS,
      siFactor: CONSTANTS.M_SUN,
      hExponent: -1,
      name: 'Solar Mass per h',
      category: 'mass'
    }));
    this.register(new UnitDefinition({
      symbol: '1e10 M_sun/h',
      dimension: DIMENSIONS.MASS,
      siFactor: 1e10 * CONSTANTS.M_SUN,
      hExponent: -1,
      name: '10^10 Solar Masses per h',
      category: 'mass'
    }));

    // ----------------------------------------------------
    // Time
    // ----------------------------------------------------
    this.register(new UnitDefinition({ symbol: 's', dimension: DIMENSIONS.TIME, siFactor: 1.0, name: 'Second', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'min', dimension: DIMENSIONS.TIME, siFactor: CONSTANTS.MINUTE, name: 'Minute', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'hr', dimension: DIMENSIONS.TIME, siFactor: CONSTANTS.HOUR, name: 'Hour', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'day', dimension: DIMENSIONS.TIME, siFactor: CONSTANTS.DAY, name: 'Day', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'yr', dimension: DIMENSIONS.TIME, siFactor: CONSTANTS.YEAR, name: 'Year', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'Myr', dimension: DIMENSIONS.TIME, siFactor: CONSTANTS.MYR, name: 'Megayear', category: 'time' }));
    this.register(new UnitDefinition({ symbol: 'Gyr', dimension: DIMENSIONS.TIME, siFactor: CONSTANTS.GYR, name: 'Gigayear', category: 'time' }));

    // ----------------------------------------------------
    // Velocity
    // ----------------------------------------------------
    this.register(new UnitDefinition({ symbol: 'm/s', dimension: DIMENSIONS.VELOCITY, siFactor: 1.0, name: 'Meter per second', category: 'velocity' }));
    this.register(new UnitDefinition({ symbol: 'km/s', dimension: DIMENSIONS.VELOCITY, siFactor: 1e3, name: 'Kilometer per second', category: 'velocity' }));
    this.register(new UnitDefinition({ symbol: 'cm/s', dimension: DIMENSIONS.VELOCITY, siFactor: 1e-2, name: 'Centimeter per second', category: 'velocity' }));
    this.register(new UnitDefinition({ symbol: 'c', dimension: DIMENSIONS.VELOCITY, siFactor: CONSTANTS.C, name: 'Speed of Light', category: 'velocity' }));

    // ----------------------------------------------------
    // Expansion Rate / Hubble Parameter (1 / Time)
    // ----------------------------------------------------
    this.register(new UnitDefinition({
      symbol: '1/s',
      dimension: DIMENSIONS.EXPANSION_RATE,
      siFactor: 1.0,
      name: 'Inverse Second',
      category: 'expansion_rate'
    }));
    this.register(new UnitDefinition({
      symbol: 'km/s/Mpc',
      dimension: DIMENSIONS.EXPANSION_RATE,
      siFactor: CONSTANTS.KMS_MPC_TO_SEC_INV,
      hExponent: 0,
      name: 'Kilometer per second per Megaparsec',
      category: 'expansion_rate'
    }));
    this.register(new UnitDefinition({
      symbol: 'km/s/(Mpc/h)',
      dimension: DIMENSIONS.EXPANSION_RATE,
      siFactor: CONSTANTS.KMS_MPC_TO_SEC_INV,
      hExponent: 1, // (km/s)/(Mpc/h) = h * (km/s/Mpc)
      name: 'Kilometer per second per (Mpc/h)',
      category: 'expansion_rate'
    }));
    this.register(new UnitDefinition({
      symbol: '100 km/s/Mpc',
      dimension: DIMENSIONS.EXPANSION_RATE,
      siFactor: 100.0 * CONSTANTS.KMS_MPC_TO_SEC_INV,
      name: '100 km/s/Mpc Hubble Benchmark',
      category: 'expansion_rate'
    }));

    // ----------------------------------------------------
    // Density
    // ----------------------------------------------------
    this.register(new UnitDefinition({
      symbol: 'kg/m^3',
      dimension: DIMENSIONS.DENSITY,
      siFactor: 1.0,
      name: 'Kilogram per cubic meter',
      category: 'density'
    }));
    this.register(new UnitDefinition({
      symbol: 'g/cm^3',
      dimension: DIMENSIONS.DENSITY,
      siFactor: 1e3,
      name: 'Gram per cubic centimeter',
      category: 'density'
    }));
    // M_sun / Mpc^3: mass / length^3
    const msunMpc3ToSi = CONSTANTS.M_SUN / Math.pow(CONSTANTS.MPC, 3);
    this.register(new UnitDefinition({
      symbol: 'M_sun/Mpc^3',
      dimension: DIMENSIONS.DENSITY,
      siFactor: msunMpc3ToSi,
      hExponent: 0,
      name: 'Solar Mass per cubic Megaparsec',
      category: 'density'
    }));
    // M_sun/(Mpc/h)^3: M_sun / (Mpc/h)^3 = M_sun * h^3 / Mpc^3 -> hExponent = +3 (or with M_sun/h -> +2)
    this.register(new UnitDefinition({
      symbol: 'M_sun/(Mpc/h)^3',
      dimension: DIMENSIONS.DENSITY,
      siFactor: msunMpc3ToSi,
      hExponent: 3,
      name: 'Solar Mass per cubic (Mpc/h)',
      category: 'density'
    }));
    this.register(new UnitDefinition({
      symbol: '(M_sun/h)/(Mpc/h)^3',
      dimension: DIMENSIONS.DENSITY,
      siFactor: msunMpc3ToSi,
      hExponent: 2, // (M_sun/h) / (Mpc/h)^3 = h^2 * M_sun / Mpc^3
      name: 'Solar Mass per h per cubic (Mpc/h)',
      category: 'density'
    }));
    // rho_crit unit: critical cosmological density at z=0 (2.77536627e11 h^2 M_sun/Mpc^3)
    this.register(new UnitDefinition({
      symbol: 'rho_crit',
      dimension: DIMENSIONS.DENSITY,
      siFactor: CONSTANTS.RHO_CRIT_COEFF_MSUN_MPC3 * msunMpc3ToSi,
      hExponent: 2,
      name: 'Critical Cosmological Density rho_crit,0',
      category: 'density'
    }));

    // ----------------------------------------------------
    // Gravitational Potential
    // ----------------------------------------------------
    this.register(new UnitDefinition({
      symbol: 'm^2/s^2',
      dimension: DIMENSIONS.GRAVITATIONAL_POTENTIAL,
      siFactor: 1.0,
      name: 'Square meter per square second',
      category: 'potential'
    }));
    this.register(new UnitDefinition({
      symbol: '(km/s)^2',
      dimension: DIMENSIONS.GRAVITATIONAL_POTENTIAL,
      siFactor: 1.0e6,
      name: 'Square kilometer per square second',
      category: 'potential'
    }));
    this.register(new UnitDefinition({
      symbol: 'c^2',
      dimension: DIMENSIONS.GRAVITATIONAL_POTENTIAL,
      siFactor: CONSTANTS.C * CONSTANTS.C,
      name: 'Speed of Light Squared',
      category: 'potential'
    }));

    // ----------------------------------------------------
    // Temperature & Energy
    // ----------------------------------------------------
    this.register(new UnitDefinition({ symbol: 'K', dimension: DIMENSIONS.TEMPERATURE, siFactor: 1.0, name: 'Kelvin', category: 'temperature' }));
    this.register(new UnitDefinition({ symbol: 'J', dimension: DIMENSIONS.ENERGY, siFactor: 1.0, name: 'Joule', category: 'energy' }));
    this.register(new UnitDefinition({ symbol: 'erg', dimension: DIMENSIONS.ENERGY, siFactor: 1e-7, name: 'Erg', category: 'energy' }));
    this.register(new UnitDefinition({ symbol: 'eV', dimension: DIMENSIONS.ENERGY, siFactor: CONSTANTS.EV, name: 'Electronvolt', category: 'energy' }));
    this.register(new UnitDefinition({ symbol: 'keV', dimension: DIMENSIONS.ENERGY, siFactor: CONSTANTS.KEV, name: 'Kiloelectronvolt', category: 'energy' }));

    // ----------------------------------------------------
    // Angle
    // ----------------------------------------------------
    this.register(new UnitDefinition({ symbol: 'rad', dimension: DIMENSIONS.ANGLE, siFactor: 1.0, name: 'Radian', category: 'angle' }));
    this.register(new UnitDefinition({ symbol: 'deg', dimension: DIMENSIONS.ANGLE, siFactor: Math.PI / 180.0, name: 'Degree', category: 'angle' }));
    this.register(new UnitDefinition({ symbol: 'arcmin', dimension: DIMENSIONS.ANGLE, siFactor: Math.PI / (180.0 * 60.0), name: 'Arcminute', category: 'angle' }));
    this.register(new UnitDefinition({ symbol: 'arcsec', dimension: DIMENSIONS.ANGLE, siFactor: Math.PI / (180.0 * 3600.0), name: 'Arcsecond', category: 'angle' }));
  }

  /**
   * Registers a unit definition.
   * @param {UnitDefinition} unitDef
   */
  register(unitDef) {
    this.units.set(unitDef.symbol, unitDef);
    // Normalize aliases (e.g. without spaces or standard variants)
    const norm = unitDef.symbol.replace(/\s+/g, '');
    if (norm !== unitDef.symbol && !this.units.has(norm)) {
      this.units.set(norm, unitDef);
    }
  }

  /**
   * Looks up a unit definition by symbol or name.
   * @param {string} symbol
   * @returns {UnitDefinition}
   */
  get(symbol) {
    const trimmed = (symbol || '').trim();
    if (this.units.has(trimmed)) {
      return this.units.get(trimmed);
    }
    const noSpace = trimmed.replace(/\s+/g, '');
    if (this.units.has(noSpace)) {
      return this.units.get(noSpace);
    }
    throw new Error(`Unknown unit symbol: "${symbol}". Registered units: ${Array.from(this.units.keys()).join(', ')}`);
  }

  /**
   * Checks whether a unit symbol is registered.
   * @param {string} symbol
   * @returns {boolean}
   */
  has(symbol) {
    const trimmed = (symbol || '').trim();
    return this.units.has(trimmed) || this.units.has(trimmed.replace(/\s+/g, ''));
  }

  /**
   * Converts a numeric value from one unit to another, verifying dimensional consistency.
   * 
   * @param {number} value - Numeric magnitude in fromUnit.
   * @param {string} fromUnit - Source unit symbol.
   * @param {string} toUnit - Target unit symbol.
   * @param {object} [context] - Cosmological context.
   * @param {number} [context.h=0.746] - Dimensionless Hubble parameter (default CF4 value 0.746).
   * @returns {number} Converted value.
   */
  convert(value, fromUnit, toUnit, context = {}) {
    if (fromUnit === toUnit) return value;

    const fromDef = this.get(fromUnit);
    const toDef = this.get(toUnit);

    if (!fromDef.dimension.equals(toDef.dimension)) {
      throw new DimensionMismatchError(fromUnit, toUnit, fromDef.dimension, toDef.dimension);
    }

    const h = typeof context.h === 'number' ? context.h : 0.746;
    const fromSi = fromDef.getEffectiveSiFactor(h);
    const toSi = toDef.getEffectiveSiFactor(h);

    const valSi = value * fromSi;
    return valSi / toSi;
  }
}

/** Global Unit Registry Singleton */
export const defaultRegistry = new UnitRegistry();

/**
 * Strongly typed Physical Quantity encapsulating magnitude, unit, and dimensional integrity.
 */
export class Quantity {
  /**
   * @param {number} value - Numeric magnitude.
   * @param {string|UnitDefinition} unit - Unit symbol or definition.
   * @param {UnitRegistry} [registry=defaultRegistry] - Unit registry to use.
   */
  constructor(value, unit = '1', registry = defaultRegistry) {
    if (typeof value !== 'number' || isNaN(value)) {
      throw new TypeError(`Quantity magnitude must be a valid number. Received: ${value}`);
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
  }

  /**
   * Converts this quantity to target unit.
   * 
   * @param {string} targetUnit - Destination unit symbol.
   * @param {object} [context] - Cosmological parameters (e.g. { h: 0.746 }).
   * @returns {Quantity} New converted Quantity.
   */
  to(targetUnit, context = {}) {
    const convertedVal = this.registry.convert(this.value, this.unit, targetUnit, context);
    return new Quantity(convertedVal, targetUnit, this.registry);
  }

  /**
   * Converts to canonical SI unit.
   * @param {object} [context]
   * @returns {Quantity}
   */
  toSI(context = {}) {
    const h = typeof context.h === 'number' ? context.h : 0.746;
    const siVal = this.value * this.unitDef.getEffectiveSiFactor(h);
    // Find matching base SI unit
    for (const u of this.registry.units.values()) {
      if (u.dimension.equals(this.dimension) && u.hExponent === 0 && u.siFactor === 1.0) {
        return new Quantity(siVal, u.symbol, this.registry);
      }
    }
    return new Quantity(siVal, '1', this.registry);
  }

  /**
   * Adds another quantity (dimensionally checked).
   * 
   * @param {Quantity} other
   * @param {object} [context]
   * @returns {Quantity}
   */
  add(other, context = {}) {
    if (!(other instanceof Quantity)) {
      throw new TypeError(`Cannot add non-Quantity to Quantity: ${other}`);
    }
    const otherInThisUnit = other.to(this.unit, context);
    return new Quantity(this.value + otherInThisUnit.value, this.unit, this.registry);
  }

  /**
   * Subtracts another quantity (dimensionally checked).
   * 
   * @param {Quantity} other
   * @param {object} [context]
   * @returns {Quantity}
   */
  subtract(other, context = {}) {
    if (!(other instanceof Quantity)) {
      throw new TypeError(`Cannot subtract non-Quantity from Quantity: ${other}`);
    }
    const otherInThisUnit = other.to(this.unit, context);
    return new Quantity(this.value - otherInThisUnit.value, this.unit, this.registry);
  }

  /**
   * Multiplies by a scalar or another Quantity.
   * 
   * @param {number|Quantity} other
   * @param {object} [context]
   * @returns {Quantity}
   */
  multiply(other, context = {}) {
    if (typeof other === 'number') {
      return new Quantity(this.value * other, this.unit, this.registry);
    }
    if (other instanceof Quantity) {
      const h = typeof context.h === 'number' ? context.h : 0.746;
      const thisSi = this.value * this.unitDef.getEffectiveSiFactor(h);
      const otherSi = other.value * other.unitDef.getEffectiveSiFactor(h);
      const newDim = this.dimension.multiply(other.dimension);
      const newValSi = thisSi * otherSi;

      // Check if exact composite unit exists
      for (const u of this.registry.units.values()) {
        if (u.dimension.equals(newDim) && u.hExponent === 0 && u.siFactor === 1.0) {
          return new Quantity(newValSi, u.symbol, this.registry);
        }
      }
      return new Quantity(newValSi, new UnitDefinition({
        symbol: `(${this.unit}*${other.unit})`,
        dimension: newDim,
        siFactor: 1.0
      }), this.registry);
    }
    throw new TypeError(`Multiply argument must be number or Quantity. Received: ${other}`);
  }

  /**
   * Divides by a scalar or another Quantity.
   * 
   * @param {number|Quantity} other
   * @param {object} [context]
   * @returns {Quantity}
   */
  divide(other, context = {}) {
    if (typeof other === 'number') {
      if (other === 0) throw new RangeError('Division by zero in Quantity');
      return new Quantity(this.value / other, this.unit, this.registry);
    }
    if (other instanceof Quantity) {
      if (other.value === 0) throw new RangeError('Division by zero in Quantity');
      const h = typeof context.h === 'number' ? context.h : 0.746;
      const thisSi = this.value * this.unitDef.getEffectiveSiFactor(h);
      const otherSi = other.value * other.unitDef.getEffectiveSiFactor(h);
      const newDim = this.dimension.divide(other.dimension);
      const newValSi = thisSi / otherSi;

      for (const u of this.registry.units.values()) {
        if (u.dimension.equals(newDim) && u.hExponent === 0 && u.siFactor === 1.0) {
          return new Quantity(newValSi, u.symbol, this.registry);
        }
      }
      return new Quantity(newValSi, new UnitDefinition({
        symbol: `(${this.unit}/${other.unit})`,
        dimension: newDim,
        siFactor: 1.0
      }), this.registry);
    }
    throw new TypeError(`Divide argument must be number or Quantity. Received: ${other}`);
  }

  /**
   * Raises quantity to a real power.
   * 
   * @param {number} p - Power exponent
   * @param {object} [context]
   * @returns {Quantity}
   */
  pow(p, context = {}) {
    const h = typeof context.h === 'number' ? context.h : 0.746;
    const thisSi = this.value * this.unitDef.getEffectiveSiFactor(h);
    const newDim = this.dimension.pow(p);
    const newValSi = Math.pow(thisSi, p);

    for (const u of this.registry.units.values()) {
      if (u.dimension.equals(newDim) && u.hExponent === 0 && u.siFactor === 1.0) {
        return new Quantity(newValSi, u.symbol, this.registry);
      }
    }
    return new Quantity(newValSi, new UnitDefinition({
      symbol: `(${this.unit})^${p}`,
      dimension: newDim,
      siFactor: 1.0
    }), this.registry);
  }

  /**
   * Equality test within tolerance.
   * 
   * @param {Quantity} other
   * @param {number} [relTol=1e-7]
   * @param {object} [context]
   * @returns {boolean}
   */
  equals(other, relTol = 1e-7, context = {}) {
    if (!(other instanceof Quantity)) return false;
    if (!this.dimension.equals(other.dimension)) return false;
    const otherInThisUnit = other.to(this.unit, context);
    const diff = Math.abs(this.value - otherInThisUnit.value);
    const maxVal = Math.max(Math.abs(this.value), Math.abs(otherInThisUnit.value), 1e-15);
    return diff / maxVal <= relTol;
  }

  /**
   * Formats quantity for display.
   * 
   * @param {number} [precision=4] - Significant digits
   * @returns {string} Formatted string (e.g. "100.0 km/s")
   */
  format(precision = 4) {
    const numStr = Math.abs(this.value) >= 1e5 || (Math.abs(this.value) < 1e-3 && this.value !== 0)
      ? this.value.toExponential(precision - 1)
      : this.value.toPrecision(precision);
    return `${parseFloat(numStr)} ${this.unit}`;
  }

  /**
   * Numeric coercion.
   * @returns {number}
   */
  valueOf() {
    return this.value;
  }

  /**
   * String representation.
   * @returns {string}
   */
  toString() {
    return `${this.value} ${this.unit}`;
  }
}

/**
 * Shorthand helper to construct a Quantity.
 * @param {number} value
 * @param {string} unit
 * @returns {Quantity}
 */
export function qty(value, unit) {
  return new Quantity(value, unit, defaultRegistry);
}
