/**
 * @file unit_conversion.js
 * @description High-precision dimensional conversion engine and dimensionally typed physical quantities
 * for the CosmicFlows-4 Research Workbench.
 * 
 * Supports:
 * - Length: Mpc/h, Mpc, kpc/h, kpc, pc, AU, km, m, ly
 * - Velocity: km/s, m/s, c, pc/Myr, AU/yr
 * - Mass: M_sun/h, M_sun, 10^10 M_sun/h, kg, g
 * - Density: (M_sun/h)/(Mpc/h)^3, M_sun/Mpc^3, rho_crit, g/cm^3, kg/m^3
 * - Potential/Energy per mass: (km/s)^2, m^2/s^2, J/kg, erg/g
 * - Acceleration: (km/s)^2/(Mpc/h), m/s^2, cm/s^2
 * - Time/Frequency: s, yr, Myr, Gyr, Hubble time, km/s/(Mpc/h), km/s/Mpc, s^-1
 * 
 * @module units/unit_conversion
 */

import { PHYSICAL_CONSTANTS, COSMOLOGICAL_MODELS } from './cosmological_constants.js';

/**
 * Base physical dimension exponents [Length, Mass, Time].
 */
export class Dimensions {
  /**
   * @param {number} [length=0] Length dimension exponent L
   * @param {number} [mass=0] Mass dimension exponent M
   * @param {number} [time=0] Time dimension exponent T
   */
  constructor(length = 0, mass = 0, time = 0) {
    this.L = length;
    this.M = mass;
    this.T = time;
  }

  equals(other) {
    if (!(other instanceof Dimensions)) return false;
    return this.L === other.L && this.M === other.M && this.T === other.T;
  }

  isDimensionless() {
    return this.L === 0 && this.M === 0 && this.T === 0;
  }

  multiply(other) {
    return new Dimensions(this.L + other.L, this.M + other.M, this.T + other.T);
  }

  divide(other) {
    return new Dimensions(this.L - other.L, this.M - other.M, this.T - other.T);
  }

  pow(p) {
    return new Dimensions(this.L * p, this.M * p, this.T * p);
  }

  toString() {
    if (this.isDimensionless()) return '1';
    const parts = [];
    if (this.L !== 0) parts.push(this.L === 1 ? 'L' : `L^${this.L}`);
    if (this.M !== 0) parts.push(this.M === 1 ? 'M' : `M^${this.M}`);
    if (this.T !== 0) parts.push(this.T === 1 ? 'T' : `T^${this.T}`);
    return parts.join(' ');
  }
}

/**
 * Standard cosmological dimensions.
 * @readonly
 */
export const DIMENSION_TYPES = Object.freeze({
  DIMENSIONLESS: Object.freeze(new Dimensions(0, 0, 0)),
  LENGTH: Object.freeze(new Dimensions(1, 0, 0)),
  MASS: Object.freeze(new Dimensions(0, 1, 0)),
  TIME: Object.freeze(new Dimensions(0, 0, 1)),
  VELOCITY: Object.freeze(new Dimensions(1, 0, -1)),
  ACCELERATION: Object.freeze(new Dimensions(1, 0, -2)),
  DENSITY: Object.freeze(new Dimensions(-3, 1, 0)),
  SPECIFIC_ENERGY: Object.freeze(new Dimensions(2, 0, -2)), // (km/s)^2 or J/kg
  FREQUENCY: Object.freeze(new Dimensions(0, 0, -1)),       // Hubble constant H0, s^-1, km/s/Mpc
  FORCE: Object.freeze(new Dimensions(1, 1, -2)),
  ENERGY: Object.freeze(new Dimensions(2, 1, -2)),
  PRESSURE: Object.freeze(new Dimensions(-1, 1, -2))
});

/**
 * UnitDefinition stores conversion factor to SI base units and h-scaling exponent.
 */
export class UnitDefinition {
  /**
   * @param {string} name Unit symbol/name.
   * @param {Dimensions} dimensions Base physical dimensions.
   * @param {number} toSI Conversion multiplier to SI base units (m, kg, s) assuming h=1.
   * @param {number} [hExponent=0] Exponent of dimensionless Hubble parameter h (e.g. Mpc/h has hExponent = -1).
   * @param {string} [description='']
   */
  constructor(name, dimensions, toSI, hExponent = 0, description = '') {
    this.name = name;
    this.dimensions = dimensions;
    this.toSI = toSI;
    this.hExponent = hExponent;
    this.description = description;
  }

  /**
   * Computes effective multiplier to SI base units given a specific h parameter.
   * @param {number} [h=0.746]
   * @returns {number}
   */
  getMultiplierToSI(h = 0.746) {
    if (this.hExponent === 0) return this.toSI;
    return this.toSI * Math.pow(h, this.hExponent);
  }
}

/**
 * Registry of standard cosmological units.
 * @readonly
 */
export const UNIT_TABLE = Object.freeze({
  // Dimensionless
  'dimensionless': new UnitDefinition('dimensionless', DIMENSION_TYPES.DIMENSIONLESS, 1.0, 0, 'Pure dimensionless scalar'),
  '1': new UnitDefinition('1', DIMENSION_TYPES.DIMENSIONLESS, 1.0, 0, 'Pure dimensionless scalar'),

  // Length
  'm': new UnitDefinition('m', DIMENSION_TYPES.LENGTH, 1.0, 0, 'Meter (SI)'),
  'km': new UnitDefinition('km', DIMENSION_TYPES.LENGTH, 1.0e3, 0, 'Kilometer'),
  'au': new UnitDefinition('au', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_AU, 0, 'Astronomical Unit'),
  'ly': new UnitDefinition('ly', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_LY, 0, 'Light-year'),
  'pc': new UnitDefinition('pc', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_PC, 0, 'Parsec'),
  'kpc': new UnitDefinition('kpc', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_KPC, 0, 'Kiloparsec'),
  'Mpc': new UnitDefinition('Mpc', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_MPC, 0, 'Megaparsec'),
  'Gpc': new UnitDefinition('Gpc', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_MPC * 1.0e3, 0, 'Gigaparsec'),
  'kpc/h': new UnitDefinition('kpc/h', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_KPC, -1, 'Kiloparsecs over h'),
  'Mpc/h': new UnitDefinition('Mpc/h', DIMENSION_TYPES.LENGTH, PHYSICAL_CONSTANTS.METERS_PER_MPC, -1, 'Megaparsecs over h (CF4 spatial coordinate)'),

  // Velocity
  'm/s': new UnitDefinition('m/s', DIMENSION_TYPES.VELOCITY, 1.0, 0, 'Meters per second'),
  'km/s': new UnitDefinition('km/s', DIMENSION_TYPES.VELOCITY, 1.0e3, 0, 'Kilometers per second (CF4 peculiar velocity)'),
  'c': new UnitDefinition('c', DIMENSION_TYPES.VELOCITY, PHYSICAL_CONSTANTS.SPEED_OF_LIGHT_SI, 0, 'Speed of light'),
  'pc/Myr': new UnitDefinition('pc/Myr', DIMENSION_TYPES.VELOCITY, PHYSICAL_CONSTANTS.METERS_PER_PC / PHYSICAL_CONSTANTS.SECONDS_PER_MYR, 0, 'Parsecs per Megayear (~0.978 km/s)'),

  // Time
  's': new UnitDefinition('s', DIMENSION_TYPES.TIME, 1.0, 0, 'Second (SI)'),
  'yr': new UnitDefinition('yr', DIMENSION_TYPES.TIME, PHYSICAL_CONSTANTS.SECONDS_PER_JULIAN_YEAR, 0, 'Julian year'),
  'Myr': new UnitDefinition('Myr', DIMENSION_TYPES.TIME, PHYSICAL_CONSTANTS.SECONDS_PER_MYR, 0, 'Megayear'),
  'Gyr': new UnitDefinition('Gyr', DIMENSION_TYPES.TIME, PHYSICAL_CONSTANTS.SECONDS_PER_GYR, 0, 'Gigayear'),

  // Mass
  'kg': new UnitDefinition('kg', DIMENSION_TYPES.MASS, 1.0, 0, 'Kilogram (SI)'),
  'g': new UnitDefinition('g', DIMENSION_TYPES.MASS, 1.0e-3, 0, 'Gram (CGS)'),
  'M_sun': new UnitDefinition('M_sun', DIMENSION_TYPES.MASS, PHYSICAL_CONSTANTS.SOLAR_MASS_KG, 0, 'Solar Mass'),
  'M_sun/h': new UnitDefinition('M_sun/h', DIMENSION_TYPES.MASS, PHYSICAL_CONSTANTS.SOLAR_MASS_KG, -1, 'Solar Mass over h'),
  '10^10 M_sun/h': new UnitDefinition('10^10 M_sun/h', DIMENSION_TYPES.MASS, PHYSICAL_CONSTANTS.SOLAR_MASS_KG * 1.0e10, -1, '10^10 Solar Masses over h'),

  // Density
  'kg/m^3': new UnitDefinition('kg/m^3', DIMENSION_TYPES.DENSITY, 1.0, 0, 'Kilograms per cubic meter'),
  'g/cm^3': new UnitDefinition('g/cm^3', DIMENSION_TYPES.DENSITY, 1.0e3, 0, 'Grams per cubic centimeter'),
  'M_sun/Mpc^3': new UnitDefinition('M_sun/Mpc^3', DIMENSION_TYPES.DENSITY, PHYSICAL_CONSTANTS.SOLAR_MASS_KG / Math.pow(PHYSICAL_CONSTANTS.METERS_PER_MPC, 3), 0, 'Solar masses per cubic Megaparsec'),
  '(M_sun/h)/(Mpc/h)^3': new UnitDefinition('(M_sun/h)/(Mpc/h)^3', DIMENSION_TYPES.DENSITY, PHYSICAL_CONSTANTS.SOLAR_MASS_KG / Math.pow(PHYSICAL_CONSTANTS.METERS_PER_MPC, 3), 2, 'Solar masses over h per (Mpc/h)^3'),

  // Specific Energy / Gravitational Potential
  'm^2/s^2': new UnitDefinition('m^2/s^2', DIMENSION_TYPES.SPECIFIC_ENERGY, 1.0, 0, 'Joules per kilogram'),
  '(km/s)^2': new UnitDefinition('(km/s)^2', DIMENSION_TYPES.SPECIFIC_ENERGY, 1.0e6, 0, '(km/s)^2 (Astrophysical potential)'),
  'J/kg': new UnitDefinition('J/kg', DIMENSION_TYPES.SPECIFIC_ENERGY, 1.0, 0, 'Joules per kilogram'),
  'erg/g': new UnitDefinition('erg/g', DIMENSION_TYPES.SPECIFIC_ENERGY, 1.0e-4, 0, 'Ergs per gram'),

  // Frequency / Hubble expansion rate
  's^-1': new UnitDefinition('s^-1', DIMENSION_TYPES.FREQUENCY, 1.0, 0, 'Per second'),
  'km/s/Mpc': new UnitDefinition('km/s/Mpc', DIMENSION_TYPES.FREQUENCY, 1.0e3 / PHYSICAL_CONSTANTS.METERS_PER_MPC, 0, 'km/s per Mpc'),
  'km/s/(Mpc/h)': new UnitDefinition('km/s/(Mpc/h)', DIMENSION_TYPES.FREQUENCY, 1.0e3 / PHYSICAL_CONSTANTS.METERS_PER_MPC, 1, 'km/s per (Mpc/h) [100 km/s/Mpc]')
});

/**
 * UnitConverter handles exact conversions across physical and astrophysical units.
 */
export class UnitConverter {
  /**
   * @param {Object} [options]
   * @param {number} [options.h=0.746] Dimensionless Hubble parameter.
   */
  constructor(options = {}) {
    this.h = options.h ?? COSMOLOGICAL_MODELS.CF4.h;
  }

  /**
   * Looks up a unit definition.
   * @param {string} unitName
   * @returns {UnitDefinition}
   */
  getUnit(unitName) {
    const def = UNIT_TABLE[unitName];
    if (!def) {
      throw new Error(`UnitConverter: Unknown unit '${unitName}'.`);
    }
    return def;
  }

  /**
   * Converts a numeric value from sourceUnit to targetUnit.
   * @param {number} value
   * @param {string} fromUnit
   * @param {string} toUnit
   * @param {number} [customH] Optional override for h.
   * @returns {number}
   */
  convert(value, fromUnit, toUnit, customH) {
    if (fromUnit === toUnit) return value;
    const h = customH ?? this.h;
    const src = this.getUnit(fromUnit);
    const dst = this.getUnit(toUnit);

    if (!src.dimensions.equals(dst.dimensions)) {
      throw new TypeError(`UnitConverter: Incompatible dimensions for conversion from '${fromUnit}' (${src.dimensions}) to '${toUnit}' (${dst.dimensions}).`);
    }

    const valueSI = value * src.getMultiplierToSI(h);
    const result = valueSI / dst.getMultiplierToSI(h);
    return result;
  }

  /**
   * Specific helper: converts spatial coordinate from Mpc/h to Mpc.
   * @param {number} valMpcOverH
   * @returns {number} valMpc = valMpcOverH / h
   */
  mpcOverH_to_mpc(valMpcOverH) {
    return valMpcOverH / this.h;
  }

  /**
   * Specific helper: converts spatial coordinate from Mpc to Mpc/h.
   * @param {number} valMpc
   * @returns {number} valMpcOverH = valMpc * h
   */
  mpc_to_mpcOverH(valMpc) {
    return valMpc * this.h;
  }

  /**
   * Specific helper: converts peculiar velocity in km/s to m/s.
   * @param {number} vKmS
   * @returns {number}
   */
  kms_to_ms(vKmS) {
    return vKmS * 1.0e3;
  }

  /**
   * Converts potential from (km/s)^2 to SI (m^2/s^2 or J/kg).
   * @param {number} potKms2
   * @returns {number}
   */
  potential_to_SI(potKms2) {
    return potKms2 * 1.0e6;
  }
}

/**
 * PhysicalQuantity wraps a numerical magnitude with a typed unit and dimensions.
 */
export class PhysicalQuantity {
  /**
   * @param {number} value Numerical magnitude.
   * @param {string} unit Unit identifier.
   * @param {number} [h=0.746] Dimensionless Hubble parameter.
   */
  constructor(value, unit, h = 0.746) {
    if (typeof value !== 'number' || isNaN(value)) {
      throw new TypeError(`PhysicalQuantity: value must be a valid number, received ${value}`);
    }
    const unitDef = UNIT_TABLE[unit];
    if (!unitDef) {
      throw new Error(`PhysicalQuantity: Unknown unit '${unit}'`);
    }

    this.value = value;
    this.unit = unit;
    this.h = h;
    this.dimensions = unitDef.dimensions;
    this.unitDef = unitDef;
  }

  /**
   * Converts this quantity to a target unit.
   * @param {string} targetUnit
   * @returns {PhysicalQuantity}
   */
  to(targetUnit) {
    const converter = new UnitConverter({ h: this.h });
    const convertedVal = converter.convert(this.value, this.unit, targetUnit, this.h);
    return new PhysicalQuantity(convertedVal, targetUnit, this.h);
  }

  /**
   * Adds another PhysicalQuantity.
   * @param {PhysicalQuantity} other
   * @returns {PhysicalQuantity}
   */
  add(other) {
    if (!(other instanceof PhysicalQuantity)) {
      throw new TypeError('PhysicalQuantity.add: argument must be a PhysicalQuantity.');
    }
    if (!this.dimensions.equals(other.dimensions)) {
      throw new TypeError(`PhysicalQuantity.add: Cannot add quantities with incompatible dimensions: ${this.unit} (${this.dimensions}) vs ${other.unit} (${other.dimensions})`);
    }
    const otherInThisUnit = other.to(this.unit).value;
    return new PhysicalQuantity(this.value + otherInThisUnit, this.unit, this.h);
  }

  /**
   * Subtracts another PhysicalQuantity.
   * @param {PhysicalQuantity} other
   * @returns {PhysicalQuantity}
   */
  subtract(other) {
    if (!(other instanceof PhysicalQuantity)) {
      throw new TypeError('PhysicalQuantity.subtract: argument must be a PhysicalQuantity.');
    }
    if (!this.dimensions.equals(other.dimensions)) {
      throw new TypeError(`PhysicalQuantity.subtract: Cannot subtract quantities with incompatible dimensions: ${this.unit} (${this.dimensions}) vs ${other.unit} (${other.dimensions})`);
    }
    const otherInThisUnit = other.to(this.unit).value;
    return new PhysicalQuantity(this.value - otherInThisUnit, this.unit, this.h);
  }

  /**
   * Multiplies by a scalar number or another PhysicalQuantity.
   * @param {number|PhysicalQuantity} factor
   * @returns {PhysicalQuantity}
   */
  multiply(factor) {
    if (typeof factor === 'number') {
      return new PhysicalQuantity(this.value * factor, this.unit, this.h);
    }
    if (factor instanceof PhysicalQuantity) {
      const valSI1 = this.value * this.unitDef.getMultiplierToSI(this.h);
      const valSI2 = factor.value * factor.unitDef.getMultiplierToSI(this.h);
      const newSI = valSI1 * valSI2;
      const newDim = this.dimensions.multiply(factor.dimensions);
      // Construct a generic SI unit representation
      return new PhysicalQuantity(newSI, 'dimensionless', this.h);
    }
    throw new TypeError('PhysicalQuantity.multiply: argument must be number or PhysicalQuantity.');
  }

  toString() {
    return `${this.value} ${this.unit}`;
  }
}
