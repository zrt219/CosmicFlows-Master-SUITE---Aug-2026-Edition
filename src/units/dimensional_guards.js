/**
 * @file dimensional_guards.js
 * @description Strict runtime scientific dimension guards, invariant assertions,
 * and cosmological physics boundary checks for the CosmicFlows-4 Research Workbench.
 * 
 * Enforces:
 * 1. Physical Dimensions Are Strict:
 *    - Supergalactic Cartesian spatial positions must be in Mpc/h.
 *    - Peculiar velocities must be in km/s.
 *    - Density contrast (delta) is dimensionless.
 *    - Never add, subtract, or combine incompatible physical dimensions (e.g. pos + vel) without explicit dt.
 *    - Incompatible operations must throw a TypeError immediately.
 * 2. Velocity Scale Factor:
 *    - x52.0 applied exactly once to CF4 velocity and error grids.
 *    - Double application must throw an Error.
 * 3. Linear Theory Continuity Diagnostics:
 *    - Validate div(v) approx -H0 * f * delta in the linear regime (|delta| < delta_linear_max).
 * 4. Coordinate frame boundaries and numerical stability guards.
 * 
 * @module units/dimensional_guards
 */

import {
  SupergalacticPosition,
  VelocityVector,
  DensityContrast,
  VelocityError,
  CF4_PUBLIC_VELOCITY_SCALE
} from '../coordinates/scientific_types.js';

/**
 * Custom error types for scientific computing failures.
 */
export class DimensionalityError extends TypeError {
  constructor(message) {
    super(`[DIMENSIONALITY ERROR] ${message}`);
    this.name = 'DimensionalityError';
  }
}

export class InvariantViolationError extends Error {
  constructor(message) {
    super(`[SCIENTIFIC INVARIANT VIOLATION] ${message}`);
    this.name = 'InvariantViolationError';
  }
}

export class BoundaryExtrapolationError extends RangeError {
  constructor(message) {
    super(`[BOUNDARY ERROR] ${message}`);
    this.name = 'BoundaryExtrapolationError';
  }
}

/**
 * Asserts that an object is a valid SupergalacticPosition in Mpc/h.
 * @param {*} obj
 * @param {string} [context='']
 */
export function assertSupergalacticPosition(obj, context = '') {
  if (!obj || typeof obj !== 'object' || typeof obj.isSupergalacticPosition !== 'function' || !obj.isSupergalacticPosition()) {
    const typeName = obj?.constructor?.name || typeof obj;
    throw new DimensionalityError(
      `${context ? context + ': ' : ''}Expected SupergalacticPosition (Mpc/h), received '${typeName}'.`
    );
  }
}

/**
 * Asserts that an object is a valid VelocityVector in km/s.
 * @param {*} obj
 * @param {string} [context='']
 */
export function assertVelocityVector(obj, context = '') {
  if (!obj || typeof obj !== 'object' || typeof obj.isVelocityVector !== 'function' || !obj.isVelocityVector()) {
    const typeName = obj?.constructor?.name || typeof obj;
    throw new DimensionalityError(
      `${context ? context + ': ' : ''}Expected VelocityVector (km/s), received '${typeName}'.`
    );
  }
}

/**
 * Asserts that an object is a valid DensityContrast (dimensionless delta).
 * @param {*} obj
 * @param {string} [context='']
 */
export function assertDensityContrast(obj, context = '') {
  if (typeof obj === 'number') {
    if (isNaN(obj) || !isFinite(obj)) {
      throw new InvariantViolationError(`${context ? context + ': ' : ''}Density contrast is non-finite (${obj}).`);
    }
    if (obj < -1.0) {
      throw new InvariantViolationError(
        `${context ? context + ': ' : ''}Unphysical negative density contrast: delta = ${obj} < -1.0 (density rho < 0).`
      );
    }
    return;
  }
  if (!obj || typeof obj !== 'object' || typeof obj.isDensityContrast !== 'function' || !obj.isDensityContrast()) {
    const typeName = obj?.constructor?.name || typeof obj;
    throw new DimensionalityError(
      `${context ? context + ': ' : ''}Expected DensityContrast, received '${typeName}'.`
    );
  }
}

/**
 * Runtime guard against direct addition/subtraction of positions and velocities.
 * @param {*} a First operand
 * @param {*} b Second operand
 * @param {string} op Operation name ('add', 'subtract', etc.)
 */
export function guardDimensionalArithmetic(a, b, op = 'arithmetic') {
  const isPosA = a instanceof SupergalacticPosition;
  const isPosB = b instanceof SupergalacticPosition;
  const isVelA = a instanceof VelocityVector;
  const isVelB = b instanceof VelocityVector;

  if ((isPosA && isVelB) || (isVelA && isPosB)) {
    throw new DimensionalityError(
      `Illegal ${op}: Cannot combine SupergalacticPosition (Mpc/h) directly with VelocityVector (km/s). ` +
      `Propagating position by velocity requires an explicit time interval dt with proper unit conversion.`
    );
  }
}

/**
 * Asserts that the official CF4 x52.0 velocity scale factor has not been applied multiple times.
 * @param {Object} field Velocity field or data container
 * @param {boolean} currentAppliedState Whether x52.0 is already recorded as applied
 */
export function assertSingleVelocityScaling(field, currentAppliedState) {
  if (currentAppliedState) {
    throw new InvariantViolationError(
      'Attempted double application of the official CosmicFlows-4 x52.0 velocity scale factor! ' +
      'Rule 3 Invariant: Scale factor must be applied exactly once.'
    );
  }
}

/**
 * Asserts that a spatial position is within the specified bounding box.
 * @param {Array<number>|SupergalacticPosition} pos [x, y, z] in Mpc/h
 * @param {Array<number>} minBound [xMin, yMin, zMin]
 * @param {Array<number>} maxBound [xMax, yMax, zMax]
 * @param {number} [tolerance=1e-6]
 */
export function assertInBoundingBox(pos, minBound, maxBound, tolerance = 1e-6) {
  const x = pos.x ?? pos[0];
  const y = pos.y ?? pos[1];
  const z = pos.z ?? pos[2];

  if (x < minBound[0] - tolerance || x > maxBound[0] + tolerance ||
      y < minBound[1] - tolerance || y > maxBound[1] + tolerance ||
      z < minBound[2] - tolerance || z > maxBound[2] + tolerance) {
    throw new BoundaryExtrapolationError(
      `Position (${x.toFixed(2)}, ${y.toFixed(2)}, ${z.toFixed(2)}) is outside bounding box ` +
      `[${minBound.map(v => v.toFixed(1))}] to [${maxBound.map(v => v.toFixed(1))}].`
    );
  }
}

/**
 * Checks whether a given cosmological region satisfies linear perturbation theory validity.
 * In the linear regime, |delta| << 1 (typically |delta| <= 0.5) and continuity equation holds.
 * @param {number} delta Density contrast
 * @param {number} [maxLinearDelta=0.5] Maximum overdensity threshold for linear regime
 * @returns {boolean} True if region is within linear regime
 */
export function isLinearPerturbationRegime(delta, maxLinearDelta = 0.5) {
  return Math.abs(delta) <= maxLinearDelta;
}

/**
 * Verifies linear continuity diagnostic: div(v) approx -H0 * f * delta.
 * @param {number} divV Velocity divergence in km/s / (Mpc/h)
 * @param {number} delta Density contrast
 * @param {number} H0_f Continuity scale factor (approx 52.0 km/s / (Mpc/h))
 * @param {number} [tolerance=0.35] Relative L2/fractional tolerance allowed in linear regime
 * @returns {{ satisfiesLinear: boolean, expectedDivV: number, actualDivV: number, relativeResidual: number }}
 */
export function evaluateLinearContinuityDiagnostic(divV, delta, H0_f = CF4_PUBLIC_VELOCITY_SCALE, tolerance = 0.35) {
  const expectedDivV = -H0_f * delta;
  const absDiff = Math.abs(divV - expectedDivV);
  const norm = Math.max(1.0, Math.abs(expectedDivV), Math.abs(divV));
  const relativeResidual = absDiff / norm;

  return {
    satisfiesLinear: relativeResidual <= tolerance,
    expectedDivV,
    actualDivV: divV,
    relativeResidual
  };
}
