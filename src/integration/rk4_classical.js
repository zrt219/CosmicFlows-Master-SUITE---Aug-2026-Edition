/**
 * @file rk4_classical.js
 * @description Classical 4th-Order Runge-Kutta (RK4) Integrator for fixed-step streamline tracing.
 * Implements the exact Butcher tableau for 4th-order ODE integration with support for both
 * forward (cosmic flow / attraction) and backward (cosmic repulsion / source tracking) modes.
 *
 * Butcher Tableau:
 *   c | A
 *  ---+-------------
 *   0 |
 * 1/2 | 1/2
 * 1/2 | 0    1/2
 *   1 | 0    0    1
 *  ---+-------------
 *     | 1/6  1/3  1/3  1/6
 *
 * @module integration/rk4_classical
 */

/**
 * Standard termination reason constants for streamline tracing.
 * @readonly
 * @enum {string}
 */
export const TerminationReason = Object.freeze({
  CONVERGED_ENDPOINT: 'CONVERGED_ENDPOINT',
  DOMAIN_EXIT: 'DOMAIN_EXIT',
  MAX_STEPS: 'MAX_STEPS',
  LOW_SPEED_STALL: 'LOW_SPEED_STALL',
  NUMERICAL_FAILURE: 'NUMERICAL_FAILURE',
  IN_PROGRESS: 'IN_PROGRESS'
});

/**
 * Integration direction mode.
 * @readonly
 * @enum {string}
 */
export const IntegrationDirection = Object.freeze({
  FORWARD: 'forward',
  BACKWARD: 'backward'
});

/**
 * Evaluates a velocity field function handling both (t, pos) and (pos) signatures.
 *
 * @param {Function} fieldFn - Vector field evaluator f(t, pos) or f(pos) returning [vx, vy, vz].
 * @param {number} t - Current parameter/time.
 * @param {Array<number>|Float32Array|Float64Array} pos - Current 3D position [x, y, z].
 * @returns {Array<number>} 3D velocity vector [vx, vy, vz].
 * @throws {Error} If evaluation yields non-finite values or invalid dimensions.
 */
export function evaluateField(fieldFn, t, pos) {
  if (typeof fieldFn !== 'function') {
    throw new TypeError('fieldFn must be a callable function');
  }

  let v;
  if (fieldFn.length >= 2) {
    v = fieldFn(t, pos);
  } else {
    v = fieldFn(pos);
  }

  if (!v || v.length < 3) {
    throw new Error(`Field evaluation returned invalid vector: ${JSON.stringify(v)}`);
  }

  const vx = Number(v[0]);
  const vy = Number(v[1]);
  const vz = Number(v[2]);

  if (!Number.isFinite(vx) || !Number.isFinite(vy) || !Number.isFinite(vz)) {
    throw new Error(`Non-finite velocity encountered at pos=[${pos[0]}, ${pos[1]}, ${pos[2]}], t=${t}: [${vx}, ${vy}, ${vz}]`);
  }

  return [vx, vy, vz];
}

/**
 * Validates if a 3D point is within specified domain bounds.
 *
 * @param {Array<number>|Float32Array|Float64Array} pos - [x, y, z] coordinate.
 * @param {Object|Array<Array<number>>|null} bounds - Domain bounds: { xMin, xMax, yMin, yMax, zMin, zMax } or [[xMin, xMax], [yMin, yMax], [zMin, zMax]].
 * @returns {boolean} True if within bounds or bounds is null.
 */
export function isWithinBounds(pos, bounds) {
  if (!bounds) return true;

  let xMin, xMax, yMin, yMax, zMin, zMax;
  if (Array.isArray(bounds) && bounds.length >= 3) {
    [[xMin, xMax], [yMin, yMax], [zMin, zMax]] = bounds;
  } else {
    ({ xMin = -Infinity, xMax = Infinity, yMin = -Infinity, yMax = Infinity, zMin = -Infinity, zMax = Infinity } = bounds);
  }

  const x = pos[0];
  const y = pos[1];
  const z = pos[2];

  return x >= xMin && x <= xMax && y >= yMin && y <= yMax && z >= zMin && z <= zMax;
}

/**
 * Computes a single classical 4th-order Runge-Kutta step.
 *
 * @param {Function} fieldFn - Vector field evaluator f(t, pos) or f(pos).
 * @param {Array<number>|Float32Array|Float64Array} pos - Current 3D position [x, y, z].
 * @param {number} t - Current parameter/time.
 * @param {number} dt - Step size (positive number).
 * @param {string|number} [direction='forward'] - Tracing direction ('forward' | 'backward' | +1 | -1).
 * @returns {{ posNext: Array<number>, tNext: number, velocity: Array<number>, kStages: Array<Array<number>> }} Next state.
 */
export function rk4Step(fieldFn, pos, t, dt, direction = IntegrationDirection.FORWARD) {
  if (dt <= 0 || !Number.isFinite(dt)) {
    throw new RangeError(`Step size dt must be a positive finite number, got ${dt}`);
  }

  const dirSign = (direction === IntegrationDirection.BACKWARD || direction === -1 || direction === '-1') ? -1.0 : 1.0;
  const h = dirSign * dt;

  const x0 = Number(pos[0]);
  const y0 = Number(pos[1]);
  const z0 = Number(pos[2]);

  // Stage 1: k1 = f(t, y0)
  const k1 = evaluateField(fieldFn, t, [x0, y0, z0]);

  // Stage 2: k2 = f(t + 0.5*h, y0 + 0.5*h*k1)
  const pos1 = [
    x0 + 0.5 * h * k1[0],
    y0 + 0.5 * h * k1[1],
    z0 + 0.5 * h * k1[2]
  ];
  const k2 = evaluateField(fieldFn, t + 0.5 * h, pos1);

  // Stage 3: k3 = f(t + 0.5*h, y0 + 0.5*h*k2)
  const pos2 = [
    x0 + 0.5 * h * k2[0],
    y0 + 0.5 * h * k2[1],
    z0 + 0.5 * h * k2[2]
  ];
  const k3 = evaluateField(fieldFn, t + 0.5 * h, pos2);

  // Stage 4: k4 = f(t + h, y0 + h*k3)
  const pos3 = [
    x0 + h * k3[0],
    y0 + h * k3[1],
    z0 + h * k3[2]
  ];
  const k4 = evaluateField(fieldFn, t + h, pos3);

  // Final 4th-order weighted combination: yNext = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
  const hSixth = h / 6.0;
  const xNext = x0 + hSixth * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
  const yNext = y0 + hSixth * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
  const zNext = z0 + hSixth * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);

  return {
    posNext: [xNext, yNext, zNext],
    tNext: t + h,
    velocity: k1,
    kStages: [k1, k2, k3, k4]
  };
}

/**
 * Traces a streamline using classical fixed-step RK4 integration.
 *
 * @param {Function} fieldFn - Vector field evaluator f(t, pos) or f(pos).
 * @param {Array<number>|Float32Array|Float64Array} seedPos - Seed coordinate [x0, y0, z0].
 * @param {Object} [options={}] - Configuration options.
 * @param {number} [options.dt=0.5] - Fixed integration step size.
 * @param {string|number} [options.direction='forward'] - 'forward' (attraction) or 'backward' (repulsion).
 * @param {number} [options.maxSteps=2000] - Maximum integration step limit.
 * @param {number} [options.tStart=0.0] - Starting time/parameter.
 * @param {Object|Array<Array<number>>|null} [options.bounds=null] - Bounding domain.
 * @param {number} [options.minVelocity=1e-6] - Minimum velocity threshold for endpoint convergence.
 * @param {number} [options.stallWindow=10] - Window size for detecting low-speed stagnation.
 * @param {number} [options.minArcLengthStep=1e-5] - Minimum displacement per step to avoid stall.
 * @returns {Object} Structured streamline trajectory result.
 */
export function traceRK4Streamline(fieldFn, seedPos, options = {}) {
  const {
    dt = 0.5,
    direction = IntegrationDirection.FORWARD,
    maxSteps = 2000,
    tStart = 0.0,
    bounds = null,
    minVelocity = 1e-6,
    stallWindow = 10,
    minArcLengthStep = 1e-5
  } = options;

  if (!seedPos || seedPos.length < 3) {
    throw new Error(`Invalid seed position: ${JSON.stringify(seedPos)}`);
  }

  const positions = [];
  const velocities = [];
  const speeds = [];
  const arcLengths = [0.0];
  const times = [tStart];

  let currentPos = [Number(seedPos[0]), Number(seedPos[1]), Number(seedPos[2])];
  let currentT = tStart;
  let totalArcLength = 0.0;
  let terminationReason = TerminationReason.MAX_STEPS;
  let stepCount = 0;

  if (!isWithinBounds(currentPos, bounds)) {
    return {
      seedPoint: [...currentPos],
      direction: (direction === IntegrationDirection.BACKWARD || direction === -1) ? 'backward' : 'forward',
      points: [currentPos],
      velocities: [],
      speeds: [],
      arcLengths: [0],
      times: [tStart],
      terminationReason: TerminationReason.DOMAIN_EXIT,
      stepCount: 0,
      stats: { totalArcLength: 0, acceptedSteps: 0, rejectedSteps: 0, minSpeed: 0, maxSpeed: 0, meanSpeed: 0 }
    };
  }

  positions.push([...currentPos]);

  // Initial field evaluation at seed
  try {
    const v0 = evaluateField(fieldFn, currentT, currentPos);
    const s0 = Math.hypot(v0[0], v0[1], v0[2]);
    velocities.push(v0);
    speeds.push(s0);

    if (s0 < minVelocity) {
      return {
        seedPoint: [...seedPos],
        direction: (direction === IntegrationDirection.BACKWARD || direction === -1) ? 'backward' : 'forward',
        points: positions,
        velocities,
        speeds,
        arcLengths,
        times,
        terminationReason: TerminationReason.CONVERGED_ENDPOINT,
        stepCount: 0,
        stats: { totalArcLength: 0, acceptedSteps: 0, rejectedSteps: 0, minSpeed: s0, maxSpeed: s0, meanSpeed: s0 }
      };
    }
  } catch (err) {
    return {
      seedPoint: [...seedPos],
      direction: (direction === IntegrationDirection.BACKWARD || direction === -1) ? 'backward' : 'forward',
      points: positions,
      velocities,
      speeds,
      arcLengths,
      times,
      terminationReason: TerminationReason.NUMERICAL_FAILURE,
      stepCount: 0,
      stats: { totalArcLength: 0, acceptedSteps: 0, rejectedSteps: 0, minSpeed: 0, maxSpeed: 0, meanSpeed: 0, error: err.message }
    };
  }

  // Integration loop
  for (let step = 0; step < maxSteps; step++) {
    let stepResult;
    try {
      stepResult = rk4Step(fieldFn, currentPos, currentT, dt, direction);
    } catch (err) {
      terminationReason = TerminationReason.NUMERICAL_FAILURE;
      break;
    }

    const nextPos = stepResult.posNext;
    const nextT = stepResult.tNext;
    const v = stepResult.velocity;
    const speed = Math.hypot(v[0], v[1], v[2]);

    const dx = nextPos[0] - currentPos[0];
    const dy = nextPos[1] - currentPos[1];
    const dz = nextPos[2] - currentPos[2];
    const stepDist = Math.hypot(dx, dy, dz);

    stepCount++;

    // Check Domain Exit
    if (!isWithinBounds(nextPos, bounds)) {
      positions.push(nextPos);
      velocities.push(v);
      speeds.push(speed);
      totalArcLength += stepDist;
      arcLengths.push(totalArcLength);
      times.push(nextT);
      terminationReason = TerminationReason.DOMAIN_EXIT;
      break;
    }

    currentPos = nextPos;
    currentT = nextT;
    totalArcLength += stepDist;

    positions.push([...currentPos]);
    velocities.push(v);
    speeds.push(speed);
    arcLengths.push(totalArcLength);
    times.push(currentT);

    // Check Convergence / Endpoint Sink
    if (speed < minVelocity) {
      terminationReason = TerminationReason.CONVERGED_ENDPOINT;
      break;
    }

    // Check Stall / Stagnation
    if (stepDist < minArcLengthStep) {
      terminationReason = TerminationReason.LOW_SPEED_STALL;
      break;
    }

    if (step >= stallWindow) {
      const recentDist = totalArcLength - arcLengths[arcLengths.length - 1 - stallWindow];
      if (recentDist < minArcLengthStep * stallWindow) {
        terminationReason = TerminationReason.LOW_SPEED_STALL;
        break;
      }
    }
  }

  // Calculate statistics
  let minSpeed = Infinity;
  let maxSpeed = -Infinity;
  let sumSpeed = 0.0;
  for (let i = 0; i < speeds.length; i++) {
    const s = speeds[i];
    if (s < minSpeed) minSpeed = s;
    if (s > maxSpeed) maxSpeed = s;
    sumSpeed += s;
  }
  const meanSpeed = speeds.length > 0 ? sumSpeed / speeds.length : 0.0;

  return {
    seedPoint: [Number(seedPos[0]), Number(seedPos[1]), Number(seedPos[2])],
    direction: (direction === IntegrationDirection.BACKWARD || direction === -1) ? 'backward' : 'forward',
    points: positions,
    velocities,
    speeds,
    arcLengths,
    times,
    terminationReason,
    stepCount,
    stats: {
      totalArcLength,
      acceptedSteps: stepCount,
      rejectedSteps: 0,
      minSpeed: Number.isFinite(minSpeed) ? minSpeed : 0.0,
      maxSpeed: Number.isFinite(maxSpeed) ? maxSpeed : 0.0,
      meanSpeed
    }
  };
}

/**
 * Classical RK4 Integrator Class with configurable defaults.
 */
export class RK4Integrator {
  /**
   * @param {Object} [config={}]
   * @param {number} [config.dt=0.5] - Default step size.
   * @param {string} [config.direction='forward'] - Default direction.
   * @param {number} [config.maxSteps=2000] - Default max step count.
   * @param {Object|null} [config.bounds=null] - Default domain bounds.
   * @param {number} [config.minVelocity=1e-6] - Convergence velocity threshold.
   */
  constructor(config = {}) {
    this.dt = config.dt ?? 0.5;
    this.direction = config.direction ?? IntegrationDirection.FORWARD;
    this.maxSteps = config.maxSteps ?? 2000;
    this.bounds = config.bounds ?? null;
    this.minVelocity = config.minVelocity ?? 1e-6;
    this.stallWindow = config.stallWindow ?? 10;
    this.minArcLengthStep = config.minArcLengthStep ?? 1e-5;
  }

  /**
   * Performs a single RK4 step.
   * @param {Function} fieldFn
   * @param {Array<number>} pos
   * @param {number} t
   * @param {number} [dt]
   * @param {string} [direction]
   */
  step(fieldFn, pos, t, dt = this.dt, direction = this.direction) {
    return rk4Step(fieldFn, pos, t, dt, direction);
  }

  /**
   * Traces a streamline starting from seed position.
   * @param {Function} fieldFn
   * @param {Array<number>} seedPos
   * @param {Object} [overrideOptions={}]
   */
  integrate(fieldFn, seedPos, overrideOptions = {}) {
    return this.trace(fieldFn, seedPos, overrideOptions);
  }

  trace(fieldFn, seedPos, overrideOptions = {}) {
    const opts = {
      dt: this.dt,
      direction: this.direction,
      maxSteps: this.maxSteps,
      bounds: this.bounds,
      minVelocity: this.minVelocity,
      stallWindow: this.stallWindow,
      minArcLengthStep: this.minArcLengthStep,
      ...overrideOptions
    };
    return traceRK4Streamline(fieldFn, seedPos, opts);
  }
}

export default RK4Integrator;
