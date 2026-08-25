/**
 * @file rk45_cash_karp.js
 * @description Adaptive Cash-Karp 6-stage embedded 4(5) Runge-Kutta integrator.
 * Implements genuine local truncation error estimation |y5 - y4|, automatic step size
 * adaptation, step rejection/retry logic, and tolerance checks for streamline tracing.
 *
 * References:
 * - Cash, J. R., & Karp, A. H. (1990). "A variable order Runge-Kutta method for initial
 *   value problems with rapidly varying right-hand sides". ACM TOMS 16(3), 201-222.
 * - Press et al., Numerical Recipes in C/C++, Section 16.2.
 *
 * @module integration/rk45_cash_karp
 */

import { TerminationReason, IntegrationDirection, evaluateField, isWithinBounds } from './rk4_classical.js';

/**
 * Cash-Karp Butcher Tableau Constants
 */
export const CASH_KARP_COEFFS = Object.freeze({
  // Nodes c_i
  c: [0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0],

  // Runge-Kutta matrix a_ij
  a21: 1.0 / 5.0,

  a31: 3.0 / 40.0,
  a32: 9.0 / 40.0,

  a41: 3.0 / 10.0,
  a42: -9.0 / 10.0,
  a43: 6.0 / 5.0,

  a51: -11.0 / 54.0,
  a52: 5.0 / 2.0,
  a53: -70.0 / 27.0,
  a54: 35.0 / 27.0,

  a61: 1631.0 / 55296.0,
  a62: 175.0 / 512.0,
  a63: 575.0 / 13824.0,
  a64: 44275.0 / 110592.0,
  a65: 253.0 / 4096.0,

  // 5th-order weights b_i
  b1: 37.0 / 378.0,
  b2: 0.0,
  b3: 250.0 / 621.0,
  b4: 125.0 / 594.0,
  b5: 0.0,
  b6: 512.0 / 1771.0,

  // 4th-order weights b*_i
  bStar1: 2825.0 / 27648.0,
  bStar2: 0.0,
  bStar3: 18575.0 / 48384.0,
  bStar4: 13525.0 / 55296.0,
  bStar5: 277.0 / 14336.0,
  bStar6: 1.0 / 4.0,

  // Exact difference weights e_i = b_i - bStar_i for local error estimate
  e1: 37.0 / 378.0 - 2825.0 / 27648.0,       // -277 / 64512
  e2: 0.0,
  e3: 250.0 / 621.0 - 18575.0 / 48384.0,     // 2825 / 72576
  e4: 125.0 / 594.0 - 13525.0 / 55296.0,     // -18575 / 770048
  e5: -277.0 / 14336.0,                      // -277 / 14336
  e6: 512.0 / 1771.0 - 1.0 / 4.0             // 277 / 7084
});

/**
 * Computes a single Cash-Karp RK4(5) step with local truncation error estimation.
 *
 * @param {Function} fieldFn - Vector field evaluator f(t, pos).
 * @param {Array<number>} pos - Current 3D position [x, y, z].
 * @param {number} t - Current parameter/time.
 * @param {number} hSigned - Signed step size (positive for forward, negative for backward).
 * @returns {{ pos5: Array<number>, pos4: Array<number>, errVec: Array<number>, kStages: Array<Array<number>> }} Step outputs.
 */
export function cashKarpStepRaw(fieldFn, pos, t, hSigned) {
  const C = CASH_KARP_COEFFS;
  const x0 = Number(pos[0]);
  const y0 = Number(pos[1]);
  const z0 = Number(pos[2]);

  // Stage 1: k1 = f(t, y0)
  const k1 = evaluateField(fieldFn, t, [x0, y0, z0]);

  // Stage 2: k2 = f(t + c2*h, y0 + h*a21*k1)
  const p1 = [
    x0 + hSigned * C.a21 * k1[0],
    y0 + hSigned * C.a21 * k1[1],
    z0 + hSigned * C.a21 * k1[2]
  ];
  const k2 = evaluateField(fieldFn, t + C.c[1] * hSigned, p1);

  // Stage 3: k3 = f(t + c3*h, y0 + h*(a31*k1 + a32*k2))
  const p2 = [
    x0 + hSigned * (C.a31 * k1[0] + C.a32 * k2[0]),
    y0 + hSigned * (C.a31 * k1[1] + C.a32 * k2[1]),
    z0 + hSigned * (C.a31 * k1[2] + C.a32 * k2[2])
  ];
  const k3 = evaluateField(fieldFn, t + C.c[2] * hSigned, p2);

  // Stage 4: k4 = f(t + c4*h, y0 + h*(a41*k1 + a42*k2 + a43*k3))
  const p3 = [
    x0 + hSigned * (C.a41 * k1[0] + C.a42 * k2[0] + C.a43 * k3[0]),
    y0 + hSigned * (C.a41 * k1[1] + C.a42 * k2[1] + C.a43 * k3[1]),
    z0 + hSigned * (C.a41 * k1[2] + C.a42 * k2[2] + C.a43 * k3[2])
  ];
  const k4 = evaluateField(fieldFn, t + C.c[3] * hSigned, p3);

  // Stage 5: k5 = f(t + c5*h, y0 + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
  const p4 = [
    x0 + hSigned * (C.a51 * k1[0] + C.a52 * k2[0] + C.a53 * k3[0] + C.a54 * k4[0]),
    y0 + hSigned * (C.a51 * k1[1] + C.a52 * k2[1] + C.a53 * k3[1] + C.a54 * k4[1]),
    z0 + hSigned * (C.a51 * k1[2] + C.a52 * k2[2] + C.a53 * k3[2] + C.a54 * k4[2])
  ];
  const k5 = evaluateField(fieldFn, t + C.c[4] * hSigned, p4);

  // Stage 6: k6 = f(t + c6*h, y0 + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))
  const p5 = [
    x0 + hSigned * (C.a61 * k1[0] + C.a62 * k2[0] + C.a63 * k3[0] + C.a64 * k4[0] + C.a65 * k5[0]),
    y0 + hSigned * (C.a61 * k1[1] + C.a62 * k2[1] + C.a63 * k3[1] + C.a64 * k4[1] + C.a65 * k5[1]),
    z0 + hSigned * (C.a61 * k1[2] + C.a62 * k2[2] + C.a63 * k3[2] + C.a64 * k4[2] + C.a65 * k5[2])
  ];
  const k6 = evaluateField(fieldFn, t + C.c[5] * hSigned, p5);

  // 5th-order candidate position y5
  const x5 = x0 + hSigned * (C.b1 * k1[0] + C.b3 * k3[0] + C.b4 * k4[0] + C.b6 * k6[0]);
  const y5 = y0 + hSigned * (C.b1 * k1[1] + C.b3 * k3[1] + C.b4 * k4[1] + C.b6 * k6[1]);
  const z5 = z0 + hSigned * (C.b1 * k1[2] + C.b3 * k3[2] + C.b4 * k4[2] + C.b6 * k6[2]);

  // 4th-order candidate position y4
  const x4 = x0 + hSigned * (C.bStar1 * k1[0] + C.bStar3 * k3[0] + C.bStar4 * k4[0] + C.bStar5 * k5[0] + C.bStar6 * k6[0]);
  const y4 = y0 + hSigned * (C.bStar1 * k1[1] + C.bStar3 * k3[1] + C.bStar4 * k4[1] + C.bStar5 * k5[1] + C.bStar6 * k6[1]);
  const z4 = z0 + hSigned * (C.bStar1 * k1[2] + C.bStar3 * k3[2] + C.bStar4 * k4[2] + C.bStar5 * k5[2] + C.bStar6 * k6[2]);

  // Local truncation error vector E = y5 - y4
  const errX = hSigned * (C.e1 * k1[0] + C.e3 * k3[0] + C.e4 * k4[0] + C.e5 * k5[0] + C.e6 * k6[0]);
  const errY = hSigned * (C.e1 * k1[1] + C.e3 * k3[1] + C.e4 * k4[1] + C.e5 * k5[1] + C.e6 * k6[1]);
  const errZ = hSigned * (C.e1 * k1[2] + C.e3 * k3[2] + C.e4 * k4[2] + C.e5 * k5[2] + C.e6 * k6[2]);

  return {
    pos5: [x5, y5, z5],
    pos4: [x4, y4, z4],
    errVec: [errX, errY, errZ],
    kStages: [k1, k2, k3, k4, k5, k6]
  };
}

/**
 * Computes scaled error norm for 3D vector ODE.
 *
 * @param {Array<number>} y0 - Starting position.
 * @param {Array<number>} y1 - Higher order candidate position.
 * @param {Array<number>} errVec - Local truncation error vector.
 * @param {number} atol - Absolute tolerance.
 * @param {number} rtol - Relative tolerance.
 * @returns {number} Normalized error norm (<= 1.0 means step meets tolerance).
 */
export function computeScaledErrorNorm(y0, y1, errVec, atol = 1e-6, rtol = 1e-6) {
  let sumSq = 0.0;
  for (let i = 0; i < 3; i++) {
    const scale = atol + Math.max(Math.abs(y0[i]), Math.abs(y1[i])) * rtol;
    const ratio = errVec[i] / scale;
    sumSq += ratio * ratio;
  }
  return Math.sqrt(sumSq / 3.0);
}

/**
 * Computes an adaptive Cash-Karp step with rejection/retry loop and next step suggestion.
 *
 * @param {Function} fieldFn - Vector field evaluator.
 * @param {Array<number>} pos - Current 3D position [x, y, z].
 * @param {number} t - Current parameter/time.
 * @param {number} dt - Proposed step size magnitude.
 * @param {Object} [options={}] - Adaptive parameters.
 * @returns {Object} Step execution result with acceptance status, next step, and error norm.
 */
export function rk45CashKarpStep(fieldFn, pos, t, dt, options = {}) {
  const {
    direction = IntegrationDirection.FORWARD,
    atol = 1e-6,
    rtol = 1e-6,
    minStep = 1e-6,
    maxStep = 10.0,
    safetyFactor = 0.90,
    facMin = 0.2,
    facMax = 5.0,
    maxRetries = 12
  } = options;

  const dirSign = (direction === IntegrationDirection.BACKWARD || direction === -1 || direction === '-1') ? -1.0 : 1.0;
  let currentDt = Math.max(minStep, Math.min(maxStep, Math.abs(dt)));
  let attempts = 0;

  while (attempts < maxRetries) {
    attempts++;
    const hSigned = dirSign * currentDt;
    const raw = cashKarpStepRaw(fieldFn, pos, t, hSigned);
    const errNorm = computeScaledErrorNorm(pos, raw.pos5, raw.errVec, atol, rtol);

    if (errNorm <= 1.0) {
      // Step accepted!
      let nextDt;
      if (errNorm < 1e-10) {
        nextDt = currentDt * facMax;
      } else {
        const factor = safetyFactor * Math.pow(1.0 / errNorm, 0.2);
        nextDt = currentDt * Math.min(facMax, Math.max(facMin, factor));
      }
      nextDt = Math.max(minStep, Math.min(maxStep, nextDt));

      return {
        accepted: true,
        posNext: raw.pos5,
        tNext: t + hSigned,
        dtUsed: currentDt,
        dtNext: nextDt,
        errNorm,
        errVec: raw.errVec,
        velocity: raw.kStages[0],
        retries: attempts - 1
      };
    }

    // Step rejected - shrink step size
    const shrinkFactor = Math.max(facMin, safetyFactor * Math.pow(1.0 / errNorm, 0.25));
    const retryDt = currentDt * shrinkFactor;

    if (retryDt < minStep) {
      return {
        accepted: false,
        posNext: pos,
        tNext: t,
        dtUsed: currentDt,
        dtNext: minStep,
        errNorm,
        errVec: raw.errVec,
        velocity: raw.kStages[0],
        retries: attempts,
        reason: 'STEP_UNDERFLOW'
      };
    }

    currentDt = retryDt;
  }

  // Max retries exceeded
  return {
    accepted: false,
    posNext: pos,
    tNext: t,
    dtUsed: currentDt,
    dtNext: minStep,
    errNorm: Infinity,
    retries: attempts,
    reason: 'MAX_RETRIES_EXCEEDED'
  };
}

/**
 * Traces a streamline using adaptive Cash-Karp RK4(5) integration.
 *
 * @param {Function} fieldFn - Vector field evaluator f(t, pos).
 * @param {Array<number>} seedPos - Seed coordinate [x0, y0, z0].
 * @param {Object} [options={}] - Configuration options.
 * @returns {Object} Structured streamline trajectory result.
 */
export function traceCashKarpStreamline(fieldFn, seedPos, options = {}) {
  const {
    initialStep = 0.5,
    direction = IntegrationDirection.FORWARD,
    maxSteps = 2000,
    tStart = 0.0,
    bounds = null,
    atol = 1e-6,
    rtol = 1e-6,
    minStep = 1e-6,
    maxStep = 10.0,
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
  let currentDt = initialStep;
  let totalArcLength = 0.0;
  let terminationReason = TerminationReason.MAX_STEPS;
  let acceptedSteps = 0;
  let rejectedSteps = 0;

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

  for (let step = 0; step < maxSteps; step++) {
    let result;
    try {
      result = rk45CashKarpStep(fieldFn, currentPos, currentT, currentDt, {
        direction,
        atol,
        rtol,
        minStep,
        maxStep
      });
    } catch (err) {
      terminationReason = TerminationReason.NUMERICAL_FAILURE;
      break;
    }

    rejectedSteps += result.retries;

    if (!result.accepted) {
      terminationReason = TerminationReason.NUMERICAL_FAILURE;
      break;
    }

    acceptedSteps++;
    const nextPos = result.posNext;
    const nextT = result.tNext;
    const v = result.velocity;
    const speed = Math.hypot(v[0], v[1], v[2]);

    const dx = nextPos[0] - currentPos[0];
    const dy = nextPos[1] - currentPos[1];
    const dz = nextPos[2] - currentPos[2];
    const stepDist = Math.hypot(dx, dy, dz);

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
    currentDt = result.dtNext;
    totalArcLength += stepDist;

    positions.push([...currentPos]);
    velocities.push(v);
    speeds.push(speed);
    arcLengths.push(totalArcLength);
    times.push(currentT);

    if (speed < minVelocity) {
      terminationReason = TerminationReason.CONVERGED_ENDPOINT;
      break;
    }

    if (stepDist < minArcLengthStep) {
      terminationReason = TerminationReason.LOW_SPEED_STALL;
      break;
    }

    if (positions.length > stallWindow) {
      const recentDist = totalArcLength - arcLengths[arcLengths.length - 1 - stallWindow];
      if (recentDist < minArcLengthStep * stallWindow) {
        terminationReason = TerminationReason.LOW_SPEED_STALL;
        break;
      }
    }
  }

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
    stepCount: acceptedSteps,
    stats: {
      totalArcLength,
      acceptedSteps,
      rejectedSteps,
      minSpeed: Number.isFinite(minSpeed) ? minSpeed : 0.0,
      maxSpeed: Number.isFinite(maxSpeed) ? maxSpeed : 0.0,
      meanSpeed
    }
  };
}

/**
 * Adaptive Cash-Karp 4(5) Integrator Class.
 */
export class RK45CashKarpIntegrator {
  /**
   * @param {Object} [config={}]
   */
  constructor(config = {}) {
    this.initialStep = config.initialStep ?? 0.5;
    this.direction = config.direction ?? IntegrationDirection.FORWARD;
    this.maxSteps = config.maxSteps ?? 2000;
    this.bounds = config.bounds ?? null;
    this.atol = config.atol ?? 1e-6;
    this.rtol = config.rtol ?? 1e-6;
    this.minStep = config.minStep ?? 1e-6;
    this.maxStep = config.maxStep ?? 10.0;
    this.minVelocity = config.minVelocity ?? 1e-6;
    this.stallWindow = config.stallWindow ?? 10;
    this.minArcLengthStep = config.minArcLengthStep ?? 1e-5;
  }

  step(fieldFn, pos, t, dt = this.initialStep, options = {}) {
    return rk45CashKarpStep(fieldFn, pos, t, dt, {
      direction: this.direction,
      atol: this.atol,
      rtol: this.rtol,
      minStep: this.minStep,
      maxStep: this.maxStep,
      ...options
    });
  }

  integrate(fieldFn, seedPos, overrideOptions = {}) {
    return this.trace(fieldFn, seedPos, overrideOptions);
  }

  trace(fieldFn, seedPos, overrideOptions = {}) {
    const opts = {
      initialStep: this.initialStep,
      direction: this.direction,
      maxSteps: this.maxSteps,
      bounds: this.bounds,
      atol: this.atol,
      rtol: this.rtol,
      minStep: this.minStep,
      maxStep: this.maxStep,
      minVelocity: this.minVelocity,
      stallWindow: this.stallWindow,
      minArcLengthStep: this.minArcLengthStep,
      ...overrideOptions
    };
    return traceCashKarpStreamline(fieldFn, seedPos, opts);
  }
}

export default RK45CashKarpIntegrator;
