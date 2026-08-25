/**
 * @file dormand_prince.js
 * @description Dormand-Prince 5(4) (DOPRI5) Adaptive Integrator with FSAL
 * (First Same As Last) property and Hairer-Wanner 4th-order continuous dense output.
 *
 * References:
 * - Dormand, J. R., & Prince, P. J. (1980). "A family of embedded Runge-Kutta formulae".
 *   Journal of Computational and Applied Mathematics, 6(1), 19-26.
 * - Hairer, E., Nørsett, S. P., & Wanner, G. (1993). "Solving Ordinary Differential Equations I:
 *   Nonstiff Problems". Springer Series in Computational Mathematics, Vol. 8.
 *
 * @module integration/dormand_prince
 */

import { TerminationReason, IntegrationDirection, evaluateField, isWithinBounds } from './rk4_classical.js';
import { computeScaledErrorNorm } from './rk45_cash_karp.js';

/**
 * Dormand-Prince 5(4) Butcher Tableau and Continuous Output Constants
 */
export const DOPRI5_COEFFS = Object.freeze({
  // Nodes c_i
  c: [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0],

  // Runge-Kutta matrix a_ij
  a21: 1.0 / 5.0,

  a31: 3.0 / 40.0,
  a32: 9.0 / 40.0,

  a41: 44.0 / 45.0,
  a42: -56.0 / 15.0,
  a43: 32.0 / 9.0,

  a51: 19372.0 / 6561.0,
  a52: -25360.0 / 2187.0,
  a53: 64448.0 / 6561.0,
  a54: -212.0 / 729.0,

  a61: 9017.0 / 3168.0,
  a62: -355.0 / 33.0,
  a63: 46732.0 / 5247.0,
  a64: 49.0 / 176.0,
  a65: -5103.0 / 18656.0,

  a71: 35.0 / 384.0,
  a72: 0.0,
  a73: 500.0 / 1113.0,
  a74: 125.0 / 192.0,
  a75: -2187.0 / 6784.0,
  a76: 11.0 / 84.0,

  // 5th-order weights b_i
  b1: 35.0 / 384.0,
  b2: 0.0,
  b3: 500.0 / 1113.0,
  b4: 125.0 / 192.0,
  b5: -2187.0 / 6784.0,
  b6: 11.0 / 84.0,
  b7: 0.0,

  // 4th-order weights bHat_i
  bHat1: 5179.0 / 57600.0,
  bHat2: 0.0,
  bHat3: 7571.0 / 16695.0,
  bHat4: 393.0 / 640.0,
  bHat5: -92097.0 / 339200.0,
  bHat6: 187.0 / 2100.0,
  bHat7: 1.0 / 40.0,

  // Exact difference weights e_i = b_i - bHat_i
  e1: 71.0 / 57600.0,
  e2: 0.0,
  e3: -71.0 / 16695.0,
  e4: 71.0 / 1920.0,
  e5: -17253.0 / 339200.0,
  e6: 22.0 / 525.0,
  e7: -1.0 / 40.0,

  // Dense continuous output coefficients (Hairer-Wanner r5 formula)
  d1: -12715105075.0 / 11282082432.0,
  d3: 87487479700.0 / 32700010779.0,
  d4: -10690763385.0 / 1880347072.0,
  d5: 701980252875.0 / 199316789632.0,
  d6: -1453857185.0 / 822651844.0,
  d7: 69997945.0 / 29380423.0
});

/**
 * Computes a single Dormand-Prince step with FSAL acceleration.
 *
 * @param {Function} fieldFn - Vector field evaluator f(t, pos).
 * @param {Array<number>} pos - Current 3D position [x, y, z].
 * @param {number} t - Current parameter/time.
 * @param {number} hSigned - Signed step size.
 * @param {Array<number>|null} [k1Cached=null] - Pre-evaluated k1 stage from previous accepted step's k7.
 * @returns {{ pos5: Array<number>, pos4: Array<number>, errVec: Array<number>, kStages: Array<Array<number>>, fsalStage: Array<number> }}
 */
export function dormandPrinceStepRaw(fieldFn, pos, t, hSigned, k1Cached = null) {
  const D = DOPRI5_COEFFS;
  const x0 = Number(pos[0]);
  const y0 = Number(pos[1]);
  const z0 = Number(pos[2]);

  // Stage 1: FSAL reuse if valid, else evaluate f(t, y0)
  const k1 = (k1Cached && k1Cached.length >= 3)
    ? [Number(k1Cached[0]), Number(k1Cached[1]), Number(k1Cached[2])]
    : evaluateField(fieldFn, t, [x0, y0, z0]);

  // Stage 2: k2 = f(t + c2*h, y0 + h*a21*k1)
  const p1 = [
    x0 + hSigned * D.a21 * k1[0],
    y0 + hSigned * D.a21 * k1[1],
    z0 + hSigned * D.a21 * k1[2]
  ];
  const k2 = evaluateField(fieldFn, t + D.c[1] * hSigned, p1);

  // Stage 3: k3 = f(t + c3*h, y0 + h*(a31*k1 + a32*k2))
  const p2 = [
    x0 + hSigned * (D.a31 * k1[0] + D.a32 * k2[0]),
    y0 + hSigned * (D.a31 * k1[1] + D.a32 * k2[1]),
    z0 + hSigned * (D.a31 * k1[2] + D.a32 * k2[2])
  ];
  const k3 = evaluateField(fieldFn, t + D.c[2] * hSigned, p2);

  // Stage 4: k4 = f(t + c4*h, y0 + h*(a41*k1 + a42*k2 + a43*k3))
  const p3 = [
    x0 + hSigned * (D.a41 * k1[0] + D.a42 * k2[0] + D.a43 * k3[0]),
    y0 + hSigned * (D.a41 * k1[1] + D.a42 * k2[1] + D.a43 * k3[1]),
    z0 + hSigned * (D.a41 * k1[2] + D.a42 * k2[2] + D.a43 * k3[2])
  ];
  const k4 = evaluateField(fieldFn, t + D.c[3] * hSigned, p3);

  // Stage 5: k5 = f(t + c5*h, y0 + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
  const p4 = [
    x0 + hSigned * (D.a51 * k1[0] + D.a52 * k2[0] + D.a53 * k3[0] + D.a54 * k4[0]),
    y0 + hSigned * (D.a51 * k1[1] + D.a52 * k2[1] + D.a53 * k3[1] + D.a54 * k4[1]),
    z0 + hSigned * (D.a51 * k1[2] + D.a52 * k2[2] + D.a53 * k3[2] + D.a54 * k4[2])
  ];
  const k5 = evaluateField(fieldFn, t + D.c[4] * hSigned, p4);

  // Stage 6: k6 = f(t + c6*h, y0 + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))
  const p5 = [
    x0 + hSigned * (D.a61 * k1[0] + D.a62 * k2[0] + D.a63 * k3[0] + D.a64 * k4[0] + D.a65 * k5[0]),
    y0 + hSigned * (D.a61 * k1[1] + D.a62 * k2[1] + D.a63 * k3[1] + D.a64 * k4[1] + D.a65 * k5[1]),
    z0 + hSigned * (D.a61 * k1[2] + D.a62 * k2[2] + D.a63 * k3[2] + D.a64 * k4[2] + D.a65 * k5[2])
  ];
  const k6 = evaluateField(fieldFn, t + D.c[5] * hSigned, p5);

  // 5th-order solution: note that a7j == bj for j=1..6
  const x5 = x0 + hSigned * (D.b1 * k1[0] + D.b3 * k3[0] + D.b4 * k4[0] + D.b5 * k5[0] + D.b6 * k6[0]);
  const y5 = y0 + hSigned * (D.b1 * k1[1] + D.b3 * k3[1] + D.b4 * k4[1] + D.b5 * k5[1] + D.b6 * k6[1]);
  const z5 = z0 + hSigned * (D.b1 * k1[2] + D.b3 * k3[2] + D.b4 * k4[2] + D.b5 * k5[2] + D.b6 * k6[2]);
  const pos5 = [x5, y5, z5];

  // Stage 7 (FSAL): evaluate f(t + h, pos5)
  const k7 = evaluateField(fieldFn, t + hSigned, pos5);

  // 4th-order candidate position yHat
  const x4 = x0 + hSigned * (D.bHat1 * k1[0] + D.bHat3 * k3[0] + D.bHat4 * k4[0] + D.bHat5 * k5[0] + D.bHat6 * k6[0] + D.bHat7 * k7[0]);
  const y4 = y0 + hSigned * (D.bHat1 * k1[1] + D.bHat3 * k3[1] + D.bHat4 * k4[1] + D.bHat5 * k5[1] + D.bHat6 * k6[1] + D.bHat7 * k7[1]);
  const z4 = z0 + hSigned * (D.bHat1 * k1[2] + D.bHat3 * k3[2] + D.bHat4 * k4[2] + D.bHat5 * k5[2] + D.bHat6 * k6[2] + D.bHat7 * k7[2]);
  const pos4 = [x4, y4, z4];

  // Local truncation error vector E = pos5 - pos4
  const errX = hSigned * (D.e1 * k1[0] + D.e3 * k3[0] + D.e4 * k4[0] + D.e5 * k5[0] + D.e6 * k6[0] + D.e7 * k7[0]);
  const errY = hSigned * (D.e1 * k1[1] + D.e3 * k3[1] + D.e4 * k4[1] + D.e5 * k5[1] + D.e6 * k6[1] + D.e7 * k7[1]);
  const errZ = hSigned * (D.e1 * k1[2] + D.e3 * k3[2] + D.e4 * k4[2] + D.e5 * k5[2] + D.e6 * k6[2] + D.e7 * k7[2]);

  return {
    pos5,
    pos4,
    errVec: [errX, errY, errZ],
    kStages: [k1, k2, k3, k4, k5, k6, k7],
    fsalStage: k7
  };
}

/**
 * Computes 4th-order dense continuous output between t and t+h for DOPRI5.
 *
 * @param {Array<number>} y0 - Step start position.
 * @param {Array<number>} y1 - Step end position (pos5).
 * @param {number} hSigned - Signed step size.
 * @param {Array<Array<number>>} kStages - Stages k1 through k7 from DOPRI5 step.
 * @param {number} theta - Normalized parameter along step, theta in [0, 1].
 * @returns {Array<number>} Interpolated 3D position [x, y, z].
 */
export function denseInterpolateDOPRI5(y0, y1, hSigned, kStages, theta) {
  if (theta <= 0.0) return [y0[0], y0[1], y0[2]];
  if (theta >= 1.0) return [y1[0], y1[1], y1[2]];

  const D = DOPRI5_COEFFS;
  const [k1, , k3, k4, k5, k6, k7] = kStages;

  const out = [0, 0, 0];
  const thetaComp = 1.0 - theta;

  for (let i = 0; i < 3; i++) {
    const y0i = y0[i];
    const y1i = y1[i];
    const dy = y1i - y0i;

    const r1 = y0i;
    const r2 = dy;
    const r3 = hSigned * k1[i] - dy;
    const r4 = dy - hSigned * k7[i] - r3;

    const r5 = hSigned * (
      D.d1 * k1[i] +
      D.d3 * k3[i] +
      D.d4 * k4[i] +
      D.d5 * k5[i] +
      D.d6 * k6[i] +
      D.d7 * k7[i]
    );

    // Continuous Hermite-Birkhoff polynomial
    out[i] = r1 + theta * (r2 + thetaComp * (r3 + theta * (r4 + thetaComp * r5)));
  }

  return out;
}

/**
 * Computes an adaptive Dormand-Prince step with FSAL reuse and step adaptation.
 *
 * @param {Function} fieldFn
 * @param {Array<number>} pos
 * @param {number} t
 * @param {number} dt
 * @param {Object} [options={}]
 * @param {Array<number>|null} [k1Cached=null]
 * @returns {Object}
 */
export function dormandPrinceStep(fieldFn, pos, t, dt, options = {}, k1Cached = null) {
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
  let currentK1 = k1Cached;

  while (attempts < maxRetries) {
    attempts++;
    const hSigned = dirSign * currentDt;
    const raw = dormandPrinceStepRaw(fieldFn, pos, t, hSigned, currentK1);
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
        fsalNext: raw.fsalStage,
        kStages: raw.kStages,
        retries: attempts - 1
      };
    }

    // Step rejected: recalculate step size; k1 remains valid for retry from same pos
    currentK1 = raw.kStages[0];
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
        fsalNext: null,
        kStages: raw.kStages,
        retries: attempts,
        reason: 'STEP_UNDERFLOW'
      };
    }

    currentDt = retryDt;
  }

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
 * Traces a streamline using Dormand-Prince 5(4) (DOPRI5) adaptive integrator with optional dense output.
 *
 * @param {Function} fieldFn - Vector field evaluator.
 * @param {Array<number>} seedPos - Seed coordinate.
 * @param {Object} [options={}] - Tracing parameters.
 * @returns {Object} Structured streamline trajectory result.
 */
export function traceDormandPrinceStreamline(fieldFn, seedPos, options = {}) {
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
    minArcLengthStep = 1e-5,
    denseOutput = false,
    denseSamplesPerStep = 5
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
  let k1Cached = null;

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
    k1Cached = v0;

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
      result = dormandPrinceStep(fieldFn, currentPos, currentT, currentDt, {
        direction,
        atol,
        rtol,
        minStep,
        maxStep
      }, k1Cached);
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
    const prevPos = [...currentPos];
    const nextPos = result.posNext;
    const nextT = result.tNext;
    const v = result.velocity;
    const speed = Math.hypot(v[0], v[1], v[2]);
    const dirSign = (direction === IntegrationDirection.BACKWARD || direction === -1) ? -1.0 : 1.0;
    const hSigned = dirSign * result.dtUsed;

    // Optional Dense Output Sub-Sampling
    if (denseOutput && denseSamplesPerStep > 1) {
      for (let s = 1; s < denseSamplesPerStep; s++) {
        const theta = s / denseSamplesPerStep;
        const subPos = denseInterpolateDOPRI5(prevPos, nextPos, hSigned, result.kStages, theta);
        if (!isWithinBounds(subPos, bounds)) {
          break;
        }
        const prevP = positions[positions.length - 1];
        const subDist = Math.hypot(subPos[0] - prevP[0], subPos[1] - prevP[1], subPos[2] - prevP[2]);
        totalArcLength += subDist;
        positions.push(subPos);
        velocities.push(v);
        speeds.push(speed);
        arcLengths.push(totalArcLength);
        times.push(currentT + theta * hSigned);
      }
    }

    const lastPoint = positions[positions.length - 1];
    const stepDist = Math.hypot(nextPos[0] - lastPoint[0], nextPos[1] - lastPoint[1], nextPos[2] - lastPoint[2]);

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
    k1Cached = result.fsalNext; // FSAL Acceleration!
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
 * Dormand-Prince 5(4) (DOPRI5) Integrator Class.
 */
export class DormandPrinceIntegrator {
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
    this.denseOutput = config.denseOutput ?? false;
    this.denseSamplesPerStep = config.denseSamplesPerStep ?? 5;
  }

  step(fieldFn, pos, t, dt = this.initialStep, options = {}, k1Cached = null) {
    return dormandPrinceStep(fieldFn, pos, t, dt, {
      direction: this.direction,
      atol: this.atol,
      rtol: this.rtol,
      minStep: this.minStep,
      maxStep: this.maxStep,
      ...options
    }, k1Cached);
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
      denseOutput: this.denseOutput,
      denseSamplesPerStep: this.denseSamplesPerStep,
      ...overrideOptions
    };
    return traceDormandPrinceStreamline(fieldFn, seedPos, opts);
  }

  denseInterpolate(y0, y1, hSigned, kStages, theta) {
    return denseInterpolateDOPRI5(y0, y1, hSigned, kStages, theta);
  }
}

export default DormandPrinceIntegrator;
