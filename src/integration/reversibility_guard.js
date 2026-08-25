/**
 * @file reversibility_guard.js
 * @description Invariant and Time-Reversibility Monitoring for Cosmological Flows.
 * Performs forward-backward integration audits to verify Hamiltonian/flow reversibility,
 * computes Finite-Time Lyapunov Exponents (FTLE), monitors phase-space divergence,
 * and classifies trajectories into laminar, critical, or chaotic flow regimes.
 *
 * @module integration/reversibility_guard
 */

import { evaluateField, TerminationReason, IntegrationDirection } from './rk4_classical.js';
import { rk45CashKarpStep } from './rk45_cash_karp.js';

export const FlowRegime = Object.freeze({
  LAMINAR: 'LAMINAR',
  REGULAR: 'REGULAR',
  SEPARATRIX: 'SEPARATRIX',
  WEAKLY_CHAOTIC: 'WEAKLY_CHAOTIC',
  STRONGLY_CHAOTIC: 'STRONGLY_CHAOTIC',
  STAGNATION: 'STAGNATION'
});

export function testTrajectoryReversibility(fieldFn, seedPos, options = {}) {
  const {
    tDuration = 5.0,
    dt = 0.1,
    atol = 1e-6,
    rtol = 1e-6,
    maxSteps = 2000
  } = options;

  const x0 = [Number(seedPos[0]), Number(seedPos[1]), Number(seedPos[2])];

  let curPos = [...x0];
  let curT = 0.0;
  let curH = dt;
  let forwardSteps = 0;
  let forwardArc = 0.0;

  while (curT < tDuration && forwardSteps < maxSteps) {
    forwardSteps++;
    const stepH = Math.min(curH, tDuration - curT);
    const res = rk45CashKarpStep(fieldFn, curPos, curT, stepH, {
      direction: IntegrationDirection.FORWARD,
      atol,
      rtol
    });

    const dx = res.posNext[0] - curPos[0];
    const dy = res.posNext[1] - curPos[1];
    const dz = res.posNext[2] - curPos[2];
    forwardArc += Math.sqrt(dx * dx + dy * dy + dz * dz);

    curPos = res.posNext;
    curT = res.tNext;
    curH = res.hNext;
  }

  const xForwardEnd = [...curPos];
  const tForwardEnd = curT;

  let revPos = [...xForwardEnd];
  let revT = tForwardEnd;
  let revH = dt;
  let backwardSteps = 0;
  let backwardArc = 0.0;

  while (revT > 0.0 && backwardSteps < maxSteps) {
    backwardSteps++;
    const stepH = Math.min(revH, revT);
    const res = rk45CashKarpStep(fieldFn, revPos, revT, stepH, {
      direction: IntegrationDirection.BACKWARD,
      atol,
      rtol
    });

    const dx = res.posNext[0] - revPos[0];
    const dy = res.posNext[1] - revPos[1];
    const dz = res.posNext[2] - revPos[2];
    backwardArc += Math.sqrt(dx * dx + dy * dy + dz * dz);

    revPos = res.posNext;
    revT = res.tNext;
    revH = res.hNext;
  }

  const xRevEnd = [...revPos];

  const errX = xRevEnd[0] - x0[0];
  const errY = xRevEnd[1] - x0[1];
  const errZ = xRevEnd[2] - x0[2];
  const absReversalError = Math.sqrt(errX * errX + errY * errY + errZ * errZ);
  const posNorm = Math.sqrt(x0[0] * x0[0] + x0[1] * x0[1] + x0[2] * x0[2]);
  const relReversalError = absReversalError / Math.max(posNorm, 1.0);

  const arcRatio = backwardArc > 0 ? forwardArc / backwardArc : 1.0;

  return {
    seedPos: x0,
    forwardEndPoint: xForwardEnd,
    reversedEndPoint: xRevEnd,
    tDuration: tForwardEnd,
    forwardSteps,
    backwardSteps,
    forwardArcLength: forwardArc,
    backwardArcLength: backwardArc,
    arcLengthRatio: arcRatio,
    absReversalError,
    relReversalError,
    isReversible: relReversalError < 1e-3
  };
}

export function computeFTLE(fieldFn, seedPos, options = {}) {
  const {
    tDuration = 5.0,
    dt = 0.05,
    perturbation = 1e-5,
    maxSteps = 2000
  } = options;

  let x = [Number(seedPos[0]), Number(seedPos[1]), Number(seedPos[2])];
  const invSqrt3 = 1.0 / Math.sqrt(3.0);
  let dx = [perturbation * invSqrt3, perturbation * invSqrt3, perturbation * invSqrt3];
  let t = 0.0;
  let steps = 0;

  const lyapunovHistory = [];

  while (t < tDuration && steps < maxSteps) {
    steps++;
    const v = evaluateField(fieldFn, t, x);

    const eps = 1e-5;
    const J = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

    for (let c = 0; c < 3; c++) {
      const xP = [...x]; xP[c] += eps;
      const xM = [...x]; xM[c] -= eps;
      const vP = evaluateField(fieldFn, t, xP);
      const vM = evaluateField(fieldFn, t, xM);
      for (let r = 0; r < 3; r++) {
        J[r][c] = (vP[r] - vM[r]) / (2.0 * eps);
      }
    }

    const d_dx = [
      J[0][0] * dx[0] + J[0][1] * dx[1] + J[0][2] * dx[2],
      J[1][0] * dx[0] + J[1][1] * dx[1] + J[1][2] * dx[2],
      J[2][0] * dx[0] + J[2][1] * dx[1] + J[2][2] * dx[2]
    ];

    const xNext = [x[0] + dt * v[0], x[1] + dt * v[1], x[2] + dt * v[2]];
    const dxNext = [dx[0] + dt * d_dx[0], dx[1] + dt * d_dx[1], dx[2] + dt * d_dx[2]];

    x = xNext;
    dx = dxNext;
    t += dt;

    const norm = Math.sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    const lambda_t = (1.0 / t) * Math.log(norm / perturbation);
    lyapunovHistory.push({ t, lambda: lambda_t, norm });
  }

  const finalNorm = Math.sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const ftle = t > 0 ? (1.0 / t) * Math.log(finalNorm / perturbation) : 0.0;

  let regime = FlowRegime.REGULAR;
  if (ftle > 0.5) {
    regime = FlowRegime.STRONGLY_CHAOTIC;
  } else if (ftle > 0.1) {
    regime = FlowRegime.WEAKLY_CHAOTIC;
  } else if (ftle < -0.1) {
    regime = FlowRegime.STAGNATION;
  } else if (Math.abs(ftle) <= 0.1) {
    regime = FlowRegime.LAMINAR;
  }

  return {
    seedPos,
    tDuration: t,
    steps,
    initialPerturbation: perturbation,
    finalPerturbation: finalNorm,
    ftle,
    regime,
    history: lyapunovHistory
  };
}

export class ReversibilityGuard {
  constructor(config = {}) {
    this.reversibilityTol = config.reversibilityTol ?? 1e-3;
    this.ftleThreshold = config.ftleThreshold ?? 0.2;
  }

  auditTrajectory(fieldFn, seedPos, options = {}) {
    const revRes = testTrajectoryReversibility(fieldFn, seedPos, options);
    const ftleRes = computeFTLE(fieldFn, seedPos, options);

    return {
      reversibility: revRes,
      lyapunov: ftleRes,
      passedAudit: revRes.isReversible && ftleRes.ftle < this.ftleThreshold,
      regime: ftleRes.regime
    };
  }
}

export default ReversibilityGuard;
