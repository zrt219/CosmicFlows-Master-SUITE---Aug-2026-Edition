/**
 * @file rk4_parameter_sweep.js
 * @description Automated Runge-Kutta Parameter Sweeper evaluating step size h,
 * maximum integration steps, velocity termination thresholds, and quantifying
 * basin boundary stability and trajectory divergence metrics.
 * 
 * Mathematical Formulations:
 * 1. Step Size Convergence:
 *    Let x_h(t) be the numerical trajectory with step size h.
 *    Global truncation error e(h) = ||x_h(t) - x_exact(t)|| ~ O(h^p).
 *    For embedded pairs (Cash-Karp 4(5), DOPRI 5(4)), local error est e_loc = ||y_5 - y_4||.
 * 
 * 2. Basin Boundary Stability Index (BBSI):
 *    BBSI(h1, h2) = 1 - (1/N) * sum_{i=1}^N delta(B_h1(x_i), B_h2(x_i))
 *    where B_h(x_i) is the terminal basin ID assigned to seed x_i under step size h.
 * 
 * 3. Trajectory Hausdorff Distance:
 *    d_H(P, Q) = max( sup_{p in P} inf_{q in Q} ||p - q||, sup_{q in Q} inf_{p in P} ||p - q|| )
 * 
 * @module integration/rk4_parameter_sweep
 */

import { traceCashKarpStreamline } from './rk45_cash_karp.js';
import { DormandPrince } from './dormand_prince.js';
import { ClassicalRK4 } from './rk4_classical.js';
import { SupergalacticPosition, VelocityVector, sha256Hex } from '../coordinates/scientific_types.js';

export class RK4ParameterSweep {
  /**
   * @param {object} velocityField Velocity field engine.
   * @param {object} [options]
   */
  constructor(velocityField, options = {}) {
    if (!velocityField) {
      throw new TypeError('RK4ParameterSweep requires a valid velocityField.');
    }
    this.field = velocityField;
    this.grid = velocityField.grid;
  }

  /**
   * Evaluates trajectory convergence and error scaling across a sequence of step sizes h.
   * 
   * @param {Array<number>} seedPos Initial position [x, y, z] in Mpc/h.
   * @param {Array<number>} stepSizes Array of test step sizes [h0, h1, ...].
   * @param {number} [tMax=50.0] Maximum integration time.
   * @returns {object} Sweep diagnostics.
   */
  sweepStepSizes(seedPos, stepSizes = [2.0, 1.0, 0.5, 0.25, 0.125], tMax = 50.0) {
    const trajectories = [];
    const fieldFn = (t, pos) => this.field.sampleVelocity(pos[0], pos[1], pos[2]);

    for (const h of stepSizes) {
      const maxSteps = Math.ceil(tMax / h) + 10;
      const res = traceCashKarpStreamline(fieldFn, seedPos, {
        initialStep: h,
        minStep: h * 0.1,
        maxStep: h,
        maxSteps
      });
      trajectories.push({
        stepSize: h,
        stepCount: res.stepCount,
        arcLength: res.stats.totalArcLength,
        terminalPos: res.points[res.points.length - 1],
        points: res.points
      });
    }

    // Compare with finest step size as reference
    const ref = trajectories[trajectories.length - 1];
    const errors = [];

    for (let i = 0; i < trajectories.length - 1; i++) {
      const traj = trajectories[i];
      const pTraj = traj.terminalPos;
      const pRef = ref.terminalPos;
      const dx = pTraj[0] - pRef[0];
      const dy = pTraj[1] - pRef[1];
      const dz = pTraj[2] - pRef[2];
      const terminalDist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const hausdorff = this.computeHausdorffDistance(traj.points, ref.points);

      errors.push({
        stepSize: traj.stepSize,
        terminalDisplacementError: terminalDist,
        hausdorffDistance: hausdorff
      });
    }

    return {
      seed: seedPos,
      referenceStepSize: ref.stepSize,
      trajectories,
      errors
    };
  }

  /**
   * Computes discrete Hausdorff distance between two 3D trajectories P and Q.
   * @param {Array<Array<number>>} P
   * @param {Array<Array<number>>} Q
   * @returns {number}
   */
  computeHausdorffDistance(P, Q) {
    const dForward = this._directedHausdorff(P, Q);
    const dBackward = this._directedHausdorff(Q, P);
    return Math.max(dForward, dBackward);
  }

  _directedHausdorff(A, B) {
    let maxMinDist = 0.0;
    for (let i = 0; i < A.length; i++) {
      let minDist = Infinity;
      const ax = A[i][0], ay = A[i][1], az = A[i][2];

      for (let j = 0; j < B.length; j++) {
        const bx = B[j][0], by = B[j][1], bz = B[j][2];
        const distSq = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2;
        if (distSq < minDist) {
          minDist = distSq;
        }
      }
      const dist = Math.sqrt(minDist);
      if (dist > maxMinDist) {
        maxMinDist = dist;
      }
    }
    return maxMinDist;
  }
}
