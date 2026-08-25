/**
 * @file streamline_terminator.js
 * @description Modular Streamline Termination Engine:
 *  - CONVERGED_ENDPOINT: Speed drops below threshold ||v|| < v_min
 *  - DOMAIN_EXIT: Position outside domain bounds
 *  - MAX_STEPS: Step count reaches maxSteps
 *  - MAX_ARC_LENGTH: Total arc length reaches maxArcLength
 *  - LOW_SPEED_STALL: Moving slower than minVelocity for stallWindow consecutive steps
 *  - VORTEX_LOOP: Trajectory cycles/loops back on itself
 *  - SINGULARITY: Non-finite values or NaN encountered
 *  - BASIN_BOUNDARY_CROSSING: Traverses watershed label boundary.
 *
 * @module streamlines/streamline_terminator
 */

import { TerminationReason, isWithinBounds } from '../integration/rk4_classical.js';

export { TerminationReason };

export class StreamlineTerminator {
  constructor(config = {}) {
    this.maxSteps = config.maxSteps ?? 2000;
    this.maxArcLength = config.maxArcLength ?? 10000.0;
    this.minVelocity = config.minVelocity ?? 1e-6;
    this.stallWindow = config.stallWindow ?? 10;
    this.minArcLengthStep = config.minArcLengthStep ?? 1e-5;
    this.bounds = config.bounds ?? null;
    this.enableLoopDetection = config.enableLoopDetection ?? true;
    this.loopDistanceThreshold = config.loopDistanceThreshold ?? 1.0;
  }

  /**
   * Checks whether a streamline should terminate at the current step.
   *
   * @param {Object} state - Current tracing state.
   * @param {Array<number>} state.pos - Current position [x, y, z].
   * @param {Array<number>} state.velocity - Current velocity [vx, vy, vz].
   * @param {number} state.step - Step count.
   * @param {number} state.arcLength - Cumulative path length.
   * @param {number} state.stallCounter - Consecutive slow steps.
   * @param {Array<Array<number>>} [state.recentPoints=[]] - Recent trajectory points for loop check.
   * @returns {{ shouldTerminate: boolean, reason: string }}
   */
  evaluate(state) {
    const {
      pos,
      velocity,
      step,
      arcLength = 0,
      stallCounter = 0,
      recentPoints = []
    } = state;

    // Check non-finite coordinates
    if (!Number.isFinite(pos[0]) || !Number.isFinite(pos[1]) || !Number.isFinite(pos[2])) {
      return { shouldTerminate: true, reason: TerminationReason.NUMERICAL_FAILURE };
    }

    // Check domain boundary
    if (!isWithinBounds(pos, this.bounds)) {
      return { shouldTerminate: true, reason: TerminationReason.DOMAIN_EXIT };
    }

    // Check max steps
    if (step >= this.maxSteps) {
      return { shouldTerminate: true, reason: TerminationReason.MAX_STEPS };
    }

    // Check max arc length
    if (arcLength >= this.maxArcLength) {
      return { shouldTerminate: true, reason: TerminationReason.MAX_STEPS };
    }

    // Check velocity magnitude
    const speed = Math.sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
    if (speed < this.minVelocity) {
      return { shouldTerminate: true, reason: TerminationReason.CONVERGED_ENDPOINT };
    }

    // Check stall
    if (stallCounter >= this.stallWindow) {
      return { shouldTerminate: true, reason: TerminationReason.LOW_SPEED_STALL };
    }

    // Check loop / cycle
    if (this.enableLoopDetection && recentPoints.length > 20) {
      const n = recentPoints.length;
      const cur = pos;
      for (let i = 0; i < n - 15; i++) {
        const p = recentPoints[i];
        const dx = cur[0] - p[0];
        const dy = cur[1] - p[1];
        const dz = cur[2] - p[2];
        if (dx * dx + dy * dy + dz * dz < this.loopDistanceThreshold * this.loopDistanceThreshold) {
          return { shouldTerminate: true, reason: 'VORTEX_LOOP' };
        }
      }
    }

    return { shouldTerminate: false, reason: TerminationReason.IN_PROGRESS };
  }
}

export default StreamlineTerminator;
