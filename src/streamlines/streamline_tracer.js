/**
 * @file streamline_tracer.js
 * @description High-level batch streamline tracing engine for cosmological vector fields.
 * Coordinates RK4, Cash-Karp RK45, and Dormand-Prince integrators with explicit termination
 * condition handling (CONVERGED_ENDPOINT, DOMAIN_EXIT, MAX_STEPS, LOW_SPEED_STALL, NUMERICAL_FAILURE),
 * bidirectional tracing, and typed array memory buffering.
 *
 * @module streamlines/streamline_tracer
 */

import { RK4Integrator, traceRK4Streamline, TerminationReason, IntegrationDirection } from '../integration/rk4_classical.js';
import { RK45CashKarpIntegrator, traceCashKarpStreamline } from '../integration/rk45_cash_karp.js';
import { DormandPrinceIntegrator, traceDormandPrinceStreamline } from '../integration/dormand_prince.js';
import { SeedManager, StreamlineSeed } from './seed_manager.js';

export { TerminationReason, IntegrationDirection };

/**
 * Supported integrator algorithm identifiers.
 * @readonly
 * @enum {string}
 */
export const IntegratorType = Object.freeze({
  RK4: 'RK4',
  RK45_CASH_KARP: 'RK45_CASH_KARP',
  DORMAND_PRINCE: 'DORMAND_PRINCE'
});

/**
 * Tracing mode.
 * @readonly
 * @enum {string}
 */
export const TracingMode = Object.freeze({
  FORWARD: 'forward',
  BACKWARD: 'backward',
  BIDIRECTIONAL: 'bidirectional'
});

/**
 * Structured container for an integrated streamline trajectory.
 */
export class Streamline {
  /**
   * @param {Object} data
   */
  constructor(data = {}) {
    this.id = data.id || `streamline_${Math.random().toString(36).substr(2, 9)}`;
    this.seedPoint = data.seedPoint ? [Number(data.seedPoint[0]), Number(data.seedPoint[1]), Number(data.seedPoint[2])] : [0, 0, 0];
    this.seedMetadata = data.seedMetadata || {};
    this.direction = data.direction || TracingMode.FORWARD;
    this.integratorType = data.integratorType || IntegratorType.DORMAND_PRINCE;

    // Trajectory arrays
    this.points = data.points || [];
    this.velocities = data.velocities || [];
    this.speeds = data.speeds || [];
    this.arcLengths = data.arcLengths || [0];
    this.times = data.times || [0];

    // Lifecycle / Termination
    this.terminationReason = data.terminationReason || TerminationReason.MAX_STEPS;
    this.stepCount = data.stepCount || 0;
    this.stats = data.stats || {
      totalArcLength: 0,
      acceptedSteps: 0,
      rejectedSteps: 0,
      minSpeed: 0,
      maxSpeed: 0,
      meanSpeed: 0
    };
  }

  get vertexCount() {
    return this.points.length;
  }

  get totalArcLength() {
    return this.stats.totalArcLength || (this.arcLengths.length > 0 ? this.arcLengths[this.arcLengths.length - 1] : 0);
  }

  /**
   * Converts trajectory vertices to a flat Float32Array suitable for WebGL line buffers.
   * @returns {Float32Array} [x0, y0, z0, x1, y1, z1, ...]
   */
  toPositionBuffer() {
    const n = this.points.length;
    const buf = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) {
      buf[i * 3 + 0] = this.points[i][0];
      buf[i * 3 + 1] = this.points[i][1];
      buf[i * 3 + 2] = this.points[i][2];
    }
    return buf;
  }

  /**
   * Converts speed scalars to a flat Float32Array suitable for WebGL color mapping.
   * @returns {Float32Array}
   */
  toSpeedBuffer() {
    return new Float32Array(this.speeds);
  }

  /**
   * Returns serializable JSON object.
   */
  toJSON() {
    return {
      id: this.id,
      seedPoint: this.seedPoint,
      seedMetadata: this.seedMetadata,
      direction: this.direction,
      integratorType: this.integratorType,
      points: this.points,
      speeds: this.speeds,
      arcLengths: this.arcLengths,
      times: this.times,
      terminationReason: this.terminationReason,
      stepCount: this.stepCount,
      stats: this.stats
    };
  }
}

/**
 * High-Level Streamline Tracing Engine.
 */
export class StreamlineTracer {
  /**
   * @param {Object} [config={}]
   * @param {string} [config.integrator='DORMAND_PRINCE'] - Integrator algorithm ('RK4', 'RK45_CASH_KARP', 'DORMAND_PRINCE').
   * @param {string} [config.mode='forward'] - Tracing mode ('forward', 'backward', 'bidirectional').
   * @param {number} [config.initialStep=0.5] - Initial step size.
   * @param {number} [config.minStep=1e-6] - Minimum step size.
   * @param {number} [config.maxStep=10.0] - Maximum step size.
   * @param {number} [config.atol=1e-6] - Absolute error tolerance.
   * @param {number} [config.rtol=1e-6] - Relative error tolerance.
   * @param {number} [config.maxSteps=2000] - Maximum integration steps.
   * @param {Object|Array<Array<number>>|null} [config.domain=null] - Bounding domain box.
   * @param {number} [config.minVelocity=1e-6] - Endpoint convergence velocity threshold.
   * @param {number} [config.stallWindow=10] - Window size for low-speed stall detection.
   * @param {number} [config.minArcLengthStep=1e-5] - Minimum displacement per step to prevent stall.
   * @param {boolean} [config.denseOutput=false] - DOPRI5 continuous dense sub-sampling.
   * @param {number} [config.denseSamplesPerStep=5] - Dense samples per DOPRI5 step.
   */
  constructor(config = {}) {
    this.integratorType = config.integrator || config.integratorType || IntegratorType.DORMAND_PRINCE;
    this.mode = config.mode || TracingMode.FORWARD;
    this.initialStep = config.initialStep ?? config.dt ?? 0.5;
    this.minStep = config.minStep ?? 1e-6;
    this.maxStep = config.maxStep ?? 10.0;
    this.atol = config.atol ?? 1e-6;
    this.rtol = config.rtol ?? 1e-6;
    this.maxSteps = config.maxSteps ?? 2000;
    this.domain = config.domain ?? config.bounds ?? null;
    this.minVelocity = config.minVelocity ?? 1e-6;
    this.stallWindow = config.stallWindow ?? 10;
    this.minArcLengthStep = config.minArcLengthStep ?? 1e-5;
    this.denseOutput = config.denseOutput ?? false;
    this.denseSamplesPerStep = config.denseSamplesPerStep ?? 5;
  }

  /**
   * Internal runner dispatching to specified single-direction integrator.
   * @private
   */
  _traceSingleDirection(fieldFn, seedPos, dir, overrideOptions = {}) {
    const opts = {
      dt: overrideOptions.initialStep ?? this.initialStep,
      initialStep: overrideOptions.initialStep ?? this.initialStep,
      direction: dir,
      maxSteps: overrideOptions.maxSteps ?? this.maxSteps,
      tStart: overrideOptions.tStart ?? 0.0,
      bounds: overrideOptions.domain ?? this.domain,
      atol: overrideOptions.atol ?? this.atol,
      rtol: overrideOptions.rtol ?? this.rtol,
      minStep: overrideOptions.minStep ?? this.minStep,
      maxStep: overrideOptions.maxStep ?? this.maxStep,
      minVelocity: overrideOptions.minVelocity ?? this.minVelocity,
      stallWindow: overrideOptions.stallWindow ?? this.stallWindow,
      minArcLengthStep: overrideOptions.minArcLengthStep ?? this.minArcLengthStep,
      denseOutput: overrideOptions.denseOutput ?? this.denseOutput,
      denseSamplesPerStep: overrideOptions.denseSamplesPerStep ?? this.denseSamplesPerStep
    };

    const type = (overrideOptions.integratorType || this.integratorType).toUpperCase();

    if (type === IntegratorType.RK4) {
      return traceRK4Streamline(fieldFn, seedPos, opts);
    } else if (type === IntegratorType.RK45_CASH_KARP) {
      return traceCashKarpStreamline(fieldFn, seedPos, opts);
    } else {
      return traceDormandPrinceStreamline(fieldFn, seedPos, opts);
    }
  }

  /**
   * Traces a single streamline from a seed position or StreamlineSeed.
   *
   * @param {Function} fieldFn - Vector field evaluator f(t, pos) or f(pos).
   * @param {StreamlineSeed|Array<number>} seedOrPos - Starting seed.
   * @param {Object} [overrideOptions={}] - Instance config overrides.
   * @returns {Streamline}
   */
  traceStreamline(fieldFn, seedOrPos, overrideOptions = {}) {
    let seedPoint;
    let seedMetadata = {};
    let seedId;

    if (seedOrPos instanceof StreamlineSeed) {
      seedPoint = seedOrPos.position;
      seedMetadata = seedOrPos.metadata;
      seedId = seedOrPos.id;
    } else if (Array.isArray(seedOrPos)) {
      seedPoint = [Number(seedOrPos[0]), Number(seedOrPos[1]), Number(seedOrPos[2])];
      seedId = overrideOptions.id || `streamline_${Math.random().toString(36).substr(2, 9)}`;
    } else if (seedOrPos && typeof seedOrPos === 'object') {
      seedPoint = seedOrPos.pos || [seedOrPos.x, seedOrPos.y, seedOrPos.z];
      seedMetadata = seedOrPos.metadata || {};
      seedId = seedOrPos.id;
    }

    const mode = overrideOptions.mode || this.mode;
    const type = overrideOptions.integratorType || this.integratorType;

    if (mode === TracingMode.BIDIRECTIONAL) {
      // 1. Trace Backward
      const backResult = this._traceSingleDirection(fieldFn, seedPoint, IntegrationDirection.BACKWARD, overrideOptions);
      // 2. Trace Forward
      const fwdResult = this._traceSingleDirection(fieldFn, seedPoint, IntegrationDirection.FORWARD, overrideOptions);

      // Combine backward (reversed, excluding seed) + forward (including seed)
      const combinedPoints = [];
      const combinedVelocities = [];
      const combinedSpeeds = [];
      const combinedTimes = [];
      const combinedArcLengths = [0.0];

      // Add reversed backward branch (excluding index 0 which is seed)
      for (let i = backResult.points.length - 1; i >= 1; i--) {
        combinedPoints.push(backResult.points[i]);
        combinedVelocities.push(backResult.velocities[i] || [0, 0, 0]);
        combinedSpeeds.push(backResult.speeds[i] || 0);
        combinedTimes.push(backResult.times[i] || 0);
      }

      // Add forward branch (starts with seed at index 0)
      for (let i = 0; i < fwdResult.points.length; i++) {
        combinedPoints.push(fwdResult.points[i]);
        combinedVelocities.push(fwdResult.velocities[i] || [0, 0, 0]);
        combinedSpeeds.push(fwdResult.speeds[i] || 0);
        combinedTimes.push(fwdResult.times[i] || 0);
      }

      // Compute cumulative arc length along combined curve
      let totalDist = 0.0;
      for (let i = 1; i < combinedPoints.length; i++) {
        const p0 = combinedPoints[i - 1];
        const p1 = combinedPoints[i];
        totalDist += Math.hypot(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);
        combinedArcLengths.push(totalDist);
      }

      let minSpeed = Infinity;
      let maxSpeed = -Infinity;
      let sumSpeed = 0.0;
      for (let i = 0; i < combinedSpeeds.length; i++) {
        const s = combinedSpeeds[i];
        if (s < minSpeed) minSpeed = s;
        if (s > maxSpeed) maxSpeed = s;
        sumSpeed += s;
      }

      // Primary termination condition prioritizes forward, then backward
      const termReason = fwdResult.terminationReason !== TerminationReason.MAX_STEPS
        ? fwdResult.terminationReason
        : backResult.terminationReason;

      return new Streamline({
        id: seedId,
        seedPoint,
        seedMetadata,
        direction: TracingMode.BIDIRECTIONAL,
        integratorType: type,
        points: combinedPoints,
        velocities: combinedVelocities,
        speeds: combinedSpeeds,
        arcLengths: combinedArcLengths,
        times: combinedTimes,
        terminationReason: termReason,
        stepCount: backResult.stepCount + fwdResult.stepCount,
        stats: {
          totalArcLength: totalDist,
          acceptedSteps: backResult.stats.acceptedSteps + fwdResult.stats.acceptedSteps,
          rejectedSteps: backResult.stats.rejectedSteps + fwdResult.stats.rejectedSteps,
          minSpeed: Number.isFinite(minSpeed) ? minSpeed : 0,
          maxSpeed: Number.isFinite(maxSpeed) ? maxSpeed : 0,
          meanSpeed: combinedSpeeds.length > 0 ? sumSpeed / combinedSpeeds.length : 0
        }
      });
    }

    // Single direction trace (forward or backward)
    const dir = mode === TracingMode.BACKWARD ? IntegrationDirection.BACKWARD : IntegrationDirection.FORWARD;
    const res = this._traceSingleDirection(fieldFn, seedPoint, dir, overrideOptions);

    return new Streamline({
      id: seedId,
      seedPoint,
      seedMetadata,
      direction: mode,
      integratorType: type,
      points: res.points,
      velocities: res.velocities,
      speeds: res.speeds,
      arcLengths: res.arcLengths,
      times: res.times,
      terminationReason: res.terminationReason,
      stepCount: res.stepCount,
      stats: res.stats
    });
  }

  /**
   * Traces a batch of seeds synchronously or with progress notifications.
   *
   * @param {Function} fieldFn - Vector field evaluator.
   * @param {SeedManager|Array<StreamlineSeed|Array<number>>} seedCollection - Seeds to trace.
   * @param {Object} [overrideOptions={}] - Tracing options.
   * @param {Function|null} [onProgress=null] - Progress callback ({ completed, total, percentage, activeStreamline }).
   * @returns {{ streamlines: Array<Streamline>, stats: Object }}
   */
  traceBatch(fieldFn, seedCollection, overrideOptions = {}, onProgress = null) {
    const startTime = performance ? performance.now() : Date.now();

    let seeds = [];
    if (seedCollection instanceof SeedManager) {
      seeds = seedCollection.getAllSeeds();
    } else if (Array.isArray(seedCollection)) {
      seeds = seedCollection;
    }

    const total = seeds.length;
    const streamlines = [];
    const terminationCounts = {
      [TerminationReason.CONVERGED_ENDPOINT]: 0,
      [TerminationReason.DOMAIN_EXIT]: 0,
      [TerminationReason.MAX_STEPS]: 0,
      [TerminationReason.LOW_SPEED_STALL]: 0,
      [TerminationReason.NUMERICAL_FAILURE]: 0
    };

    let totalVertices = 0;
    let totalArcLengthSum = 0;

    for (let i = 0; i < total; i++) {
      const seed = seeds[i];
      const sl = this.traceStreamline(fieldFn, seed, overrideOptions);
      streamlines.push(sl);

      // Accumulate stats
      if (terminationCounts[sl.terminationReason] !== undefined) {
        terminationCounts[sl.terminationReason]++;
      }
      totalVertices += sl.vertexCount;
      totalArcLengthSum += sl.totalArcLength;

      if (onProgress && typeof onProgress === 'function') {
        const completed = i + 1;
        const percentage = Math.round((completed / total) * 100);
        onProgress({
          completed,
          total,
          percentage,
          activeStreamline: sl
        });
      }
    }

    const endTime = performance ? performance.now() : Date.now();
    const durationMs = endTime - startTime;

    return {
      streamlines,
      stats: {
        totalStreamlines: total,
        totalVertices,
        averageVerticesPerStreamline: total > 0 ? totalVertices / total : 0,
        averageArcLength: total > 0 ? totalArcLengthSum / total : 0,
        durationMs,
        terminationCounts
      }
    };
  }
}

export default StreamlineTracer;
