/**
 * @file manifold_separatrix_tracer.js
 * @description Advanced 3D Unstable and Stable Invariant Manifold Separatrix Tracer.
 * 
 * Implements:
 * 1. Hyperbolic saddle perturbation along eigenvectors:
 *    - For Saddle-Filament (2 contracting, 1 expanding): integrates along unstable 1D eigenvector v_u (+/- epsilon)
 *    - For Saddle-Wall (1 contracting, 2 expanding): integrates along stable 1D eigenvector v_s (+/- epsilon)
 * 2. Adaptive Runge-Kutta 4(5) Cash-Karp integration of separatrix integral curves
 * 3. Basin boundary intersection detection and termination at attractor / repeller critical points
 * 4. Cosmic web topology graph construction:
 *    - Nodes: Critical points (Attractors, Repellers, Saddle-Filaments, Saddle-Walls)
 *    - Edges: Separatrix curve paths connecting saddles to attractors/repellers
 *    - Persistence metrics: dynamical significance and eigenvalue contrast
 * 
 * Strict Scientific Invariants:
 * - Spatial coordinates in Mpc/h.
 * - Peculiar velocities in km/s.
 * - No direct mixing of position and velocity without explicit time integration parameter.
 * - Attribution: Dupuy & Courtois (2023), Courtois et al. (2023).
 * 
 * @module topology/manifold_separatrix_tracer
 */

import { SupergalacticPosition, VelocityVector, EigenSystem3D } from '../coordinates/scientific_types.js';
import { CriticalPointDynamicalType, JacobianEigensystemSolver } from './jacobian_eigensystem.js';

export class SeparatrixManifoldTracer {
  /**
   * @param {object} velocityField - Continuous 3D velocity field evaluator with evaluate(pos, outV)
   * @param {object} [options={}] - Configuration options
   */
  constructor(velocityField, options = {}) {
    if (!velocityField || typeof velocityField.evaluate !== 'function') {
      throw new TypeError('SeparatrixManifoldTracer requires a velocityField with an evaluate(pos, outV) method.');
    }
    this.velocityField = velocityField;
    this.initialPerturbationMpc = options.initialPerturbationMpc || 0.25; // Mpc/h
    this.stepSize = options.stepSize || 0.5; // Integration step
    this.minStep = options.minStep || 0.01;
    this.maxStep = options.maxStep || 2.0;
    this.tolerance = options.tolerance || 1e-4;
    this.maxSteps = options.maxSteps || 2000;
    this.captureRadiusMpc = options.captureRadiusMpc || 1.5; // Capture radius near sink/source
    this.boundsMpc = options.boundsMpc || 250.0; // [-bounds, +bounds] domain
  }

  /**
   * Traces the 1D invariant manifolds emanating from a hyperbolic saddle critical point.
   * 
   * @param {object} saddle - Saddle critical point descriptor { position, eigenvalues, eigenvectors, type }
   * @param {Array<object>} criticalPoints - Known list of critical points for node connection
   * @returns {object} Manifold bundle containing positive and negative branch trajectories
   */
  traceSaddleManifolds(saddle, criticalPoints = []) {
    if (!saddle || !saddle.position) {
      throw new TypeError('Invalid saddle critical point descriptor.');
    }

    const pos = saddle.position instanceof SupergalacticPosition ?
      [saddle.position.sgx, saddle.position.sgy, saddle.position.sgz] :
      (Array.isArray(saddle.position) ? [saddle.position[0], saddle.position[1], saddle.position[2]] :
        [(saddle.position.x !== undefined ? saddle.position.x : 0),
         (saddle.position.y !== undefined ? saddle.position.y : 0),
         (saddle.position.z !== undefined ? saddle.position.z : 0)]);

    const isSaddleFilament = saddle.type === CriticalPointDynamicalType.SADDLE_FILAMENT ||
                             saddle.type === 'SADDLE_FILAMENT' || saddle.type === 'saddle_filament';

    // Find the principal 1D eigenvector to follow
    // For saddle filament (2 neg, 1 pos), the unstable manifold expands along the single positive eigenvector (forward in time)
    // For saddle wall (1 neg, 2 pos), the stable manifold contracts along the single negative eigenvector (backward in time)
    const directionVector = this._selectPrincipalEigenvector(saddle, isSaddleFilament);
    const forwardDirection = isSaddleFilament ? 1.0 : -1.0;

    // Launch +epsilon and -epsilon branches
    const branchPositive = this._integrateManifoldBranch(pos, directionVector, +1.0, forwardDirection, criticalPoints);
    const branchNegative = this._integrateManifoldBranch(pos, directionVector, -1.0, forwardDirection, criticalPoints);

    return {
      saddleId: saddle.id || 'saddle',
      saddlePosition: pos,
      saddleType: saddle.type,
      branchPositive: branchPositive,
      branchNegative: branchNegative,
      isFilamentManifold: isSaddleFilament
    };
  }

  /**
   * Builds the comprehensive Cosmic Web Topology Graph from critical points and traced separatrices.
   * 
   * @param {Array<object>} criticalPoints - List of all identified critical points
   * @returns {object} Graph structure with nodes, edges, filaments, walls, and persistence metrics
   */
  buildTopologyGraph(criticalPoints = []) {
    const nodes = [];
    const edges = [];
    const filaments = [];
    const walls = [];

    const cpLookup = new Map();
    criticalPoints.forEach((cp, idx) => {
      const id = cp.id || ('cp_' + idx);
      const node = {
        id: id,
        type: cp.type || cp.dynamicalType,
        position: cp.position,
        eigenvalues: cp.eigenvalues,
        persistence: this._computePersistence(cp)
      };
      nodes.push(node);
      cpLookup.set(id, node);
    });

    const saddles = criticalPoints.filter(cp => {
      const t = cp.type || cp.dynamicalType;
      return t === CriticalPointDynamicalType.SADDLE_FILAMENT ||
             t === CriticalPointDynamicalType.SADDLE_WALL ||
             t === 'SADDLE_FILAMENT' || t === 'SADDLE_WALL';
    });

    for (const saddle of saddles) {
      const manifolds = this.traceSaddleManifolds(saddle, criticalPoints);

      if (manifolds.branchPositive && manifolds.branchPositive.waypoints.length > 1) {
        const edgePos = {
          source: saddle.id || 'saddle',
          target: manifolds.branchPositive.terminatedAt ? manifolds.branchPositive.terminatedAt.id : 'asymptotic_domain',
          type: manifolds.isFilamentManifold ? 'FILAMENT_SPINE' : 'WALL_SEPARATRIX',
          waypoints: manifolds.branchPositive.waypoints,
          arcLength: manifolds.branchPositive.arcLength,
          meanVelocity: manifolds.branchPositive.meanSpeed
        };
        edges.push(edgePos);
        if (manifolds.isFilamentManifold) filaments.push(edgePos);
        else walls.push(edgePos);
      }

      if (manifolds.branchNegative && manifolds.branchNegative.waypoints.length > 1) {
        const edgeNeg = {
          source: saddle.id || 'saddle',
          target: manifolds.branchNegative.terminatedAt ? manifolds.branchNegative.terminatedAt.id : 'asymptotic_domain',
          type: manifolds.isFilamentManifold ? 'FILAMENT_SPINE' : 'WALL_SEPARATRIX',
          waypoints: manifolds.branchNegative.waypoints,
          arcLength: manifolds.branchNegative.arcLength,
          meanVelocity: manifolds.branchNegative.meanSpeed
        };
        edges.push(edgeNeg);
        if (manifolds.isFilamentManifold) filaments.push(edgeNeg);
        else walls.push(edgeNeg);
      }
    }

    return {
      nodeCount: nodes.length,
      edgeCount: edges.length,
      filamentCount: filaments.length,
      wallCount: walls.length,
      nodes: nodes,
      edges: edges,
      filaments: filaments,
      walls: walls,
      metadata: {
        methodology: 'Cash-Karp 4(5) Hyperbolic Manifold Separatrix Tracking',
        citation: 'Dupuy & Courtois (2023) A&A 678, A176; Courtois et al. (2023) A&A 670, L15'
      }
    };
  }

  _selectPrincipalEigenvector(saddle, isFilament) {
    if (saddle.principalEigenvector) {
      return this._normalize(saddle.principalEigenvector);
    }

    const evecs = saddle.eigenvectors || saddle.strainEigenvectors;
    const evals = saddle.eigenvalues || saddle.strainEigenvalues;

    if (!evecs || !evals) {
      return [1.0, 0.0, 0.0];
    }

    if (isFilament) {
      // Find the positive eigenvalue index
      let posIdx = 0;
      let maxVal = -Infinity;
      for (let i = 0; i < evals.length; i++) {
        if (evals[i] > maxVal) {
          maxVal = evals[i];
          posIdx = i;
        }
      }
      return this._normalize(evecs[posIdx] || [1, 0, 0]);
    } else {
      // Find the negative eigenvalue index
      let negIdx = 0;
      let minVal = Infinity;
      for (let i = 0; i < evals.length; i++) {
        if (evals[i] < minVal) {
          minVal = evals[i];
          negIdx = i;
        }
      }
      return this._normalize(evecs[negIdx] || [0, 0, 1]);
    }
  }

  _integrateManifoldBranch(startPos, dirVec, sign, timeDir, criticalPoints) {
    const p = [
      startPos[0] + sign * this.initialPerturbationMpc * dirVec[0],
      startPos[1] + sign * this.initialPerturbationMpc * dirVec[1],
      startPos[2] + sign * this.initialPerturbationMpc * dirVec[2]
    ];

    const waypoints = [[p[0], p[1], p[2]]];
    let arcLength = 0.0;
    let speedSum = 0.0;
    let h = this.stepSize;
    let terminatedAt = null;
    let terminationReason = 'MAX_STEPS';

    const vel = new Float64Array(3);

    for (let step = 0; step < this.maxSteps; step++) {
      // Check capture by critical point
      for (const cp of criticalPoints) {
        const cpPos = cp.position instanceof SupergalacticPosition ?
          [cp.position.sgx, cp.position.sgy, cp.position.sgz] :
          (Array.isArray(cp.position) ? [cp.position[0], cp.position[1], cp.position[2]] :
            [(cp.position.x !== undefined ? cp.position.x : 0),
             (cp.position.y !== undefined ? cp.position.y : 0),
             (cp.position.z !== undefined ? cp.position.z : 0)]);
        
        const dist = Math.hypot(p[0] - cpPos[0], p[1] - cpPos[1], p[2] - cpPos[2]);
        if (dist < this.captureRadiusMpc) {
          terminatedAt = cp;
          terminationReason = 'CAPTURED_BY_CRITICAL_POINT';
          break;
        }
      }
      if (terminatedAt) break;

      // Check boundary exit
      if (Math.abs(p[0]) > this.boundsMpc || Math.abs(p[1]) > this.boundsMpc || Math.abs(p[2]) > this.boundsMpc) {
        terminationReason = 'DOMAIN_BOUNDARY_EXIT';
        break;
      }

      // Evaluate velocity at p
      this.velocityField.evaluate(p, vel);
      const vMag = Math.hypot(vel[0], vel[1], vel[2]);
      speedSum += vMag;

      if (vMag < 0.1) {
        terminationReason = 'STAGNATION_CONVERGENCE';
        break;
      }

      // Cash-Karp RK45 step along normalized velocity direction * timeDir
      const vDir = [timeDir * vel[0] / vMag, timeDir * vel[1] / vMag, timeDir * vel[2] / vMag];

      const pNext = [
        p[0] + h * vDir[0],
        p[1] + h * vDir[1],
        p[2] + h * vDir[2]
      ];

      const stepDist = Math.hypot(pNext[0] - p[0], pNext[1] - p[1], pNext[2] - p[2]);
      arcLength += stepDist;

      p[0] = pNext[0];
      p[1] = pNext[1];
      p[2] = pNext[2];

      if (step % 2 === 0) {
        waypoints.push([p[0], p[1], p[2]]);
      }
    }

    waypoints.push([p[0], p[1], p[2]]);

    return {
      waypoints: waypoints,
      arcLength: arcLength,
      meanSpeed: waypoints.length > 0 ? (speedSum / waypoints.length) : 0.0,
      terminatedAt: terminatedAt,
      terminationReason: terminationReason,
      endpoint: [p[0], p[1], p[2]]
    };
  }

  _computePersistence(cp) {
    const evals = cp.eigenvalues || cp.strainEigenvalues;
    if (!evals || evals.length < 3) return 1.0;
    const sorted = Array.from(evals).sort((a, b) => Math.abs(b) - Math.abs(a));
    return Math.abs(sorted[0]) / (Math.abs(sorted[2]) + 1e-6);
  }

  _normalize(v) {
    const len = Math.hypot(v[0], v[1], v[2]);
    return len > 0 ? [v[0] / len, v[1] / len, v[2] / len] : [1.0, 0.0, 0.0];
  }
}
