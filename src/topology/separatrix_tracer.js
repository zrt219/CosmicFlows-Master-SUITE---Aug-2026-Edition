/**
 * @file separatrix_tracer.js
 * @description Numerical Manifold Integration along Hyperbolic Saddle Eigenvectors.
 * Traces 1D unstable/stable separatrices to construct the Cosmic Web Topology Graph:
 * - Filament Spines connecting Type-2 Saddles to Supercluster Attractors
 * - Void Drains / Wall Sheets connecting Type-1 Saddles to Repellers
 * - Graph construction, shortest path traversal, and geometric metric computation.
 */

import { sampleVelocityField } from '../watershed/endpoint_segmentation.js';
import { analyzeVelocityTensor } from './eigen_topology.js';
import { solveNewtonRaphson3D } from './root_finder.js';

/**
 * @typedef {Object} TopologyNode
 * @property {string} id - Unique identifier (e.g. 'ATT_1', 'SAD_FIL_3')
 * @property {'ATTRACTOR'|'SADDLE_FILAMENT'|'SADDLE_WALL'|'REPELLER'|'DEGENERATE'} type - Node classification
 * @property {number} morseIndex - Morse index (0 to 3)
 * @property {[number, number, number]} position - Coordinates [x, y, z] in Mpc/h
 * @property {number[]} eigenvalues - Jacobian eigenvalues [λ1, λ2, λ3]
 * @property {number[][]} eigenvectors - Eigenvectors [e1, e2, e3]
 * @property {string} label - Display label
 */

/**
 * @typedef {Object} TopologyEdge
 * @property {string} id - Edge identifier
 * @property {string} sourceId - Source node ID (usually saddle)
 * @property {string} targetId - Target node ID (attractor, repeller, or 'BOUNDARY')
 * @property {'FILAMENT_SPINE'|'WALL_BOUNDARY'|'VOID_DRAIN'} edgeType - Manifold type
 * @property {Array<[number, number, number]>} points - 3D polyline vertices
 * @property {number} arcLength - Total curve arc length in Mpc/h
 * @property {number} meanVelocity - Average velocity norm along curve (km/s)
 * @property {number} meanCurvature - Average geometric curvature along polyline
 * @property {'converged'|'boundary'|'max_steps'} status - Termination reason
 */

/**
 * @typedef {Object} CosmicWebGraph
 * @property {TopologyNode[]} nodes - All topological critical points
 * @property {TopologyEdge[]} edges - All traced separatrices / filament spines
 * @property {Record<string, string[]>} adjacency - Adjacency map: nodeId -> [connectedNodeIds]
 * @property {Record<string, TopologyNode>} nodeMap - Fast lookup map: nodeId -> node
 */

/**
 * Integrates a single 1D separatrix manifold curve starting from a displaced saddle point.
 *
 * @param {Float32Array} velocityGrid - Velocity field buffer
 * @param {number} nx
 * @param {number} ny
 * @param {number} nz
 * @param {number} halfExtent - Box half-extent in Mpc/h
 * @param {[number, number, number]} startPos - Initial displaced position x0
 * @param {1|-1} direction - Time direction: +1 = forward (unstable manifold), -1 = backward (stable manifold)
 * @param {TopologyNode[]} targetCandidates - Candidate nodes to snap/converge into
 * @param {Object} [options={}]
 * @param {number} [options.dt=0.05] - Time step
 * @param {number} [options.maxSteps=800] - Max integration steps
 * @param {number} [options.captureRadiusMpc=5.0] - Target node capture radius (Mpc/h)
 * @param {number} [options.velTol=1e-3] - Velocity stagnation threshold
 * @returns {{
 *   points: Array<[number, number, number]>,
 *   targetNodeId: string|null,
 *   arcLength: number,
 *   meanVelocity: number,
 *   meanCurvature: number,
 *   status: 'converged'|'boundary'|'max_steps'
 * }}
 */
export function traceManifoldCurve(
  velocityGrid,
  nx,
  ny,
  nz,
  halfExtent,
  startPos,
  direction,
  targetCandidates,
  options = {}
) {
  const dt = options.dt ?? 0.05;
  const maxSteps = options.maxSteps ?? 800;
  const captureRadius = options.captureRadiusMpc ?? 5.0;
  const captureRadiusSq = captureRadius * captureRadius;
  const vTol = options.velTol ?? 1e-3;

  const points = [[...startPos]];
  let [x, y, z] = startPos;
  let arcLength = 0;
  let sumVel = 0;
  let targetNodeId = null;
  let status = 'max_steps';

  for (let s = 0; s < maxSteps; s++) {
    // Check boundary
    if (Math.abs(x) > halfExtent || Math.abs(y) > halfExtent || Math.abs(z) > halfExtent) {
      status = 'boundary';
      break;
    }

    // Check proximity to target nodes
    for (const target of targetCandidates) {
      const dx = x - target.position[0];
      const dy = y - target.position[1];
      const dz = z - target.position[2];
      if (dx * dx + dy * dy + dz * dz < captureRadiusSq) {
        targetNodeId = target.id;
        status = 'converged';
        points.push([...target.position]);
        arcLength += Math.hypot(dx, dy, dz);
        break;
      }
    }
    if (targetNodeId) break;

    // Sample velocity
    const [vx, vy, vz] = sampleVelocityField(velocityGrid, nx, ny, nz, halfExtent, x, y, z);
    const speed = Math.hypot(vx, vy, vz);
    sumVel += speed;

    if (speed < vTol) {
      status = 'converged';
      break;
    }

    // RK4 integration step along (direction * v)
    const sign = direction;
    const k1x = sign * vx, k1y = sign * vy, k1z = sign * vz;

    const [k2vx, k2vy, k2vz] = sampleVelocityField(
      velocityGrid, nx, ny, nz, halfExtent,
      x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, z + 0.5 * dt * k1z
    );
    const k2x = sign * k2vx, k2y = sign * k2vy, k2z = sign * k2vz;

    const [k3vx, k3vy, k3vz] = sampleVelocityField(
      velocityGrid, nx, ny, nz, halfExtent,
      x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, z + 0.5 * dt * k2z
    );
    const k3x = sign * k3vx, k3y = sign * k3vy, k3z = sign * k3vz;

    const [k4vx, k4vy, k4vz] = sampleVelocityField(
      velocityGrid, nx, ny, nz, halfExtent,
      x + dt * k3x, y + dt * k3y, z + dt * k3z
    );
    const k4x = sign * k4vx, k4y = sign * k4vy, k4z = sign * k4vz;

    const stepX = (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    const stepY = (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
    const stepZ = (dt / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);

    const stepDist = Math.hypot(stepX, stepY, stepZ);
    arcLength += stepDist;

    x += stepX;
    y += stepY;
    z += stepZ;

    points.push([x, y, z]);
  }

  // Compute mean curvature of polyline
  let sumCurvature = 0;
  let curvatureSegments = 0;
  for (let i = 1; i < points.length - 1; i++) {
    const p0 = points[i - 1];
    const p1 = points[i];
    const p2 = points[i + 1];

    const v1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    const v2 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    const d1 = Math.hypot(v1[0], v1[1], v1[2]);
    const d2 = Math.hypot(v2[0], v2[1], v2[2]);

    if (d1 > 1e-6 && d2 > 1e-6) {
      const dot = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) / (d1 * d2);
      const clampedDot = Math.max(-1.0, Math.min(1.0, dot));
      const angle = Math.acos(clampedDot);
      sumCurvature += angle / (0.5 * (d1 + d2));
      curvatureSegments++;
    }
  }

  return {
    points,
    targetNodeId,
    arcLength,
    meanVelocity: points.length > 0 ? sumVel / points.length : 0,
    meanCurvature: curvatureSegments > 0 ? sumCurvature / curvatureSegments : 0,
    status
  };
}

/**
 * Builds the complete Cosmic Web Topology Graph from velocity field and critical points.
 *
 * @param {Float32Array} velocityGrid
 * @param {number} nx
 * @param {number} ny
 * @param {number} nz
 * @param {number} halfExtent
 * @param {Array<{ root: [number, number, number], jacobian: number[][] }>} criticalRoots
 * @param {Object} [options={}]
 * @param {number} [options.eigenDisplacementMpc=0.5] - Initial displacement epsilon along eigenvector (Mpc/h)
 * @param {number} [options.captureRadiusMpc=8.0] - Attractor capture radius (Mpc/h)
 * @returns {CosmicWebGraph}
 */
export function buildCosmicWebGraph(
  velocityGrid,
  nx,
  ny,
  nz,
  halfExtent,
  criticalRoots,
  options = {}
) {
  const eps = options.eigenDisplacementMpc ?? 0.5;

  // Step 1: Classify all critical roots into topology nodes
  const nodes = [];
  const nodeMap = {};
  const attractors = [];
  const repellers = [];
  const filamentSaddles = [];
  const wallSaddles = [];

  for (let i = 0; i < criticalRoots.length; i++) {
    const cand = criticalRoots[i];
    const profile = analyzeVelocityTensor(cand.jacobian);
    const id = `${profile.type.substring(0, 3)}_${i + 1}`;

    const node = {
      id,
      type: profile.type,
      morseIndex: profile.morseIndex,
      position: [...cand.root],
      eigenvalues: profile.eigenvalues,
      eigenvectors: profile.eigenvectors,
      label: `${profile.type} #${i + 1}`
    };

    nodes.push(node);
    nodeMap[id] = node;

    if (profile.type === 'ATTRACTOR') attractors.push(node);
    else if (profile.type === 'REPELLER') repellers.push(node);
    else if (profile.type === 'SADDLE_FILAMENT') filamentSaddles.push(node);
    else if (profile.type === 'SADDLE_WALL') wallSaddles.push(node);
  }

  // Step 2: Trace Separatrices from Saddles
  const edges = [];
  const adjacency = {};
  for (const n of nodes) {
    adjacency[n.id] = [];
  }

  // Trace 1D Filament Spines from SADDLE_FILAMENT (unstable eigenvector e3 with lambda3 > 0)
  for (const saddle of filamentSaddles) {
    const e3 = saddle.eigenvectors[2]; // expanding direction

    // Forward direction +e3
    const startPlus = [
      saddle.position[0] + eps * e3[0],
      saddle.position[1] + eps * e3[1],
      saddle.position[2] + eps * e3[2]
    ];
    const tracePlus = traceManifoldCurve(
      velocityGrid, nx, ny, nz, halfExtent,
      startPlus, 1, attractors, options
    );
    const edgeId1 = `EDGE_${saddle.id}_PLUS`;
    const target1 = tracePlus.targetNodeId || 'BOUNDARY';
    edges.push({
      id: edgeId1,
      sourceId: saddle.id,
      targetId: target1,
      edgeType: 'FILAMENT_SPINE',
      points: [[...saddle.position], ...tracePlus.points],
      arcLength: tracePlus.arcLength,
      meanVelocity: tracePlus.meanVelocity,
      meanCurvature: tracePlus.meanCurvature,
      status: tracePlus.status
    });
    if (tracePlus.targetNodeId) {
      adjacency[saddle.id].push(tracePlus.targetNodeId);
      adjacency[tracePlus.targetNodeId].push(saddle.id);
    }

    // Forward direction -e3
    const startMinus = [
      saddle.position[0] - eps * e3[0],
      saddle.position[1] - eps * e3[1],
      saddle.position[2] - eps * e3[2]
    ];
    const traceMinus = traceManifoldCurve(
      velocityGrid, nx, ny, nz, halfExtent,
      startMinus, 1, attractors, options
    );
    const edgeId2 = `EDGE_${saddle.id}_MINUS`;
    const target2 = traceMinus.targetNodeId || 'BOUNDARY';
    edges.push({
      id: edgeId2,
      sourceId: saddle.id,
      targetId: target2,
      edgeType: 'FILAMENT_SPINE',
      points: [[...saddle.position], ...traceMinus.points],
      arcLength: traceMinus.arcLength,
      meanVelocity: traceMinus.meanVelocity,
      meanCurvature: traceMinus.meanCurvature,
      status: traceMinus.status
    });
    if (traceMinus.targetNodeId) {
      adjacency[saddle.id].push(traceMinus.targetNodeId);
      adjacency[traceMinus.targetNodeId].push(saddle.id);
    }
  }

  // Trace 1D Void Drains / Sheet normal from SADDLE_WALL (stable eigenvector e1 with lambda1 < 0)
  for (const saddle of wallSaddles) {
    const e1 = saddle.eigenvectors[0]; // contracting direction

    // Backward direction +e1
    const startPlus = [
      saddle.position[0] + eps * e1[0],
      saddle.position[1] + eps * e1[1],
      saddle.position[2] + eps * e1[2]
    ];
    const tracePlus = traceManifoldCurve(
      velocityGrid, nx, ny, nz, halfExtent,
      startPlus, -1, repellers, options
    );
    const edgeId1 = `EDGE_${saddle.id}_DRAIN_PLUS`;
    const target1 = tracePlus.targetNodeId || 'BOUNDARY';
    edges.push({
      id: edgeId1,
      sourceId: saddle.id,
      targetId: target1,
      edgeType: 'VOID_DRAIN',
      points: [[...saddle.position], ...tracePlus.points],
      arcLength: tracePlus.arcLength,
      meanVelocity: tracePlus.meanVelocity,
      meanCurvature: tracePlus.meanCurvature,
      status: tracePlus.status
    });
    if (tracePlus.targetNodeId) {
      adjacency[saddle.id].push(tracePlus.targetNodeId);
      adjacency[tracePlus.targetNodeId].push(saddle.id);
    }
  }

  return {
    nodes,
    edges,
    adjacency,
    nodeMap
  };
}

/**
 * Finds shortest topological filament path between two critical nodes using Dijkstra's algorithm.
 *
 * @param {CosmicWebGraph} graph
 * @param {string} startNodeId
 * @param {string} endNodeId
 * @returns {{ path: string[], totalDistanceMpc: number, edges: TopologyEdge[] }|null}
 */
export function findShortestFilamentPath(graph, startNodeId, endNodeId) {
  const { nodeMap, adjacency, edges } = graph;
  if (!nodeMap[startNodeId] || !nodeMap[endNodeId]) return null;

  // Build edge distance lookup
  const edgeDistMap = {};
  for (const e of edges) {
    if (e.targetId !== 'BOUNDARY') {
      const key1 = `${e.sourceId}->${e.targetId}`;
      const key2 = `${e.targetId}->${e.sourceId}`;
      edgeDistMap[key1] = e;
      edgeDistMap[key2] = e;
    }
  }

  const dist = {};
  const prev = {};
  const unvisited = new Set(Object.keys(nodeMap));

  for (const id of unvisited) {
    dist[id] = Infinity;
  }
  dist[startNodeId] = 0;

  while (unvisited.size > 0) {
    // Find min dist node
    let current = null;
    let minDist = Infinity;
    for (const id of unvisited) {
      if (dist[id] < minDist) {
        minDist = dist[id];
        current = id;
      }
    }

    if (current === null || minDist === Infinity) break;
    if (current === endNodeId) break;

    unvisited.delete(current);

    const neighbors = adjacency[current] || [];
    for (const neighbor of neighbors) {
      if (!unvisited.has(neighbor)) continue;

      const edgeObj = edgeDistMap[`${current}->${neighbor}`];
      const weight = edgeObj ? edgeObj.arcLength : 1.0;
      const alt = dist[current] + weight;

      if (alt < dist[neighbor]) {
        dist[neighbor] = alt;
        prev[neighbor] = current;
      }
    }
  }

  if (dist[endNodeId] === Infinity) return null;

  // Reconstruct path
  const path = [];
  let curr = endNodeId;
  while (curr) {
    path.unshift(curr);
    curr = prev[curr];
  }

  const pathEdges = [];
  for (let i = 0; i < path.length - 1; i++) {
    const e = edgeDistMap[`${path[i]}->${path[i+1]}`];
    if (e) pathEdges.push(e);
  }

  return {
    path,
    totalDistanceMpc: dist[endNodeId],
    edges: pathEdges
  };
}
