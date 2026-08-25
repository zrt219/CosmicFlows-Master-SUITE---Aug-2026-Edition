/**
 * @file bulk_flow_estimator.js
 * @description Cosmic Bulk Flow Estimator Suite.
 * Implements triple estimator algorithms:
 * 1. FIELD_VOLUME_WEIGHTED: Uniform voxel volume integration across radial spheres/shells
 * 2. CATALOG_WEIGHTED: Galaxy catalog point tracer projection tensor estimator
 * 3. INVERSE_VARIANCE_WEIGHTED: Minimum-variance / Maximum Likelihood estimator with cosmic dispersion
 *
 * Provides radial shell decomposition, apex dipole angles (to CMB & Shapley), and covariance propagation.
 */

import { invert3x3, matVec3 } from '../topology/root_finder.js';

/**
 * Standard astronomical reference directions for dipole apex comparisons.
 */
export const ASTRONOMICAL_APEX_REFERENCES = Object.freeze({
  CMB_DIPOLE: Object.freeze({
    name: 'CMB Dipole (Planck 2018)',
    galacticL: 264.021,
    galacticB: 48.253,
    velocityKms: 369.82,
    // Supergalactic coordinates in km/s (approximate conversion)
    supergalacticKms: [-245.0, 269.0, -68.0]
  }),
  SHAPLEY_CORE: Object.freeze({
    name: 'Shapley Supercluster Core (A3558)',
    supergalacticMpc: [-170.7, 118.7, -21.3],
    supergalacticKms: [-12734.0, 8855.0, -1589.0]
  }),
  GREAT_ATTRACTOR: Object.freeze({
    name: 'Great Attractor / Norma (A3627)',
    supergalacticMpc: [-46.0, 22.9, -6.8],
    supergalacticKms: [-3432.0, 1708.0, -507.0]
  })
});

/**
 * Transforms Cartesian Supergalactic [x, y, z] to Supergalactic Longitude (SGL) and Latitude (SGB) in degrees.
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @returns {{ r: number, sgl: number, sgb: number }}
 */
export function cartesianToSupergalacticAngles(x, y, z) {
  const r = Math.hypot(x, y, z);
  if (r < 1e-9) {
    return { r: 0, sgl: 0, sgb: 0 };
  }

  // SGB = asin(z / r)
  const sgb = Math.asin(Math.max(-1.0, Math.min(1.0, z / r))) * (180.0 / Math.PI);
  // SGL = atan2(y, x)
  let sgl = Math.atan2(y, x) * (180.0 / Math.PI);
  if (sgl < 0) sgl += 360.0;

  return { r, sgl, sgb };
}

/**
 * Transforms Supergalactic [SGL, SGB] in degrees and distance r to Cartesian [x, y, z].
 * @param {number} r
 * @param {number} sglDeg
 * @param {number} sgbDeg
 * @returns {[number, number, number]}
 */
export function supergalacticAnglesToCartesian(r, sglDeg, sgbDeg) {
  const sglRad = (sglDeg * Math.PI) / 180.0;
  const sgbRad = (sgbDeg * Math.PI) / 180.0;

  const cosB = Math.cos(sgbRad);
  const x = r * cosB * Math.cos(sglRad);
  const y = r * cosB * Math.sin(sglRad);
  const z = r * Math.sin(sgbRad);

  return [x, y, z];
}

/**
 * Computes angular separation in degrees between two 3D vectors.
 * @param {number[]} v1 - [x1, y1, z1]
 * @param {number[]} v2 - [x2, y2, z2]
 * @returns {number} Angle in degrees [0, 180]
 */
export function vectorAngleDeg(v1, v2) {
  const m1 = Math.hypot(v1[0], v1[1], v1[2]);
  const m2 = Math.hypot(v2[0], v2[1], v2[2]);
  if (m1 < 1e-9 || m2 < 1e-9) return 0;

  const dot = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (m1 * m2);
  const clamped = Math.max(-1.0, Math.min(1.0, dot));
  return Math.acos(clamped) * (180.0 / Math.PI);
}

/**
 * @typedef {Object} BulkFlowResult
 * @property {'FIELD_VOLUME_WEIGHTED'|'CATALOG_WEIGHTED'|'INVERSE_VARIANCE_WEIGHTED'} method
 * @property {[number, number, number]} bulkFlowVector - [Vx, Vy, Vz] in km/s
 * @property {number} magnitude - ||V_bulk|| in km/s
 * @property {number} sglDeg - Supergalactic Longitude of apex (degrees)
 * @property {number} sgbDeg - Supergalactic Latitude of apex (degrees)
 * @property {number} angleToShapleyDeg - Angle between bulk flow apex and Shapley Core
 * @property {number} angleToCMBDipoleDeg - Angle between bulk flow apex and CMB Dipole
 * @property {number} sampleCount - Number of voxels or galaxies in estimator
 * @property {number[][]} [covarianceMatrix] - 3x3 error covariance matrix C_ij (for weighted estimators)
 * @property {[number, number, number]} [uncertainty1Sigma] - [σ_Vx, σ_Vy, σ_Vz] (km/s)
 * @property {number} [radialMinMpc] - Inner shell radius in Mpc/h
 * @property {number} [radialMaxMpc] - Outer shell radius in Mpc/h
 */

/**
 * Estimator 1: FIELD_VOLUME_WEIGHTED
 * Computes bulk flow by direct volume-weighted integration of velocity field on a regular 3D grid.
 *
 * @param {Float32Array} velocityGrid - Flattened velocity grid [nx * ny * nz * 3]
 * @param {number} nx
 * @param {number} ny
 * @param {number} nz
 * @param {number} halfExtent - Box half-extent in Mpc/h
 * @param {number} rMinMpc - Minimum radial shell distance in Mpc/h (0 for sphere)
 * @param {number} rMaxMpc - Maximum radial shell distance in Mpc/h
 * @param {number} [hubbleH0=74.6]
 * @returns {BulkFlowResult}
 */
export function estimateBulkFlowFieldVolumeWeighted(
  velocityGrid,
  nx,
  ny,
  nz,
  halfExtent,
  rMinMpc,
  rMaxMpc,
  hubbleH0 = 74.6
) {
  const boxSize = 2.0 * halfExtent;
  let sumVx = 0;
  let sumVy = 0;
  let sumVz = 0;
  let count = 0;

  const rMinSq = rMinMpc * rMinMpc;
  const rMaxSq = rMaxMpc * rMaxMpc;

  for (let ix = 0; ix < nx; ix++) {
    const x = ((ix + 0.5) / nx) * boxSize - halfExtent;
    for (let iy = 0; iy < ny; iy++) {
      const y = ((iy + 0.5) / ny) * boxSize - halfExtent;
      for (let iz = 0; iz < nz; iz++) {
        const z = ((iz + 0.5) / nz) * boxSize - halfExtent;
        const rSq = x * x + y * y + z * z;

        if (rSq >= rMinSq && rSq <= rMaxSq) {
          const idx = ((ix * ny + iy) * nz + iz) * 3;
          sumVx += velocityGrid[idx];
          sumVy += velocityGrid[idx + 1];
          sumVz += velocityGrid[idx + 2];
          count++;
        }
      }
    }
  }

  const vX = count > 0 ? sumVx / count : 0;
  const vY = count > 0 ? sumVy / count : 0;
  const vZ = count > 0 ? sumVz / count : 0;
  const mag = Math.hypot(vX, vY, vZ);
  const angles = cartesianToSupergalacticAngles(vX, vY, vZ);

  const angleToShapley = vectorAngleDeg([vX, vY, vZ], ASTRONOMICAL_APEX_REFERENCES.SHAPLEY_CORE.supergalacticKms);
  const angleToCMB = vectorAngleDeg([vX, vY, vZ], ASTRONOMICAL_APEX_REFERENCES.CMB_DIPOLE.supergalacticKms);

  return {
    method: 'FIELD_VOLUME_WEIGHTED',
    bulkFlowVector: [vX, vY, vZ],
    magnitude: mag,
    sglDeg: angles.sgl,
    sgbDeg: angles.sgb,
    angleToShapleyDeg: angleToShapley,
    angleToCMBDipoleDeg: angleToCMB,
    sampleCount: count,
    radialMinMpc: rMinMpc,
    radialMaxMpc: rMaxMpc
  };
}

/**
 * Estimator 2: CATALOG_WEIGHTED
 * Computes bulk flow from galaxy survey point tracers with line-of-sight peculiar velocities u_n.
 * Solves (sum w_n r_hat_n r_hat_n^T) V_bulk = sum w_n u_n r_hat_n.
 *
 * @param {Array<{ x: number, y: number, z: number, u: number, weight?: number }>} galaxies
 * @param {number} [rMinMpc=0]
 * @param {number} [rMaxMpc=Infinity]
 * @returns {BulkFlowResult}
 */
export function estimateBulkFlowCatalogWeighted(galaxies, rMinMpc = 0, rMaxMpc = Infinity) {
  let A = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
  ];
  let b = [0, 0, 0];
  let count = 0;

  const rMinSq = rMinMpc * rMinMpc;
  const rMaxSq = rMaxMpc * rMaxMpc;

  for (const g of galaxies) {
    const rSq = g.x * g.x + g.y * g.y + g.z * g.z;
    if (rSq < rMinSq || rSq > rMaxSq) continue;

    const r = Math.sqrt(rSq);
    if (r < 1e-6) continue;

    const nx = g.x / r;
    const ny = g.y / r;
    const nz = g.z / r;
    const w = g.weight ?? 1.0;
    const u = g.u; // line-of-sight peculiar velocity

    A[0][0] += w * nx * nx;
    A[0][1] += w * nx * ny;
    A[0][2] += w * nx * nz;

    A[1][0] += w * ny * nx;
    A[1][1] += w * ny * ny;
    A[1][2] += w * ny * nz;

    A[2][0] += w * nz * nx;
    A[2][1] += w * nz * ny;
    A[2][2] += w * nz * nz;

    b[0] += w * u * nx;
    b[1] += w * u * ny;
    b[2] += w * u * nz;

    count++;
  }

  if (count < 3) {
    return {
      method: 'CATALOG_WEIGHTED',
      bulkFlowVector: [0, 0, 0],
      magnitude: 0,
      sglDeg: 0,
      sgbDeg: 0,
      angleToShapleyDeg: 0,
      angleToCMBDipoleDeg: 0,
      sampleCount: count,
      radialMinMpc: rMinMpc,
      radialMaxMpc: rMaxMpc
    };
  }

  const invA = invert3x3(A);
  if (!invA) {
    throw new Error('Singular tensor in catalog bulk flow estimator: insufficient geometric sky coverage.');
  }

  const V = matVec3(invA, b);
  const mag = Math.hypot(V[0], V[1], V[2]);
  const angles = cartesianToSupergalacticAngles(V[0], V[1], V[2]);

  const angleToShapley = vectorAngleDeg(V, ASTRONOMICAL_APEX_REFERENCES.SHAPLEY_CORE.supergalacticKms);
  const angleToCMB = vectorAngleDeg(V, ASTRONOMICAL_APEX_REFERENCES.CMB_DIPOLE.supergalacticKms);

  return {
    method: 'CATALOG_WEIGHTED',
    bulkFlowVector: V,
    magnitude: mag,
    sglDeg: angles.sgl,
    sgbDeg: angles.sgb,
    angleToShapleyDeg: angleToShapley,
    angleToCMBDipoleDeg: angleToCMB,
    sampleCount: count,
    radialMinMpc: rMinMpc,
    radialMaxMpc: rMaxMpc
  };
}

/**
 * Estimator 3: INVERSE_VARIANCE_WEIGHTED (Maximum Likelihood / Minimum Variance)
 * Weights galaxies by 1 / (sigma_u^2 + sigma_star^2) where sigma_star is cosmic thermal dispersion.
 *
 * @param {Array<{ x: number, y: number, z: number, u: number, error: number }>} galaxies
 * @param {number} [rMinMpc=0]
 * @param {number} [rMaxMpc=Infinity]
 * @param {number} [cosmicDispersionKms=187.0] - Thermal velocity dispersion sigma_* (km/s)
 * @returns {BulkFlowResult}
 */
export function estimateBulkFlowInverseVarianceWeighted(
  galaxies,
  rMinMpc = 0,
  rMaxMpc = Infinity,
  cosmicDispersionKms = 187.0
) {
  let A = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
  ];
  let b = [0, 0, 0];
  let count = 0;

  const rMinSq = rMinMpc * rMinMpc;
  const rMaxSq = rMaxMpc * rMaxMpc;
  const sigStarSq = cosmicDispersionKms * cosmicDispersionKms;

  for (const g of galaxies) {
    const rSq = g.x * g.x + g.y * g.y + g.z * g.z;
    if (rSq < rMinSq || rSq > rMaxSq) continue;

    const r = Math.sqrt(rSq);
    if (r < 1e-6) continue;

    const nx = g.x / r;
    const ny = g.y / r;
    const nz = g.z / r;
    const sigMeas = g.error ?? 150.0;
    const totalVar = sigMeas * sigMeas + sigStarSq;
    const w = 1.0 / totalVar;
    const u = g.u;

    A[0][0] += w * nx * nx;
    A[0][1] += w * nx * ny;
    A[0][2] += w * nx * nz;

    A[1][0] += w * ny * nx;
    A[1][1] += w * ny * ny;
    A[1][2] += w * ny * nz;

    A[2][0] += w * nz * nx;
    A[2][1] += w * nz * ny;
    A[2][2] += w * nz * nz;

    b[0] += w * u * nx;
    b[1] += w * u * ny;
    b[2] += w * u * nz;

    count++;
  }

  if (count < 3) {
    return {
      method: 'INVERSE_VARIANCE_WEIGHTED',
      bulkFlowVector: [0, 0, 0],
      magnitude: 0,
      sglDeg: 0,
      sgbDeg: 0,
      angleToShapleyDeg: 0,
      angleToCMBDipoleDeg: 0,
      sampleCount: count,
      radialMinMpc: rMinMpc,
      radialMaxMpc: rMaxMpc,
      covarianceMatrix: [[0,0,0],[0,0,0],[0,0,0]],
      uncertainty1Sigma: [0, 0, 0]
    };
  }

  const cov = invert3x3(A);
  if (!cov) {
    throw new Error('Singular covariance matrix in ML bulk flow estimator.');
  }

  const V = matVec3(cov, b);
  const mag = Math.hypot(V[0], V[1], V[2]);
  const angles = cartesianToSupergalacticAngles(V[0], V[1], V[2]);

  const sigmaVx = Math.sqrt(Math.max(0, cov[0][0]));
  const sigmaVy = Math.sqrt(Math.max(0, cov[1][1]));
  const sigmaVz = Math.sqrt(Math.max(0, cov[2][2]));

  const angleToShapley = vectorAngleDeg(V, ASTRONOMICAL_APEX_REFERENCES.SHAPLEY_CORE.supergalacticKms);
  const angleToCMB = vectorAngleDeg(V, ASTRONOMICAL_APEX_REFERENCES.CMB_DIPOLE.supergalacticKms);

  return {
    method: 'INVERSE_VARIANCE_WEIGHTED',
    bulkFlowVector: V,
    magnitude: mag,
    sglDeg: angles.sgl,
    sgbDeg: angles.sgb,
    angleToShapleyDeg: angleToShapley,
    angleToCMBDipoleDeg: angleToCMB,
    sampleCount: count,
    radialMinMpc: rMinMpc,
    radialMaxMpc: rMaxMpc,
    covarianceMatrix: cov,
    uncertainty1Sigma: [sigmaVx, sigmaVy, sigmaVz]
  };
}

/**
 * Computes radial shell profile decomposition of bulk flow across consecutive radial shells.
 *
 * @param {Array<number>} shellEdgesMpc - Radial boundary edges in Mpc/h, e.g. [0, 20, 40, 60, 80, 100, 150, 200]
 * @param {'FIELD'|'CATALOG'} inputType
 * @param {Object} data - Input grid or galaxy catalog
 * @param {Object} [options={}]
 * @returns {Array<BulkFlowResult & { shellLabel: string, rMidMpc: number }>}
 */
export function computeRadialShellProfile(shellEdgesMpc, inputType, data, options = {}) {
  const profile = [];

  for (let i = 0; i < shellEdgesMpc.length - 1; i++) {
    const rMin = shellEdgesMpc[i];
    const rMax = shellEdgesMpc[i + 1];
    const rMid = 0.5 * (rMin + rMax);
    const shellLabel = `${rMin}-${rMax} Mpc/h`;

    let res;
    if (inputType === 'FIELD') {
      const { velocityGrid, nx, ny, nz, halfExtent, hubbleH0 } = data;
      res = estimateBulkFlowFieldVolumeWeighted(velocityGrid, nx, ny, nz, halfExtent, rMin, rMax, hubbleH0);
    } else {
      const { galaxies, cosmicDispersionKms, useInverseVariance } = data;
      if (useInverseVariance) {
        res = estimateBulkFlowInverseVarianceWeighted(galaxies, rMin, rMax, cosmicDispersionKms);
      } else {
        res = estimateBulkFlowCatalogWeighted(galaxies, rMin, rMax);
      }
    }

    profile.push({
      ...res,
      shellLabel,
      rMidMpc: rMid
    });
  }

  return profile;
}
