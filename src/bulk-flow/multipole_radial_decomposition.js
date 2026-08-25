/**
 * @file multipole_radial_decomposition.js
 * @description Radial Multipolar Expansion Suite: Monopole, Dipole, Quadrupole (Shear), and Octupole Tensors out to R = 300 Mpc/h.
 *
 * Implements:
 * 1. Monopole expansion M(R): Mean radial inflow/outflow and Hubble bubble fluctuation.
 * 2. Dipole vector V_bulk(R): Bulk flow velocity vector and direction angles (SGL, SGB, l, b).
 * 3. Quadrupole cosmic shear tensor Q_ij(R): 3x3 symmetric trace-free tensor, principal shear axes, eigenvalues.
 * 4. Octupole 3rd-order tensor O_ijk(R): 7 independent degrees of freedom representing octupolar flow patterns.
 * 5. Continuous radial profile evaluation from R = 10 to R = 300 Mpc/h in spherical top-hat or differential shells.
 *
 * @module bulk-flow/multipole_radial_decomposition
 */

import { jacobiDiagonalize3x3 } from '../topology/eigen_topology.js';
import { cartesianToSupergalacticAngles, vectorAngleDeg, ASTRONOMICAL_APEX_REFERENCES } from './bulk_flow_estimator.js';

/**
 * Computes full multipolar radial decomposition (Monopole, Dipole, Quadrupole, Octupole)
 * on a 3D regular Cartesian velocity grid.
 */
export class MultipoleRadialDecomposer {
  /**
   * @param {Object} grid - Grid indexer with { nx, ny, nz, boxSize, origin, gridIndexToCoord }
   * @param {Float32Array|Float64Array} velocityGrid - Flattened velocity grid [nx * ny * nz * 3] in km/s
   */
  constructor(grid, velocityGrid) {
    this.grid = grid;
    this.velocityGrid = velocityGrid;
  }

  /**
   * Evaluates multipole expansion inside a sphere of radius R (Mpc/h) centered at center.
   *
   * @param {number} R - Sphere radius in Mpc/h (up to 300 Mpc/h)
   * @param {number[]} [center=[0,0,0]] - Sphere center in Supergalactic coordinates
   * @returns {{
   *   radiusMpc: number,
   *   sampleCount: number,
   *   monopole: { radialVelocityKms: number, divergenceKmsMpc: number },
   *   dipole: { vector: [number, number, number], magnitude: number, sglDeg: number, sgbDeg: number, angleToShapleyDeg: number, angleToCMBDipoleDeg: number },
   *   quadrupole: { tensor: number[][], shearMagnitude: number, eigenvalues: number[], principalAxes: number[][] },
   *   octupole: { tensor: number[][][], octupoleMagnitude: number }
   * }}
   */
  decomposeSphere(R, center = [0.0, 0.0, 0.0]) {
    const { nx, ny, nz } = this.grid;
    const rSqMax = R * R;

    let count = 0;
    let sumVr = 0.0;
    let sumVx = 0.0, sumVy = 0.0, sumVz = 0.0;

    // Quadrupole accumulators Q_ij (3x3)
    const Q = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];

    // Octupole accumulators O_ijk (3x3x3)
    const O = Array.from({ length: 3 }, () =>
      Array.from({ length: 3 }, () => [0, 0, 0])
    );

    for (let ix = 0; ix < nx; ix++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let iz = 0; iz < nz; iz++) {
          const [x, y, z] = this.grid.gridIndexToCoord(ix, iy, iz);
          const dx = x - center[0];
          const dy = y - center[1];
          const dz = z - center[2];
          const distSq = dx * dx + dy * dy + dz * dz;

          if (distSq <= rSqMax && distSq > 1e-6) {
            const dist = Math.sqrt(distSq);
            const nHat = [dx / dist, dy / dist, dz / dist];

            const idx = ((ix * ny + iy) * nz + iz) * 3;
            const vx = this.velocityGrid[idx];
            const vy = this.velocityGrid[idx + 1];
            const vz = this.velocityGrid[idx + 2];
            const v = [vx, vy, vz];

            // Radial velocity v_r = v . nHat
            const v_r = vx * nHat[0] + vy * nHat[1] + vz * nHat[2];
            sumVr += v_r;

            // Dipole accumulation
            sumVx += vx;
            sumVy += vy;
            sumVz += vz;

            // Quadrupole accumulation: Q_ij = 0.5 * (v_i n_j + v_j n_i) - (1/3) delta_ij (v . n)
            for (let i = 0; i < 3; i++) {
              for (let j = 0; j < 3; j++) {
                const delta_ij = i === j ? 1.0 : 0.0;
                const q_val = 0.5 * (v[i] * nHat[j] + v[j] * nHat[i]) - (1.0 / 3.0) * delta_ij * v_r;
                Q[i][j] += q_val;
              }
            }

            // Octupole accumulation: symmetric trace-free 3rd order tensor
            for (let i = 0; i < 3; i++) {
              for (let j = 0; j < 3; j++) {
                for (let k = 0; k < 3; k++) {
                  // Symmetrized product of v and nHat
                  const symTerm = (
                    v[i] * nHat[j] * nHat[k] +
                    v[j] * nHat[k] * nHat[i] +
                    v[k] * nHat[i] * nHat[j]
                  ) / 3.0;

                  // Trace subtractions
                  const delta_ij = i === j ? 1.0 : 0.0;
                  const delta_jk = j === k ? 1.0 : 0.0;
                  const delta_ki = k === i ? 1.0 : 0.0;

                  const traceSub = (1.0 / 5.0) * (
                    delta_ij * v_r * nHat[k] +
                    delta_jk * v_r * nHat[i] +
                    delta_ki * v_r * nHat[j]
                  );

                  O[i][j][k] += (symTerm - traceSub);
                }
              }
            }

            count++;
          }
        }
      }
    }

    if (count === 0) {
      return {
        radiusMpc: R,
        sampleCount: 0,
        monopole: { radialVelocityKms: 0, divergenceKmsMpc: 0 },
        dipole: { vector: [0, 0, 0], magnitude: 0, sglDeg: 0, sgbDeg: 0, angleToShapleyDeg: 0, angleToCMBDipoleDeg: 0 },
        quadrupole: { tensor: [[0,0,0],[0,0,0],[0,0,0]], shearMagnitude: 0, eigenvalues: [0,0,0], principalAxes: [[1,0,0],[0,1,0],[0,0,1]] },
        octupole: { tensor: O, octupoleMagnitude: 0 }
      };
    }

    // 1. Monopole
    const v_r_mean = sumVr / count;
    // delta H = 3 * v_r_mean / R
    const deltaH = (3.0 * v_r_mean) / R;

    // 2. Dipole
    const vX = sumVx / count;
    const vY = sumVy / count;
    const vZ = sumVz / count;
    const dipoleMag = Math.hypot(vX, vY, vZ);
    const angles = cartesianToSupergalacticAngles(vX, vY, vZ);
    const angleToShapley = vectorAngleDeg([vX, vY, vZ], ASTRONOMICAL_APEX_REFERENCES.SHAPLEY_CORE.supergalacticKms);
    const angleToCMB = vectorAngleDeg([vX, vY, vZ], ASTRONOMICAL_APEX_REFERENCES.CMB_DIPOLE.supergalacticKms);

    // 3. Quadrupole normalization & diagonalization
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        Q[i][j] /= count;
      }
    }
    const { eigenvalues, eigenvectors } = jacobiDiagonalize3x3(Q);
    let shearFrobSq = 0.0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        shearFrobSq += Q[i][j] * Q[i][j];
      }
    }
    const shearMag = Math.sqrt(shearFrobSq);

    // 4. Octupole normalization & Frobenius magnitude
    let octupoleFrobSq = 0.0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        for (let k = 0; k < 3; k++) {
          O[i][j][k] /= count;
          octupoleFrobSq += O[i][j][k] * O[i][j][k];
        }
      }
    }
    const octupoleMag = Math.sqrt(octupoleFrobSq);

    return {
      radiusMpc: R,
      sampleCount: count,
      monopole: {
        radialVelocityKms: v_r_mean,
        divergenceKmsMpc: deltaH
      },
      dipole: {
        vector: [vX, vY, vZ],
        magnitude: dipoleMag,
        sglDeg: angles.sgl,
        sgbDeg: angles.sgb,
        angleToShapleyDeg: angleToShapley,
        angleToCMBDipoleDeg: angleToCMB
      },
      quadrupole: {
        tensor: Q,
        shearMagnitude: shearMag,
        eigenvalues,
        principalAxes: eigenvectors
      },
      octupole: {
        tensor: O,
        octupoleMagnitude: octupoleMag
      }
    };
  }

  /**
   * Computes full radial multipolar profiles across array of radial shell/sphere radii.
   *
   * @param {number[]} [radii=[20, 40, 60, 80, 100, 150, 200, 250, 300]]
   * @param {number[]} [center=[0,0,0]]
   * @returns {Array<Object>}
   */
  computeRadialMultipoleProfile(radii = [20, 40, 60, 80, 100, 150, 200, 250, 300], center = [0, 0, 0]) {
    return radii.map(r => this.decomposeSphere(r, center));
  }
}
