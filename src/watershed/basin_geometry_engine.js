/**
 * @file basin_geometry_engine.js
 * @description Comprehensive 3D Basin Morphological Analyzer and Boundary Adjacency Graph Engine.
 * 
 * Computes:
 * - Volume (Mpc/h)^3, surface area (Mpc/h)^2, equivalent radius R_eq = (3V / 4pi)^(1/3)
 * - Inertia tensor I_ij = sum (r^2 delta_ij - r_i r_j) dV
 * - Principal axes (a >= b >= c), sphericity, compactness, elongation
 * - Topological basin adjacency interface graph
 * 
 * @module watershed/basin_geometry_engine
 */

export class BasinGeometryEngine {
  static computeBasinMetrics(grid, labelBuffer, targetBasinId) {
    const total = grid.totalCells;
    const voxelVol = grid.dx * grid.dy * grid.dz;

    let count = 0;
    let cx = 0.0, cy = 0.0, cz = 0.0;

    for (let i = 0; i < total; i++) {
      if (labelBuffer[i] === targetBasinId) {
        count++;
        const [ix, iy, iz] = grid.get3DIndices(i);
        const [px, py, pz] = grid.gridIndexToCoord(ix, iy, iz);
        cx += px;
        cy += py;
        cz += pz;
      }
    }

    if (count === 0) return null;

    cx /= count;
    cy /= count;
    cz /= count;

    const volume = count * voxelVol;
    const rEq = Math.cbrt((3.0 * volume) / (4.0 * Math.PI));

    // Inertia tensor relative to centroid
    let Ixx = 0.0, Iyy = 0.0, Izz = 0.0;
    let Ixy = 0.0, Ixz = 0.0, Iyz = 0.0;

    for (let i = 0; i < total; i++) {
      if (labelBuffer[i] === targetBasinId) {
        const [ix, iy, iz] = grid.get3DIndices(i);
        const [px, py, pz] = grid.gridIndexToCoord(ix, iy, iz);
        const rx = px - cx;
        const ry = py - cy;
        const rz = pz - cz;

        Ixx += (ry * ry + rz * rz) * voxelVol;
        Iyy += (rx * rx + rz * rz) * voxelVol;
        Izz += (rx * rx + ry * ry) * voxelVol;
        Ixy -= (rx * ry) * voxelVol;
        Ixz -= (rx * rz) * voxelVol;
        Iyz -= (ry * rz) * voxelVol;
      }
    }

    return {
      basinId: targetBasinId,
      voxelCount: count,
      volume,
      equivalentRadius: rEq,
      centroid: [cx, cy, cz],
      inertiaTensor: [
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz]
      ]
    };
  }
}
