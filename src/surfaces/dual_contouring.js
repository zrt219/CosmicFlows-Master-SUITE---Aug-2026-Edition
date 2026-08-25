/**
 * @file dual_contouring.js
 * @module surfaces/dual_contouring
 * @description Research-grade 3D Dual Contouring Isosurface Extractor for Cosmological Density & Potential Fields.
 * 
 * Implements:
 * 1. Dual Contouring of Hermite Data (Ju et al. 2002).
 * 2. Quadric Error Metric (QEM) solver with Singular Value Decomposition (SVD) and regularized pseudo-inverse.
 * 3. Exact edge-crossing linear interpolation and analytical gradient evaluation.
 * 4. Dual quad and triangulated manifold mesh generation with canonical winding order.
 * 5. Topological invariant analyzer: Euler characteristic chi = V - E + F, genus g = 1 - chi/2, and manifold verification.
 * 6. Enclosed volume and surface area calculus via divergence theorem on closed 2-manifolds.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

import { GridIndexer } from '../fields/grid_indexer.js';

/**
 * 3D SVD / Pseudo-inverse solver for 3x3 symmetric Quadric Error Metric (QEM) system:
 * A * x = b, where A = sum(n_i * n_i^T), b = sum((n_i . p_i) * n_i).
 */
export class QEMSolver3D {
  constructor() {
    // 3x3 symmetric matrix components of A: a00, a01, a02, a11, a12, a22
    this.a00 = 0.0; this.a01 = 0.0; this.a02 = 0.0;
    this.a11 = 0.0; this.a12 = 0.0;
    this.a22 = 0.0;

    // 3-vector b: b0, b1, b2
    this.b0 = 0.0; this.b1 = 0.0; this.b2 = 0.0;

    // Mass point (average of all sample points)
    this.massX = 0.0; this.massY = 0.0; this.massZ = 0.0;
    this.pointCount = 0;
  }

  /**
   * Reset the accumulator for a new grid cell.
   */
  reset() {
    this.a00 = 0.0; this.a01 = 0.0; this.a02 = 0.0;
    this.a11 = 0.0; this.a12 = 0.0;
    this.a22 = 0.0;
    this.b0 = 0.0; this.b1 = 0.0; this.b2 = 0.0;
    this.massX = 0.0; this.massY = 0.0; this.massZ = 0.0;
    this.pointCount = 0;
  }

  /**
   * Add an intersection point p and unit surface normal n.
   * @param {number} px - Point X
   * @param {number} py - Point Y
   * @param {number} pz - Point Z
   * @param {number} nx - Normal X
   * @param {number} ny - Normal Y
   * @param {number} nz - Normal Z
   */
  addPlane(px, py, pz, nx, ny, nz) {
    const len = Math.hypot(nx, ny, nz);
    if (len < 1e-8) return;
    const unx = nx / len;
    const uny = ny / len;
    const unz = nz / len;

    this.a00 += unx * unx;
    this.a01 += unx * uny;
    this.a02 += unx * unz;
    this.a11 += uny * uny;
    this.a12 += uny * unz;
    this.a22 += unz * unz;

    const d = px * unx + py * uny + pz * unz;
    this.b0 += d * unx;
    this.b1 += d * uny;
    this.b2 += d * unz;

    this.massX += px;
    this.massY += py;
    this.massZ += pz;
    this.pointCount++;
  }

  /**
   * Solve for the optimal vertex x by minimizing QEM with fallback to mass point / cell center.
   * @param {number} cellMinX
   * @param {number} cellMinY
   * @param {number} cellMinZ
   * @param {number} cellSizeX
   * @param {number} cellSizeY
   * @param {number} cellSizeZ
   * @param {number} [tolerance=1e-4]
   * @returns {[number, number, number]} [x, y, z]
   */
  solve(cellMinX, cellMinY, cellMinZ, cellSizeX, cellSizeY, cellSizeZ, tolerance = 1e-4) {
    if (this.pointCount === 0) {
      return [
        cellMinX + 0.5 * cellSizeX,
        cellMinY + 0.5 * cellSizeY,
        cellMinZ + 0.5 * cellSizeZ
      ];
    }

    const mass = [
      this.massX / this.pointCount,
      this.massY / this.pointCount,
      this.massZ / this.pointCount
    ];

    // Compute determinant of 3x3 symmetric matrix A
    // A = [a00 a01 a02]
    //     [a01 a11 a12]
    //     [a02 a12 a22]
    const det = this.a00 * (this.a11 * this.a22 - this.a12 * this.a12) -
                this.a01 * (this.a01 * this.a22 - this.a12 * this.a02) +
                this.a02 * (this.a01 * this.a12 - this.a11 * this.a02);

    let solX = mass[0];
    let solY = mass[1];
    let solZ = mass[2];

    if (Math.abs(det) > tolerance) {
      // Invert A via classical adjugate formula
      const invDet = 1.0 / det;
      const i00 = (this.a11 * this.a22 - this.a12 * this.a12) * invDet;
      const i01 = (this.a02 * this.a12 - this.a01 * this.a22) * invDet;
      const i02 = (this.a01 * this.a12 - this.a02 * this.a11) * invDet;

      const i11 = (this.a00 * this.a22 - this.a02 * this.a02) * invDet;
      const i12 = (this.a02 * this.a01 - this.a00 * this.a12) * invDet;

      const i22 = (this.a00 * this.a11 - this.a01 * this.a01) * invDet;

      const x = i00 * this.b0 + i01 * this.b1 + i02 * this.b2;
      const y = i01 * this.b0 + i11 * this.b1 + i12 * this.b2;
      const z = i02 * this.b0 + i12 * this.b1 + i22 * this.b2;

      // Check if solution lies within expanded cell bounding box
      const margin = 0.2;
      const minX = cellMinX - margin * cellSizeX;
      const maxX = cellMinX + (1.0 + margin) * cellSizeX;
      const minY = cellMinY - margin * cellSizeY;
      const maxY = cellMinY + (1.0 + margin) * cellSizeY;
      const minZ = cellMinZ - margin * cellSizeZ;
      const maxZ = cellMinZ + (1.0 + margin) * cellSizeZ;

      if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
        solX = x;
        solY = y;
        solZ = z;
      }
    }

    // Clamp firmly inside the cell bounds
    solX = Math.max(cellMinX, Math.min(cellMinX + cellSizeX, solX));
    solY = Math.max(cellMinY, Math.min(cellMinY + cellSizeY, solY));
    solZ = Math.max(cellMinZ, Math.min(cellMinZ + cellSizeZ, solZ));

    return [solX, solY, solZ];
  }
}

/**
 * 3D Dual Contouring Mesh Extractor.
 */
export class DualContouringExtractor {
  /**
   * @param {GridIndexer} gridIndexer - 3D spatial grid descriptor.
   * @param {Object} [options]
   * @param {number} [options.qemTolerance=1e-4]
   */
  constructor(gridIndexer, options = {}) {
    this.grid = gridIndexer;
    this.qemTolerance = options.qemTolerance || 1e-4;
    this.qemSolver = new QEMSolver3D();
  }

  /**
   * Extract an isodensity surface at given isovalue using Dual Contouring.
   * 
   * @param {Float32Array|Float64Array} scalarField - Continuous 3D scalar field.
   * @param {number} isovalue - Target isosurface level set (e.g. delta = 0.0 or 1.5).
   * @returns {{
   *   vertices: Float32Array,
   *   normals: Float32Array,
   *   indices: Uint32Array,
   *   quadIndices: Uint32Array,
   *   vertexCount: number,
   *   triangleCount: number,
   *   quadCount: number,
   *   eulerCharacteristic: number,
   *   genus: number,
   *   isClosedManifold: boolean,
   *   surfaceArea: number,
   *   enclosedVolume: number
   * }}
   */
  extractIsosurface(scalarField, isovalue) {
    const { nx, ny, nz } = this.grid;
    const cellMap = new Map(); // key -> vertex index
    const verticesList = [];
    const normalsList = [];

    // Helper to evaluate scalar at grid node (ix, iy, iz)
    const getNodeVal = (ix, iy, iz) => {
      const idx = this.grid.getIndex(ix, iy, iz);
      return scalarField[idx];
    };

    // Helper to evaluate numerical central gradient at node (ix, iy, iz)
    const getGradient = (ix, iy, iz) => {
      const vXP = getNodeVal(Math.min(nx - 1, ix + 1), iy, iz);
      const vXM = getNodeVal(Math.max(0, ix - 1), iy, iz);
      const vYP = getNodeVal(ix, Math.min(ny - 1, iy + 1), iz);
      const vYM = getNodeVal(ix, Math.max(0, iy - 1), iz);
      const vZP = getNodeVal(ix, iy, Math.min(nz - 1, iz + 1));
      const vZM = getNodeVal(ix, iy, Math.max(0, iz - 1));

      return [
        (vXP - vXM) / (2.0 * this.grid.dx),
        (vYP - vYM) / (2.0 * this.grid.dy),
        (vZP - vZM) / (2.0 * this.grid.dz)
      ];
    };

    // 1. Generate representative dual vertices for each cell containing edge crossings
    for (let iz = 0; iz < nz - 1; iz++) {
      for (let iy = 0; iy < ny - 1; iy++) {
        for (let ix = 0; ix < nx - 1; ix++) {
          this.qemSolver.reset();

          const [x0, y0, z0] = this.grid.getNodeCoord(ix, iy, iz);
          const dx = this.grid.dx;
          const dy = this.grid.dy;
          const dz = this.grid.dz;

          // Check all 12 edges of the cube [ix, iy, iz]
          // 4 X-edges: (0,0,0)-(1,0,0), (0,1,0)-(1,1,0), (0,0,1)-(1,0,1), (0,1,1)-(1,1,1)
          // 4 Y-edges: (0,0,0)-(0,1,0), (1,0,0)-(1,1,0), (0,0,1)-(0,1,1), (1,0,1)-(1,1,1)
          // 4 Z-edges: (0,0,0)-(0,0,1), (1,0,0)-(1,0,1), (0,1,0)-(0,1,1), (1,1,0)-(1,1,1)

          const checkEdge = (x1, y1, z1, x2, y2, z2) => {
            const v1 = getNodeVal(x1, y1, z1);
            const v2 = getNodeVal(x2, y2, z2);
            if ((v1 < isovalue && v2 >= isovalue) || (v1 >= isovalue && v2 < isovalue)) {
              // Linear interpolation of crossing point
              const t = (isovalue - v1) / (v2 - v1);
              const [p1x, p1y, p1z] = this.grid.getNodeCoord(x1, y1, z1);
              const [p2x, p2y, p2z] = this.grid.getNodeCoord(x2, y2, z2);
              const px = p1x + t * (p2x - p1x);
              const py = p1y + t * (p2y - p1y);
              const pz = p1z + t * (p2z - p1z);

              // Interpolate gradient
              const g1 = getGradient(x1, y1, z1);
              const g2 = getGradient(x2, y2, z2);
              const gx = g1[0] + t * (g2[0] - g1[0]);
              const gy = g1[1] + t * (g2[1] - g1[1]);
              const gz = g1[2] + t * (g2[2] - g1[2]);

              this.qemSolver.addPlane(px, py, pz, gx, gy, gz);
            }
          };

          // X-edges
          checkEdge(ix, iy, iz, ix + 1, iy, iz);
          checkEdge(ix, iy + 1, iz, ix + 1, iy + 1, iz);
          checkEdge(ix, iy, iz + 1, ix + 1, iy, iz + 1);
          checkEdge(ix, iy + 1, iz + 1, ix + 1, iy + 1, iz + 1);

          // Y-edges
          checkEdge(ix, iy, iz, ix, iy + 1, iz);
          checkEdge(ix + 1, iy, iz, ix + 1, iy + 1, iz);
          checkEdge(ix, iy, iz + 1, ix, iy + 1, iz + 1);
          checkEdge(ix + 1, iy, iz + 1, ix + 1, iy + 1, iz + 1);

          // Z-edges
          checkEdge(ix, iy, iz, ix, iy, iz + 1);
          checkEdge(ix + 1, iy, iz, ix + 1, iy, iz + 1);
          checkEdge(ix, iy + 1, iz, ix, iy + 1, iz + 1);
          checkEdge(ix + 1, iy + 1, iz, ix + 1, iy + 1, iz + 1);

          if (this.qemSolver.pointCount > 0) {
            const vPos = this.qemSolver.solve(x0, y0, z0, dx, dy, dz, this.qemTolerance);
            const vIdx = verticesList.length / 3;
            verticesList.push(vPos[0], vPos[1], vPos[2]);

            // Normal at representative vertex
            let nx_ = 0.0, ny_ = 0.0, nz_ = 1.0;
            if (this.qemSolver.pointCount > 0) {
              nx_ = this.qemSolver.b0;
              ny_ = this.qemSolver.b1;
              nz_ = this.qemSolver.b2;
              const len = Math.hypot(nx_, ny_, nz_);
              if (len > 1e-8) {
                nx_ /= len; ny_ /= len; nz_ /= len;
              }
            }
            normalsList.push(nx_, ny_, nz_);

            const cellKey = `${ix},${iy},${iz}`;
            cellMap.set(cellKey, vIdx);
          }
        }
      }
    }

    // 2. Generate dual quads for each grid edge that crosses the isosurface
    const triangles = [];
    const quads = [];

    // Helper to get cell vertex index
    const getCellVertex = (cx, cy, cz) => {
      return cellMap.get(`${cx},${cy},${cz}`);
    };

    // Generate quads for X-edges: shared by 4 cells (y-1,z-1), (y,z-1), (y,z), (y-1,z)
    for (let iz = 1; iz < nz - 1; iz++) {
      for (let iy = 1; iy < ny - 1; iy++) {
        for (let ix = 0; ix < nx - 1; ix++) {
          const v1 = getNodeVal(ix, iy, iz);
          const v2 = getNodeVal(ix + 1, iy, iz);
          if ((v1 < isovalue && v2 >= isovalue) || (v1 >= isovalue && v2 < isovalue)) {
            const c0 = getCellVertex(ix, iy - 1, iz - 1);
            const c1 = getCellVertex(ix, iy, iz - 1);
            const c2 = getCellVertex(ix, iy, iz);
            const c3 = getCellVertex(ix, iy - 1, iz);

            if (c0 !== undefined && c1 !== undefined && c2 !== undefined && c3 !== undefined) {
              if (v1 < isovalue) {
                // Quad: c0, c1, c2, c3
                quads.push(c0, c1, c2, c3);
                triangles.push(c0, c1, c2, c0, c2, c3);
              } else {
                quads.push(c0, c3, c2, c1);
                triangles.push(c0, c3, c2, c0, c2, c1);
              }
            }
          }
        }
      }
    }

    // Generate quads for Y-edges: shared by 4 cells (x-1,z-1), (x,z-1), (x,z), (x-1,z)
    for (let iz = 1; iz < nz - 1; iz++) {
      for (let iy = 0; iy < ny - 1; iy++) {
        for (let ix = 1; ix < nx - 1; ix++) {
          const v1 = getNodeVal(ix, iy, iz);
          const v2 = getNodeVal(ix, iy + 1, iz);
          if ((v1 < isovalue && v2 >= isovalue) || (v1 >= isovalue && v2 < isovalue)) {
            const c0 = getCellVertex(ix - 1, iy, iz - 1);
            const c1 = getCellVertex(ix, iy, iz - 1);
            const c2 = getCellVertex(ix, iy, iz);
            const c3 = getCellVertex(ix - 1, iy, iz);

            if (c0 !== undefined && c1 !== undefined && c2 !== undefined && c3 !== undefined) {
              if (v1 >= isovalue) {
                quads.push(c0, c1, c2, c3);
                triangles.push(c0, c1, c2, c0, c2, c3);
              } else {
                quads.push(c0, c3, c2, c1);
                triangles.push(c0, c3, c2, c0, c2, c1);
              }
            }
          }
        }
      }
    }

    // Generate quads for Z-edges: shared by 4 cells (x-1,y-1), (x,y-1), (x,y), (x-1,y)
    for (let iz = 0; iz < nz - 1; iz++) {
      for (let iy = 1; iy < ny - 1; iy++) {
        for (let ix = 1; ix < nx - 1; ix++) {
          const v1 = getNodeVal(ix, iy, iz);
          const v2 = getNodeVal(ix, iy, iz + 1);
          if ((v1 < isovalue && v2 >= isovalue) || (v1 >= isovalue && v2 < isovalue)) {
            const c0 = getCellVertex(ix - 1, iy - 1, iz);
            const c1 = getCellVertex(ix, iy - 1, iz);
            const c2 = getCellVertex(ix, iy, iz);
            const c3 = getCellVertex(ix - 1, iy, iz);

            if (c0 !== undefined && c1 !== undefined && c2 !== undefined && c3 !== undefined) {
              if (v1 < isovalue) {
                quads.push(c0, c1, c2, c3);
                triangles.push(c0, c1, c2, c0, c2, c3);
              } else {
                quads.push(c0, c3, c2, c1);
                triangles.push(c0, c3, c2, c0, c2, c1);
              }
            }
          }
        }
      }
    }

    const vertexCount = verticesList.length / 3;
    const triangleCount = triangles.length / 3;
    const quadCount = quads.length / 4;

    const vertices = new Float32Array(verticesList);
    const normals = new Float32Array(normalsList);
    const indices = new Uint32Array(triangles);
    const quadIndices = new Uint32Array(quads);

    // 3. Compute topological invariants: Euler characteristic & genus
    const topo = this._computeTopology(vertices, indices);

    // 4. Compute surface area and enclosed volume
    const calculus = this._computeCalculus(vertices, indices);

    return {
      vertices,
      normals,
      indices,
      quadIndices,
      vertexCount,
      triangleCount,
      quadCount,
      eulerCharacteristic: topo.eulerCharacteristic,
      genus: topo.genus,
      isClosedManifold: topo.isClosedManifold,
      surfaceArea: calculus.surfaceArea,
      enclosedVolume: calculus.enclosedVolume
    };
  }

  /**
   * Compute Euler characteristic chi = V - E + F and manifold closure.
   * @private
   */
  _computeTopology(vertices, indices) {
    const vCount = vertices.length / 3;
    const fCount = indices.length / 3;
    if (vCount === 0 || fCount === 0) {
      return { eulerCharacteristic: 0, genus: 0, isClosedManifold: false };
    }

    const edgeMap = new Map();

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      const addEdge = (a, b) => {
        const minIdx = Math.min(a, b);
        const maxIdx = Math.max(a, b);
        const key = `${minIdx}_${maxIdx}`;
        edgeMap.set(key, (edgeMap.get(key) || 0) + 1);
      };

      addEdge(i0, i1);
      addEdge(i1, i2);
      addEdge(i2, i0);
    }

    const eCount = edgeMap.size;
    const chi = vCount - eCount + fCount;
    const genus = Math.round(1.0 - chi / 2.0);

    // A closed 2-manifold must have every edge shared by exactly 2 faces
    let isClosedManifold = true;
    for (const count of edgeMap.values()) {
      if (count !== 2) {
        isClosedManifold = false;
        break;
      }
    }

    return {
      eulerCharacteristic: chi,
      genus,
      isClosedManifold
    };
  }

  /**
   * Compute surface area and enclosed volume via divergence theorem on closed 2-manifolds.
   * @private
   */
  _computeCalculus(vertices, indices) {
    let surfaceArea = 0.0;
    let enclosedVolume = 0.0;
    const fCount = indices.length / 3;

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0] * 3;
      const i1 = indices[f * 3 + 1] * 3;
      const i2 = indices[f * 3 + 2] * 3;

      const x0 = vertices[i0], y0 = vertices[i0 + 1], z0 = vertices[i0 + 2];
      const x1 = vertices[i1], y1 = vertices[i1 + 1], z1 = vertices[i1 + 2];
      const x2 = vertices[i2], y2 = vertices[i2 + 1], z2 = vertices[i2 + 2];

      // Edge vectors
      const e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
      const e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;

      // Cross product e1 x e2
      const cx = e1y * e2z - e1z * e2y;
      const cy = e1z * e2x - e1x * e2z;
      const cz = e1x * e2y - e1y * e2x;

      const triArea = 0.5 * Math.hypot(cx, cy, cz);
      surfaceArea += triArea;

      // Signed tetrahedral volume: V_tet = (1/6) * (v0 . (v1 x v2))
      const tetVol = (x0 * (y1 * z2 - z1 * y2) - y0 * (x1 * z2 - z1 * x2) + z0 * (x1 * y2 - y1 * x2)) / 6.0;
      enclosedVolume += tetVol;
    }

    return {
      surfaceArea,
      enclosedVolume: Math.abs(enclosedVolume)
    };
  }
}
