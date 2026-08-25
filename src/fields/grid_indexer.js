/**
 * @file grid_indexer.js
 * @description Exact FITS (SGZ, SGY, SGX) <-> Canonical (SGX, SGY, SGZ) stride mapping,
 * voxel indexing, continuous physical coordinate transformations, and boundary clamping
 * for cosmological 3D scalar and vector grids.
 * 
 * Cosmicflows / FITS Coordinate Conventions:
 * - FITS arrays in astronomical formats typically order axes as NAXIS1 (SGX), NAXIS2 (SGY), NAXIS3 (SGZ).
 * - In row-major (C-order) 3D memory layouts, flat index = ix + nx * (iy + ny * iz).
 * - Physical coordinates: Supergalactic coordinates (SGX, SGY, SGZ) in h^-1 Mpc or km/s.
 * 
 * @module fields/grid_indexer
 */

/**
 * Supported boundary handling modes.
 * @readonly
 * @enum {string}
 */
export const BoundaryMode = Object.freeze({
  CLAMP: 'clamp',
  PERIODIC: 'periodic',
  REFLECT: 'reflect',
  ZERO: 'zero'
});

/**
 * Supported axis storage ordering conventions.
 * @readonly
 * @enum {string}
 */
export const StrideOrder = Object.freeze({
  CANONICAL_XYZ: 'canonical_xyz', // X fastest (stride 1), Y middle (stride nx), Z slowest (stride nx*ny)
  FITS_ZYX: 'fits_zyx',           // Z slowest (axis 3), Y middle (axis 2), X fastest (axis 1) - standard FITS flat buffer
  FORTRAN_XYZ: 'fortran_xyz',     // Same flat order as FITS_ZYX (1-based column major equivalent)
  CUSTOM: 'custom'
});

/**
 * GridIndexer provides high-performance 3D voxel indexing, coordinate projections,
 * and boundary condition management.
 */
export class GridIndexer {
  /**
   * Constructs a 3D GridIndexer.
   * 
   * @param {Object} options Configuration parameters.
   * @param {number} options.nx Grid dimensions along SGX (integer >= 2).
   * @param {number} options.ny Grid dimensions along SGY (integer >= 2).
   * @param {number} options.nz Grid dimensions along SGZ (integer >= 2).
   * @param {Array<number>} [options.origin=[-100.0, -100.0, -100.0]] Physical coordinates of grid minimum [xMin, yMin, zMin].
   * @param {Array<number>} [options.boxSize=[200.0, 200.0, 200.0]] Physical domain dimensions [Lx, Ly, Lz].
   * @param {Array<number>} [options.spacing] Optional explicit grid spacing [dx, dy, dz].
   * @param {string} [options.strideOrder='canonical_xyz'] Memory layout stride order.
   * @param {string} [options.boundaryMode='clamp'] Boundary handling strategy ('clamp', 'periodic', 'reflect', 'zero').
   * @param {boolean} [options.isCellCentered=false] True if voxels represent cell centers, false if vertex/node grid.
   */
  constructor(options) {
    if (!options) {
      throw new Error('GridIndexer: options object must be provided.');
    }

    const {
      nx,
      ny,
      nz,
      origin = [-100.0, -100.0, -100.0],
      boxSize = [200.0, 200.0, 200.0],
      spacing = null,
      strideOrder = StrideOrder.CANONICAL_XYZ,
      boundaryMode = BoundaryMode.CLAMP,
      isCellCentered = false
    } = options;

    if (!Number.isInteger(nx) || nx < 2 || !Number.isInteger(ny) || ny < 2 || !Number.isInteger(nz) || nz < 2) {
      throw new RangeError(`GridIndexer: Grid dimensions (nx, ny, nz) must be integers >= 2. Received (${nx}, ${ny}, ${nz})`);
    }

    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.totalCells = nx * ny * nz;

    if (!Array.isArray(origin) || origin.length !== 3 || origin.some(v => typeof v !== 'number' || !Number.isFinite(v))) {
      throw new TypeError('GridIndexer: origin must be a 3-element array of finite numbers [x0, y0, z0].');
    }
    this.origin = Float64Array.from(origin);

    if (spacing) {
      if (!Array.isArray(spacing) || spacing.length !== 3 || spacing.some(v => typeof v !== 'number' || v <= 0 || !Number.isFinite(v))) {
        throw new TypeError('GridIndexer: spacing must be a 3-element array of positive numbers [dx, dy, dz].');
      }
      this.dx = spacing[0];
      this.dy = spacing[1];
      this.dz = spacing[2];
      this.boxSize = Float64Array.of(
        isCellCentered ? this.dx * nx : this.dx * (nx - 1),
        isCellCentered ? this.dy * ny : this.dy * (ny - 1),
        isCellCentered ? this.dz * nz : this.dz * (nz - 1)
      );
    } else {
      if (!Array.isArray(boxSize) || boxSize.length !== 3 || boxSize.some(v => typeof v !== 'number' || v <= 0 || !Number.isFinite(v))) {
        throw new TypeError('GridIndexer: boxSize must be a 3-element array of positive numbers [Lx, Ly, Lz].');
      }
      this.boxSize = Float64Array.from(boxSize);
      this.dx = isCellCentered ? this.boxSize[0] / nx : this.boxSize[0] / (nx - 1);
      this.dy = isCellCentered ? this.boxSize[1] / ny : this.boxSize[1] / (ny - 1);
      this.dz = isCellCentered ? this.boxSize[2] / nz : this.boxSize[2] / (nz - 1);
    }

    this.invDx = 1.0 / this.dx;
    this.invDy = 1.0 / this.dy;
    this.invDz = 1.0 / this.dz;

    this.strideOrder = strideOrder;
    this.boundaryMode = boundaryMode;
    this.isCellCentered = Boolean(isCellCentered);

    // Compute strides
    this.strideX = 1;
    this.strideY = nx;
    this.strideZ = nx * ny;

    // Physical domain bounds
    this.xMin = this.origin[0];
    this.yMin = this.origin[1];
    this.zMin = this.origin[2];
    this.xMax = this.xMin + this.boxSize[0];
    this.yMax = this.yMin + this.boxSize[1];
    this.zMax = this.zMin + this.boxSize[2];
  }

  /**
   * Computes the 1D flat buffer index from 3D grid cell indices (ix, iy, iz).
   * 
   * @param {number} ix Integer index along SGX [0, nx-1].
   * @param {number} iy Integer index along SGY [0, ny-1].
   * @param {number} iz Integer index along SGZ [0, nz-1].
   * @returns {number} 1D flat index.
   */
  getLinearIndex(ix, iy, iz) {
    return ix * this.strideX + iy * this.strideY + iz * this.strideZ;
  }

  /**
   * Alias for getLinearIndex.
   * 
   * @param {number} ix Integer index along SGX [0, nx-1].
   * @param {number} iy Integer index along SGY [0, ny-1].
   * @param {number} iz Integer index along SGZ [0, nz-1].
   * @returns {number} 1D flat index.
   */
  index(ix, iy, iz) {
    return this.getLinearIndex(ix, iy, iz);
  }

  /**
   * Converts a 1D flat buffer index back into 3D integer indices [ix, iy, iz].
   * 
   * @param {number} linearIndex 1D buffer index [0, totalCells - 1].
   * @returns {[number, number, number]} [ix, iy, iz] grid indices.
   */
  get3DIndices(linearIndex) {
    if (linearIndex < 0 || linearIndex >= this.totalCells || !Number.isInteger(linearIndex)) {
      throw new RangeError(`GridIndexer: linearIndex ${linearIndex} out of range [0, ${this.totalCells - 1}].`);
    }
    const iz = Math.floor(linearIndex / (this.nx * this.ny));
    const rem = linearIndex - iz * (this.nx * this.ny);
    const iy = Math.floor(rem / this.nx);
    const ix = rem - iy * this.nx;
    return [ix, iy, iz];
  }

  /**
   * Maps FITS primary HDU coordinate indices (SGZ_fits, SGY_fits, SGX_fits)
   * where SGX is NAXIS1, SGY is NAXIS2, SGZ is NAXIS3 into canonical flat index.
   * 
   * @param {number} ix NAXIS1 index (SGX).
   * @param {number} iy NAXIS2 index (SGY).
   * @param {number} iz NAXIS3 index (SGZ).
   * @returns {number} Canonical linear index.
   */
  fitsToCanonicalIndex(ix, iy, iz) {
    return this.getLinearIndex(ix, iy, iz);
  }

  /**
   * Clamps or wraps integer 3D indices according to the configured BoundaryMode.
   * 
   * @param {number} ix Integer index along SGX.
   * @param {number} iy Integer index along SGY.
   * @param {number} iz Integer index along SGZ.
   * @param {string} [modeOverride] Optional boundary mode override.
   * @returns {[number, number, number]|null} Clamped/wrapped [ix, iy, iz], or null if ZERO mode is out of bounds.
   */
  handleBoundaryIndex(ix, iy, iz, modeOverride = null) {
    const mode = modeOverride || this.boundaryMode;

    switch (mode) {
      case BoundaryMode.CLAMP: {
        const cx = Math.max(0, Math.min(this.nx - 1, ix));
        const cy = Math.max(0, Math.min(this.ny - 1, iy));
        const cz = Math.max(0, Math.min(this.nz - 1, iz));
        return [cx, cy, cz];
      }

      case BoundaryMode.PERIODIC: {
        const px = ((ix % this.nx) + this.nx) % this.nx;
        const py = ((iy % this.ny) + this.ny) % this.ny;
        const pz = ((iz % this.nz) + this.nz) % this.nz;
        return [px, py, pz];
      }

      case BoundaryMode.REFLECT: {
        const reflect1D = (i, n) => {
          if (n <= 1) return 0;
          let p = ((i % (2 * n - 2)) + (2 * n - 2)) % (2 * n - 2);
          if (p >= n) {
            p = 2 * n - 2 - p;
          }
          return p;
        };
        return [reflect1D(ix, this.nx), reflect1D(iy, this.ny), reflect1D(iz, this.nz)];
      }

      case BoundaryMode.ZERO: {
        if (ix < 0 || ix >= this.nx || iy < 0 || iy >= this.ny || iz < 0 || iz >= this.nz) {
          return null;
        }
        return [ix, iy, iz];
      }

      default:
        throw new Error(`GridIndexer: Unknown boundary mode "${mode}".`);
    }
  }

  /**
   * Checks whether a continuous physical point (x, y, z) lies within the domain bounding box.
   * 
   * @param {number} x SGX coordinate.
   * @param {number} y SGY coordinate.
   * @param {number} z SGZ coordinate.
   * @returns {boolean} True if inside [xMin, xMax] x [yMin, yMax] x [zMin, zMax].
   */
  isInsideDomain(x, y, z) {
    return (
      x >= this.xMin && x <= this.xMax &&
      y >= this.yMin && y <= this.yMax &&
      z >= this.zMin && z <= this.zMax
    );
  }

  /**
   * Clamps continuous physical coordinates (x, y, z) to the domain bounding box.
   * 
   * @param {number} x SGX coordinate.
   * @param {number} y SGY coordinate.
   * @param {number} z SGZ coordinate.
   * @returns {[number, number, number]} Clamped [x, y, z].
   */
  clampCoordinates(x, y, z) {
    const cx = Math.max(this.xMin, Math.min(this.xMax, x));
    const cy = Math.max(this.yMin, Math.min(this.yMax, y));
    const cz = Math.max(this.zMin, Math.min(this.zMax, z));
    return [cx, cy, cz];
  }

  /**
   * Converts continuous physical coordinates (x, y, z) into continuous fractional grid indices (gx, gy, gz).
   * 
   * @param {number} x SGX coordinate.
   * @param {number} y SGY coordinate.
   * @param {number} z SGZ coordinate.
   * @returns {[number, number, number]} Continuous fractional grid coordinates.
   */
  coordToGridIndex(x, y, z) {
    let gx = (x - this.xMin) * this.invDx;
    let gy = (y - this.yMin) * this.invDy;
    let gz = (z - this.zMin) * this.invDz;

    if (this.isCellCentered) {
      gx -= 0.5;
      gy -= 0.5;
      gz -= 0.5;
    }

    return [gx, gy, gz];
  }

  /**
   * Converts continuous fractional grid coordinates (gx, gy, gz) into continuous physical coordinates (x, y, z).
   * 
   * @param {number} gx Grid index along SGX.
   * @param {number} gy Grid index along SGY.
   * @param {number} gz Grid index along SGZ.
   * @returns {[number, number, number]} Physical coordinates [x, y, z].
   */
  gridIndexToCoord(gx, gy, gz) {
    let x, y, z;
    if (this.isCellCentered) {
      x = this.xMin + (gx + 0.5) * this.dx;
      y = this.yMin + (gy + 0.5) * this.dy;
      z = this.zMin + (gz + 0.5) * this.dz;
    } else {
      x = this.xMin + gx * this.dx;
      y = this.yMin + gy * this.dy;
      z = this.zMin + gz * this.dz;
    }
    return [x, y, z];
  }

  /**
   * Returns the physical coordinates of the node or cell center at integer index (ix, iy, iz).
   * 
   * @param {number} ix Integer index along SGX.
   * @param {number} iy Integer index along SGY.
   * @param {number} iz Integer index along SGZ.
   * @returns {[number, number, number]} Physical coordinates [x, y, z].
   */
  getNodeCoord(ix, iy, iz) {
    return this.gridIndexToCoord(ix, iy, iz);
  }

  /**
   * Returns the 8 corner coordinates bounding the cell containing continuous point (x, y, z).
   * 
   * @param {number} x SGX coordinate.
   * @param {number} y SGY coordinate.
   * @param {number} z SGZ coordinate.
   * @returns {Array<[number, number, number]>} 8 physical corner coordinates.
   */
  getCellCorners(x, y, z) {
    const [gx, gy, gz] = this.coordToGridIndex(x, y, z);
    const i0 = Math.floor(gx);
    const j0 = Math.floor(gy);
    const k0 = Math.floor(gz);

    const corners = [];
    for (let k = 0; k <= 1; k++) {
      for (let j = 0; j <= 1; j++) {
        for (let i = 0; i <= 1; i++) {
          corners.push(this.getNodeCoord(i0 + i, j0 + j, k0 + k));
        }
      }
    }
    return corners;
  }

  /**
   * Retrieves indices of the 6 face-sharing nearest neighbors for a cell (ix, iy, iz).
   * 
   * @param {number} ix Integer index along SGX.
   * @param {number} iy Integer index along SGY.
   * @param {number} iz Integer index along SGZ.
   * @returns {Array<[number, number, number]|null>} Array of 6 neighbor indices (+x, -x, +y, -y, +z, -z).
   */
  get6Neighbors(ix, iy, iz) {
    const offsets = [
      [1, 0, 0], [-1, 0, 0],
      [0, 1, 0], [0, -1, 0],
      [0, 0, 1], [0, 0, -1]
    ];
    return offsets.map(([ox, oy, oz]) => this.handleBoundaryIndex(ix + ox, iy + oy, iz + oz));
  }
}
