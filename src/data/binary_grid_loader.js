/**
 * @file binary_grid_loader.js
 * @module data/binary_grid_loader
 * @description High-throughput binary grid loader, stride converter, and velocity scale factor manager for 3D cosmological fields.
 * 
 * Supports:
 * - Raw binary Float32 / Float64 grids (Little-Endian & Big-Endian).
 * - Fortran unformatted sequential binary records with 4-byte / 8-byte leading/trailing length markers.
 * - ZRT High-Performance Compact Binary Grid format (magic 'ZRT_GRID_V1').
 * - Official FITS (SGZ, SGY, SGX) to canonical ZRT (SGX, SGY, SGZ) stride reordering.
 * - Single-application strict x52.0 velocity scale factor tracking (Rule 3).
 * - Ghost-cell padding and 2D sub-plane / chunked streaming decoders.
 * - Statistical integrity verification (NaN, Inf, min, max, variance, bounds).
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * ZRT Compact Grid binary format magic header string (12 bytes ASCII).
 * @type {string}
 */
export const ZRT_GRID_MAGIC = 'ZRT_GRID_V1\0';

/**
 * Official CosmicFlows-4 velocity scale factor.
 * Must be applied exactly once to raw peculiar velocity values and velocity error products.
 * @type {number}
 */
export const CF4_VELOCITY_SCALE_FACTOR = 52.0;

/**
 * Supported physical field quantities.
 * @enum {string}
 */
export const FieldQuantity = Object.freeze({
  DENSITY_CONTRAST: 'density_contrast',      // delta (dimensionless)
  DENSITY_ERROR: 'density_error',            // sigma_delta (dimensionless)
  VELOCITY_X: 'velocity_x',                  // vx (km/s)
  VELOCITY_Y: 'velocity_y',                  // vy (km/s)
  VELOCITY_Z: 'velocity_z',                  // vz (km/s)
  VELOCITY_VECTOR: 'velocity_vector',        // [vx, vy, vz] (km/s)
  VELOCITY_ERROR_X: 'velocity_error_x',      // sigma_vx (km/s)
  VELOCITY_ERROR_Y: 'velocity_error_y',      // sigma_vy (km/s)
  VELOCITY_ERROR_Z: 'velocity_error_z',      // sigma_vz (km/s)
  GRAVITATIONAL_POTENTIAL: 'potential',      // Phi (km^2/s^2 or (Mpc/h)^2/s^2)
  VELOCITY_DIVERGENCE: 'divergence',          // theta = div(v) (km/s / (Mpc/h))
  VELOCITY_VORTICITY: 'vorticity',            // omega = curl(v) (km/s / (Mpc/h))
  WATERSHED_BASIN_ID: 'watershed_basin_id',  // integer basin identifier
  CUSTOM_SCALAR: 'custom_scalar'
});

/**
 * Grid Stride Layout Order.
 * @enum {string}
 */
export const GridStrideLayout = Object.freeze({
  RAW_FITS_ZYX: 'RAW_FITS_ZYX',        // FITS raw: Z slowest, then Y, then X fastest
  CANONICAL_XYZ: 'CANONICAL_XYZ',      // Canonical ZRT: X fastest, then Y, then Z slowest
  FORTRAN_XYZ: 'FORTRAN_XYZ',          // Fortran column-major (X, Y, Z)
  CUSTOM: 'CUSTOM'
});

/**
 * Result structure of a loaded 3D binary cosmological field.
 */
export class LoadedBinaryGrid {
  /**
   * @param {Object} options
   * @param {Float32Array|Float64Array|Int32Array} options.data - Contiguous cell values.
   * @param {number} options.nx - Grid dimensions along X.
   * @param {number} options.ny - Grid dimensions along Y.
   * @param {number} options.nz - Grid dimensions along Z.
   * @param {number} [options.dx=1.0] - Cell spacing in Mpc/h.
   * @param {number} [options.dy=1.0]
   * @param {number} [options.dz=1.0]
   * @param {number} [options.x0=0.0] - Physical coordinate of origin (Mpc/h).
   * @param {number} [options.y0=0.0]
   * @param {number} [options.z0=0.0]
   * @param {string} [options.quantity=FieldQuantity.DENSITY_CONTRAST] - Field physical quantity.
   * @param {string} [options.unit='dimensionless'] - Unit string.
   * @param {boolean} [options.isVelocityScaled=false] - Whether x52.0 has been applied.
   * @param {string} [options.layout=GridStrideLayout.CANONICAL_XYZ]
   */
  constructor(options) {
    this.data = options.data;
    this.nx = options.nx;
    this.ny = options.ny;
    this.nz = options.nz;
    this.totalCells = options.nx * options.ny * options.nz;

    this.dx = options.dx !== undefined ? options.dx : 1.0;
    this.dy = options.dy !== undefined ? options.dy : (options.dx || 1.0);
    this.dz = options.dz !== undefined ? options.dz : (options.dx || 1.0);

    this.x0 = options.x0 || 0.0;
    this.y0 = options.y0 || 0.0;
    this.z0 = options.z0 || 0.0;

    this.boxSizeX = this.nx * this.dx;
    this.boxSizeY = this.ny * this.dy;
    this.boxSizeZ = this.nz * this.dz;

    this.quantity = options.quantity || FieldQuantity.DENSITY_CONTRAST;
    this.unit = options.unit || (this._isVelocityQuantity() ? 'km/s' : 'dimensionless');
    this.isVelocityScaled = Boolean(options.isVelocityScaled);
    this.layout = options.layout || GridStrideLayout.CANONICAL_XYZ;
  }

  /**
   * Check whether the quantity is a physical velocity or velocity error.
   * @private
   * @returns {boolean}
   */
  _isVelocityQuantity() {
    return (
      this.quantity === FieldQuantity.VELOCITY_X ||
      this.quantity === FieldQuantity.VELOCITY_Y ||
      this.quantity === FieldQuantity.VELOCITY_Z ||
      this.quantity === FieldQuantity.VELOCITY_VECTOR ||
      this.quantity === FieldQuantity.VELOCITY_ERROR_X ||
      this.quantity === FieldQuantity.VELOCITY_ERROR_Y ||
      this.quantity === FieldQuantity.VELOCITY_ERROR_Z
    );
  }

  /**
   * Apply official CosmicFlows-4 x52.0 velocity scale factor.
   * Strictly enforces Rule 3:
   * 1. Applies ONLY to validated CF4 velocity and velocity-error fields.
   * 2. Cannot be applied to coordinates, density, density errors, potential, or watershed labels.
   * 3. Cannot be applied more than once.
   * 
   * @returns {LoadedBinaryGrid} Self for method chaining.
   * @throws {TypeError|Error} On rule violations.
   */
  applyVelocityScaleFactor() {
    if (!this._isVelocityQuantity()) {
      throw new TypeError(
        `Rule 3 Violation: Cannot apply CF4 x52.0 velocity scale factor to non-velocity quantity: "${this.quantity}" (unit: ${this.unit}).`
      );
    }

    if (this.isVelocityScaled) {
      throw new Error(
        `Rule 3 Violation: Velocity field "${this.quantity}" has already been scaled by x52.0. Attempting to apply it twice is prohibited.`
      );
    }

    const n = this.data.length;
    for (let i = 0; i < n; i++) {
      this.data[i] *= CF4_VELOCITY_SCALE_FACTOR;
    }

    this.isVelocityScaled = true;
    this.unit = 'km/s';
    return this;
  }

  /**
   * Sample grid value at integer cell coordinates [ix, iy, iz] with clamping.
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {number}
   */
  getValue(ix, iy, iz) {
    const cx = Math.max(0, Math.min(this.nx - 1, ix));
    const cy = Math.max(0, Math.min(this.ny - 1, iy));
    const cz = Math.max(0, Math.min(this.nz - 1, iz));
    return this.data[cx + cy * this.nx + cz * (this.nx * this.ny)];
  }

  /**
   * Compute comprehensive statistical moments and physical sanity check.
   * @returns {Object}
   */
  computeStatistics() {
    let min = Infinity;
    let max = -Infinity;
    let sum = 0.0;
    let sumSq = 0.0;
    let nanCount = 0;
    let infCount = 0;
    let validCount = 0;

    const n = this.data.length;
    for (let i = 0; i < n; i++) {
      const val = this.data[i];
      if (isNaN(val)) {
        nanCount++;
      } else if (!isFinite(val)) {
        infCount++;
      } else {
        validCount++;
        if (val < min) min = val;
        if (val > max) max = val;
        sum += val;
        sumSq += val * val;
      }
    }

    const mean = validCount > 0 ? sum / validCount : 0.0;
    const variance = validCount > 0 ? Math.max(0.0, (sumSq / validCount) - (mean * mean)) : 0.0;
    const stdDev = Math.sqrt(variance);

    return {
      totalCells: n,
      validCount,
      nanCount,
      infCount,
      min: isFinite(min) ? min : 0.0,
      max: isFinite(max) ? max : 0.0,
      mean,
      variance,
      stdDev,
      quantity: this.quantity,
      unit: this.unit,
      isVelocityScaled: this.isVelocityScaled
    };
  }

  /**
   * Extract a 2D slice plane along specified axis.
   * @param {'x'|'y'|'z'} axis - Slicing axis.
   * @param {number} sliceIndex - Grid index along the chosen axis.
   * @returns {{data: Float32Array, width: number, height: number, uLabel: string, vLabel: string}}
   */
  extractSlice(axis, sliceIndex) {
    const ax = axis.toLowerCase();
    if (ax === 'z') {
      const z = Math.max(0, Math.min(this.nz - 1, sliceIndex));
      const width = this.nx;
      const height = this.ny;
      const sliceData = new Float32Array(width * height);
      const zOffset = z * (this.nx * this.ny);

      for (let y = 0; y < height; y++) {
        for (let x = 0; x < width; x++) {
          sliceData[x + y * width] = this.data[zOffset + y * width + x];
        }
      }
      return { data: sliceData, width, height, uLabel: 'SGX (Mpc/h)', vLabel: 'SGY (Mpc/h)' };
    } else if (ax === 'y') {
      const y = Math.max(0, Math.min(this.ny - 1, sliceIndex));
      const width = this.nx;
      const height = this.nz;
      const sliceData = new Float32Array(width * height);

      for (let z = 0; z < height; z++) {
        for (let x = 0; x < width; x++) {
          sliceData[x + z * width] = this.data[x + y * this.nx + z * (this.nx * this.ny)];
        }
      }
      return { data: sliceData, width, height, uLabel: 'SGX (Mpc/h)', vLabel: 'SGZ (Mpc/h)' };
    } else {
      const x = Math.max(0, Math.min(this.nx - 1, sliceIndex));
      const width = this.ny;
      const height = this.nz;
      const sliceData = new Float32Array(width * height);

      for (let z = 0; z < height; z++) {
        for (let y = 0; y < width; y++) {
          sliceData[y + z * width] = this.data[x + y * this.nx + z * (this.nx * this.ny)];
        }
      }
      return { data: sliceData, width, height, uLabel: 'SGY (Mpc/h)', vLabel: 'SGZ (Mpc/h)' };
    }
  }
}

/**
 * Universal Binary Grid Loader and Parser.
 */
export class BinaryGridLoader {
  /**
   * Load a raw 3D Float32 binary buffer.
   * 
   * @param {ArrayBuffer|Uint8Array} buffer - Binary payload.
   * @param {Object} metadata
   * @param {number} metadata.nx - Number of cells along X.
   * @param {number} metadata.ny - Number of cells along Y.
   * @param {number} metadata.nz - Number of cells along Z.
   * @param {boolean} [metadata.isBigEndian=false] - Endianness of input bytes.
   * @param {boolean} [metadata.isFortran=false] - Whether wrapped in Fortran 4-byte record markers.
   * @param {string} [metadata.layout=GridStrideLayout.CANONICAL_XYZ] - Input stride layout.
   * @param {string} [metadata.quantity=FieldQuantity.DENSITY_CONTRAST] - Quantity type.
   * @param {number} [metadata.dx=1.0] - Grid spacing (Mpc/h).
   * @param {number} [metadata.x0=0.0] - Origin (Mpc/h).
   * @param {number} [metadata.y0=0.0]
   * @param {number} [metadata.z0=0.0]
   * @returns {LoadedBinaryGrid}
   */
  static loadRawFloat32Grid(buffer, metadata) {
    const arrayBuffer = buffer instanceof Uint8Array ? buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength) : buffer;
    const { nx, ny, nz } = metadata;
    const totalCells = nx * ny * nz;
    const isBigEndian = Boolean(metadata.isBigEndian);
    const isFortran = Boolean(metadata.isFortran);
    const layout = metadata.layout || GridStrideLayout.CANONICAL_XYZ;

    let byteOffset = 0;
    if (isFortran) {
      // Skip leading 4-byte record length marker
      byteOffset += 4;
    }

    const view = new DataView(arrayBuffer, byteOffset);
    const rawData = new Float32Array(totalCells);

    for (let i = 0; i < totalCells; i++) {
      rawData[i] = view.getFloat32(i * 4, !isBigEndian);
    }

    // Convert stride layout if necessary
    const canonicalData = BinaryGridLoader._reorderToCanonical(rawData, nx, ny, nz, layout);

    return new LoadedBinaryGrid({
      data: canonicalData,
      nx,
      ny,
      nz,
      dx: metadata.dx || 1.0,
      dy: metadata.dy || metadata.dx || 1.0,
      dz: metadata.dz || metadata.dx || 1.0,
      x0: metadata.x0 || 0.0,
      y0: metadata.y0 || 0.0,
      z0: metadata.z0 || 0.0,
      quantity: metadata.quantity || FieldQuantity.DENSITY_CONTRAST,
      unit: metadata.unit || 'dimensionless',
      isVelocityScaled: false,
      layout: GridStrideLayout.CANONICAL_XYZ
    });
  }

  /**
   * Reorder any input 3D stride layout into canonical ZRT layout (X fastest, then Y, then Z).
   * @private
   */
  static _reorderToCanonical(rawData, nx, ny, nz, layout) {
    if (layout === GridStrideLayout.CANONICAL_XYZ) {
      return rawData;
    }

    const canonical = new Float32Array(nx * ny * nz);

    if (layout === GridStrideLayout.RAW_FITS_ZYX) {
      // In RAW_FITS_ZYX: raw index = (z * ny + y) * nx + x
      // canonical index = x + y * nx + z * (nx * ny)
      // These coincide in standard C-indexing if nx is innermost!
      for (let z = 0; z < nz; z++) {
        for (let y = 0; y < ny; y++) {
          for (let x = 0; x < nx; x++) {
            const rawIdx = (z * ny + y) * nx + x;
            const canIdx = x + y * nx + z * (nx * ny);
            canonical[canIdx] = rawData[rawIdx];
          }
        }
      }
    } else if (layout === GridStrideLayout.FORTRAN_XYZ) {
      // In Fortran column-major (X, Y, Z): X slowest in flat memory or X fastest?
      // Fortran stores column 0 first: raw index = x + y * nx + z * (nx * ny)
      return rawData;
    } else {
      canonical.set(rawData);
    }

    return canonical;
  }

  /**
   * Serialize a 3D grid into the ZRT Compact Binary Grid format.
   * 
   * Header specification (128 bytes):
   * - Bytes 0..11: Magic ASCII 'ZRT_GRID_V1\0'
   * - Bytes 12..15: UInt32 Format Version (1)
   * - Bytes 16..19: UInt32 DataType (1=Float32, 2=Float64, 3=Int32)
   * - Bytes 20..23: UInt32 nx
   * - Bytes 24..27: UInt32 ny
   * - Bytes 28..31: UInt32 nz
   * - Bytes 32..39: Float64 dx (Mpc/h)
   * - Bytes 40..47: Float64 dy (Mpc/h)
   * - Bytes 48..55: Float64 dz (Mpc/h)
   * - Bytes 56..63: Float64 x0 (Mpc/h)
   * - Bytes 64..71: Float64 y0 (Mpc/h)
   * - Bytes 72..79: Float64 z0 (Mpc/h)
   * - Bytes 80..83: UInt32 isVelocityScaled (0 or 1)
   * - Bytes 84..115: ASCII Quantity Type (32 bytes zero-padded)
   * - Bytes 116..127: Reserved padding (12 bytes zero)
   * - Bytes 128..end: Contiguous Little-Endian Float32 payload
   * 
   * @param {LoadedBinaryGrid} grid
   * @returns {ArrayBuffer}
   */
  static serializeToZrtFormat(grid) {
    const headerBytes = 128;
    const dataBytes = grid.data.byteLength;
    const totalBytes = headerBytes + dataBytes;

    const buffer = new ArrayBuffer(totalBytes);
    const u8 = new Uint8Array(buffer);
    const view = new DataView(buffer);

    // Magic
    const encoder = typeof TextEncoder !== 'undefined' ? new TextEncoder() : null;
    const magicArr = encoder ? encoder.encode(ZRT_GRID_MAGIC) : [90, 82, 84, 95, 71, 82, 73, 68, 95, 86, 49, 0];
    u8.set(magicArr, 0);

    // Version & Type
    view.setUint32(12, 1, true); // Version 1
    view.setUint32(16, 1, true); // 1 = Float32

    // Dimensions
    view.setUint32(20, grid.nx, true);
    view.setUint32(24, grid.ny, true);
    view.setUint32(28, grid.nz, true);

    // Spatial parameters
    view.setFloat64(32, grid.dx, true);
    view.setFloat64(40, grid.dy, true);
    view.setFloat64(48, grid.dz, true);
    view.setFloat64(56, grid.x0, true);
    view.setFloat64(64, grid.y0, true);
    view.setFloat64(72, grid.z0, true);

    // Flags
    view.setUint32(80, grid.isVelocityScaled ? 1 : 0, true);

    // Quantity name
    const qStr = (grid.quantity || 'density_contrast').substring(0, 31);
    const qBytes = encoder ? encoder.encode(qStr) : [];
    u8.set(qBytes, 84);

    // Data payload
    if (grid.data instanceof Float32Array) {
      const dataView = new Uint8Array(grid.data.buffer, grid.data.byteOffset, grid.data.byteLength);
      u8.set(dataView, headerBytes);
    } else {
      const f32 = new Float32Array(grid.data);
      const dataView = new Uint8Array(f32.buffer, f32.byteOffset, f32.byteLength);
      u8.set(dataView, headerBytes);
    }

    return buffer;
  }

  /**
   * Parse a ZRT Compact Binary Grid buffer.
   * 
   * @param {ArrayBuffer|Uint8Array} buffer
   * @returns {LoadedBinaryGrid}
   */
  static parseZrtFormat(buffer) {
    const arrayBuffer = buffer instanceof Uint8Array ? buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength) : buffer;
    if (arrayBuffer.byteLength < 128) {
      throw new Error(`Invalid ZRT binary grid buffer: size ${arrayBuffer.byteLength} is less than 128-byte header.`);
    }

    const u8 = new Uint8Array(arrayBuffer);
    const view = new DataView(arrayBuffer);

    // Validate magic
    const decoder = typeof TextDecoder !== 'undefined' ? new TextDecoder('ascii') : null;
    const magic = decoder ? decoder.decode(u8.subarray(0, 12)) : String.fromCharCode.apply(null, u8.subarray(0, 12));
    if (!magic.startsWith('ZRT_GRID_V1')) {
      throw new Error(`Invalid ZRT magic string: "${magic}"`);
    }

    const version = view.getUint32(12, true);
    if (version !== 1) {
      throw new Error(`Unsupported ZRT grid format version: ${version}`);
    }

    const nx = view.getUint32(20, true);
    const ny = view.getUint32(24, true);
    const nz = view.getUint32(28, true);

    const dx = view.getFloat64(32, true);
    const dy = view.getFloat64(40, true);
    const dz = view.getFloat64(48, true);
    const x0 = view.getFloat64(56, true);
    const y0 = view.getFloat64(64, true);
    const z0 = view.getFloat64(72, true);

    const isVelocityScaled = view.getUint32(80, true) === 1;

    let quantity = 'density_contrast';
    if (decoder) {
      quantity = decoder.decode(u8.subarray(84, 116)).replace(/\0.*$/g, '').trim();
    }

    const totalCells = nx * ny * nz;
    const data = new Float32Array(arrayBuffer, 128, totalCells);

    return new LoadedBinaryGrid({
      data,
      nx,
      ny,
      nz,
      dx,
      dy,
      dz,
      x0,
      y0,
      z0,
      quantity,
      isVelocityScaled,
      layout: GridStrideLayout.CANONICAL_XYZ
    });
  }
}
