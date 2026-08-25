/**
 * @file watershed_classifier.js
 * @description Official Dupuy & Courtois (2023) Table A.1 Basin Taxonomy and 128³ FITS/Binary Segmentation Loader.
 * Provides complete astronomical basin metadata, coordinate transformations (voxel <-> Supergalactic Cartesian),
 * FITS/binary voxel buffer decoders, and geometric moments for cosmic velocity watershed basins.
 *
 * References:
 * - Dupuy, A., & Courtois, H. M. (2023). "Cosmicflows-4: The Watershed Basins of Attraction". MNRAS / A&A. Table A.1.
 * - Tully, R. B., Courtois, H., Hoffman, Y., & Pomarède, D. (2014). "The Laniakea supercluster of galaxies". Nature, 513(7516), 71-73.
 */

/**
 * @typedef {Object} BasinDefinition
 * @property {number} id - Unique numeric basin identifier (0 = unsegmented)
 * @property {string} name - Primary astronomical basin name
 * @property {string} abbreviation - Short identification code
 * @property {string} description - Astrometric description and associated superclusters
 * @property {number[]} attractorSGKms - Canonical attractor location [SGX, SGY, SGZ] in km/s
 * @property {number[]} attractorSGMpc - Canonical attractor location [SGX, SGY, SGZ] in Mpc/h (h=0.75 / h=1.0)
 * @property {number} volumeCF4Ind1e6 - Published volume in individual CF4 in 10^6 (Mpc/h)^3
 * @property {number} volumeCF4Grp1e6 - Published volume in grouped CF4 in 10^6 (Mpc/h)^3
 * @property {number} charVelocityKms - Characteristic basin flow velocity in km/s
 * @property {number} distToShapleyMpc - Characteristic separation from Shapley Core in Mpc/h
 * @property {string} hexColor - Publication color palette hexadecimal code
 * @property {number[]} rgbColor - Normalized RGB components [r, g, b] in range [0, 1]
 */

/**
 * Dupuy & Courtois (2023) Table A.1 Official Basin Taxonomy (Basins 0 to 8 + sub-basins).
 * @type {Readonly<Record<number, BasinDefinition>>}
 */
export const BASIN_TAXONOMY = Object.freeze({
  0: Object.freeze({
    id: 0,
    name: 'Unsegmented/Background',
    abbreviation: 'UNSEG',
    description: 'Unsegmented background domain, void regions, or trajectories escaping the computational box',
    attractorSGKms: [0, 0, 0],
    attractorSGMpc: [0, 0, 0],
    volumeCF4Ind1e6: 483.38,
    volumeCF4Grp1e6: 512.40,
    charVelocityKms: 326.6,
    distToShapleyMpc: 284.8,
    hexColor: '#4b5563',
    rgbColor: [0.294, 0.333, 0.388]
  }),
  1: Object.freeze({
    id: 1,
    name: 'Laniakea',
    abbreviation: 'LAN',
    description: 'Home supercluster basin containing the Milky Way, Virgo, Hydra-Centaurus, and Great Attractor (Norma)',
    attractorSGKms: [-3450, 1720, -510],
    attractorSGMpc: [-46.0, 22.9, -6.8],
    volumeCF4Ind1e6: 2.07,
    volumeCF4Grp1e6: 1.95,
    charVelocityKms: 444.2,
    distToShapleyMpc: 179.9,
    hexColor: '#3b82f6',
    rgbColor: [0.231, 0.510, 0.965]
  }),
  2: Object.freeze({
    id: 2,
    name: 'Apus',
    abbreviation: 'APU',
    description: 'Southern hemispheric basin encompassing the South Pole Wall and Pavo-Indus concentration',
    attractorSGKms: [-4820, -7100, -2340],
    attractorSGMpc: [-64.3, -94.7, -31.2],
    volumeCF4Ind1e6: 9.51,
    volumeCF4Grp1e6: 8.84,
    charVelocityKms: 643.8,
    distToShapleyMpc: 205.1,
    hexColor: '#10b981',
    rgbColor: [0.063, 0.725, 0.506]
  }),
  3: Object.freeze({
    id: 3,
    name: 'Hercules',
    abbreviation: 'HER',
    description: 'Hercules Supercluster basin extending towards the Great Wall (CfA2) and Abell 2199/2147 complexes',
    attractorSGKms: [2980, 9450, 4810],
    attractorSGMpc: [39.7, 126.0, 64.1],
    volumeCF4Ind1e6: 3.26,
    volumeCF4Grp1e6: 3.82,
    charVelocityKms: 402.9,
    distToShapleyMpc: 170.2,
    hexColor: '#f59e0b',
    rgbColor: [0.961, 0.620, 0.043]
  }),
  4: Object.freeze({
    id: 4,
    name: 'Lepus',
    abbreviation: 'LEP',
    description: 'Lepus/Columba region basin situated between the Local Void and southern gravitational sinks',
    attractorSGKms: [-9800, -3200, 1450],
    attractorSGMpc: [-130.7, -42.7, 19.3],
    volumeCF4Ind1e6: 8.02,
    volumeCF4Grp1e6: 7.45,
    charVelocityKms: 665.8,
    distToShapleyMpc: 203.0,
    hexColor: '#ec4899',
    rgbColor: [0.925, 0.282, 0.600]
  }),
  5: Object.freeze({
    id: 5,
    name: 'Perseus-Pisces',
    abbreviation: 'PSP',
    description: 'Prominent Perseus-Pisces Supercluster chain and Pisces-Cetus supercluster complex filament',
    attractorSGKms: [5400, -2100, -1800],
    attractorSGMpc: [72.0, -28.0, -24.0],
    volumeCF4Ind1e6: 4.84,
    volumeCF4Grp1e6: 5.12,
    charVelocityKms: 574.0,
    distToShapleyMpc: 260.5,
    hexColor: '#8b5cf6',
    rgbColor: [0.545, 0.361, 0.965]
  }),
  6: Object.freeze({
    id: 6,
    name: 'Shapley',
    abbreviation: 'SHA',
    description: 'The dominant cosmic attractor basin centered on the Shapley Supercluster Concentration (A3558/A3528)',
    attractorSGKms: [-12800, 8900, -1600],
    attractorSGMpc: [-170.7, 118.7, -21.3],
    volumeCF4Ind1e6: 7.82,
    volumeCF4Grp1e6: 26.95,
    charVelocityKms: 485.4,
    distToShapleyMpc: 26.5,
    hexColor: '#ef4444',
    rgbColor: [0.937, 0.267, 0.267]
  }),
  7: Object.freeze({
    id: 7,
    name: 'SDSS-1a',
    abbreviation: 'SD1',
    description: 'Northern SDSS Great Wall basin component alpha (high-redshift Northern Galactic Cap)',
    attractorSGKms: [14200, 18500, 9200],
    attractorSGMpc: [189.3, 246.7, 122.7],
    volumeCF4Ind1e6: 183.86,
    volumeCF4Grp1e6: 165.20,
    charVelocityKms: 396.8,
    distToShapleyMpc: 407.6,
    hexColor: '#06b6d4',
    rgbColor: [0.024, 0.714, 0.831]
  }),
  8: Object.freeze({
    id: 8,
    name: 'SDSS-2a',
    abbreviation: 'SD2',
    description: 'Southern SDSS Great Wall basin component alpha (high-redshift Southern Galactic Cap)',
    attractorSGKms: [-16400, 15100, -11200],
    attractorSGMpc: [-218.7, 201.3, -149.3],
    volumeCF4Ind1e6: 168.79,
    volumeCF4Grp1e6: 154.10,
    charVelocityKms: 472.5,
    distToShapleyMpc: 225.0,
    hexColor: '#14b8a6',
    rgbColor: [0.078, 0.722, 0.651]
  })
});

/**
 * Extended sub-basins mapping for detailed 128³ cosmological research.
 * @type {Readonly<Record<number, BasinDefinition>>}
 */
export const EXTENDED_BASIN_TAXONOMY = Object.freeze({
  ...BASIN_TAXONOMY,
  9: Object.freeze({
    id: 9,
    name: 'SDSS-2b',
    abbreviation: 'SD2b',
    description: 'SDSS southern peripheral basin branch beta',
    attractorSGKms: [-22000, 11000, -15000],
    attractorSGMpc: [-293.3, 146.7, -200.0],
    volumeCF4Ind1e6: 128.45,
    volumeCF4Grp1e6: 110.30,
    charVelocityKms: 341.1,
    distToShapleyMpc: 451.4,
    hexColor: '#6366f1',
    rgbColor: [0.388, 0.400, 0.945]
  })
});

/**
 * Grid Metadata Interface
 * @typedef {Object} GridMetadata
 * @property {number} nx - Grid dimension in X (default: 128)
 * @property {number} ny - Grid dimension in Y (default: 128)
 * @property {number} nz - Grid dimension in Z (default: 128)
 * @property {number} boxSizeMpc - Total box width in Mpc/h (e.g. 1000.0)
 * @property {number} boxSizeKms - Total box width in km/s (e.g. 74600.0)
 * @property {number} halfExtentMpc - Grid half-extent in Mpc/h (e.g. 500.0)
 * @property {number} halfExtentKms - Grid half-extent in km/s (e.g. 37300.0)
 * @property {number} dxMpc - Voxel cell resolution in Mpc/h (e.g. 500*2/128 = 7.8125)
 * @property {number} dxKms - Voxel cell resolution in km/s
 * @property {number} hubbleH0 - Hubble parameter H0 (km/s/Mpc) (default: 74.6)
 * @property {string} coordinateSystem - Name of coordinate system (Supergalactic Cartesian)
 */

/**
 * Creates default grid metadata for standard 128³ or 64³ cosmological boxes.
 * @param {Partial<GridMetadata>} [overrides={}]
 * @returns {GridMetadata}
 */
export function createGridMetadata(overrides = {}) {
  const nx = overrides.nx ?? 128;
  const ny = overrides.ny ?? 128;
  const nz = overrides.nz ?? 128;
  const halfExtentMpc = overrides.halfExtentMpc ?? 500.0;
  const hubbleH0 = overrides.hubbleH0 ?? 74.6;
  const halfExtentKms = overrides.halfExtentKms ?? (halfExtentMpc * hubbleH0);
  const boxSizeMpc = overrides.boxSizeMpc ?? (2.0 * halfExtentMpc);
  const boxSizeKms = overrides.boxSizeKms ?? (2.0 * halfExtentKms);
  const dxMpc = overrides.dxMpc ?? (boxSizeMpc / nx);
  const dxKms = overrides.dxKms ?? (boxSizeKms / nx);
  const coordinateSystem = overrides.coordinateSystem ?? 'Supergalactic Cartesian';

  return {
    nx,
    ny,
    nz,
    boxSizeMpc,
    boxSizeKms,
    halfExtentMpc,
    halfExtentKms,
    dxMpc,
    dxKms,
    hubbleH0,
    coordinateSystem
  };
}

/**
 * 3D Voxel Grid container holding integer basin labels or float field values.
 */
export class VoxelGrid {
  /**
   * @param {number} nx - Dimension X
   * @param {number} ny - Dimension Y
   * @param {number} nz - Dimension Z
   * @param {TypedArray|ArrayBuffer} [data] - Optional backing typed array or ArrayBuffer
   * @param {GridMetadata} [metadata] - Coordinate metadata
   * @param {string} [dataType='Uint16Array'] - Type of array if allocating ('Uint8Array'|'Uint16Array'|'Int32Array'|'Float32Array')
   */
  constructor(nx, ny, nz, data = null, metadata = null, dataType = 'Uint16Array') {
    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.totalVoxels = nx * ny * nz;
    this.metadata = metadata || createGridMetadata({ nx, ny, nz });

    if (data instanceof ArrayBuffer) {
      if (dataType === 'Uint8Array') this.data = new Uint8Array(data);
      else if (dataType === 'Int32Array') this.data = new Int32Array(data);
      else if (dataType === 'Float32Array') this.data = new Float32Array(data);
      else this.data = new Uint16Array(data);
    } else if (ArrayBuffer.isView(data)) {
      this.data = data;
    } else {
      if (dataType === 'Uint8Array') this.data = new Uint8Array(this.totalVoxels);
      else if (dataType === 'Int32Array') this.data = new Int32Array(this.totalVoxels);
      else if (dataType === 'Float32Array') this.data = new Float32Array(this.totalVoxels);
      else this.data = new Uint16Array(this.totalVoxels);
    }
  }

  /**
   * Computes flat array index for (ix, iy, iz).
   * Row-major indexing: index = (ix * ny + iy) * nz + iz.
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {number}
   */
  getIndex(ix, iy, iz) {
    return (ix * this.ny + iy) * this.nz + iz;
  }

  /**
   * Unflattens a flat index into [ix, iy, iz] grid indices.
   * @param {number} idx
   * @returns {[number, number, number]}
   */
  getGridCoords(idx) {
    const iz = idx % this.nz;
    const rem = Math.floor(idx / this.nz);
    const iy = rem % this.ny;
    const ix = Math.floor(rem / this.ny);
    return [ix, iy, iz];
  }

  /**
   * Retrieves voxel value at (ix, iy, iz).
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {number}
   */
  get(ix, iy, iz) {
    if (ix < 0 || ix >= this.nx || iy < 0 || iy >= this.ny || iz < 0 || iz >= this.nz) {
      return 0;
    }
    return this.data[(ix * this.ny + iy) * this.nz + iz];
  }

  /**
   * Sets voxel value at (ix, iy, iz).
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @param {number} val
   */
  set(ix, iy, iz, val) {
    if (ix >= 0 && ix < this.nx && iy >= 0 && iy < this.ny && iz >= 0 && iz < this.nz) {
      this.data[(ix * this.ny + iy) * this.nz + iz] = val;
    }
  }

  /**
   * Transforms grid indices [ix, iy, iz] to physical Supergalactic coordinates.
   * Center of voxel (ix + 0.5) is mapped to physical space.
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @param {'Mpc'|'kms'} [unit='Mpc']
   * @returns {{x: number, y: number, z: number}}
   */
  gridToWorld(ix, iy, iz, unit = 'Mpc') {
    const halfExtent = unit === 'kms' ? this.metadata.halfExtentKms : this.metadata.halfExtentMpc;
    const boxSize = unit === 'kms' ? this.metadata.boxSizeKms : this.metadata.boxSizeMpc;

    const x = ((ix + 0.5) / this.nx) * boxSize - halfExtent;
    const y = ((iy + 0.5) / this.ny) * boxSize - halfExtent;
    const z = ((iz + 0.5) / this.nz) * boxSize - halfExtent;

    return { x, y, z };
  }

  /**
   * Transforms physical Supergalactic coordinates to continuous and discrete grid indices.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {'Mpc'|'kms'} [unit='Mpc']
   * @returns {{ix: number, iy: number, iz: number, fx: number, fy: number, fz: number, inside: boolean}}
   */
  worldToGrid(x, y, z, unit = 'Mpc') {
    const halfExtent = unit === 'kms' ? this.metadata.halfExtentKms : this.metadata.halfExtentMpc;
    const boxSize = unit === 'kms' ? this.metadata.boxSizeKms : this.metadata.boxSizeMpc;

    const fx = ((x + halfExtent) / boxSize) * this.nx - 0.5;
    const fy = ((y + halfExtent) / boxSize) * this.ny - 0.5;
    const fz = ((z + halfExtent) / boxSize) * this.nz - 0.5;

    const ix = Math.round(fx);
    const iy = Math.round(fy);
    const iz = Math.round(fz);

    const inside = ix >= 0 && ix < this.nx && iy >= 0 && iy < this.ny && iz >= 0 && iz < this.nz;

    return {
      ix: Math.max(0, Math.min(this.nx - 1, ix)),
      iy: Math.max(0, Math.min(this.ny - 1, iy)),
      iz: Math.max(0, Math.min(this.nz - 1, iz)),
      fx,
      fy,
      fz,
      inside
    };
  }

  /**
   * Clones this VoxelGrid with independent data buffer.
   * @returns {VoxelGrid}
   */
  clone() {
    const copyData = new this.data.constructor(this.data.length);
    copyData.set(this.data);
    return new VoxelGrid(this.nx, this.ny, this.nz, copyData, { ...this.metadata });
  }
}

/**
 * Parses raw FITS 3D image buffer or binary segmentation array.
 * Reads FITS header cards or standard binary format into a VoxelGrid.
 *
 * @param {ArrayBuffer|Uint8Array} buffer - Binary buffer containing FITS or raw voxel data
 * @param {Object} [options={}] - Parsing configuration
 * @param {number} [options.nx=128]
 * @param {number} [options.ny=128]
 * @param {number} [options.nz=128]
 * @param {number} [options.halfExtentMpc=500.0]
 * @param {number} [options.hubbleH0=74.6]
 * @returns {VoxelGrid}
 */
export function loadSegmentationBuffer(buffer, options = {}) {
  const byteBuf = buffer instanceof Uint8Array ? buffer.buffer : buffer;
  const view = new DataView(byteBuf);

  // Check if buffer starts with FITS standard ASCII keyword "SIMPLE  ="
  const isFits = byteBuf.byteLength >= 80 &&
    String.fromCharCode(view.getUint8(0), view.getUint8(1), view.getUint8(2), view.getUint8(3), view.getUint8(4), view.getUint8(5)) === 'SIMPLE';

  if (isFits) {
    return parseFitsSegmentation(byteBuf, options);
  }

  // Raw binary buffer: determine dimensions from byte length or options
  let nx = options.nx ?? 128;
  let ny = options.ny ?? 128;
  let nz = options.nz ?? 128;

  const totalVoxels = nx * ny * nz;
  const byteLen = byteBuf.byteLength;

  let typedData;
  if (byteLen === totalVoxels * 1) {
    typedData = new Uint8Array(byteBuf);
  } else if (byteLen === totalVoxels * 2) {
    typedData = new Uint16Array(byteBuf);
  } else if (byteLen === totalVoxels * 4) {
    // Determine if int or float
    if (options.isFloat) {
      typedData = new Float32Array(byteBuf);
    } else {
      typedData = new Int32Array(byteBuf);
    }
  } else {
    // Attempt dimension inference for cubic volume
    const elementSizes = [1, 2, 4];
    let matched = false;
    for (const sz of elementSizes) {
      const numElem = byteLen / sz;
      const cubeRoot = Math.round(Math.cbrt(numElem));
      if (cubeRoot * cubeRoot * cubeRoot === numElem) {
        nx = ny = nz = cubeRoot;
        if (sz === 1) typedData = new Uint8Array(byteBuf);
        else if (sz === 2) typedData = new Uint16Array(byteBuf);
        else typedData = new Int32Array(byteBuf);
        matched = true;
        break;
      }
    }
    if (!matched) {
      typedData = new Uint16Array(byteBuf);
    }
  }

  const meta = createGridMetadata({
    nx,
    ny,
    nz,
    halfExtentMpc: options.halfExtentMpc ?? 500.0,
    hubbleH0: options.hubbleH0 ?? 74.6
  });

  return new VoxelGrid(nx, ny, nz, typedData, meta);
}

/**
 * Parses standard FITS 3D primary HDU containing 3D image arrays.
 * Handles big-endian header cards and binary payload.
 *
 * @param {ArrayBuffer} buffer
 * @param {Object} options
 * @returns {VoxelGrid}
 */
export function parseFitsSegmentation(buffer, options = {}) {
  const bytes = new Uint8Array(buffer);
  let headerEnd = -1;
  const cards = {};

  // Scan 2880-byte FITS header blocks
  for (let blockStart = 0; blockStart < bytes.length; blockStart += 2880) {
    for (let cardStart = blockStart; cardStart < blockStart + 2880; cardStart += 80) {
      const cardStr = String.fromCharCode(...bytes.subarray(cardStart, cardStart + 80));
      const key = cardStr.substring(0, 8).trim();
      if (key === 'END') {
        headerEnd = blockStart + 2880;
        break;
      }
      if (cardStr.includes('=')) {
        const valStr = cardStr.substring(10, 30).trim();
        cards[key] = isNaN(Number(valStr)) ? valStr.replace(/'/g, '').trim() : Number(valStr);
      }
    }
    if (headerEnd !== -1) break;
  }

  if (headerEnd === -1) {
    throw new Error('Invalid FITS file: END keyword not found in header blocks.');
  }

  const naxis = cards['NAXIS'] || 3;
  const nx = cards['NAXIS1'] || options.nx || 128;
  const ny = cards['NAXIS2'] || options.ny || 128;
  const nz = cards['NAXIS3'] || options.nz || 128;
  const bitpix = cards['BITPIX'] || 16;
  const bzero = cards['BZERO'] ?? 0;
  const bscale = cards['BSCALE'] ?? 1;

  const dataView = new DataView(buffer, headerEnd);
  const total = nx * ny * nz;
  const result = new Uint16Array(total);

  if (bitpix === 8) {
    for (let i = 0; i < total; i++) {
      result[i] = Math.round(dataView.getUint8(i) * bscale + bzero);
    }
  } else if (bitpix === 16) {
    for (let i = 0; i < total; i++) {
      // FITS standard is big-endian
      result[i] = Math.round(dataView.getInt16(i * 2, false) * bscale + bzero);
    }
  } else if (bitpix === 32) {
    for (let i = 0; i < total; i++) {
      result[i] = Math.round(dataView.getInt32(i * 4, false) * bscale + bzero);
    }
  } else if (bitpix === -32) {
    for (let i = 0; i < total; i++) {
      result[i] = Math.round(dataView.getFloat32(i * 4, false) * bscale + bzero);
    }
  } else {
    throw new Error(`Unsupported FITS BITPIX: ${bitpix}`);
  }

  const meta = createGridMetadata({
    nx,
    ny,
    nz,
    halfExtentMpc: options.halfExtentMpc ?? 500.0,
    hubbleH0: options.hubbleH0 ?? 74.6
  });

  return new VoxelGrid(nx, ny, nz, result, meta);
}

/**
 * Looks up basin metadata by numeric ID.
 * @param {number} basinId
 * @returns {BasinDefinition}
 */
export function getBasinById(basinId) {
  if (EXTENDED_BASIN_TAXONOMY[basinId]) {
    return EXTENDED_BASIN_TAXONOMY[basinId];
  }
  return {
    id: basinId,
    name: `Basin-${basinId}`,
    abbreviation: `B${basinId}`,
    description: `Dynamic Watershed Basin #${basinId}`,
    attractorSGKms: [0, 0, 0],
    attractorSGMpc: [0, 0, 0],
    volumeCF4Ind1e6: 0,
    volumeCF4Grp1e6: 0,
    charVelocityKms: 0,
    distToShapleyMpc: 0,
    hexColor: '#94a3b8',
    rgbColor: [0.580, 0.639, 0.722]
  };
}

/**
 * Looks up basin metadata by name (case-insensitive fuzzy match).
 * @param {string} name
 * @returns {BasinDefinition|null}
 */
export function getBasinByName(name) {
  if (!name) return null;
  const search = name.trim().toLowerCase();
  for (const basin of Object.values(EXTENDED_BASIN_TAXONOMY)) {
    if (basin.name.toLowerCase() === search || basin.abbreviation.toLowerCase() === search) {
      return basin;
    }
  }
  return null;
}

/**
 * Returns list of all official Table A.1 basins.
 * @returns {BasinDefinition[]}
 */
export function getAllBasins() {
  return Object.values(BASIN_TAXONOMY);
}

/**
 * Computes spatial statistics (voxel count, volume, centroid, bounding box, inertia tensor)
 * for each basin in a VoxelGrid segmentation.
 *
 * @param {VoxelGrid} grid - Input segmentation grid
 * @returns {Record<number, {
 *   basinId: number,
 *   name: string,
 *   voxelCount: number,
 *   volumeMpc3: number,
 *   volume1e6Mpc3: number,
 *   centroidMpc: [number, number, number],
 *   centroidKms: [number, number, number],
 *   boundingBoxMpc: { minX: number, maxX: number, minY: number, maxY: number, minZ: number, maxZ: number },
 *   inertiaTensor: number[][],
 *   principalSemiAxesMpc: [number, number, number]
 * }>}
 */
export function computeBasinStatistics(grid) {
  const stats = {};
  const voxelVolMpc3 = grid.metadata.dxMpc * grid.metadata.dxMpc * grid.metadata.dxMpc;

  const total = grid.totalVoxels;
  for (let idx = 0; idx < total; idx++) {
    const label = grid.data[idx];
    if (!stats[label]) {
      const basinDef = getBasinById(label);
      stats[label] = {
        basinId: label,
        name: basinDef.name,
        voxelCount: 0,
        sumX: 0,
        sumY: 0,
        sumZ: 0,
        minIx: grid.nx,
        maxIx: -1,
        minIy: grid.ny,
        maxIy: -1,
        minIz: grid.nz,
        maxIz: -1,
        // Second moments for inertia
        sumXX: 0,
        sumYY: 0,
        sumZZ: 0,
        sumXY: 0,
        sumXZ: 0,
        sumYZ: 0
      };
    }

    const s = stats[label];
    s.voxelCount++;

    const [ix, iy, iz] = grid.getGridCoords(idx);
    const { x, y, z } = grid.gridToWorld(ix, iy, iz, 'Mpc');

    s.sumX += x;
    s.sumY += y;
    s.sumZ += z;

    s.sumXX += x * x;
    s.sumYY += y * y;
    s.sumZZ += z * z;
    s.sumXY += x * y;
    s.sumXZ += x * z;
    s.sumYZ += y * z;

    if (ix < s.minIx) s.minIx = ix;
    if (ix > s.maxIx) s.maxIx = ix;
    if (iy < s.minIy) s.minIy = iy;
    if (iy > s.maxIy) s.maxIy = iy;
    if (iz < s.minIz) s.minIz = iz;
    if (iz > s.maxIz) s.maxIz = iz;
  }

  const result = {};
  for (const [key, s] of Object.entries(stats)) {
    const n = s.voxelCount;
    const cx = s.sumX / n;
    const cy = s.sumY / n;
    const cz = s.sumZ / n;

    const minWorld = grid.gridToWorld(s.minIx, s.minIy, s.minIz, 'Mpc');
    const maxWorld = grid.gridToWorld(s.maxIx, s.maxIy, s.maxIz, 'Mpc');

    // Central second moments
    const Ixx = (s.sumXX / n) - (cx * cx);
    const Iyy = (s.sumYY / n) - (cy * cy);
    const Izz = (s.sumZZ / n) - (cz * cz);
    const Ixy = (s.sumXY / n) - (cx * cy);
    const Ixz = (s.sumXZ / n) - (cx * cz);
    const Iyz = (s.sumYZ / n) - (cy * cz);

    const volumeMpc3 = n * voxelVolMpc3;
    const h0 = grid.metadata.hubbleH0;

    result[key] = {
      basinId: s.basinId,
      name: s.name,
      voxelCount: n,
      volumeMpc3,
      volume1e6Mpc3: volumeMpc3 / 1e6,
      centroidMpc: [cx, cy, cz],
      centroidKms: [cx * h0, cy * h0, cz * h0],
      boundingBoxMpc: {
        minX: minWorld.x - grid.metadata.dxMpc * 0.5,
        maxX: maxWorld.x + grid.metadata.dxMpc * 0.5,
        minY: minWorld.y - grid.metadata.dxMpc * 0.5,
        maxY: maxWorld.y + grid.metadata.dxMpc * 0.5,
        minZ: minWorld.z - grid.metadata.dxMpc * 0.5,
        maxZ: maxWorld.z + grid.metadata.dxMpc * 0.5
      },
      inertiaTensor: [
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz]
      ],
      principalSemiAxesMpc: [
        Math.sqrt(Math.max(0, Ixx)),
        Math.sqrt(Math.max(0, Iyy)),
        Math.sqrt(Math.max(0, Izz))
      ]
    };
  }

  return result;
}
