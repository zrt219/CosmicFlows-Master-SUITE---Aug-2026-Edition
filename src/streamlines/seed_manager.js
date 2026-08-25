/**
 * @file seed_manager.js
 * @description Multi-class streamline seeding manager supporting:
 *  - GRID_VOXEL: Regular 64^3 / 128^3 / custom 3D grid lattices with strides, jitter, and density masking.
 *  - OBSERVED_OBJECT: Real CF4 (Cosmicflows-4) grouped and ungrouped galaxy catalogs with full astrometric metadata.
 *  - VISUALIZATION: Fibonacci spherical shells, Supergalactic planar slices, Poisson/random distributions, and interactive cursor seeds.
 *
 * @module streamlines/seed_manager
 */

/**
 * Seeding generation modes.
 * @readonly
 * @enum {string}
 */
export const SeedMode = Object.freeze({
  GRID_VOXEL: 'GRID_VOXEL',
  OBSERVED_OBJECT: 'OBSERVED_OBJECT',
  VISUALIZATION: 'VISUALIZATION'
});

/**
 * Well-known Cosmographic Attractor and Repeller Landmarks (Coordinates in SGX, SGY, SGZ [h^-1 Mpc]).
 */
export const COSMIC_LANDMARKS = Object.freeze({
  LOCAL_GROUP: { id: 'local_group', name: 'Local Group', pos: [0.0, 0.0, 0.0], radius: 5.0 },
  VIRGO_CLUSTER: { id: 'virgo_cl', name: 'Virgo Cluster', pos: [-2.1, 11.8, -0.6], radius: 10.0 },
  COMA_CLUSTER: { id: 'coma_cl', name: 'Coma Cluster', pos: [7.2, 69.1, 6.8], radius: 15.0 },
  CENTAURUS_CLUSTER: { id: 'centaurus_cl', name: 'Centaurus Cluster', pos: [-31.2, 14.5, 3.2], radius: 12.0 },
  NORMA_GREAT_ATTRACTOR: { id: 'ga_norma', name: 'Great Attractor / Norma', pos: [-44.0, 11.0, -9.0], radius: 20.0 },
  SHAPLEY_CORE: { id: 'shapley_core', name: 'Shapley Supercluster Core', pos: [-105.0, 52.0, 15.0], radius: 30.0 },
  DIPOLE_REPELLER: { id: 'dipole_repeller', name: 'Dipole Repeller', pos: [85.0, -110.0, -45.0], radius: 25.0 },
  PERSEUS_PISCES: { id: 'perseus_pisces', name: 'Perseus-Pisces Supercluster', pos: [48.0, -18.0, -15.0], radius: 25.0 },
  VELA_SUPERCLUSTER: { id: 'vela_scl', name: 'Vela Supercluster', pos: [-130.0, -45.0, -60.0], radius: 35.0 },
  FORNAX_CLUSTER: { id: 'fornax_cl', name: 'Fornax Cluster', pos: [-2.5, -14.2, -12.8], radius: 8.0 }
});

/**
 * Represents an individual seed with 3D coordinates, classification, and arbitrary metadata.
 */
export class StreamlineSeed {
  /**
   * @param {Array<number>} position - [x, y, z] in simulation units (typically h^-1 Mpc).
   * @param {string} [mode=SeedMode.VISUALIZATION] - Generation mode.
   * @param {Object} [metadata={}] - Custom astrometric or visualization metadata.
   */
  constructor(position, mode = SeedMode.VISUALIZATION, metadata = {}) {
    if (!position || position.length < 3) {
      throw new Error(`Invalid seed position: ${JSON.stringify(position)}`);
    }
    this.id = metadata.id || `seed_${Math.random().toString(36).substr(2, 9)}`;
    this.x = Number(position[0]);
    this.y = Number(position[1]);
    this.z = Number(position[2]);
    this.mode = mode;
    this.metadata = { ...metadata };
  }

  get position() {
    return [this.x, this.y, this.z];
  }

  toJSON() {
    return {
      id: this.id,
      pos: [this.x, this.y, this.z],
      mode: this.mode,
      metadata: this.metadata
    };
  }
}

/**
 * Multi-Class Streamline Seeding Engine.
 */
export class SeedManager {
  constructor() {
    /** @type {Map<string, StreamlineSeed>} */
    this.seeds = new Map();
  }

  /**
   * Clears all stored seeds.
   */
  clear() {
    this.seeds.clear();
  }

  /**
   * Returns current seed count.
   * @returns {number}
   */
  get size() {
    return this.seeds.size;
  }

  /**
   * Returns array of all registered seeds.
   * @returns {Array<StreamlineSeed>}
   */
  getAllSeeds() {
    return Array.from(this.seeds.values());
  }

  /**
   * Adds a single seed object.
   * @param {StreamlineSeed|Array<number>} seedOrPos
   * @param {string} [mode]
   * @param {Object} [metadata]
   * @returns {StreamlineSeed}
   */
  addSeed(seedOrPos, mode = SeedMode.VISUALIZATION, metadata = {}) {
    let seed;
    if (seedOrPos instanceof StreamlineSeed) {
      seed = seedOrPos;
    } else {
      seed = new StreamlineSeed(seedOrPos, mode, metadata);
    }
    this.seeds.set(seed.id, seed);
    return seed;
  }

  /**
   * Adds a batch of seeds.
   * @param {Array<StreamlineSeed|Array<number>>} seedList
   * @param {string} [defaultMode]
   * @returns {number} Count of added seeds.
   */
  addBatch(seedList, defaultMode = SeedMode.VISUALIZATION) {
    if (!Array.isArray(seedList)) return 0;
    let count = 0;
    for (let i = 0; i < seedList.length; i++) {
      const item = seedList[i];
      if (item instanceof StreamlineSeed) {
        this.seeds.set(item.id, item);
      } else if (Array.isArray(item) && item.length >= 3) {
        const s = new StreamlineSeed(item, defaultMode, { batchIndex: count });
        this.seeds.set(s.id, s);
      }
      count++;
    }
    return count;
  }

  // =========================================================================
  // 1. GRID_VOXEL SEED GENERATORS
  // =========================================================================

  /**
   * Generates regular 3D voxel grid seeds across a specified bounding box.
   *
   * @param {Object} [options={}]
   * @param {Array<number>} [options.resolution=[64, 64, 64]] - [Nx, Ny, Nz] grid dimensions.
   * @param {Object|Array<Array<number>>} [options.domain=[[-100, 100], [-100, 100], [-100, 100]]] - Bounding box.
   * @param {Array<number>} [options.strides=[1, 1, 1]] - Stride step in voxel cells [sx, sy, sz].
   * @param {boolean} [options.jitter=false] - Apply random cell jitter.
   * @param {number} [options.jitterFraction=0.5] - Jitter fraction relative to cell size.
   * @param {Function|null} [options.maskFn=null] - Filtering callback (x, y, z, i, j, k) => boolean.
   * @param {boolean} [options.addToStore=true] - Whether to store seeds in instance.
   * @returns {Array<StreamlineSeed>} Generated seeds.
   */
  generateGridSeeds(options = {}) {
    const {
      resolution = [64, 64, 64],
      domain = [[-100, 100], [-100, 100], [-100, 100]],
      strides = [1, 1, 1],
      jitter = false,
      jitterFraction = 0.5,
      maskFn = null,
      addToStore = true
    } = options;

    const [nx, ny, nz] = resolution;
    const [sx, sy, sz] = strides;

    let xMin = -100, xMax = 100, yMin = -100, yMax = 100, zMin = -100, zMax = 100;
    if (Array.isArray(domain)) {
      [[xMin, xMax], [yMin, yMax], [zMin, zMax]] = domain;
    } else if (domain) {
      ({ xMin = -100, xMax = 100, yMin = -100, yMax = 100, zMin = -100, zMax = 100 } = domain);
    }

    const dx = nx > 1 ? (xMax - xMin) / (nx - 1) : 0;
    const dy = ny > 1 ? (yMax - yMin) / (ny - 1) : 0;
    const dz = nz > 1 ? (zMax - zMin) / (nz - 1) : 0;

    const generated = [];

    for (let i = 0; i < nx; i += sx) {
      const baseX = xMin + i * dx;
      for (let j = 0; j < ny; j += sy) {
        const baseY = yMin + j * dy;
        for (let k = 0; k < nz; k += sz) {
          const baseZ = zMin + k * dz;

          let posX = baseX;
          let posY = baseY;
          let posZ = baseZ;

          if (jitter) {
            posX += (Math.random() - 0.5) * dx * jitterFraction * 2.0;
            posY += (Math.random() - 0.5) * dy * jitterFraction * 2.0;
            posZ += (Math.random() - 0.5) * dz * jitterFraction * 2.0;
          }

          if (maskFn && !maskFn(posX, posY, posZ, i, j, k)) {
            continue;
          }

          const seed = new StreamlineSeed([posX, posY, posZ], SeedMode.GRID_VOXEL, {
            gridIndex: [i, j, k],
            gridResolution: [nx, ny, nz]
          });

          generated.push(seed);
          if (addToStore) {
            this.seeds.set(seed.id, seed);
          }
        }
      }
    }

    return generated;
  }

  // =========================================================================
  // 2. OBSERVED_OBJECT SEED GENERATORS (CF4 Catalog)
  // =========================================================================

  /**
   * Ingests Cosmicflows-4 (CF4) galaxy catalog records.
   *
   * @param {Array<Object>} catalogEntries - Array of galaxy entries with pos [SGX, SGY, SGZ] or {sgx, sgy, sgz}.
   * @param {Object} [filterOptions={}] - Filtering criteria.
   * @param {boolean} [addToStore=true]
   * @returns {Array<StreamlineSeed>}
   */
  loadObservedCatalog(catalogEntries, filterOptions = {}, addToStore = true) {
    if (!Array.isArray(catalogEntries)) return [];

    const {
      maxDistance = Infinity,
      minDistance = 0,
      center = [0, 0, 0],
      groupOnly = false,
      ungroupedOnly = false,
      minMassWeight = -Infinity,
      bounds = null,
      maxEntries = Infinity
    } = filterOptions;

    const generated = [];

    for (let i = 0; i < catalogEntries.length; i++) {
      if (generated.length >= maxEntries) break;
      const entry = catalogEntries[i];

      let pos;
      if (Array.isArray(entry.pos)) {
        pos = entry.pos;
      } else if (Array.isArray(entry.position)) {
        pos = entry.position;
      } else if (entry.sgx !== undefined && entry.sgy !== undefined && entry.sgz !== undefined) {
        pos = [entry.sgx, entry.sgy, entry.sgz];
      } else if (entry.x !== undefined && entry.y !== undefined && entry.z !== undefined) {
        pos = [entry.x, entry.y, entry.z];
      } else {
        continue;
      }

      const x = Number(pos[0]);
      const y = Number(pos[1]);
      const z = Number(pos[2]);

      const dist = Math.hypot(x - center[0], y - center[1], z - center[2]);
      if (dist < minDistance || dist > maxDistance) continue;

      if (bounds) {
        let xMin, xMax, yMin, yMax, zMin, zMax;
        if (Array.isArray(bounds)) {
          [[xMin, xMax], [yMin, yMax], [zMin, zMax]] = bounds;
        } else {
          ({ xMin = -Infinity, xMax = Infinity, yMin = -Infinity, yMax = Infinity, zMin = -Infinity, zMax = Infinity } = bounds);
        }
        if (x < xMin || x > xMax || y < yMin || y > yMax || z < zMin || z > zMax) continue;
      }

      const isGroup = !!(entry.isGroupMember || (entry.groupId && entry.groupId !== 0 && entry.groupId !== -1));
      if (groupOnly && !isGroup) continue;
      if (ungroupedOnly && isGroup) continue;

      const massWeight = entry.massWeight ?? entry.luminosity ?? 1.0;
      if (massWeight < minMassWeight) continue;

      const seed = new StreamlineSeed([x, y, z], SeedMode.OBSERVED_OBJECT, {
        catalogId: entry.id || entry.pgc || `gal_${i}`,
        name: entry.name || `CF4_${entry.pgc || i}`,
        cz: entry.cz ?? 0,
        dMpc: entry.dMpc ?? dist,
        vpec: entry.vpec ?? [0, 0, 0],
        groupId: entry.groupId ?? null,
        isGroupMember: isGroup,
        uncertainty: entry.uncertainty ?? entry.dErr ?? 0.0,
        massWeight
      });

      generated.push(seed);
      if (addToStore) {
        this.seeds.set(seed.id, seed);
      }
    }

    return generated;
  }

  // =========================================================================
  // 3. VISUALIZATION SEED GENERATORS (Shells, Slices, Poisson, Landmarks)
  // =========================================================================

  /**
   * Generates spherical shell seeds using golden spiral / Fibonacci sphere lattice.
   *
   * @param {Array<number>} center - [cx, cy, cz].
   * @param {number} radius - Shell radius in simulation units.
   * @param {number} [count=100] - Number of points on sphere.
   * @param {Object} [metadata={}]
   * @param {boolean} [addToStore=true]
   * @returns {Array<StreamlineSeed>}
   */
  generateFibonacciSphereSeeds(center, radius, count = 100, metadata = {}, addToStore = true) {
    const generated = [];
    const goldenRatio = (1.0 + Math.sqrt(5.0)) / 2.0;
    const goldenAngle = 2.0 * Math.PI * (1.0 - 1.0 / goldenRatio);

    for (let i = 0; i < count; i++) {
      const y = 1.0 - (i / (count - 1)) * 2.0; // from 1 to -1
      const radiusAtY = Math.sqrt(Math.max(0.0, 1.0 - y * y));
      const theta = goldenAngle * i;

      const x = Math.cos(theta) * radiusAtY;
      const z = Math.sin(theta) * radiusAtY;

      const px = center[0] + radius * x;
      const py = center[1] + radius * y;
      const pz = center[2] + radius * z;

      const seed = new StreamlineSeed([px, py, pz], SeedMode.VISUALIZATION, {
        type: 'FIBONACCI_SPHERE',
        center: [...center],
        radius,
        index: i,
        ...metadata
      });

      generated.push(seed);
      if (addToStore) {
        this.seeds.set(seed.id, seed);
      }
    }

    return generated;
  }

  /**
   * Generates seeds distributed around cosmological landmark attractors or repellers.
   *
   * @param {string} landmarkKey - Key from COSMIC_LANDMARKS (e.g. 'SHAPLEY_CORE', 'VIRGO_CLUSTER').
   * @param {number} [shellCount=3] - Number of concentric shells.
   * @param {number} [pointsPerShell=50] - Number of seeds per shell.
   * @param {boolean} [addToStore=true]
   * @returns {Array<StreamlineSeed>}
   */
  generateLandmarkSeeds(landmarkKey, shellCount = 3, pointsPerShell = 50, addToStore = true) {
    const landmark = COSMIC_LANDMARKS[landmarkKey.toUpperCase()] ||
      Object.values(COSMIC_LANDMARKS).find(lm => lm.id === landmarkKey || lm.name.toLowerCase() === landmarkKey.toLowerCase());

    if (!landmark) {
      throw new Error(`Cosmic landmark '${landmarkKey}' not found in known registry.`);
    }

    const generated = [];
    const baseRadius = landmark.radius || 15.0;

    for (let s = 1; s <= shellCount; s++) {
      const r = (s / shellCount) * baseRadius;
      const shellSeeds = this.generateFibonacciSphereSeeds(landmark.pos, r, pointsPerShell, {
        landmarkId: landmark.id,
        landmarkName: landmark.name,
        shellIndex: s
      }, false);

      for (let i = 0; i < shellSeeds.length; i++) {
        generated.push(shellSeeds[i]);
        if (addToStore) {
          this.seeds.set(shellSeeds[i].id, shellSeeds[i]);
        }
      }
    }

    return generated;
  }

  /**
   * Generates seeds across a 2D planar slice in 3D space (e.g. Supergalactic plane SGZ = 0).
   *
   * @param {Object} [options={}]
   * @param {string} [options.plane='SGZ'] - 'SGX', 'SGY', or 'SGZ'.
   * @param {number} [options.sliceOffset=0.0] - Coordinate offset on normal axis.
   * @param {Array<number>} [options.range1=[-100, 100]] - Span of first in-plane axis.
   * @param {Array<number>} [options.range2=[-100, 100]] - Span of second in-plane axis.
   * @param {Array<number>} [options.gridSize=[32, 32]] - [N1, N2] resolution.
   * @param {boolean} [addToStore=true]
   * @returns {Array<StreamlineSeed>}
   */
  generatePlanarSliceSeeds(options = {}, addToStore = true) {
    const {
      plane = 'SGZ',
      sliceOffset = 0.0,
      range1 = [-100, 100],
      range2 = [-100, 100],
      gridSize = [32, 32]
    } = options;

    const [n1, n2] = gridSize;
    const [min1, max1] = range1;
    const [min2, max2] = range2;

    const d1 = n1 > 1 ? (max1 - min1) / (n1 - 1) : 0;
    const d2 = n2 > 1 ? (max2 - min2) / (n2 - 1) : 0;

    const generated = [];

    for (let i = 0; i < n1; i++) {
      const u = min1 + i * d1;
      for (let j = 0; j < n2; j++) {
        const v = min2 + j * d2;
        let pos;
        if (plane.toUpperCase() === 'SGZ' || plane.toUpperCase() === 'XY') {
          pos = [u, v, sliceOffset];
        } else if (plane.toUpperCase() === 'SGY' || plane.toUpperCase() === 'XZ') {
          pos = [u, sliceOffset, v];
        } else {
          pos = [sliceOffset, u, v];
        }

        const seed = new StreamlineSeed(pos, SeedMode.VISUALIZATION, {
          plane,
          sliceOffset,
          planeCoords: [u, v],
          gridIndex: [i, j]
        });

        generated.push(seed);
        if (addToStore) {
          this.seeds.set(seed.id, seed);
        }
      }
    }

    return generated;
  }

  /**
   * Generates uniform random or pseudo-Poisson disk sampled points within a bounding sphere.
   *
   * @param {Array<number>} center
   * @param {number} radius
   * @param {number} count
   * @param {number} [minDistance=0] - Optional minimum spacing for Poisson-like distribution.
   * @param {boolean} [addToStore=true]
   * @returns {Array<StreamlineSeed>}
   */
  generateSphericalDistribution(center, radius, count = 100, minDistance = 0, addToStore = true) {
    const generated = [];
    const maxAttempts = count * 50;
    let attempts = 0;

    while (generated.length < count && attempts < maxAttempts) {
      attempts++;
      // Uniform distribution in sphere volume
      const u = Math.random();
      const r = radius * Math.cbrt(u);
      const theta = Math.random() * 2.0 * Math.PI;
      const phi = Math.acos(2.0 * Math.random() - 1.0);

      const px = center[0] + r * Math.sin(phi) * Math.cos(theta);
      const py = center[1] + r * Math.sin(phi) * Math.sin(theta);
      const pz = center[2] + r * Math.cos(phi);

      if (minDistance > 0) {
        let tooClose = false;
        for (let k = 0; k < generated.length; k++) {
          const g = generated[k];
          const d = Math.hypot(px - g.x, py - g.y, pz - g.z);
          if (d < minDistance) {
            tooClose = true;
            break;
          }
        }
        if (tooClose) continue;
      }

      const seed = new StreamlineSeed([px, py, pz], SeedMode.VISUALIZATION, {
        type: 'SPHERE_UNIFORM',
        center: [...center],
        radius
      });

      generated.push(seed);
      if (addToStore) {
        this.seeds.set(seed.id, seed);
      }
    }

    return generated;
  }

  // =========================================================================
  // UTILITIES & EXPORTS
  // =========================================================================

  /**
   * Partitions current seeds into uniform batches for parallel/worker processing.
   *
   * @param {number} [batchSize=100]
   * @returns {Array<Array<StreamlineSeed>>}
   */
  createBatches(batchSize = 100) {
    const all = this.getAllSeeds();
    const batches = [];
    for (let i = 0; i < all.length; i += batchSize) {
      batches.push(all.slice(i, i + batchSize));
    }
    return batches;
  }

  /**
   * Exports all seed coordinates as a flat Float32Array [x0, y0, z0, x1, y1, z1, ...].
   * @returns {Float32Array}
   */
  toFloat32Array() {
    const all = this.getAllSeeds();
    const arr = new Float32Array(all.length * 3);
    for (let i = 0; i < all.length; i++) {
      arr[i * 3 + 0] = all[i].x;
      arr[i * 3 + 1] = all[i].y;
      arr[i * 3 + 2] = all[i].z;
    }
    return arr;
  }

  /**
   * Serializes seed manager contents to plain JSON structure.
   * @returns {Array<Object>}
   */
  toJSON() {
    return this.getAllSeeds().map(s => s.toJSON());
  }

  /**
   * Deserializes seeds from JSON structure.
   * @param {Array<Object>} jsonList
   * @param {boolean} [clearExisting=true]
   */
  fromJSON(jsonList, clearExisting = true) {
    if (clearExisting) this.clear();
    if (!Array.isArray(jsonList)) return;
    for (let i = 0; i < jsonList.length; i++) {
      const item = jsonList[i];
      const pos = item.pos || [item.x, item.y, item.z];
      const seed = new StreamlineSeed(pos, item.mode, { ...item.metadata, id: item.id });
      this.seeds.set(seed.id, seed);
    }
  }
}

export default SeedManager;
