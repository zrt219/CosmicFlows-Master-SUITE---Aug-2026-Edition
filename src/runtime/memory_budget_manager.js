/**
 * @file memory_budget_manager.js
 * @module runtime/memory_budget_manager
 * @description Production-grade WebGL / Web Worker VRAM & System Memory Budget Manager,
 * LRU Cache Eviction Sentinel, High-Frequency Telemetry Poller, and OOM Prevention Guard
 * for the ZRT Cosmicflows Scientific Computational Engine.
 *
 * Implements:
 * 1. Exact VRAM memory budgeting across WebGL buffers (geometries, textures 2D/3D, instanced meshes,
 *    particle arrays, render targets, materials, shader programs) and Web Worker transferable ArrayBuffers.
 * 2. Soft (512 MB default) and Hard (1024 MB default) memory thresholds with multi-level hysteresis
 *    and state transitions (NORMAL, MODERATE, CRITICAL, EXCEEDED).
 * 3. Multi-tier priority-weighted LRU eviction with automatic Three.js BufferGeometry/Material/Texture
 *    lifecycle disposal, dependency graph resolution, and ref-counting.
 * 4. High-frequency telemetry polling (10-100 Hz), moving-average smoothing (EMA/DEMA), peak tracking,
 *    allocation rate derivatives, memory leak detection via linear regression, and JS Heap integration.
 * 5. Out-of-memory prevention guard for 64^3, 128^3, and 256^3 cosmological vector/scalar grids,
 *    Runge-Kutta RK45 streamline traces, stream surfaces, tube ribbon mesh allocations, and adaptive LOD downscaling.
 * 6. FinalizationRegistry & WeakRef GC Sentinel for detecting orphaned unmanaged WebGL/TypedArray leaks.
 * 7. Multi-Worker Memory Budget Coordinator for partitioning quotas across distributed worker threads.
 * 8. Strict mathematical and memory conservation invariants verification.
 *
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Byte multiplier constants.
 */
export const BYTES_PER_KB = 1024;
export const BYTES_PER_MB = 1024 * 1024;
export const BYTES_PER_GB = 1024 * 1024 * 1024;

/**
 * Standard typed array element byte sizes.
 * @enum {number}
 */
export const DATA_TYPE_BYTES = Object.freeze({
  Int8: 1,
  Uint8: 1,
  Uint8Clamped: 1,
  Int16: 2,
  Uint16: 2,
  Int32: 4,
  Uint32: 4,
  Float32: 4,
  Float64: 8,
  BigInt64: 8,
  BigUint64: 8,
  HalfFloat: 2,
  Byte: 1,
  UnsignedByte: 1,
  Short: 2,
  UnsignedShort: 2,
  Int: 4,
  UnsignedInt: 4,
  Float: 4,
  // Three.js constant names
  1009: 1, // ByteType
  1010: 1, // UnsignedByteType
  1011: 2, // ShortType
  1012: 2, // UnsignedShortType
  1013: 4, // IntType
  1014: 4, // UnsignedIntType
  1015: 4, // FloatType
  1016: 2  // HalfFloatType
});

/**
 * Resource categories tracked in VRAM and Worker memory.
 * @enum {string}
 */
export const RESOURCE_TYPE = Object.freeze({
  GEOMETRY: 'geometry',
  TEXTURE_2D: 'texture_2d',
  TEXTURE_3D: 'texture_3d',
  INSTANCED_MESH: 'instanced_mesh',
  PARTICLE_SYSTEM: 'particle_system',
  RENDER_TARGET: 'render_target',
  WORKER_BUFFER: 'worker_buffer',
  MATERIAL: 'material',
  SHADER_PROGRAM: 'shader_program',
  CUSTOM: 'custom'
});

/**
 * Priority tiers for LRU cache eviction resistance.
 * Lower numerical value = evicted first.
 * Higher numerical value = protected longer.
 * PINNED = strictly never evicted automatically.
 * @enum {number}
 */
export const RESOURCE_PRIORITY = Object.freeze({
  TRANSIENT: 0,   // Scratch buffers, temporary intermediate passes
  LOW: 1,         // Distant LODs, background velocity stream chunks
  MEDIUM: 2,      // Standard interactive visual geometries, active textures
  HIGH: 3,        // Core UI overlays, primary cosmic field meshes, sink nodes
  PINNED: 4       // Master grid tensors, permanent coordinate frames, shaders
});

/**
 * Memory pressure state machine levels.
 * @enum {string}
 */
export const MEMORY_PRESSURE_LEVEL = Object.freeze({
  NORMAL: 'normal',       // Tracked memory < moderate threshold
  MODERATE: 'moderate',   // Moderate threshold <= tracked memory < soft limit
  CRITICAL: 'critical',   // Soft limit <= tracked memory < hard limit
  EXCEEDED: 'exceeded'    // Tracked memory >= hard limit
});

/**
 * WebGL Texture formats and their channel counts.
 * @enum {number}
 */
export const TEXTURE_FORMAT_CHANNELS = Object.freeze({
  AlphaFormat: 1,
  RedFormat: 1,
  RedIntegerFormat: 1,
  LuminanceFormat: 1,
  LuminanceAlphaFormat: 2,
  RGFormat: 2,
  RGIntegerFormat: 2,
  RGBFormat: 3,
  RGBIntegerFormat: 3,
  RGBAFormat: 4,
  RGBAIntegerFormat: 4,
  DepthFormat: 1,
  DepthStencilFormat: 2,
  // Standard string aliases
  R: 1,
  RED: 1,
  RG: 2,
  RGB: 3,
  RGBA: 4,
  DEPTH: 1,
  DEPTH_STENCIL: 2,
  // Three.js constant integer mappings
  1021: 1, // AlphaFormat
  1022: 3, // RGBFormat
  1023: 4, // RGBAFormat
  1024: 1, // LuminanceFormat
  1025: 2, // LuminanceAlphaFormat
  1026: 1, // DepthFormat
  1027: 2, // DepthStencilFormat
  1028: 1, // RedFormat
  1029: 1, // RedIntegerFormat
  1030: 2, // RGFormat
  1031: 2, // RGIntegerFormat
  1032: 3, // RGBIntegerFormat
  1033: 4  // RGBAIntegerFormat
});

/**
 * Compressed texture format bits per pixel multipliers.
 * @type {Readonly<Object>}
 */
export const COMPRESSED_FORMAT_BPP = Object.freeze({
  // S3TC / DXT
  RGB_S3TC_DXT1_Format: 0.5,       // 4 bits per pixel
  RGBA_S3TC_DXT1_Format: 0.5,
  RGBA_S3TC_DXT3_Format: 1.0,      // 8 bits per pixel
  RGBA_S3TC_DXT5_Format: 1.0,
  // ETC
  RGB_ETC1_Format: 0.5,
  RGB_ETC2_Format: 0.5,
  RGBA_ETC2_EAC_Format: 1.0,
  // ASTC
  RGBA_ASTC_4x4_Format: 1.0,
  RGBA_ASTC_8x8_Format: 0.25
});

/**
 * Default configuration parameters for the Memory Budget Manager.
 * @type {Readonly<Object>}
 */
export const DEFAULT_BUDGET_CONFIG = Object.freeze({
  softLimitBytes: 512 * BYTES_PER_MB,         // 512 MB default soft limit
  hardLimitBytes: 1024 * BYTES_PER_MB,        // 1024 MB (1.0 GB) default hard limit
  moderateThresholdRatio: 0.75,               // 75% of soft limit = 384 MB
  criticalThresholdRatio: 1.00,               // 100% of soft limit = 512 MB
  targetEvictionHeadroom: 0.15,               // Evict to 15% below soft limit on pressure
  telemetryPollIntervalMs: 100,               // 10 Hz polling rate
  maxHistorySamples: 120,                     // 120 samples history window (12 seconds at 10 Hz)
  autoEvictOnHardLimit: true,                 // Synchronously trigger LRU eviction when hard limit reached
  enforceAlignment: true,                     // Align byte offsets to 4-byte boundaries
  alignmentBytes: 4,                          // GPU VRAM buffer alignment in bytes
  enableConsoleWarnings: false,               // Emit diagnostic console warnings on pressure
  enableHeapTracking: true,                   // Query performance.memory or process.memoryUsage when available
  enableGCSentinel: true,                     // Enable FinalizationRegistry memory leak sentinel
  maxWorkerQuotaRatio: 0.60                   // Web Workers pool max memory quota (60% of total budget)
});

// ============================================================================
// Custom Error Hierarchy
// ============================================================================

/**
 * Base error class for memory budget manager exceptions.
 */
export class MemoryManagerError extends Error {
  /**
   * @param {string} message
   * @param {Object} [details={}]
   */
  constructor(message, details = {}) {
    super(message);
    this.name = this.constructor.name;
    this.details = details;
    this.timestamp = Date.now();
    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, this.constructor);
    }
  }
}

/**
 * Error thrown when an allocation exceeds hard limits or available budget.
 */
export class MemoryBudgetExceededError extends MemoryManagerError {
  /**
   * @param {number} requestedBytes
   * @param {number} currentTotalBytes
   * @param {number} hardLimitBytes
   * @param {Object} [context={}]
   */
  constructor(requestedBytes, currentTotalBytes, hardLimitBytes, context = {}) {
    const msg = `Memory budget exceeded: Requested ${(requestedBytes / BYTES_PER_MB).toFixed(2)} MB, ` +
      `Current Total: ${(currentTotalBytes / BYTES_PER_MB).toFixed(2)} MB, ` +
      `Hard Limit: ${(hardLimitBytes / BYTES_PER_MB).toFixed(2)} MB`;
    super(msg, { requestedBytes, currentTotalBytes, hardLimitBytes, ...context });
  }
}

/**
 * Error thrown when resource disposal or eviction fails.
 */
export class ResourceDisposalError extends MemoryManagerError {
  /**
   * @param {string} resourceId
   * @param {string} reason
   * @param {Object} [details={}]
   */
  constructor(resourceId, reason, details = {}) {
    super(`Failed to dispose resource '${resourceId}': ${reason}`, { resourceId, reason, ...details });
  }
}

/**
 * Error thrown when an invalid resource is registered or queried.
 */
export class InvalidResourceError extends MemoryManagerError {
  /**
   * @param {string} resourceId
   * @param {string} reason
   */
  constructor(resourceId, reason) {
    super(`Invalid resource specification for '${resourceId}': ${reason}`, { resourceId, reason });
  }
}

/**
 * Error thrown when a stream reconstruction or allocation is refused by OOM guard.
 */
export class AllocationRefusedError extends MemoryManagerError {
  /**
   * @param {string} task
   * @param {number} requiredBytes
   * @param {number} availableBytes
   * @param {Object} [recommendation={}]
   */
  constructor(task, requiredBytes, availableBytes, recommendation = {}) {
    const msg = `Allocation refused for task '${task}': Requires ${(requiredBytes / BYTES_PER_MB).toFixed(2)} MB, ` +
      `only ${(availableBytes / BYTES_PER_MB).toFixed(2)} MB safe headroom available.`;
    super(msg, { task, requiredBytes, availableBytes, recommendation });
  }
}

// ============================================================================
// Resource Entry Data Structure
// ============================================================================

/**
 * Encapsulates a tracked memory resource with its metadata, lifecycle timestamps,
 * priority tier, and dependency graph connections.
 */
export class ResourceEntry {
  /**
   * @param {Object} options
   * @param {string} options.id - Unique resource identifier
   * @param {string} [options.name] - Human-readable name / debug label
   * @param {string} options.type - RESOURCE_TYPE enum value
   * @param {number} options.sizeBytes - Exact memory footprint in bytes
   * @param {number} [options.priority=RESOURCE_PRIORITY.MEDIUM] - Eviction priority
   * @param {*} [options.rawRef=null] - Direct reference or object to manage/dispose
   * @param {Function} [options.disposeCallback=null] - Custom explicit disposal hook
   * @param {Object} [options.metadata={}] - Arbitrary scientific metadata
   */
  constructor({
    id,
    name = '',
    type = RESOURCE_TYPE.CUSTOM,
    sizeBytes = 0,
    priority = RESOURCE_PRIORITY.MEDIUM,
    rawRef = null,
    disposeCallback = null,
    metadata = {}
  }) {
    if (!id || typeof id !== 'string') {
      throw new InvalidResourceError(String(id), 'Resource ID must be a non-empty string.');
    }
    if (typeof sizeBytes !== 'number' || isNaN(sizeBytes) || sizeBytes < 0) {
      throw new InvalidResourceError(id, `sizeBytes must be a non-negative number, got ${sizeBytes}`);
    }

    this.id = id;
    this.name = name || id;
    this.type = type;
    this.sizeBytes = Math.floor(sizeBytes);
    this.priority = priority;
    this.rawRef = rawRef;
    this.disposeCallback = disposeCallback;
    this.metadata = { ...metadata };

    const now = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    this.createdAt = now;
    this.lastAccessedAt = now;
    this.accessCount = 0;
    this.isDisposed = false;
    this.isPinned = (priority === RESOURCE_PRIORITY.PINNED);

    /**
     * IDs of other resources this resource depends upon (outgoing refs).
     * @type {Set<string>}
     */
    this.dependencies = new Set();

    /**
     * IDs of other resources that depend on this resource (incoming refs).
     * @type {Set<string>}
     */
    this.dependents = new Set();
  }

  /**
   * Updates last-accessed timestamp and increments access counter.
   * @returns {ResourceEntry}
   */
  touch() {
    this.lastAccessedAt = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    this.accessCount++;
    return this;
  }

  /**
   * Sets priority and updates pinned status.
   * @param {number} priority
   * @returns {ResourceEntry}
   */
  setPriority(priority) {
    this.priority = priority;
    this.isPinned = (priority === RESOURCE_PRIORITY.PINNED);
    return this;
  }

  /**
   * Pins resource to prevent automatic LRU eviction.
   * @returns {ResourceEntry}
   */
  pin() {
    this.priority = RESOURCE_PRIORITY.PINNED;
    this.isPinned = true;
    return this;
  }

  /**
   * Unpins resource and sets priority to MEDIUM or specified level.
   * @param {number} [fallbackPriority=RESOURCE_PRIORITY.MEDIUM]
   * @returns {ResourceEntry}
   */
  unpin(fallbackPriority = RESOURCE_PRIORITY.MEDIUM) {
    this.priority = fallbackPriority;
    this.isPinned = false;
    return this;
  }

  /**
   * Links a dependency to this resource.
   * @param {string} resourceId
   */
  addDependency(resourceId) {
    this.dependencies.add(resourceId);
  }

  /**
   * Removes a dependency.
   * @param {string} resourceId
   */
  removeDependency(resourceId) {
    this.dependencies.delete(resourceId);
  }

  /**
   * Links a dependent resource.
   * @param {string} resourceId
   */
  addDependent(resourceId) {
    this.dependents.add(resourceId);
  }

  /**
   * Removes a dependent resource.
   * @param {string} resourceId
   */
  removeDependent(resourceId) {
    this.dependents.delete(resourceId);
  }

  /**
   * Serializes metadata and state for telemetry / IPC.
   * @returns {Object}
   */
  toJSON() {
    return {
      id: this.id,
      name: this.name,
      type: this.type,
      sizeBytes: this.sizeBytes,
      sizeMB: +(this.sizeBytes / BYTES_PER_MB).toFixed(4),
      priority: this.priority,
      isPinned: this.isPinned,
      isDisposed: this.isDisposed,
      createdAt: this.createdAt,
      lastAccessedAt: this.lastAccessedAt,
      accessCount: this.accessCount,
      dependenciesCount: this.dependencies.size,
      dependentsCount: this.dependents.size,
      metadata: this.metadata
    };
  }
}

// ============================================================================
// Memory Estimator Utility
// ============================================================================

/**
 * High-precision memory footprint estimator for Three.js objects, typed buffers,
 * cosmological grids, and particle streamline systems.
 */
export class MemoryEstimator {
  /**
   * Aligns byte count to a specified boundary (default 4 bytes).
   * @param {number} bytes
   * @param {number} [alignment=4]
   * @returns {number}
   */
  static alignBytes(bytes, alignment = 4) {
    if (alignment <= 1) return Math.ceil(bytes);
    return Math.ceil(bytes / alignment) * alignment;
  }

  /**
   * Calculates element byte size from TypedArray instance or constructor name.
   * @param {TypedArray|string|Function} arrayOrType
   * @returns {number}
   */
  static getElementByteSize(arrayOrType) {
    if (!arrayOrType) return 4;

    if (arrayOrType.BYTES_PER_ELEMENT) {
      return arrayOrType.BYTES_PER_ELEMENT;
    }

    const typeStr = typeof arrayOrType === 'function' ? arrayOrType.name : String(arrayOrType);

    for (const [key, size] of Object.entries(DATA_TYPE_BYTES)) {
      if (typeStr.includes(key)) {
        return size;
      }
    }

    return 4; // Default float32 / int32
  }

  /**
   * Estimates memory for a single BufferAttribute or InterleavedBufferAttribute.
   * @param {Object} attr - Three.js BufferAttribute
   * @returns {number}
   */
  static estimateBufferAttributeBytes(attr) {
    if (!attr) return 0;

    // Direct typed array allocation
    if (attr.array && attr.array.byteLength !== undefined) {
      return attr.array.byteLength;
    }

    // Interleaved buffer
    if (attr.data && attr.data.array && attr.data.array.byteLength !== undefined) {
      return attr.data.array.byteLength;
    }

    // Array or count * itemSize fallback
    const count = attr.count || (attr.array ? attr.array.length : 0);
    const itemSize = attr.itemSize || 1;
    const bytesPerElem = attr.array ? (attr.array.BYTES_PER_ELEMENT || 4) : 4;

    return count * itemSize * bytesPerElem;
  }

  /**
   * Estimates memory for a Three.js BufferGeometry including all attributes and index.
   * @param {Object} geometry - Three.js BufferGeometry
   * @returns {number}
   */
  static estimateBufferGeometryBytes(geometry) {
    if (!geometry) return 0;

    let totalBytes = 0;

    // Attributes (position, normal, uv, color, velocities, tensors, etc.)
    if (geometry.attributes) {
      for (const key of Object.keys(geometry.attributes)) {
        const attr = geometry.attributes[key];
        totalBytes += MemoryEstimator.estimateBufferAttributeBytes(attr);
      }
    }

    // Index buffer
    if (geometry.index) {
      totalBytes += MemoryEstimator.estimateBufferAttributeBytes(geometry.index);
    }

    // Morph attributes
    if (geometry.morphAttributes) {
      for (const targetKey of Object.keys(geometry.morphAttributes)) {
        const targetList = geometry.morphAttributes[targetKey];
        if (Array.isArray(targetList)) {
          for (const morphAttr of targetList) {
            totalBytes += MemoryEstimator.estimateBufferAttributeBytes(morphAttr);
          }
        }
      }
    }

    return totalBytes;
  }

  /**
   * Estimates memory for a 2D Texture.
   * @param {Object} texture - Three.js Texture or texture descriptor
   * @returns {number}
   */
  static estimateTextureBytes(texture) {
    if (!texture) return 0;

    const image = texture.image;
    let width = 0;
    let height = 0;

    if (image) {
      width = image.width || 0;
      height = image.height || 0;
      if (image.data && image.data.byteLength) {
        return image.data.byteLength;
      }
    } else {
      width = texture.width || 0;
      height = texture.height || 0;
    }

    if (width === 0 || height === 0) return 0;

    // Compressed texture check
    if (texture.format && COMPRESSED_FORMAT_BPP[texture.format]) {
      const bpp = COMPRESSED_FORMAT_BPP[texture.format];
      let bytes = Math.ceil(width * height * bpp);
      if (texture.generateMipmaps || (texture.mipmaps && texture.mipmaps.length > 0)) {
        bytes = Math.floor(bytes * 1.3333333333333333);
      }
      return bytes;
    }

    // Determine channels
    let channels = 4;
    if (texture.format) {
      channels = TEXTURE_FORMAT_CHANNELS[texture.format] || 4;
    }

    // Determine bytes per channel
    let bytesPerChannel = 1;
    if (texture.type) {
      const typeVal = texture.type;
      if (DATA_TYPE_BYTES[typeVal] !== undefined) {
        bytesPerChannel = DATA_TYPE_BYTES[typeVal];
      } else {
        const typeStr = String(typeVal);
        if (typeStr.includes('FloatType') && !typeStr.includes('Half')) {
          bytesPerChannel = 4;
        } else if (typeStr.includes('HalfFloatType')) {
          bytesPerChannel = 2;
        } else if (typeStr.includes('IntType') || typeStr.includes('UnsignedIntType')) {
          bytesPerChannel = 4;
        } else if (typeStr.includes('ShortType')) {
          bytesPerChannel = 2;
        }
      }
    }

    let baseBytes = width * height * channels * bytesPerChannel;

    // Mipmaps calculation (geometric series sum = 1 + 1/4 + 1/16 + ... ~= 1.3333)
    if (texture.generateMipmaps || (texture.mipmaps && texture.mipmaps.length > 0)) {
      baseBytes = Math.floor(baseBytes * 1.3333333333333333);
    }

    // Cube texture multiplier
    if (texture.isCubeTexture || (texture.images && texture.images.length === 6)) {
      baseBytes *= 6;
    }

    return baseBytes;
  }

  /**
   * Estimates memory for a 3D Data Texture (e.g. 64x64x64 or 128x128x128 cosmological field).
   * @param {Object} texture3D - Three.js Data3DTexture or 3D grid descriptor
   * @returns {number}
   */
  static estimateData3DTextureBytes(texture3D) {
    if (!texture3D) return 0;

    const image = texture3D.image;
    let width = 0;
    let height = 0;
    let depth = 0;

    if (image) {
      width = image.width || 0;
      height = image.height || 0;
      depth = image.depth || 0;
      if (image.data && image.data.byteLength) {
        return image.data.byteLength;
      }
    } else {
      width = texture3D.width || 0;
      height = texture3D.height || 0;
      depth = texture3D.depth || 0;
    }

    if (width === 0 || height === 0 || depth === 0) return 0;

    let channels = 1;
    if (texture3D.format) {
      channels = TEXTURE_FORMAT_CHANNELS[texture3D.format] || 1;
    }

    let bytesPerChannel = 4; // Data3D is typically Float32 for scientific fields
    if (texture3D.type) {
      const typeVal = texture3D.type;
      if (DATA_TYPE_BYTES[typeVal] !== undefined) {
        bytesPerChannel = DATA_TYPE_BYTES[typeVal];
      } else {
        const typeStr = String(typeVal);
        if (typeStr.includes('HalfFloatType')) {
          bytesPerChannel = 2;
        } else if (typeStr.includes('UnsignedByteType') || typeStr.includes('ByteType')) {
          bytesPerChannel = 1;
        }
      }
    }

    let totalBytes = width * height * depth * channels * bytesPerChannel;

    if (texture3D.generateMipmaps) {
      // 3D mipmaps sum = 1 + 1/8 + 1/64 + ... ~= 1.142857
      totalBytes = Math.floor(totalBytes * 1.1428571428571428);
    }

    return totalBytes;
  }

  /**
   * Estimates memory for an InstancedMesh (geometry + instance matrices + instance colors + instance attributes).
   * @param {Object} instancedMesh - Three.js InstancedMesh
   * @returns {number}
   */
  static estimateInstancedMeshBytes(instancedMesh) {
    if (!instancedMesh) return 0;

    let totalBytes = 0;

    // Base geometry
    if (instancedMesh.geometry) {
      totalBytes += MemoryEstimator.estimateBufferGeometryBytes(instancedMesh.geometry);
    }

    // Instance matrix attribute: count * 16 floats * 4 bytes = count * 64 bytes
    const count = instancedMesh.count || 0;
    if (instancedMesh.instanceMatrix) {
      totalBytes += MemoryEstimator.estimateBufferAttributeBytes(instancedMesh.instanceMatrix);
    } else if (count > 0) {
      totalBytes += count * 16 * 4;
    }

    // Instance color attribute: count * 3 floats * 4 bytes = count * 12 bytes
    if (instancedMesh.instanceColor) {
      totalBytes += MemoryEstimator.estimateBufferAttributeBytes(instancedMesh.instanceColor);
    }

    // Custom instance attributes
    if (instancedMesh.geometry && instancedMesh.geometry.attributes) {
      for (const key of Object.keys(instancedMesh.geometry.attributes)) {
        const attr = instancedMesh.geometry.attributes[key];
        if (attr && (attr.isInstancedBufferAttribute || attr.isInstanced)) {
          totalBytes += MemoryEstimator.estimateBufferAttributeBytes(attr);
        }
      }
    }

    return totalBytes;
  }

  /**
   * Estimates memory for a WebGLRenderTarget (color textures + depth texture + stencil).
   * @param {Object} renderTarget - Three.js WebGLRenderTarget
   * @returns {number}
   */
  static estimateRenderTargetBytes(renderTarget) {
    if (!renderTarget) return 0;

    const width = renderTarget.width || 0;
    const height = renderTarget.height || 0;
    if (width === 0 || height === 0) return 0;

    let totalBytes = 0;

    // Color texture(s)
    if (Array.isArray(renderTarget.textures)) {
      for (const tex of renderTarget.textures) {
        totalBytes += MemoryEstimator.estimateTextureBytes({ ...tex, width, height });
      }
    } else if (renderTarget.texture) {
      totalBytes += MemoryEstimator.estimateTextureBytes({ ...renderTarget.texture, width, height });
    } else {
      // Default RGBA UnsignedByte
      totalBytes += width * height * 4 * 1;
    }

    // Multi-sample renderbuffers (MSAA)
    const samples = renderTarget.samples || 0;
    if (samples > 0) {
      totalBytes *= samples;
    }

    // Depth buffer
    if (renderTarget.depthBuffer) {
      // 24-bit depth = 3 bytes or 32-bit depth = 4 bytes
      totalBytes += width * height * 4;
    }

    // Stencil buffer
    if (renderTarget.stencilBuffer) {
      totalBytes += width * height * 1;
    }

    return totalBytes;
  }

  /**
   * Estimates memory for a 3D scalar or vector cosmological grid.
   * @param {Object} params
   * @param {number} params.resolution - Grid resolution (e.g. 64, 128, 256)
   * @param {number} [params.channels=1] - 1 for scalar density, 3 for velocity vector, 6 for strain tensor
   * @param {number} [params.bytesPerElement=4] - 4 for Float32Array, 8 for Float64Array
   * @param {boolean} [params.includeDerivatives=false] - If true, adds gradient and divergence tensors
   * @returns {number}
   */
  static estimateGridBytes({ resolution, channels = 1, bytesPerElement = 4, includeDerivatives = false }) {
    if (typeof resolution !== 'number' || resolution <= 0) {
      throw new Error(`Invalid grid resolution: ${resolution}`);
    }

    const cellCount = resolution * resolution * resolution;
    let totalChannels = channels;

    if (includeDerivatives) {
      // For scalar: + 3 channels (gradient); For vector: + 1 (div) + 3 (curl) + 6 (shear tensor) = + 10
      totalChannels += (channels === 1 ? 3 : 10);
    }

    return cellCount * totalChannels * bytesPerElement;
  }

  /**
   * Exact analytical memory calculation for Runge-Kutta RK45 streamline reconstruction.
   * @param {Object} params
   * @param {number} params.seedCount - Number of streamline seed trajectories
   * @param {number} [params.avgStepsPerStreamline=128] - Average integration steps
   * @param {number} [params.attributesPerVertex=10] - Position (3), Velocity (3), Density (1), Curvature (1), Time (1), Color (1)
   * @param {boolean} [params.renderAsTubes=true] - If true, generates radial tube mesh vertices
   * @param {number} [params.tubeRadialSegments=6] - Number of radial segments per tube ring
   * @returns {Object} Detailed byte breakdown and total required bytes
   */
  static estimateStreamlineReconstructionBytes({
    seedCount,
    avgStepsPerStreamline = 128,
    attributesPerVertex = 10,
    renderAsTubes = true,
    tubeRadialSegments = 6
  }) {
    if (typeof seedCount !== 'number' || seedCount <= 0) {
      throw new Error(`Invalid seedCount: ${seedCount}`);
    }

    const totalRawPoints = seedCount * avgStepsPerStreamline;

    // 1. Raw trace trajectory buffer (Position + Attributes): Float32Array
    const rawTraceBytes = totalRawPoints * attributesPerVertex * 4;

    // 2. Render geometry buffers
    let renderGeometryBytes = 0;
    let triangleCount = 0;

    if (renderAsTubes) {
      // Tube mesh: each step segment produces (tubeRadialSegments * 2) triangles
      // Vertices: (totalRawPoints) * (tubeRadialSegments + 1)
      const tubeVertices = totalRawPoints * (tubeRadialSegments + 1);
      // Vertex attributes: Position(3), Normal(3), Color(3), UV(2), Tangent(4) = 15 floats * 4 bytes = 60 bytes/vertex
      const vertexBytes = tubeVertices * 15 * 4;
      // Indices: (avgStepsPerStreamline - 1) * tubeRadialSegments * 6 indices * 4 bytes (Uint32Array)
      const indexBytes = seedCount * (avgStepsPerStreamline - 1) * tubeRadialSegments * 6 * 4;
      renderGeometryBytes = vertexBytes + indexBytes;
      triangleCount = seedCount * (avgStepsPerStreamline - 1) * tubeRadialSegments * 2;
    } else {
      // Line geometry: totalRawPoints * Position(3) + Color(3) = 6 floats * 4 bytes = 24 bytes/vertex
      renderGeometryBytes = totalRawPoints * 6 * 4;
      triangleCount = 0;
    }

    // 3. Worker transfer payload overhead (serialization and transfer header)
    const workerTransferOverheadBytes = Math.floor(rawTraceBytes * 0.05);

    const totalBytes = rawTraceBytes + renderGeometryBytes + workerTransferOverheadBytes;

    return {
      seedCount,
      avgStepsPerStreamline,
      totalRawPoints,
      rawTraceBytes,
      renderGeometryBytes,
      workerTransferOverheadBytes,
      triangleCount,
      totalBytes,
      totalMB: +(totalBytes / BYTES_PER_MB).toFixed(3)
    };
  }

  /**
   * Estimates memory for stream surfaces (2D sheet tracing through 3D vector fields).
   * @param {Object} params
   * @param {number} params.seedCurvePoints - Points along generating rake curve
   * @param {number} [params.integrationSteps=100] - Steps along flow direction
   * @param {number} [params.attributesPerVertex=12] - Position(3), Normal(3), Velocity(3), Scalar(3)
   * @returns {Object}
   */
  static estimateStreamSurfaceBytes({
    seedCurvePoints,
    integrationSteps = 100,
    attributesPerVertex = 12
  }) {
    const totalVertices = seedCurvePoints * integrationSteps;
    const vertexBytes = totalVertices * attributesPerVertex * 4;
    // Triangles: (seedCurvePoints - 1) * (integrationSteps - 1) * 2 * 3 indices * 4 bytes
    const indexCount = (seedCurvePoints - 1) * (integrationSteps - 1) * 6;
    const indexBytes = indexCount * 4;
    const totalBytes = vertexBytes + indexBytes;

    return {
      seedCurvePoints,
      integrationSteps,
      totalVertices,
      indexCount,
      triangleCount: indexCount / 3,
      vertexBytes,
      indexBytes,
      totalBytes,
      totalMB: +(totalBytes / BYTES_PER_MB).toFixed(3)
    };
  }
}

// ============================================================================
// Out-of-Memory (OOM) Prevention Guard
// ============================================================================

/**
 * Diagnostic and proactive OOM prevention guard for high-resolution 3D cosmological
 * field volumes (64^3, 128^3, 256^3) and massive RK45 streamline reconstruction tasks.
 */
export class StreamlineOOMGuard {
  /**
   * Assesses safety of a proposed streamline reconstruction against current VRAM & memory headroom.
   * @param {Object} config - Streamline reconstruction parameters
   * @param {number} availableBytes - Currently available safe budget headroom
   * @param {number} softLimitBytes - Configured soft memory limit
   * @param {number} hardLimitBytes - Configured hard memory limit
   * @returns {Object} Safety assessment report
   */
  static assessReconstructionSafety(config, availableBytes, softLimitBytes, hardLimitBytes) {
    const estimate = MemoryEstimator.estimateStreamlineReconstructionBytes(config);
    const requiredBytes = estimate.totalBytes;

    const isSafe = requiredBytes <= availableBytes;
    const wouldExceedSoftLimit = requiredBytes > availableBytes;
    const wouldExceedHardLimit = requiredBytes > (hardLimitBytes * 0.95);

    let recommendedConfig = null;
    if (!isSafe) {
      recommendedConfig = StreamlineOOMGuard.recommendAdaptiveDownscaling(config, availableBytes);
    }

    return {
      isSafe,
      wouldExceedSoftLimit,
      wouldExceedHardLimit,
      requiredBytes,
      requiredMB: estimate.totalMB,
      availableBytes,
      availableMB: +(availableBytes / BYTES_PER_MB).toFixed(3),
      estimate,
      recommendedConfig,
      diagnosticMessage: isSafe
        ? `Reconstruction safe: ${(requiredBytes / BYTES_PER_MB).toFixed(2)} MB fits in ${(availableBytes / BYTES_PER_MB).toFixed(2)} MB headroom.`
        : `OOM Guard Warning: Reconstruction requires ${(requiredBytes / BYTES_PER_MB).toFixed(2)} MB, exceeding safe headroom of ${(availableBytes / BYTES_PER_MB).toFixed(2)} MB.`
    };
  }

  /**
   * Recommends downscaled parameters (seed reduction, tube -> line downgrade, step truncation)
   * to strictly fit within available memory budget.
   * @param {Object} originalConfig
   * @param {number} targetBytes
   * @returns {Object} Optimized parameter configuration
   */
  static recommendAdaptiveDownscaling(originalConfig, targetBytes) {
    const safeTarget = Math.max(1024 * 1024, targetBytes * 0.85); // 15% safety buffer

    // Step 1: Try disabling tubes and rendering as lines
    let candidate = { ...originalConfig, renderAsTubes: false };
    let est = MemoryEstimator.estimateStreamlineReconstructionBytes(candidate);
    if (est.totalBytes <= safeTarget) {
      return {
        ...candidate,
        strategy: 'DISABLE_TUBES_USE_LINES',
        projectedMB: est.totalMB,
        reductionRatio: +(est.totalBytes / MemoryEstimator.estimateStreamlineReconstructionBytes(originalConfig).totalBytes).toFixed(3)
      };
    }

    // Step 2: Try reducing radial segments to 3
    candidate = { ...originalConfig, renderAsTubes: true, tubeRadialSegments: 3 };
    est = MemoryEstimator.estimateStreamlineReconstructionBytes(candidate);
    if (est.totalBytes <= safeTarget) {
      return {
        ...candidate,
        strategy: 'REDUCE_RADIAL_SEGMENTS_TO_3',
        projectedMB: est.totalMB,
        reductionRatio: +(est.totalBytes / MemoryEstimator.estimateStreamlineReconstructionBytes(originalConfig).totalBytes).toFixed(3)
      };
    }

    // Step 3: Decimate seed count to fit target
    const currentSeeds = originalConfig.seedCount || 1000;
    const estLine = MemoryEstimator.estimateStreamlineReconstructionBytes({ ...originalConfig, renderAsTubes: false });
    const bytesPerSeedLine = estLine.totalBytes / currentSeeds;

    const maxFeasibleSeeds = Math.max(10, Math.floor(safeTarget / bytesPerSeedLine));

    return {
      ...originalConfig,
      seedCount: maxFeasibleSeeds,
      renderAsTubes: false,
      strategy: 'DECIMATE_SEEDS_AND_USE_LINES',
      projectedMB: +(MemoryEstimator.estimateStreamlineReconstructionBytes({
        ...originalConfig,
        seedCount: maxFeasibleSeeds,
        renderAsTubes: false
      }).totalBytes / BYTES_PER_MB).toFixed(3),
      reductionRatio: +(maxFeasibleSeeds / currentSeeds).toFixed(3)
    };
  }

  /**
   * Assesses safety for 3D grid allocation (e.g. 64^3 vs 128^3 vs 256^3 field data textures).
   * @param {number} resolution - Grid dimension (e.g. 64, 128, 256)
   * @param {number} channels - Number of channels (1 for scalar, 3 for velocity vector)
   * @param {number} availableBytes - Available headroom
   * @returns {Object}
   */
  static assessGridSafety(resolution, channels, availableBytes) {
    const requiredBytes = MemoryEstimator.estimateGridBytes({ resolution, channels, bytesPerElement: 4 });
    const isSafe = requiredBytes <= availableBytes;

    let fallbackResolution = resolution;
    if (!isSafe) {
      if (resolution >= 256 && MemoryEstimator.estimateGridBytes({ resolution: 128, channels, bytesPerElement: 4 }) <= availableBytes) {
        fallbackResolution = 128;
      } else if (resolution >= 128 && MemoryEstimator.estimateGridBytes({ resolution: 64, channels, bytesPerElement: 4 }) <= availableBytes) {
        fallbackResolution = 64;
      } else {
        fallbackResolution = 32;
      }
    }

    return {
      resolution,
      channels,
      requiredBytes,
      requiredMB: +(requiredBytes / BYTES_PER_MB).toFixed(3),
      availableBytes,
      isSafe,
      fallbackResolution,
      diagnosticMessage: isSafe
        ? `Grid ${resolution}^3 (${channels}ch) safe: ${(requiredBytes / BYTES_PER_MB).toFixed(2)} MB.`
        : `Grid ${resolution}^3 (${channels}ch) exceeds headroom! Recommend fallback to ${fallbackResolution}^3.`
    };
  }
}

// ============================================================================
// Garbage Collection & Finalization Sentinel
// ============================================================================

/**
 * Weak reference and FinalizationRegistry sentinel to detect orphaned objects
 * and notify the budget manager when JS GC collects objects that were not explicitly disposed.
 */
export class GarbageCollectionSentinel {
  /**
   * @param {Function} [onOrphanCollected=null]
   */
  constructor(onOrphanCollected = null) {
    this._onOrphanCollected = onOrphanCollected;
    this._trackedTargets = new Map();
    this._collectedCount = 0;

    if (typeof FinalizationRegistry !== 'undefined') {
      this._registry = new FinalizationRegistry((heldValue) => {
        this._handleFinalization(heldValue);
      });
    } else {
      this._registry = null;
    }
  }

  /**
   * Registers a target object for GC finalization monitoring.
   * @param {Object} target
   * @param {Object} tokenInfo
   */
  track(target, tokenInfo) {
    if (!this._registry || !target || typeof target !== 'object') return;
    const id = tokenInfo.id || String(Math.random());
    this._trackedTargets.set(id, tokenInfo);
    try {
      this._registry.register(target, id, target);
    } catch {
      // Ignore primitive registration errors
    }
  }

  /**
   * Unregisters a target when explicitly disposed.
   * @param {Object} target
   * @param {string} id
   */
  untrack(target, id) {
    if (!this._registry || !target) return;
    this._trackedTargets.delete(id);
    try {
      this._registry.unregister(target);
    } catch {
      // Ignore unregistration errors
    }
  }

  /**
   * Handles GC cleanup event for orphaned object.
   * @param {string} id
   * @private
   */
  _handleFinalization(id) {
    const tokenInfo = this._trackedTargets.get(id);
    this._trackedTargets.delete(id);
    this._collectedCount++;
    if (this._onOrphanCollected && tokenInfo) {
      this._onOrphanCollected(tokenInfo);
    }
  }

  /**
   * Returns count of orphaned resources collected by engine.
   * @returns {number}
   */
  get collectedCount() {
    return this._collectedCount;
  }
}

// ============================================================================
// Multi-Worker Memory Budget Coordinator
// ============================================================================

/**
 * Manages partitioned memory allowances and monitors transfer buffers across
 * distributed Web Worker pool threads.
 */
export class WorkerMemoryCoordinator {
  /**
   * @param {MemoryBudgetManager} budgetManager
   * @param {Object} [options={}]
   */
  constructor(budgetManager, options = {}) {
    this.manager = budgetManager;
    this.maxTotalWorkerBytes = options.maxTotalWorkerBytes || (this.manager.softLimitBytes * 0.60);
    this._workers = new Map(); // workerId -> { quotaBytes, activeBytes, taskCount }
  }

  /**
   * Registers a worker thread and assigns memory quota.
   * @param {string} workerId
   * @param {number} [quotaBytes=0]
   */
  registerWorker(workerId, quotaBytes = 0) {
    const quota = quotaBytes > 0 ? quotaBytes : Math.floor(this.maxTotalWorkerBytes / 4);
    this._workers.set(workerId, {
      workerId,
      quotaBytes: quota,
      activeBytes: 0,
      completedTasks: 0,
      totalTransferredBytes: 0
    });
  }

  /**
   * Unregisters a worker thread.
   * @param {string} workerId
   */
  unregisterWorker(workerId) {
    this._workers.delete(workerId);
  }

  /**
   * Verifies if worker can accept a zero-copy buffer of given size.
   * @param {string} workerId
   * @param {number} sizeBytes
   * @returns {boolean}
   */
  canWorkerAccept(workerId, sizeBytes) {
    const w = this._workers.get(workerId);
    if (!w) return false;
    return (w.activeBytes + sizeBytes) <= w.quotaBytes;
  }

  /**
   * Records buffer transfer to worker.
   * @param {string} workerId
   * @param {number} sizeBytes
   */
  recordDispatch(workerId, sizeBytes) {
    const w = this._workers.get(workerId);
    if (w) {
      w.activeBytes += sizeBytes;
      w.totalTransferredBytes += sizeBytes;
    }
  }

  /**
   * Records buffer return or task completion from worker.
   * @param {string} workerId
   * @param {number} sizeBytes
   */
  recordCompletion(workerId, sizeBytes) {
    const w = this._workers.get(workerId);
    if (w) {
      w.activeBytes = Math.max(0, w.activeBytes - sizeBytes);
      w.completedTasks++;
    }
  }

  /**
   * Returns worker pool telemetry metrics.
   * @returns {Object}
   */
  getWorkerTelemetry() {
    const workerStats = {};
    let totalActiveWorkerBytes = 0;

    for (const [id, w] of this._workers.entries()) {
      totalActiveWorkerBytes += w.activeBytes;
      workerStats[id] = {
        ...w,
        activeMB: +(w.activeBytes / BYTES_PER_MB).toFixed(2),
        quotaMB: +(w.quotaBytes / BYTES_PER_MB).toFixed(2)
      };
    }

    return {
      workerCount: this._workers.size,
      totalActiveWorkerBytes,
      totalActiveWorkerMB: +(totalActiveWorkerBytes / BYTES_PER_MB).toFixed(2),
      maxTotalWorkerBytes: this.maxTotalWorkerBytes,
      workerStats
    };
  }
}

// ============================================================================
// Core Memory Budget Manager & Sentinel Implementation
// ============================================================================

/**
 * Main WebGL and Web Worker Memory Budget Manager, LRU Eviction Sentinel,
 * and Telemetry Engine.
 */
export class MemoryBudgetManager {
  /**
   * @param {Object} [config={}] - Configuration options
   */
  constructor(config = {}) {
    this.config = {
      ...DEFAULT_BUDGET_CONFIG,
      ...config
    };

    // Validate limit invariants
    if (this.config.softLimitBytes <= 0) {
      throw new Error(`Invalid softLimitBytes: ${this.config.softLimitBytes}. Must be > 0.`);
    }
    if (this.config.hardLimitBytes < this.config.softLimitBytes) {
      throw new Error(
        `hardLimitBytes (${this.config.hardLimitBytes}) cannot be less than softLimitBytes (${this.config.softLimitBytes}).`
      );
    }

    /**
     * Primary registry of all tracked active resources (Key: ID, Value: ResourceEntry).
     * @type {Map<string, ResourceEntry>}
     * @private
     */
    this._resources = new Map();

    /**
     * Category byte counters for fast O(1) breakdown reporting.
     * @type {Map<string, number>}
     * @private
     */
    this._categoryBytes = new Map();
    for (const type of Object.values(RESOURCE_TYPE)) {
      this._categoryBytes.set(type, 0);
    }

    /**
     * Total tracked VRAM / Worker memory bytes.
     * @type {number}
     * @private
     */
    this._totalTrackedBytes = 0;

    /**
     * Peak tracked memory watermark in bytes.
     * @type {number}
     * @private
     */
    this._peakTrackedBytes = 0;

    /**
     * Cumulative total bytes allocated since startup.
     * @type {number}
     * @private
     */
    this._cumulativeAllocatedBytes = 0;

    /**
     * Cumulative total bytes freed / disposed since startup.
     * @type {number}
     * @private
     */
    this._cumulativeFreedBytes = 0;

    /**
     * Cumulative count of resources evicted via LRU.
     * @type {number}
     * @private
     */
    this._totalEvictionCount = 0;

    /**
     * Current memory pressure level state.
     * @type {string}
     * @private
     */
    this._pressureLevel = MEMORY_PRESSURE_LEVEL.NORMAL;

    /**
     * Active Worker transferable buffers registry.
     * @type {Map<string, Object>}
     * @private
     */
    this._workerBuffers = new Map();

    /**
     * Event listeners map.
     * @type {Map<string, Set<Function>>}
     * @private
     */
    this._eventListeners = new Map();

    /**
     * Telemetry history ring buffer.
     * @type {Array<Object>}
     * @private
     */
    this._telemetryHistory = [];

    /**
     * High-frequency polling timer handle.
     * @type {number|null}
     * @private
     */
    this._pollTimer = null;

    /**
     * Timestamp of last telemetry poll.
     * @type {number}
     * @private
     */
    this._lastPollTimestamp = 0;

    /**
     * Running moving average of allocation rate (bytes/sec).
     * @type {number}
     * @private
     */
    this._smoothedAllocationRate = 0;

    /**
     * Flag indicating if manager has been destroyed.
     * @type {boolean}
     * @private
     */
    this._isDestroyed = false;

    /**
     * Garbage collection sentinel.
     * @type {GarbageCollectionSentinel}
     */
    this.gcSentinel = new GarbageCollectionSentinel((tokenInfo) => {
      this._emit('gcOrphanCollected', tokenInfo);
      if (this.hasResource(tokenInfo.id)) {
        this.unregisterResource(tokenInfo.id, false);
      }
    });

    /**
     * Multi-worker coordinator.
     * @type {WorkerMemoryCoordinator}
     */
    this.workerCoordinator = new WorkerMemoryCoordinator(this, {
      maxTotalWorkerBytes: this.config.softLimitBytes * this.config.maxWorkerQuotaRatio
    });
  }

  // ==========================================================================
  // Getters & Invariant State
  // ==========================================================================

  /**
   * Gets total tracked bytes in memory.
   * @returns {number}
   */
  get totalTrackedBytes() {
    return this._totalTrackedBytes;
  }

  /**
   * Gets peak memory watermark.
   * @returns {number}
   */
  get peakTrackedBytes() {
    return this._peakTrackedBytes;
  }

  /**
   * Gets current memory pressure level.
   * @returns {string}
   */
  get pressureLevel() {
    return this._pressureLevel;
  }

  /**
   * Gets total number of tracked resources.
   * @returns {number}
   */
  get resourceCount() {
    return this._resources.size;
  }

  /**
   * Gets soft memory limit in bytes.
   * @returns {number}
   */
  get softLimitBytes() {
    return this.config.softLimitBytes;
  }

  /**
   * Gets hard memory limit in bytes.
   * @returns {number}
   */
  get hardLimitBytes() {
    return this.config.hardLimitBytes;
  }

  /**
   * Returns safe available memory before hitting soft limit.
   * @returns {number}
   */
  get availableSoftHeadroomBytes() {
    return Math.max(0, this.config.softLimitBytes - this._totalTrackedBytes);
  }

  /**
   * Returns safe available memory before hitting hard limit.
   * @returns {number}
   */
  get availableHardHeadroomBytes() {
    return Math.max(0, this.config.hardLimitBytes - this._totalTrackedBytes);
  }

  // ==========================================================================
  // Event Subscription System
  // ==========================================================================

  /**
   * Subscribes a listener to a memory event.
   * Events: 'allocation', 'disposal', 'eviction', 'pressureChange', 'oomWarning', 'telemetryTick', 'gcOrphanCollected'
   * @param {string} event
   * @param {Function} callback
   * @returns {MemoryBudgetManager}
   */
  on(event, callback) {
    if (!this._eventListeners.has(event)) {
      this._eventListeners.set(event, new Set());
    }
    this._eventListeners.get(event).add(callback);
    return this;
  }

  /**
   * Unsubscribes a listener from a memory event.
   * @param {string} event
   * @param {Function} callback
   * @returns {MemoryBudgetManager}
   */
  off(event, callback) {
    if (this._eventListeners.has(event)) {
      this._eventListeners.get(event).delete(callback);
    }
    return this;
  }

  /**
   * Dispatches an event to all registered listeners.
   * @param {string} event
   * @param {*} data
   * @private
   */
  _emit(event, data) {
    if (this._isDestroyed) return;
    const listeners = this._eventListeners.get(event);
    if (listeners) {
      for (const listener of listeners) {
        try {
          listener(data);
        } catch (err) {
          if (this.config.enableConsoleWarnings && typeof console !== 'undefined') {
            console.warn(`[MemoryBudgetManager] Listener error on event '${event}':`, err);
          }
        }
      }
    }
  }

  // ==========================================================================
  // Resource Registration & Allocation
  // ==========================================================================

  /**
   * Registers a Three.js BufferGeometry or raw vertex buffer.
   * @param {string} id - Unique identifier
   * @param {Object} geometry - Three.js BufferGeometry or descriptor
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  registerGeometry(id, geometry, options = {}) {
    const sizeBytes = options.sizeBytes !== undefined
      ? options.sizeBytes
      : MemoryEstimator.estimateBufferGeometryBytes(geometry);

    return this.registerResource({
      id,
      name: options.name || (geometry && geometry.name) || id,
      type: RESOURCE_TYPE.GEOMETRY,
      sizeBytes,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.MEDIUM,
      rawRef: geometry,
      disposeCallback: options.disposeCallback || (() => {
        if (geometry && typeof geometry.dispose === 'function') {
          geometry.dispose();
        }
      }),
      metadata: options.metadata || {}
    });
  }

  /**
   * Registers a 2D WebGL Texture.
   * @param {string} id
   * @param {Object} texture - Three.js Texture or image descriptor
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  registerTexture(id, texture, options = {}) {
    const sizeBytes = options.sizeBytes !== undefined
      ? options.sizeBytes
      : MemoryEstimator.estimateTextureBytes(texture);

    return this.registerResource({
      id,
      name: options.name || (texture && texture.name) || id,
      type: RESOURCE_TYPE.TEXTURE_2D,
      sizeBytes,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.MEDIUM,
      rawRef: texture,
      disposeCallback: options.disposeCallback || (() => {
        if (texture && typeof texture.dispose === 'function') {
          texture.dispose();
        }
      }),
      metadata: options.metadata || {}
    });
  }

  /**
   * Registers a 3D Data Texture (e.g. cosmological field grid).
   * @param {string} id
   * @param {Object} texture3D - Three.js Data3DTexture or 3D grid
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  registerData3DTexture(id, texture3D, options = {}) {
    const sizeBytes = options.sizeBytes !== undefined
      ? options.sizeBytes
      : MemoryEstimator.estimateData3DTextureBytes(texture3D);

    return this.registerResource({
      id,
      name: options.name || (texture3D && texture3D.name) || id,
      type: RESOURCE_TYPE.TEXTURE_3D,
      sizeBytes,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.HIGH,
      rawRef: texture3D,
      disposeCallback: options.disposeCallback || (() => {
        if (texture3D && typeof texture3D.dispose === 'function') {
          texture3D.dispose();
        }
      }),
      metadata: options.metadata || {}
    });
  }

  /**
   * Registers an InstancedMesh.
   * @param {string} id
   * @param {Object} instancedMesh - Three.js InstancedMesh
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  registerInstancedMesh(id, instancedMesh, options = {}) {
    const sizeBytes = options.sizeBytes !== undefined
      ? options.sizeBytes
      : MemoryEstimator.estimateInstancedMeshBytes(instancedMesh);

    return this.registerResource({
      id,
      name: options.name || (instancedMesh && instancedMesh.name) || id,
      type: RESOURCE_TYPE.INSTANCED_MESH,
      sizeBytes,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.MEDIUM,
      rawRef: instancedMesh,
      disposeCallback: options.disposeCallback || (() => {
        if (instancedMesh) {
          if (instancedMesh.geometry && typeof instancedMesh.geometry.dispose === 'function') {
            instancedMesh.geometry.dispose();
          }
          if (instancedMesh.material && typeof instancedMesh.material.dispose === 'function') {
            instancedMesh.material.dispose();
          }
          if (typeof instancedMesh.dispose === 'function') {
            instancedMesh.dispose();
          }
        }
      }),
      metadata: options.metadata || {}
    });
  }

  /**
   * Registers a Particle System / Point Cloud.
   * @param {string} id
   * @param {Object} particleSystem - Three.js Points or custom particle descriptor
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  registerParticleSystem(id, particleSystem, options = {}) {
    let sizeBytes = 0;
    if (options.sizeBytes !== undefined) {
      sizeBytes = options.sizeBytes;
    } else if (particleSystem && particleSystem.geometry) {
      sizeBytes = MemoryEstimator.estimateBufferGeometryBytes(particleSystem.geometry);
    } else if (particleSystem && particleSystem.byteLength) {
      sizeBytes = particleSystem.byteLength;
    }

    return this.registerResource({
      id,
      name: options.name || (particleSystem && particleSystem.name) || id,
      type: RESOURCE_TYPE.PARTICLE_SYSTEM,
      sizeBytes,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.MEDIUM,
      rawRef: particleSystem,
      disposeCallback: options.disposeCallback || (() => {
        if (particleSystem) {
          if (particleSystem.geometry && typeof particleSystem.geometry.dispose === 'function') {
            particleSystem.geometry.dispose();
          }
          if (particleSystem.material && typeof particleSystem.material.dispose === 'function') {
            particleSystem.material.dispose();
          }
        }
      }),
      metadata: options.metadata || {}
    });
  }

  /**
   * Registers a WebGLRenderTarget.
   * @param {string} id
   * @param {Object} renderTarget - Three.js WebGLRenderTarget
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  registerRenderTarget(id, renderTarget, options = {}) {
    const sizeBytes = options.sizeBytes !== undefined
      ? options.sizeBytes
      : MemoryEstimator.estimateRenderTargetBytes(renderTarget);

    return this.registerResource({
      id,
      name: options.name || id,
      type: RESOURCE_TYPE.RENDER_TARGET,
      sizeBytes,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.HIGH,
      rawRef: renderTarget,
      disposeCallback: options.disposeCallback || (() => {
        if (renderTarget && typeof renderTarget.dispose === 'function') {
          renderTarget.dispose();
        }
      }),
      metadata: options.metadata || {}
    });
  }

  /**
   * Registers a Material and tracks its associated textures.
   * @param {string} id
   * @param {Object} material - Three.js Material
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  registerMaterial(id, material, options = {}) {
    const sizeBytes = options.sizeBytes || 2048; // Uniform buffer overhead approximation
    const entry = this.registerResource({
      id,
      name: options.name || (material && material.name) || id,
      type: RESOURCE_TYPE.MATERIAL,
      sizeBytes,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.MEDIUM,
      rawRef: material,
      disposeCallback: options.disposeCallback || (() => {
        if (material && typeof material.dispose === 'function') {
          material.dispose();
        }
      }),
      metadata: options.metadata || {}
    });

    // Auto-discover child textures attached to the material
    if (material) {
      const textureSlots = ['map', 'alphaMap', 'normalMap', 'bumpMap', 'roughnessMap', 'metalnessMap', 'envMap'];
      for (const slot of textureSlots) {
        if (material[slot] && typeof material[slot] === 'object') {
          const tex = material[slot];
          const texId = `${id}_${slot}_tex`;
          if (!this.hasResource(texId)) {
            const texEntry = this.registerTexture(texId, tex, {
              priority: entry.priority,
              name: `${entry.name} [${slot}]`
            });
            entry.addDependency(texId);
            texEntry.addDependent(entry.id);
          }
        }
      }
    }

    return entry;
  }

  /**
   * Generic resource registration method with byte-alignment and pressure enforcement.
   * @param {Object} resourceOptions
   * @returns {ResourceEntry}
   */
  registerResource(resourceOptions) {
    if (this._isDestroyed) {
      throw new MemoryManagerError('Cannot register resource on a destroyed MemoryBudgetManager.');
    }

    let { id, type, sizeBytes, priority, rawRef, disposeCallback, metadata, name } = resourceOptions;

    if (!id || typeof id !== 'string') {
      throw new InvalidResourceError(String(id), 'Resource ID must be a non-empty string.');
    }

    // Align bytes if enabled
    if (this.config.enforceAlignment && this.config.alignmentBytes > 1) {
      sizeBytes = MemoryEstimator.alignBytes(sizeBytes, this.config.alignmentBytes);
    }

    // Check if ID is already registered; if so, replace it cleanly
    if (this._resources.has(id)) {
      this.unregisterResource(id);
    }

    // Enforce hard budget limit if enabled
    const projectedTotal = this._totalTrackedBytes + sizeBytes;
    if (projectedTotal > this.config.hardLimitBytes) {
      if (this.config.autoEvictOnHardLimit) {
        const bytesNeeded = projectedTotal - this.config.hardLimitBytes;
        const freed = this.evictLRU(bytesNeeded + (this.config.hardLimitBytes * this.config.targetEvictionHeadroom));
        if (this._totalTrackedBytes + sizeBytes > this.config.hardLimitBytes) {
          throw new MemoryBudgetExceededError(sizeBytes, this._totalTrackedBytes, this.config.hardLimitBytes, {
            resourceId: id,
            type,
            freedBytes: freed
          });
        }
      } else {
        throw new MemoryBudgetExceededError(sizeBytes, this._totalTrackedBytes, this.config.hardLimitBytes, {
          resourceId: id,
          type
        });
      }
    }

    const entry = new ResourceEntry({
      id,
      name,
      type: type || RESOURCE_TYPE.CUSTOM,
      sizeBytes,
      priority: priority !== undefined ? priority : RESOURCE_PRIORITY.MEDIUM,
      rawRef,
      disposeCallback,
      metadata
    });

    // Track in GC sentinel if enabled
    if (this.config.enableGCSentinel && rawRef && typeof rawRef === 'object') {
      this.gcSentinel.track(rawRef, { id, type, sizeBytes });
    }

    this._resources.set(id, entry);
    this._totalTrackedBytes += sizeBytes;
    this._cumulativeAllocatedBytes += sizeBytes;

    if (this._totalTrackedBytes > this._peakTrackedBytes) {
      this._peakTrackedBytes = this._totalTrackedBytes;
    }

    // Update category bytes
    const catBytes = this._categoryBytes.get(entry.type) || 0;
    this._categoryBytes.set(entry.type, catBytes + sizeBytes);

    // Update memory pressure state
    this._updatePressureState();

    this._emit('allocation', entry.toJSON());

    return entry;
  }

  // ==========================================================================
  // Resource Access, Touch & Pinning
  // ==========================================================================

  /**
   * Retrieves a resource entry by ID.
   * @param {string} id
   * @returns {ResourceEntry|null}
   */
  getResource(id) {
    return this._resources.get(id) || null;
  }

  /**
   * Checks if a resource ID is currently registered.
   * @param {string} id
   * @returns {boolean}
   */
  hasResource(id) {
    return this._resources.has(id);
  }

  /**
   * Touches a resource, refreshing its LRU timestamp and incrementing access count.
   * @param {string} id
   * @returns {boolean}
   */
  touchResource(id) {
    const entry = this._resources.get(id);
    if (entry) {
      entry.touch();
      return true;
    }
    return false;
  }

  /**
   * Pins a resource to protect it from LRU eviction.
   * @param {string} id
   * @returns {boolean}
   */
  pinResource(id) {
    const entry = this._resources.get(id);
    if (entry) {
      entry.pin();
      return true;
    }
    return false;
  }

  /**
   * Unpins a resource.
   * @param {string} id
   * @param {number} [fallbackPriority=RESOURCE_PRIORITY.MEDIUM]
   * @returns {boolean}
   */
  unpinResource(id, fallbackPriority = RESOURCE_PRIORITY.MEDIUM) {
    const entry = this._resources.get(id);
    if (entry) {
      entry.unpin(fallbackPriority);
      return true;
    }
    return false;
  }

  /**
   * Updates priority tier of a resource.
   * @param {string} id
   * @param {number} priority
   * @returns {boolean}
   */
  setResourcePriority(id, priority) {
    const entry = this._resources.get(id);
    if (entry) {
      entry.setPriority(priority);
      return true;
    }
    return false;
  }

  // ==========================================================================
  // De-registration, Disposal & LRU Eviction
  // ==========================================================================

  /**
   * Unregisters a resource from tracking, updates byte accounting, and calls disposal hooks.
   * @param {string} id
   * @param {boolean} [disposeObject=true] - Whether to invoke object disposal
   * @returns {boolean}
   */
  unregisterResource(id, disposeObject = true) {
    const entry = this._resources.get(id);
    if (!entry) return false;

    if (entry.rawRef && typeof entry.rawRef === 'object') {
      this.gcSentinel.untrack(entry.rawRef, id);
    }

    if (disposeObject && !entry.isDisposed) {
      this._disposeEntry(entry);
    }

    // Clean up dependency connections
    for (const depId of entry.dependencies) {
      const dep = this._resources.get(depId);
      if (dep) dep.removeDependent(id);
    }
    for (const parentId of entry.dependents) {
      const parent = this._resources.get(parentId);
      if (parent) parent.removeDependency(id);
    }

    this._resources.delete(id);
    this._totalTrackedBytes = Math.max(0, this._totalTrackedBytes - entry.sizeBytes);
    this._cumulativeFreedBytes += entry.sizeBytes;

    const catBytes = this._categoryBytes.get(entry.type) || 0;
    this._categoryBytes.set(entry.type, Math.max(0, catBytes - entry.sizeBytes));

    this._updatePressureState();

    this._emit('disposal', entry.toJSON());
    return true;
  }

  /**
   * Disposes raw underlying object safely.
   * @param {ResourceEntry} entry
   * @private
   */
  _disposeEntry(entry) {
    try {
      if (typeof entry.disposeCallback === 'function') {
        entry.disposeCallback();
      } else if (entry.rawRef && typeof entry.rawRef.dispose === 'function') {
        entry.rawRef.dispose();
      }
      entry.isDisposed = true;
    } catch (err) {
      entry.isDisposed = true;
      throw new ResourceDisposalError(entry.id, err.message, { originalError: err });
    }
  }

  /**
   * Evicts least recently used resources until target bytes are freed or no more unpinned candidates exist.
   * @param {number} targetBytesToFree - Number of bytes to free
   * @param {number} [maxPriorityToEvict=RESOURCE_PRIORITY.HIGH - 1] - Evict up to this priority tier
   * @returns {number} Total bytes successfully freed
   */
  evictLRU(targetBytesToFree, maxPriorityToEvict = RESOURCE_PRIORITY.HIGH - 1) {
    if (targetBytesToFree <= 0) return 0;

    let freedBytes = 0;

    // Collect eviction candidates (unpinned, priority <= maxPriorityToEvict)
    const candidates = [];
    for (const entry of this._resources.values()) {
      if (!entry.isPinned && entry.priority <= maxPriorityToEvict) {
        candidates.push(entry);
      }
    }

    // Sort candidates: Lowest priority first; then oldest lastAccessedAt
    candidates.sort((a, b) => {
      if (a.priority !== b.priority) {
        return a.priority - b.priority;
      }
      return a.lastAccessedAt - b.lastAccessedAt;
    });

    for (const entry of candidates) {
      if (freedBytes >= targetBytesToFree) break;

      const size = entry.sizeBytes;
      const success = this.unregisterResource(entry.id, true);
      if (success) {
        freedBytes += size;
        this._totalEvictionCount++;
        this._emit('eviction', { id: entry.id, freedBytes: size, type: entry.type });
      }
    }

    return freedBytes;
  }

  /**
   * Purges all resources of a specified category.
   * @param {string} type - RESOURCE_TYPE enum value
   * @param {boolean} [forcePinned=false] - If true, evicts pinned resources as well
   * @returns {number} Total bytes freed
   */
  purgeCategory(type, forcePinned = false) {
    let freedBytes = 0;
    const toDelete = [];

    for (const entry of this._resources.values()) {
      if (entry.type === type) {
        if (forcePinned || !entry.isPinned) {
          toDelete.push(entry);
        }
      }
    }

    for (const entry of toDelete) {
      const size = entry.sizeBytes;
      if (this.unregisterResource(entry.id, true)) {
        freedBytes += size;
      }
    }

    return freedBytes;
  }

  /**
   * Disposes all registered resources and clears memory completely.
   * @param {boolean} [forcePinned=true]
   * @returns {number} Total bytes freed
   */
  disposeAll(forcePinned = true) {
    let freedBytes = 0;
    const toDelete = Array.from(this._resources.values());

    for (const entry of toDelete) {
      if (forcePinned || !entry.isPinned) {
        const size = entry.sizeBytes;
        if (this.unregisterResource(entry.id, true)) {
          freedBytes += size;
        }
      }
    }

    return freedBytes;
  }

  // ==========================================================================
  // Worker Transferable ArrayBuffer Registry
  // ==========================================================================

  /**
   * Registers a zero-copy transferable ArrayBuffer dispatched to or received from a Web Worker.
   * @param {string} id - Buffer identifier
   * @param {ArrayBuffer|SharedArrayBuffer|TypedArray} buffer - Buffer instance
   * @param {Object} [options={}]
   * @returns {ResourceEntry}
   */
  trackWorkerTransfer(id, buffer, options = {}) {
    let byteLength = 0;
    if (buffer && buffer.byteLength !== undefined) {
      byteLength = buffer.byteLength;
    } else if (buffer && buffer.buffer && buffer.buffer.byteLength !== undefined) {
      byteLength = buffer.buffer.byteLength;
    }

    const workerId = options.workerId || 'pool';

    this._workerBuffers.set(id, {
      id,
      workerId,
      direction: options.direction || 'outbound',
      transferredAt: Date.now(),
      byteLength,
      isDetached: buffer ? buffer.byteLength === 0 : false
    });

    if (this.workerCoordinator) {
      this.workerCoordinator.recordDispatch(workerId, byteLength);
    }

    return this.registerResource({
      id: `worker_buf_${id}`,
      name: options.name || `WorkerBuffer [${id}]`,
      type: RESOURCE_TYPE.WORKER_BUFFER,
      sizeBytes: byteLength,
      priority: options.priority !== undefined ? options.priority : RESOURCE_PRIORITY.TRANSIENT,
      rawRef: buffer,
      disposeCallback: () => {
        const info = this._workerBuffers.get(id);
        if (info && this.workerCoordinator) {
          this.workerCoordinator.recordCompletion(info.workerId, info.byteLength);
        }
        this._workerBuffers.delete(id);
      },
      metadata: { workerId, ...options.metadata }
    });
  }

  /**
   * Marks a worker buffer as detached (after postMessage transfer zero-copy detach).
   * @param {string} id
   * @returns {boolean}
   */
  markWorkerBufferDetached(id) {
    const regId = `worker_buf_${id}`;
    const info = this._workerBuffers.get(id);
    if (info) {
      info.isDetached = true;
    }
    return this.unregisterResource(regId, false);
  }

  // ==========================================================================
  // Allocation Guard & OOM Preflight Sentinel
  // ==========================================================================

  /**
   * Evaluates whether an allocation of the specified size is safe without breaching limits.
   * @param {number} bytesToAllocate
   * @param {Object} [options={}]
   * @returns {{canAllocate: boolean, isSoftSafe: boolean, requiredBytes: number, availableHeadroom: number, softHeadroom: number, projectedPressure: string}}
   */
  canAllocate(bytesToAllocate, options = {}) {
    const projectedTotal = this._totalTrackedBytes + bytesToAllocate;
    const isHardSafe = projectedTotal <= this.config.hardLimitBytes;
    const isSoftSafe = projectedTotal <= this.config.softLimitBytes;

    return {
      canAllocate: isHardSafe,
      isSoftSafe,
      requiredBytes: bytesToAllocate,
      availableHeadroom: Math.max(0, this.config.hardLimitBytes - this._totalTrackedBytes),
      softHeadroom: Math.max(0, this.config.softLimitBytes - this._totalTrackedBytes),
      projectedPressure: projectedTotal >= this.config.hardLimitBytes
        ? MEMORY_PRESSURE_LEVEL.EXCEEDED
        : projectedTotal >= this.config.softLimitBytes
          ? MEMORY_PRESSURE_LEVEL.CRITICAL
          : projectedTotal >= (this.config.softLimitBytes * this.config.moderateThresholdRatio)
            ? MEMORY_PRESSURE_LEVEL.MODERATE
            : MEMORY_PRESSURE_LEVEL.NORMAL
    };
  }

  /**
   * Preflight verification guard for 3D cosmological grid / streamline reconstructions.
   * Throws AllocationRefusedError if unsafe and strict mode is active.
   * @param {Object} config - Streamline reconstruction configuration
   * @param {Object} [options={}]
   * @param {boolean} [options.strict=false] - Throw error on failure if true
   * @returns {Object} Safety assessment report
   */
  guardStreamlineReconstruction(config, options = {}) {
    const assessment = StreamlineOOMGuard.assessReconstructionSafety(
      config,
      this.availableSoftHeadroomBytes,
      this.config.softLimitBytes,
      this.config.hardLimitBytes
    );

    if (!assessment.isSafe) {
      this._emit('oomWarning', assessment);
      if (options.strict) {
        throw new AllocationRefusedError(
          'StreamlineReconstruction',
          assessment.requiredBytes,
          assessment.availableBytes,
          assessment.recommendedConfig
        );
      }
    }

    return assessment;
  }

  /**
   * Preflight verification guard for 3D cosmological scalar/vector volume grid.
   * @param {number} resolution - Grid resolution (e.g. 64, 128, 256)
   * @param {number} channels - Number of channels (1 for scalar, 3 for velocity)
   * @param {Object} [options={}]
   * @returns {Object}
   */
  guardGridAllocation(resolution, channels, options = {}) {
    const assessment = StreamlineOOMGuard.assessGridSafety(
      resolution,
      channels,
      this.availableSoftHeadroomBytes
    );

    if (!assessment.isSafe) {
      this._emit('oomWarning', assessment);
      if (options.strict) {
        throw new AllocationRefusedError(
          `GridVolume_${resolution}^3`,
          assessment.requiredBytes,
          assessment.availableBytes,
          { fallbackResolution: assessment.fallbackResolution }
        );
      }
    }

    return assessment;
  }

  // ==========================================================================
  // High-Frequency Telemetry & Polling System
  // ==========================================================================

  /**
   * Updates memory pressure state and fires events on transition.
   * @private
   */
  _updatePressureState() {
    let newLevel = MEMORY_PRESSURE_LEVEL.NORMAL;

    if (this._totalTrackedBytes >= this.config.hardLimitBytes) {
      newLevel = MEMORY_PRESSURE_LEVEL.EXCEEDED;
    } else if (this._totalTrackedBytes >= this.config.softLimitBytes) {
      newLevel = MEMORY_PRESSURE_LEVEL.CRITICAL;
    } else if (this._totalTrackedBytes >= (this.config.softLimitBytes * this.config.moderateThresholdRatio)) {
      newLevel = MEMORY_PRESSURE_LEVEL.MODERATE;
    }

    if (newLevel !== this._pressureLevel) {
      const oldLevel = this._pressureLevel;
      this._pressureLevel = newLevel;
      this._emit('pressureChange', { oldLevel, newLevel, totalTrackedBytes: this._totalTrackedBytes });
    }
  }

  /**
   * Gathers JS Heap statistics from performance.memory (Chrome) or process.memoryUsage (Node).
   * @returns {Object}
   */
  getHeapMetrics() {
    // Browser Chrome performance.memory
    if (typeof performance !== 'undefined' && performance.memory) {
      return {
        usedJSHeapSize: performance.memory.usedJSHeapSize || 0,
        totalJSHeapSize: performance.memory.totalJSHeapSize || 0,
        jsHeapSizeLimit: performance.memory.jsHeapSizeLimit || 0,
        source: 'performance.memory'
      };
    }

    // Node.js process.memoryUsage
    if (typeof process !== 'undefined' && typeof process.memoryUsage === 'function') {
      const mem = process.memoryUsage();
      return {
        usedJSHeapSize: mem.heapUsed || 0,
        totalJSHeapSize: mem.heapTotal || 0,
        jsHeapSizeLimit: mem.rss || 0,
        external: mem.external || 0,
        arrayBuffers: mem.arrayBuffers || 0,
        source: 'process.memoryUsage'
      };
    }

    return {
      usedJSHeapSize: 0,
      totalJSHeapSize: 0,
      jsHeapSizeLimit: 0,
      source: 'unavailable'
    };
  }

  /**
   * Generates a comprehensive real-time telemetry snapshot.
   * @returns {Object}
   */
  getTelemetrySnapshot() {
    const now = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    const heap = this.config.enableHeapTracking ? this.getHeapMetrics() : {};

    const categoryBreakdown = {};
    for (const [cat, bytes] of this._categoryBytes.entries()) {
      categoryBreakdown[cat] = {
        bytes,
        mb: +(bytes / BYTES_PER_MB).toFixed(3),
        percent: this._totalTrackedBytes > 0 ? +((bytes / this._totalTrackedBytes) * 100).toFixed(1) : 0
      };
    }

    return {
      timestamp: now,
      totalTrackedBytes: this._totalTrackedBytes,
      totalTrackedMB: +(this._totalTrackedBytes / BYTES_PER_MB).toFixed(3),
      peakTrackedBytes: this._peakTrackedBytes,
      peakTrackedMB: +(this._peakTrackedBytes / BYTES_PER_MB).toFixed(3),
      softLimitBytes: this.config.softLimitBytes,
      softLimitMB: +(this.config.softLimitBytes / BYTES_PER_MB).toFixed(1),
      hardLimitBytes: this.config.hardLimitBytes,
      hardLimitMB: +(this.config.hardLimitBytes / BYTES_PER_MB).toFixed(1),
      availableSoftHeadroomMB: +(this.availableSoftHeadroomBytes / BYTES_PER_MB).toFixed(3),
      availableHardHeadroomMB: +(this.availableHardHeadroomBytes / BYTES_PER_MB).toFixed(3),
      pressureLevel: this._pressureLevel,
      resourceCount: this._resources.size,
      evictionCount: this._totalEvictionCount,
      cumulativeAllocatedMB: +(this._cumulativeAllocatedBytes / BYTES_PER_MB).toFixed(3),
      cumulativeFreedMB: +(this._cumulativeFreedBytes / BYTES_PER_MB).toFixed(3),
      smoothedAllocationRateMBps: +(this._smoothedAllocationRate / BYTES_PER_MB).toFixed(3),
      categories: categoryBreakdown,
      workerPool: this.workerCoordinator ? this.workerCoordinator.getWorkerTelemetry() : {},
      gcCollectedOrphans: this.gcSentinel ? this.gcSentinel.collectedCount : 0,
      heap
    };
  }

  /**
   * Executes a single telemetry polling tick.
   * @returns {Object}
   */
  pollTelemetry() {
    const now = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    const snapshot = this.getTelemetrySnapshot();

    if (this._lastPollTimestamp > 0) {
      const dt = (now - this._lastPollTimestamp) / 1000.0;
      if (dt > 0.001) {
        // Calculate instantaneous allocation rate
        const lastSnapshot = this._telemetryHistory[this._telemetryHistory.length - 1];
        if (lastSnapshot) {
          const deltaBytes = snapshot.totalTrackedBytes - lastSnapshot.totalTrackedBytes;
          const instantRate = deltaBytes / dt;
          // Exponential moving average smoothing (alpha = 0.2)
          this._smoothedAllocationRate = (0.2 * instantRate) + (0.8 * this._smoothedAllocationRate);
        }
      }
    }

    this._lastPollTimestamp = now;

    // Add to history ring buffer
    this._telemetryHistory.push(snapshot);
    if (this._telemetryHistory.length > this.config.maxHistorySamples) {
      this._telemetryHistory.shift();
    }

    this._emit('telemetryTick', snapshot);
    return snapshot;
  }

  /**
   * Detects potential memory leaks by calculating ordinary least squares (OLS) linear regression
   * slope on recent telemetry samples.
   * @param {number} [sampleCount=30]
   * @returns {{slopeBytesPerSec: number, slopeMBPerMin: number, isLeaking: boolean, confidence: number}}
   */
  detectMemoryLeakTrend(sampleCount = 30) {
    const samples = this._telemetryHistory.slice(-sampleCount);
    if (samples.length < 5) {
      return { slopeBytesPerSec: 0, slopeMBPerMin: 0, isLeaking: false, confidence: 0 };
    }

    const n = samples.length;
    let sumT = 0;
    let sumM = 0;
    let sumTM = 0;
    let sumTT = 0;

    const t0 = samples[0].timestamp;

    for (let i = 0; i < n; i++) {
      const t = (samples[i].timestamp - t0) / 1000.0; // Seconds
      const m = samples[i].totalTrackedBytes;
      sumT += t;
      sumM += m;
      sumTM += (t * m);
      sumTT += (t * t);
    }

    const denom = (n * sumTT) - (sumT * sumT);
    if (Math.abs(denom) < 1e-9) {
      return { slopeBytesPerSec: 0, slopeMBPerMin: 0, isLeaking: false, confidence: 0 };
    }

    const slope = ((n * sumTM) - (sumT * sumM)) / denom; // bytes per second
    const slopeMBPerMin = (slope * 60.0) / BYTES_PER_MB;
    const isLeaking = slopeMBPerMin > 5.0 && this._totalTrackedBytes > (this.config.softLimitBytes * 0.5);

    return {
      slopeBytesPerSec: +slope.toFixed(2),
      slopeMBPerMin: +slopeMBPerMin.toFixed(3),
      isLeaking,
      confidence: Math.min(1.0, n / 30.0)
    };
  }

  /**
   * Starts periodic high-frequency telemetry polling.
   * @param {number} [intervalMs] - Optional override interval
   * @returns {MemoryBudgetManager}
   */
  startTelemetryPolling(intervalMs) {
    if (this._pollTimer) {
      this.stopTelemetryPolling();
    }

    const interval = intervalMs || this.config.telemetryPollIntervalMs;
    this._pollTimer = setInterval(() => {
      this.pollTelemetry();
    }, interval);

    return this;
  }

  /**
   * Stops periodic telemetry polling.
   * @returns {MemoryBudgetManager}
   */
  stopTelemetryPolling() {
    if (this._pollTimer) {
      clearInterval(this._pollTimer);
      this._pollTimer = null;
    }
    return this;
  }

  /**
   * Returns recent telemetry history array.
   * @returns {Array<Object>}
   */
  getHistory() {
    return [...this._telemetryHistory];
  }

  // ==========================================================================
  // Scientific Invariants & Consistency Verification
  // ==========================================================================

  /**
   * Verifies mathematical and memory conservation invariants:
   * 1. totalTrackedBytes == sum of all individual resource sizeBytes.
   * 2. totalTrackedBytes == sum of all category bytes.
   * 3. No negative byte sizes.
   * 4. Soft limit <= hard limit.
   * @returns {{isValid: boolean, discrepancies: Array<string>}}
   */
  verifyInvariants() {
    const discrepancies = [];

    // Invariant 1: Sum of resource sizes == totalTrackedBytes
    let computedResourceSum = 0;
    for (const [id, entry] of this._resources.entries()) {
      if (entry.sizeBytes < 0) {
        discrepancies.push(`Resource '${id}' has negative sizeBytes: ${entry.sizeBytes}`);
      }
      computedResourceSum += entry.sizeBytes;
    }

    if (computedResourceSum !== this._totalTrackedBytes) {
      discrepancies.push(
        `Total bytes discrepancy: tracked = ${this._totalTrackedBytes}, computed resource sum = ${computedResourceSum}`
      );
    }

    // Invariant 2: Sum of category bytes == totalTrackedBytes
    let computedCategorySum = 0;
    for (const [cat, bytes] of this._categoryBytes.entries()) {
      if (bytes < 0) {
        discrepancies.push(`Category '${cat}' has negative bytes: ${bytes}`);
      }
      computedCategorySum += bytes;
    }

    if (computedCategorySum !== this._totalTrackedBytes) {
      discrepancies.push(
        `Category sum discrepancy: tracked = ${this._totalTrackedBytes}, computed category sum = ${computedCategorySum}`
      );
    }

    // Invariant 3: Limits
    if (this.config.softLimitBytes > this.config.hardLimitBytes) {
      discrepancies.push(`softLimitBytes (${this.config.softLimitBytes}) > hardLimitBytes (${this.config.hardLimitBytes})`);
    }

    return {
      isValid: discrepancies.length === 0,
      discrepancies
    };
  }

  // ==========================================================================
  // Reset & Teardown
  // ==========================================================================

  /**
   * Cleans up all resources and stops polling.
   */
  reset() {
    this.stopTelemetryPolling();
    this.disposeAll(true);
    this._resources.clear();
    for (const type of Object.values(RESOURCE_TYPE)) {
      this._categoryBytes.set(type, 0);
    }
    this._totalTrackedBytes = 0;
    this._peakTrackedBytes = 0;
    this._cumulativeAllocatedBytes = 0;
    this._cumulativeFreedBytes = 0;
    this._totalEvictionCount = 0;
    this._pressureLevel = MEMORY_PRESSURE_LEVEL.NORMAL;
    this._workerBuffers.clear();
    this._telemetryHistory = [];
    this._smoothedAllocationRate = 0;
    this._lastPollTimestamp = 0;
  }

  /**
   * Permanently destroys this manager instance and releases all event listeners.
   */
  destroy() {
    this.reset();
    this._eventListeners.clear();
    this._isDestroyed = true;
  }
}
