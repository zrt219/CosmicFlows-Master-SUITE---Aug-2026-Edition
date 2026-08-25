/**
 * @fileoverview Cryptographic Audit Logger, Parameter Tracking, and W3C PROV Execution Manifests for Cosmological Workflows.
 * 
 * Provides an immutable, tamper-evident hash-chained audit log for computational cosmography runs.
 * Tracks algorithm parameters, coordinate transforms, filter execution stages, data dependencies,
 * and outputs machine-readable W3C PROV-O / PROV-DM JSON-LD execution manifests.
 * 
 * @module provenance/audit_logger
 */

import { pureJsSha256 } from './manifest.js';

/**
 * Genesis root constant for cosmological audit trail hash chaining.
 * @type {string}
 */
export const GENESIS_ROOT_HASH = '0000000000000000000000000000000000000000000000000000000000000000';

/**
 * Standard log levels for audit events.
 * @readonly
 * @enum {string}
 */
export const LogLevel = Object.freeze({
  DEBUG: 'DEBUG',
  INFO: 'INFO',
  PARAM: 'PARAM',
  TRANSFORM: 'TRANSFORM',
  STAGE: 'STAGE',
  CHECKPOINT: 'CHECKPOINT',
  WARN: 'WARN',
  ERROR: 'ERROR'
});

/**
 * Deterministically serializes an arbitrary object to a canonical JSON string (sorted keys).
 * 
 * @param {*} obj - Object or value to serialize.
 * @returns {string} Canonical JSON representation.
 */
export function canonicalJson(obj) {
  if (obj === null || typeof obj !== 'object') {
    return JSON.stringify(obj);
  }
  if (Array.isArray(obj)) {
    return '[' + obj.map(item => canonicalJson(item)).join(',') + ']';
  }
  const keys = Object.keys(obj).sort();
  const pairs = keys.map(k => `${JSON.stringify(k)}:${canonicalJson(obj[k])}`);
  return '{' + pairs.join(',') + '}';
}

/**
 * Computes SHA-256 for a canonical string.
 * @param {string} str - Input text.
 * @returns {string} 64-character hex hash.
 */
function hashString(str) {
  const utf8 = unescape(encodeURIComponent(str));
  const arr = new Uint8Array(utf8.length);
  for (let i = 0; i < utf8.length; i++) {
    arr[i] = utf8.charCodeAt(i);
  }
  return pureJsSha256(arr);
}

/**
 * Represents a single verifiable event record in the audit chain.
 */
export class AuditRecord {
  /**
   * @param {object} params
   * @param {number} params.index - 0-based sequential index.
   * @param {string} params.timestamp - ISO-8601 UTC timestamp.
   * @param {string} params.runId - Execution run identifier.
   * @param {string} params.level - Event severity level.
   * @param {string} params.type - Categorical event type.
   * @param {*} params.payload - Arbitrary event payload.
   * @param {string} params.prevHash - SHA-256 hash of previous record.
   * @param {string} [params.hash] - Optional precalculated record hash.
   */
  constructor({ index, timestamp, runId, level, type, payload, prevHash, hash }) {
    this.index = index;
    this.timestamp = timestamp;
    this.runId = runId;
    this.level = level;
    this.type = type;
    this.payload = payload;
    this.prevHash = prevHash;
    this.hash = hash || this.computeHash();
  }

  /**
   * Calculates the tamper-evident hash for this record.
   * @returns {string} SHA-256 hex string.
   */
  computeHash() {
    const rawContent = `${this.prevHash}::${this.index}::${this.timestamp}::${this.runId}::${this.level}::${this.type}::${canonicalJson(this.payload)}`;
    return hashString(rawContent);
  }

  /**
   * Validates that this record's internal hash matches its contents.
   * @returns {boolean}
   */
  isValid() {
    return this.hash === this.computeHash();
  }
}

/**
 * Cryptographic Audit Logger providing tamper-evident run tracking and PROV manifest generation.
 */
export class AuditLogger {
  /**
   * @param {object} [options]
   * @param {string} [options.runId] - Unique ID for this computation session.
   * @param {string} [options.appName='ZRT-Cosmicflows-Workbench'] - Host application name.
   * @param {string} [options.appVersion='2026.2.0'] - Software semantic version.
   * @param {object} [options.environment] - Environment telemetry (Node, Browser, OS, WebGL).
   */
  constructor(options = {}) {
    this.runId = options.runId || `run_${Date.now()}_${Math.random().toString(36).slice(2, 9)}`;
    this.appName = options.appName || 'ZRT-Cosmicflows-Workbench';
    this.appVersion = options.appVersion || '2026.2.0';
    this.sessionStartTime = new Date().toISOString();
    this.environment = options.environment || this._detectEnvironment();

    /** @type {AuditRecord[]} */
    this.records = [];
    /** @type {Map<string, object>} */
    this.parameterStore = new Map();
    /** @type {Set<string>} */
    this.registeredEntities = new Set();
    /** @type {Set<string>} */
    this.registeredActivities = new Set();

    // Create Genesis Record
    this._appendRecord(LogLevel.CHECKPOINT, 'GENESIS', {
      appName: this.appName,
      appVersion: this.appVersion,
      sessionStartTime: this.sessionStartTime,
      environment: this.environment
    });
  }

  /**
   * Detects current runtime environment characteristics.
   * @private
   * @returns {object}
   */
  _detectEnvironment() {
    const env = {
      isNode: typeof process !== 'undefined' && Boolean(process.versions?.node),
      isBrowser: typeof window !== 'undefined' && typeof document !== 'undefined',
      platform: typeof process !== 'undefined' ? process.platform : (typeof navigator !== 'undefined' ? navigator.platform : 'unknown')
    };

    if (env.isNode) {
      env.nodeVersion = process.version;
      env.arch = process.arch;
    } else if (env.isBrowser) {
      env.userAgent = navigator.userAgent;
    }
    return env;
  }

  /**
   * Appends an event to the cryptographic audit chain.
   * @private
   * @param {string} level - Log level.
   * @param {string} type - Event type string.
   * @param {*} payload - Structured data payload.
   * @returns {AuditRecord} Created record.
   */
  _appendRecord(level, type, payload) {
    const index = this.records.length;
    const prevHash = index === 0 ? GENESIS_ROOT_HASH : this.records[index - 1].hash;
    const timestamp = new Date().toISOString();

    const record = new AuditRecord({
      index,
      timestamp,
      runId: this.runId,
      level,
      type,
      payload,
      prevHash
    });

    this.records.push(record);
    return record;
  }

  /**
   * Log a general information or telemetry event.
   * @param {string} type - Event type description.
   * @param {*} payload - Structured payload.
   * @returns {AuditRecord}
   */
  logInfo(type, payload) {
    return this._appendRecord(LogLevel.INFO, type, payload);
  }

  /**
   * Log an event with a custom level.
   * @param {string} level - Log level from LogLevel enum.
   * @param {string} type - Event type.
   * @param {*} payload - Event payload.
   * @returns {AuditRecord}
   */
  logEvent(level, type, payload) {
    return this._appendRecord(level, type, payload);
  }

  /**
   * Log a coordinate transformation event.
   * @param {string} fromFrame - Source frame (e.g. 'equatorial_spherical').
   * @param {string} toFrame - Destination frame (e.g. 'supergalactic_cartesian').
   * @param {object} details - Number of points, bounds, precision.
   * @returns {AuditRecord}
   */
  logCoordinateTransform(fromFrame, toFrame, details) {
    return this._appendRecord(LogLevel.TRANSFORM, 'COORDINATE_TRANSFORM', {
      fromFrame,
      toFrame,
      ...details
    });
  }

  /**
   * Records and tracks cosmological and algorithm parameters.
   * 
   * @param {string} category - Category name (e.g. 'cosmology', 'wiener_filter', 'grid').
   * @param {object} params - Key-value parameter dictionary.
   * @returns {AuditRecord}
   */
  trackParameters(category, params) {
    const existing = this.parameterStore.get(category) || {};
    const updated = { ...existing, ...params };
    this.parameterStore.set(category, updated);

    return this._appendRecord(LogLevel.PARAM, `PARAM_UPDATE_${category.toUpperCase()}`, {
      category,
      parameters: params,
      mergedParameters: updated
    });
  }

  /**
   * Retrieves active parameters for a category.
   * @param {string} category - Category name.
   * @returns {object|null} Active parameter map or null.
   */
  getParameters(category) {
    return this.parameterStore.get(category) ? { ...this.parameterStore.get(category) } : null;
  }

  /**
   * Wraps and executes a computation stage, logging duration, telemetry, and output integrity.
   * 
   * @template T
   * @param {string} stageName - Human-readable stage name.
   * @param {() => Promise<T>|T} stageFn - Computation logic.
   * @param {object} [metadata] - Input entity IDs, parameters, expected counts.
   * @returns {Promise<T>} Result of stage execution.
   */
  async logStage(stageName, stageFn, metadata = {}) {
    const stageId = `activity_${stageName.toLowerCase().replace(/[^a-z0-9_]/g, '_')}_${Date.now()}`;
    this.registeredActivities.add(stageId);

    const startTime = performance.now();
    this._appendRecord(LogLevel.STAGE, 'STAGE_START', {
      stageId,
      stageName,
      metadata
    });

    try {
      const result = await stageFn();
      const durationMs = performance.now() - startTime;

      let outputTelemetry = {};
      if (result && typeof result === 'object') {
        if (ArrayBuffer.isView(result) || result instanceof ArrayBuffer) {
          outputTelemetry.byteLength = result.byteLength;
        } else if (Array.isArray(result)) {
          outputTelemetry.itemCount = result.length;
        } else if (typeof result.length === 'number') {
          outputTelemetry.length = result.length;
        }
      }

      this._appendRecord(LogLevel.STAGE, 'STAGE_END', {
        stageId,
        stageName,
        status: 'SUCCESS',
        durationMs: Math.round(durationMs * 100) / 100,
        outputTelemetry
      });

      return result;
    } catch (error) {
      const durationMs = performance.now() - startTime;
      this._appendRecord(LogLevel.ERROR, 'STAGE_ERROR', {
        stageId,
        stageName,
        status: 'FAILED',
        durationMs: Math.round(durationMs * 100) / 100,
        errorMessage: error?.message || String(error),
        stack: error?.stack
      });
      throw error;
    }
  }

  /**
   * Verifies the cryptographic integrity of the entire audit chain from genesis to head.
   * 
   * @returns {{valid: boolean, recordCount: number, headHash: string, corruptedIndex: number|null, error?: string}}
   */
  verifyIntegrity() {
    if (this.records.length === 0) {
      return { valid: false, recordCount: 0, headHash: '', corruptedIndex: 0, error: 'Empty audit log' };
    }

    // Check genesis record
    if (this.records[0].prevHash !== GENESIS_ROOT_HASH) {
      return {
        valid: false,
        recordCount: this.records.length,
        headHash: this.records[this.records.length - 1].hash,
        corruptedIndex: 0,
        error: `Genesis record prevHash mismatch. Expected ${GENESIS_ROOT_HASH}, got ${this.records[0].prevHash}`
      };
    }

    for (let i = 0; i < this.records.length; i++) {
      const rec = this.records[i];

      // Verify sequence index
      if (rec.index !== i) {
        return {
          valid: false,
          recordCount: this.records.length,
          headHash: this.records[this.records.length - 1].hash,
          corruptedIndex: i,
          error: `Sequence index broken at position ${i}. Expected ${i}, found ${rec.index}`
        };
      }

      // Verify internal hash calculation
      const calculatedHash = rec.computeHash();
      if (rec.hash !== calculatedHash) {
        return {
          valid: false,
          recordCount: this.records.length,
          headHash: this.records[this.records.length - 1].hash,
          corruptedIndex: i,
          error: `Record content corrupted at index ${i}. Recorded hash: ${rec.hash}, computed hash: ${calculatedHash}`
        };
      }

      // Verify chain linkage to previous record
      if (i > 0) {
        const prev = this.records[i - 1];
        if (rec.prevHash !== prev.hash) {
          return {
            valid: false,
            recordCount: this.records.length,
            headHash: this.records[this.records.length - 1].hash,
            corruptedIndex: i,
            error: `Chain broken between index ${i - 1} and ${i}. Expected prevHash ${prev.hash}, got ${rec.prevHash}`
          };
        }
      }
    }

    return {
      valid: true,
      recordCount: this.records.length,
      headHash: this.records[this.records.length - 1].hash,
      corruptedIndex: null
    };
  }

  /**
   * Generates a W3C PROV-DM / PROV-O JSON-LD compliant machine-readable execution manifest.
   * 
   * @param {object} [options]
   * @param {string[]} [options.usedDatasets] - Dataset IDs consumed.
   * @param {string[]} [options.generatedArtifacts] - Artifact IDs produced.
   * @returns {object} PROV-O JSON-LD document.
   */
  generateExecutionManifest(options = {}) {
    const integrity = this.verifyIntegrity();
    const manifestId = `prov:run:${this.runId}`;

    const manifest = {
      '@context': {
        'prov': 'http://www.w3.org/ns/prov#',
        'xsd': 'http://www.w3.org/2001/XMLSchema#',
        'zrt': 'https://ip2i.in2p3.fr/cosmography/ontology#'
      },
      '@id': manifestId,
      '@type': ['prov:Bundle', 'zrt:ExecutionManifest'],
      'prov:generatedAtTime': new Date().toISOString(),
      'zrt:integrity': {
        'zrt:valid': integrity.valid,
        'zrt:recordCount': integrity.recordCount,
        'zrt:headHash': integrity.headHash,
        'zrt:runId': this.runId
      },
      'prov:agent': {
        '@id': `agent:software:${this.appName}`,
        '@type': ['prov:SoftwareAgent', 'prov:Agent'],
        'zrt:appName': this.appName,
        'zrt:appVersion': this.appVersion,
        'zrt:environment': this.environment
      },
      'zrt:parameters': Object.fromEntries(this.parameterStore.entries()),
      'prov:activity': [],
      'prov:entity': [],
      'prov:wasGeneratedBy': [],
      'prov:used': []
    };

    // Synthesize activities from stage log records
    const stageRecords = this.records.filter(r => r.type === 'STAGE_START' || r.type === 'STAGE_END');
    const stageMap = new Map();

    for (const r of stageRecords) {
      const stageId = r.payload.stageId;
      if (!stageMap.has(stageId)) {
        stageMap.set(stageId, { id: stageId, name: r.payload.stageName });
      }
      const entry = stageMap.get(stageId);
      if (r.type === 'STAGE_START') {
        entry.startTime = r.timestamp;
        entry.metadata = r.payload.metadata;
      } else if (r.type === 'STAGE_END') {
        entry.endTime = r.timestamp;
        entry.status = r.payload.status;
        entry.durationMs = r.payload.durationMs;
      }
    }

    for (const [id, act] of stageMap.entries()) {
      manifest['prov:activity'].push({
        '@id': `activity:${id}`,
        '@type': 'prov:Activity',
        'prov:label': act.name,
        'prov:startedAtTime': act.startTime,
        'prov:endedAtTime': act.endTime,
        'prov:wasAssociatedWith': `agent:software:${this.appName}`,
        'zrt:durationMs': act.durationMs,
        'zrt:status': act.status || 'UNKNOWN'
      });
    }

    // Register datasets used
    if (options.usedDatasets) {
      for (const ds of options.usedDatasets) {
        const entityId = `entity:dataset:${ds}`;
        manifest['prov:entity'].push({
          '@id': entityId,
          '@type': ['prov:Entity', 'zrt:CatalogDataset'],
          'prov:label': ds
        });
        manifest['prov:used'].push({
          'prov:activity': manifest['prov:activity'][0]?.['@id'] || `activity:main:${this.runId}`,
          'prov:entity': entityId
        });
      }
    }

    // Register generated artifacts
    if (options.generatedArtifacts) {
      for (const art of options.generatedArtifacts) {
        const entityId = `entity:artifact:${art}`;
        manifest['prov:entity'].push({
          '@id': entityId,
          '@type': ['prov:Entity', 'zrt:ReconstructionField'],
          'prov:label': art
        });
        manifest['prov:wasGeneratedBy'].push({
          'prov:entity': entityId,
          'prov:activity': manifest['prov:activity'][manifest['prov:activity'].length - 1]?.['@id'] || `activity:main:${this.runId}`
        });
      }
    }

    return manifest;
  }

  /**
   * Exports the complete audit log in JSON format.
   * @param {boolean} [pretty=true] - Format with indents.
   * @returns {string} JSON text.
   */
  toJson(pretty = true) {
    const data = {
      runId: this.runId,
      appName: this.appName,
      appVersion: this.appVersion,
      sessionStartTime: this.sessionStartTime,
      environment: this.environment,
      integrity: this.verifyIntegrity(),
      parameters: Object.fromEntries(this.parameterStore.entries()),
      records: this.records.map(r => ({
        index: r.index,
        timestamp: r.timestamp,
        runId: r.runId,
        level: r.level,
        type: r.type,
        payload: r.payload,
        prevHash: r.prevHash,
        hash: r.hash
      }))
    };
    return pretty ? JSON.stringify(data, null, 2) : JSON.stringify(data);
  }

  /**
   * Generates a readable Markdown summary report of the audit trail.
   * @returns {string} Markdown text.
   */
  toMarkdownReport() {
    const integrity = this.verifyIntegrity();
    const lines = [];
    lines.push(`# Cryptographic Execution Audit Report`);
    lines.push(`**Run ID:** \`${this.runId}\`  `);
    lines.push(`**Application:** ${this.appName} v${this.appVersion}  `);
    lines.push(`**Session Start:** ${this.sessionStartTime}  `);
    lines.push(`**Audit Chain Integrity:** ${integrity.valid ? '✅ VERIFIED' : '❌ CORRUPTED'} (${integrity.recordCount} records, Head Hash: \`${integrity.headHash.slice(0, 16)}...\`)  \n`);

    lines.push(`## Active Parameters`);
    for (const [cat, params] of this.parameterStore.entries()) {
      lines.push(`### ${cat.toUpperCase()}`);
      lines.push('```json');
      lines.push(JSON.stringify(params, null, 2));
      lines.push('```');
    }

    lines.push(`\n## Audit Event Chronology`);
    lines.push(`| # | Timestamp | Level | Type | Record Hash |`);
    lines.push(`|---|---|---|---|---|`);
    for (const r of this.records) {
      lines.push(`| ${r.index} | ${r.timestamp.slice(11, 23)} | \`${r.level}\` | ${r.type} | \`${r.hash.slice(0, 12)}...\` |`);
    }

    return lines.join('\n');
  }

  /**
   * Deserializes and verifies an audit log from JSON.
   * 
   * @param {string|object} jsonInput - Serialized JSON string or parsed object.
   * @returns {AuditLogger} Reconstructed and verified AuditLogger instance.
   */
  static fromJson(jsonInput) {
    const data = typeof jsonInput === 'string' ? JSON.parse(jsonInput) : jsonInput;
    const logger = new AuditLogger({
      runId: data.runId,
      appName: data.appName,
      appVersion: data.appVersion,
      environment: data.environment
    });

    logger.sessionStartTime = data.sessionStartTime;
    logger.records = []; // Clear default genesis, reconstruct exact chain

    for (const r of data.records) {
      const record = new AuditRecord({
        index: r.index,
        timestamp: r.timestamp,
        runId: r.runId,
        level: r.level,
        type: r.type,
        payload: r.payload,
        prevHash: r.prevHash,
        hash: r.hash
      });
      logger.records.push(record);
    }

    if (data.parameters) {
      for (const [k, v] of Object.entries(data.parameters)) {
        logger.parameterStore.set(k, v);
      }
    }

    // Verify reconstructed log
    const verification = logger.verifyIntegrity();
    if (!verification.valid) {
      throw new Error(`Imported audit log failed cryptographic verification at index ${verification.corruptedIndex}: ${verification.error}`);
    }

    return logger;
  }
}
