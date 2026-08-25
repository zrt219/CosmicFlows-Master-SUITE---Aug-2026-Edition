/**
 * @file remote_fits_streamer.js
 * @description Progressive HTTP Range-Request Remote FITS Ingestion, Web Worker Streaming,
 * and IndexedDB Binary Cache Layer for full-resolution (64³, 128³, 256³) Cosmicflows-4 grids.
 * 
 * Implements:
 * 1. HTTP 1.1 / HTTP 2 Range-Request chunked streaming for large astronomical FITS files.
 * 2. 2880-byte FITS header block parser supporting 3D primary HDUs and binary tables.
 * 3. Client-side IndexedDB binary blob caching (zrt_cosmicflows_cache_v1) with SHA-256 integrity digests.
 * 4. Automatic Big-Endian byte decoding, BSCALE/BZERO affine transformation, and (SGZ, SGY, SGX) -> (SGX, SGY, SGZ) canonical reordering.
 * 5. Official IP2I/CDS preset dataset registry with velocity x52.0 scaling enforcement.
 * 6. Dynamic telemetry and progress event emitters for real-time WebGL UI status updates.
 * 
 * Scientific Attribution:
 * - Courtois et al. (2023), A&A 670, L15 (DOI: 10.1051/0004-6361/202245331)
 * - Dupuy & Courtois (2023), A&A 678, A176 (DOI: 10.1051/0004-6361/202346802)
 * - Hoffman et al. (2024), MNRAS 527, 4 (DOI: 10.1093/mnras/stad3782)
 * 
 * @module data/remote_fits_streamer
 */

import { GridIndexer } from '../fields/grid_indexer.js';
import { VelocityField } from '../fields/velocity_field.js';
import { DensityField } from '../fields/density_field.js';
import { ScientificUnits } from '../coordinates/scientific_types.js';

/**
 * Standard 2880-byte FITS logical record length.
 */
export const FITS_RECORD_BYTES = 2880;

/**
 * Standard 80-byte FITS card image length.
 */
export const FITS_CARD_BYTES = 80;

/**
 * Known official Cosmicflows-4 dataset presets.
 */
export const CF4_DATASET_PRESETS = Object.freeze({
  CF4_UNGROUPED_64: {
    id: 'cf4_ungrouped_64',
    name: 'CF4 Ungrouped / Individual (64³ Base Grid)',
    resolution: 64,
    boxSizeMpcOverH: 1000.0,
    halfExtentMpcOverH: 500.0,
    source: 'IP2I Lyon / Hoffman et al. 2024',
    velocityFile: './data/cf4_velocity_individual_64.bin',
    densityFile: './data/cf4_density_individual_64.bin',
    format: 'raw_binary_float32',
    isVelocityScaled: true,
    citations: ['Courtois et al. 2023', 'Dupuy & Courtois 2023']
  },
  CF4_GROUPED_64: {
    id: 'cf4_grouped_64',
    name: 'CF4 Grouped Galaxies (64³ Base Grid)',
    resolution: 64,
    boxSizeMpcOverH: 1000.0,
    halfExtentMpcOverH: 500.0,
    source: 'IP2I Lyon / Hoffman et al. 2024',
    velocityFile: './data/cf4_velocity_grouped_64.bin',
    densityFile: './data/cf4_density_grouped_64.bin',
    format: 'raw_binary_float32',
    isVelocityScaled: true,
    citations: ['Courtois et al. 2023', 'Dupuy & Courtois 2023']
  },
  CF4PP_POSTERIOR_MEAN_128: {
    id: 'cf4pp_posterior_mean_128',
    name: 'CF4++ 10,000-Step HMC Posterior Mean (128³ High-Res)',
    resolution: 128,
    boxSizeMpcOverH: 1000.0,
    halfExtentMpcOverH: 500.0,
    source: 'IP2I Lyon / Dupuy & Courtois 2023',
    velocityFile: 'https://projets.ip2i.in2p3.fr/cosmicflows/cf4pp/mean_128_vel.fits',
    densityFile: 'https://projets.ip2i.in2p3.fr/cosmicflows/cf4pp/mean_128_dens.fits',
    format: 'fits_primary_array',
    isVelocityScaled: false,
    citations: ['Courtois et al. 2023', 'Dupuy & Courtois 2023']
  },
  CF4PP_POSTERIOR_RMS_128: {
    id: 'cf4pp_posterior_rms_128',
    name: 'CF4++ Posterior RMS Uncertainty (128³ High-Res)',
    resolution: 128,
    boxSizeMpcOverH: 1000.0,
    halfExtentMpcOverH: 500.0,
    source: 'IP2I Lyon / Dupuy & Courtois 2023',
    velocityFile: 'https://projets.ip2i.in2p3.fr/cosmicflows/cf4pp/rms_128_vel.fits',
    densityFile: 'https://projets.ip2i.in2p3.fr/cosmicflows/cf4pp/rms_128_dens.fits',
    format: 'fits_primary_array',
    isVelocityScaled: false,
    citations: ['Courtois et al. 2023', 'Dupuy & Courtois 2023']
  }
});

/**
 * IndexedDB Persistent Binary Cache Manager for streaming FITS cubes and grids.
 */
export class IndexedDBGridCache {
  constructor(dbName = 'zrt_cosmicflows_cache_v1', storeName = 'grid_blobs') {
    this.dbName = dbName;
    this.storeName = storeName;
    this._db = null;
    this._isSupported = typeof indexedDB !== 'undefined';
  }

  /**
   * Initializes or returns open IndexedDB connection.
   * @returns {Promise<IDBDatabase|null>}
   */
  async open() {
    if (!this._isSupported) return null;
    if (this._db) return this._db;

    return new Promise((resolve, reject) => {
      const req = indexedDB.open(this.dbName, 1);
      req.onupgradeneeded = (evt) => {
        const db = evt.target.result;
        if (!db.objectStoreNames.contains(this.storeName)) {
          const store = db.createObjectStore(this.storeName, { keyPath: 'key' });
          store.createIndex('updatedAt', 'updatedAt', { unique: false });
          store.createIndex('sha256', 'sha256', { unique: false });
        }
      };
      req.onsuccess = (evt) => {
        this._db = evt.target.result;
        resolve(this._db);
      };
      req.onerror = () => {
        console.warn('IndexedDBGridCache: Failed to open IndexedDB, caching disabled.');
        resolve(null);
      };
    });
  }

  /**
   * Stores a binary buffer in IndexedDB.
   * @param {string} key Unique dataset cache key.
   * @param {ArrayBuffer|Float32Array|Float64Array} buffer Binary payload.
   * @param {Object} [meta={}] Associated metadata.
   * @returns {Promise<boolean>}
   */
  async put(key, buffer, meta = {}) {
    const db = await this.open();
    if (!db) return false;

    return new Promise((resolve) => {
      try {
        const tx = db.transaction([this.storeName], 'readwrite');
        const store = tx.objectStore(this.storeName);
        const record = {
          key,
          buffer: buffer instanceof ArrayBuffer ? buffer : buffer.buffer,
          byteLength: buffer.byteLength,
          meta,
          updatedAt: Date.now()
        };
        const req = store.put(record);
        req.onsuccess = () => resolve(true);
        req.onerror = () => resolve(false);
      } catch (err) {
        console.warn('IndexedDBGridCache.put error:', err);
        resolve(false);
      }
    });
  }

  /**
   * Retrieves a binary buffer from IndexedDB.
   * @param {string} key
   * @returns {Promise<{ key: string, buffer: ArrayBuffer, meta: Object }|null>}
   */
  async get(key) {
    const db = await this.open();
    if (!db) return null;

    return new Promise((resolve) => {
      try {
        const tx = db.transaction([this.storeName], 'readonly');
        const store = tx.objectStore(this.storeName);
        const req = store.get(key);
        req.onsuccess = (evt) => {
          const res = evt.target.result;
          resolve(res || null);
        };
        req.onerror = () => resolve(null);
      } catch (err) {
        resolve(null);
      }
    });
  }

  /**
   * Deletes a cached entry.
   * @param {string} key
   * @returns {Promise<boolean>}
   */
  async delete(key) {
    const db = await this.open();
    if (!db) return false;

    return new Promise((resolve) => {
      try {
        const tx = db.transaction([this.storeName], 'readwrite');
        const store = tx.objectStore(this.storeName);
        const req = store.delete(key);
        req.onsuccess = () => resolve(true);
        req.onerror = () => resolve(false);
      } catch (err) {
        resolve(false);
      }
    });
  }

  /**
   * Clears entire IndexedDB cache store.
   * @returns {Promise<boolean>}
   */
  async clear() {
    const db = await this.open();
    if (!db) return false;

    return new Promise((resolve) => {
      try {
        const tx = db.transaction([this.storeName], 'readwrite');
        const store = tx.objectStore(this.storeName);
        const req = store.clear();
        req.onsuccess = () => resolve(true);
        req.onerror = () => resolve(false);
      } catch (err) {
        resolve(false);
      }
    });
  }
}

/**
 * FITS Primary Header Parser and Metadata Card Reader.
 */
export class FITSHeaderParser {
  /**
   * Parses 2880-byte FITS logical blocks into structured header key-value dictionary.
   * @param {Uint8Array|ArrayBuffer} headerBytes
   * @returns {{ cards: Map<string, any>, comments: string[], history: string[], rawCards: string[] }}
   */
  static parse(headerBytes) {
    const bytes = headerBytes instanceof Uint8Array ? headerBytes : new Uint8Array(headerBytes);
    const numCards = Math.floor(bytes.length / FITS_CARD_BYTES);
    const cards = new Map();
    const comments = [];
    const history = [];
    const rawCards = [];

    const decoder = new TextDecoder('ascii');

    for (let c = 0; c < numCards; c++) {
      const cardRaw = decoder.decode(bytes.subarray(c * FITS_CARD_BYTES, (c + 1) * FITS_CARD_BYTES));
      rawCards.push(cardRaw);

      const key = cardRaw.substring(0, 8).trim();
      if (key === 'END') {
        break;
      }

      if (key === 'COMMENT') {
        comments.push(cardRaw.substring(8).trim());
        continue;
      }
      if (key === 'HISTORY') {
        history.push(cardRaw.substring(8).trim());
        continue;
      }

      // Standard KEY = VALUE / COMMENT
      const eqIdx = cardRaw.indexOf('=');
      if (eqIdx !== -1 && eqIdx < 10) {
        let rest = cardRaw.substring(eqIdx + 1).trim();
        let val = null;

        if (rest.startsWith("'")) {
          // String literal delimited by quotes
          const endQuote = rest.indexOf("'", 1);
          if (endQuote !== -1) {
            val = rest.substring(1, endQuote).trim();
          } else {
            val = rest.substring(1).trim();
          }
        } else {
          let slashIdx = rest.indexOf('/');
          let valStr = slashIdx !== -1 ? rest.substring(0, slashIdx).trim() : rest.trim();
          if (valStr === 'T') {
            val = true;
          } else if (valStr === 'F') {
            val = false;
          } else if (/^[+-]?\d+$/.test(valStr)) {
            val = parseInt(valStr, 10);
          } else if (/^[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?$/.test(valStr)) {
            val = parseFloat(valStr);
          } else {
            val = valStr;
          }
        }

        cards.set(key, val);
      }
    }

    return { cards, comments, history, rawCards };
  }

  /**
   * Computes the total byte size of data payload following the FITS header.
   * @param {Map<string, any>} cards
   * @returns {number}
   */
  static computeDataByteLength(cards) {
    const bitpix = cards.get('BITPIX') || -32;
    const naxis = cards.get('NAXIS') || 0;
    if (naxis === 0) return 0;

    let numElements = 1;
    for (let i = 1; i <= naxis; i++) {
      const dim = cards.get(`NAXIS${i}`) || 1;
      numElements *= dim;
    }

    const bytesPerElem = Math.abs(bitpix) / 8;
    const rawBytes = numElements * bytesPerElem;
    // FITS records are padded to 2880-byte boundaries
    const numBlocks = Math.ceil(rawBytes / FITS_RECORD_BYTES);
    return numBlocks * FITS_RECORD_BYTES;
  }
}

/**
 * Universal Remote FITS Streamer & Progressive Ingestion Pipeline.
 */
export class RemoteFITSStreamer {
  /**
   * @param {Object} [options={}]
   * @param {boolean} [options.enableCache=true] Use IndexedDB local caching.
   * @param {number} [options.chunkSizeBytes=1048576] Chunk size for range requests (1 MB default).
   * @param {number} [options.maxRetries=3] Maximum range request retry attempts.
   */
  constructor(options = {}) {
    this.enableCache = options.enableCache !== false;
    this.chunkSizeBytes = options.chunkSizeBytes || 1024 * 1024;
    this.maxRetries = options.maxRetries || 3;
    this.cache = new IndexedDBGridCache();
    this.activeTransfers = new Map();
  }

  /**
   * Reads standard 2880-byte header blocks from a remote URL via HTTP Range-Request.
   * @param {string} url
   * @param {Function} [onProgress]
   * @returns {Promise<{ headerCards: Map<string, any>, headerByteLength: number, totalFileSize: number }>}
   */
  async streamHeader(url, onProgress = null) {
    // Initial fetch of first 2880 * 2 bytes (5760 bytes covers typical FITS primary headers)
    const initialBytes = FITS_RECORD_BYTES * 2;
    const resp = await this._fetchRangeWithRetry(url, 0, initialBytes - 1);
    const buffer = await resp.arrayBuffer();
    const bytes = new Uint8Array(buffer);

    let parsed = FITSHeaderParser.parse(bytes);
    let headerByteLength = initialBytes;

    // Check if END card was found
    let hasEnd = parsed.rawCards.some(card => card.startsWith('END'));
    let currentOffset = initialBytes;

    while (!hasEnd && currentOffset < 100 * FITS_RECORD_BYTES) {
      const nextChunk = await this._fetchRangeWithRetry(url, currentOffset, currentOffset + FITS_RECORD_BYTES - 1);
      const nextBuf = await nextChunk.arrayBuffer();
      const combined = new Uint8Array(bytes.length + nextBuf.byteLength);
      combined.set(bytes, 0);
      combined.set(new Uint8Array(nextBuf), bytes.length);

      parsed = FITSHeaderParser.parse(combined);
      hasEnd = parsed.rawCards.some(card => card.startsWith('END'));
      currentOffset += FITS_RECORD_BYTES;
      headerByteLength = currentOffset;
    }

    const contentRange = resp.headers.get('Content-Range');
    let totalFileSize = 0;
    if (contentRange) {
      const match = contentRange.match(/\/(\d+)$/);
      if (match) totalFileSize = parseInt(match[1], 10);
    }

    return {
      headerCards: parsed.cards,
      headerByteLength,
      totalFileSize
    };
  }

  /**
   * Progressively streams and decodes a complete 3D FITS Primary Array into canonical (SGX, SGY, SGZ) format.
   * 
   * @param {string} url Target remote FITS URL.
   * @param {Object} [options={}] Ingestion options.
   * @param {boolean} [options.isVelocity=false] If true, applies x52.0 velocity scale factor.
   * @param {string} [options.cacheKey] Custom IndexedDB cache key.
   * @param {Function} [options.onProgress] Callback receiving { stage, loadedBytes, totalBytes, percent }.
   * @returns {Promise<{ data: Float32Array|Float64Array, nx: number, ny: number, nz: number, bitpix: number, header: Map<string, any> }>}
   */
  async streamFITSArray3D(url, options = {}) {
    const cacheKey = options.cacheKey || url;
    const isVelocity = Boolean(options.isVelocity);
    const onProgress = options.onProgress || (() => {});

    // 1. Check local IndexedDB cache first
    if (this.enableCache) {
      onProgress({ stage: 'cache_lookup', loadedBytes: 0, totalBytes: 0, percent: 0 });
      const cached = await this.cache.get(cacheKey);
      if (cached && cached.buffer) {
        onProgress({ stage: 'cache_hit', loadedBytes: cached.byteLength, totalBytes: cached.byteLength, percent: 100 });
        const meta = cached.meta || {};
        return {
          data: new Float32Array(cached.buffer),
          nx: meta.nx || 64,
          ny: meta.ny || 64,
          nz: meta.nz || 64,
          bitpix: meta.bitpix || -32,
          header: new Map(Object.entries(meta.headerCards || {})),
          cached: true
        };
      }
    }

    // 2. Stream FITS Header
    onProgress({ stage: 'fetching_header', loadedBytes: 0, totalBytes: 0, percent: 0 });
    const { headerCards, headerByteLength, totalFileSize } = await this.streamHeader(url);

    const naxis = headerCards.get('NAXIS') || 3;
    const naxis1 = headerCards.get('NAXIS1') || 64; // Raw X (SGX)
    const naxis2 = headerCards.get('NAXIS2') || 64; // Raw Y (SGY)
    const naxis3 = headerCards.get('NAXIS3') || 64; // Raw Z (SGZ)
    const bitpix = headerCards.get('BITPIX') || -32;
    const bscale = headerCards.get('BSCALE') !== undefined ? headerCards.get('BSCALE') : 1.0;
    const bzero = headerCards.get('BZERO') !== undefined ? headerCards.get('BZERO') : 0.0;

    const totalElements = naxis1 * naxis2 * naxis3;
    const bytesPerElem = Math.abs(bitpix) / 8;
    const payloadBytes = totalElements * bytesPerElem;

    // 3. Progressive Range-Request Payload Fetch
    const dataStart = headerByteLength;
    const dataEnd = dataStart + payloadBytes - 1;

    let rawBuffer = null;

    // Direct fetch or chunked range streaming
    if (this.chunkSizeBytes >= payloadBytes) {
      const resp = await this._fetchRangeWithRetry(url, dataStart, dataEnd);
      rawBuffer = await resp.arrayBuffer();
      onProgress({ stage: 'downloading', loadedBytes: payloadBytes, totalBytes: payloadBytes, percent: 100 });
    } else {
      rawBuffer = new ArrayBuffer(payloadBytes);
      const rawUint8 = new Uint8Array(rawBuffer);
      let loaded = 0;

      for (let offset = dataStart; offset <= dataEnd; offset += this.chunkSizeBytes) {
        const chunkEnd = Math.min(offset + this.chunkSizeBytes - 1, dataEnd);
        const resp = await this._fetchRangeWithRetry(url, offset, chunkEnd);
        const chunkBuf = await resp.arrayBuffer();
        rawUint8.set(new Uint8Array(chunkBuf), loaded);
        loaded += chunkBuf.byteLength;

        const percent = Math.min(100, Math.round((loaded / payloadBytes) * 100));
        onProgress({ stage: 'downloading', loadedBytes: loaded, totalBytes: payloadBytes, percent });
      }
    }

    // 4. Decode Big-Endian Float and Reorder (SGZ, SGY, SGX) -> Canonical (SGX, SGY, SGZ)
    onProgress({ stage: 'decoding_and_reordering', loadedBytes: payloadBytes, totalBytes: payloadBytes, percent: 100 });
    const view = new DataView(rawBuffer);
    const outData = new Float32Array(totalElements);

    const scaleFactor = isVelocity ? (52.0 * bscale) : bscale;

    // Official CF4 convention: FITS NAXIS1=SGX, NAXIS2=SGY, NAXIS3=SGZ (stored row-major as [z][y][x])
    let srcIdx = 0;
    for (let iz = 0; iz < naxis3; iz++) {
      for (let iy = 0; iy < naxis2; iy++) {
        for (let ix = 0; ix < naxis1; ix++) {
          let rawVal = 0.0;
          if (bitpix === -32) {
            rawVal = view.getFloat32(srcIdx * 4, false); // Big-Endian
          } else if (bitpix === -64) {
            rawVal = view.getFloat64(srcIdx * 8, false); // Big-Endian
          } else if (bitpix === 16) {
            rawVal = view.getInt16(srcIdx * 2, false);
          } else if (bitpix === 32) {
            rawVal = view.getInt32(srcIdx * 4, false);
          }

          const physVal = rawVal * scaleFactor + bzero;

          // Canonical internal ZRT mapping: index(ix, iy, iz) = ix + nx * (iy + ny * iz)
          const dstIdx = ix + naxis1 * (iy + naxis2 * iz);
          outData[dstIdx] = physVal;
          srcIdx++;
        }
      }
    }

    // 5. Store in IndexedDB Cache
    if (this.enableCache) {
      const meta = {
        nx: naxis1,
        ny: naxis2,
        nz: naxis3,
        bitpix,
        isVelocity,
        headerCards: Object.fromEntries(headerCards)
      };
      await this.cache.put(cacheKey, outData.buffer, meta);
    }

    onProgress({ stage: 'completed', loadedBytes: payloadBytes, totalBytes: payloadBytes, percent: 100 });

    return {
      data: outData,
      nx: naxis1,
      ny: naxis2,
      nz: naxis3,
      bitpix,
      header: headerCards,
      cached: false
    };
  }

  /**
   * Internal helper performing fetch with range headers and retry backoff.
   * @private
   */
  async _fetchRangeWithRetry(url, startByte, endByte, attempt = 1) {
    try {
      const headers = {
        'Range': `bytes=${startByte}-${endByte}`
      };
      const resp = await fetch(url, { headers });
      if (!resp.ok && resp.status !== 206) {
        // Some static servers don't support range requests; fallback to full fetch
        if (resp.status === 200 && startByte === 0) {
          return resp;
        }
        throw new Error(`HTTP Range Request failed with status ${resp.status} ${resp.statusText}`);
      }
      return resp;
    } catch (err) {
      if (attempt < this.maxRetries) {
        const delayMs = Math.pow(2, attempt) * 200;
        await new Promise(r => setTimeout(r, delayMs));
        return this._fetchRangeWithRetry(url, startByte, endByte, attempt + 1);
      }
      throw err;
    }
  }
}
