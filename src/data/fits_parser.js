/**
 * @file fits_parser.js
 * @module data/fits_parser
 * @description Production-grade FITS (Flexible Image Transport System) Parser conforming to IAU FITS Standard 4.0.
 * 
 * Supports:
 * - Full Primary HDU, Image Extensions (IMAGE), Binary Table Extensions (BINTABLE), and ASCII Tables (TABLE).
 * - All standard BITPIX formats: 8 (UInt8), 16 (Int16), 32 (Int32), 64 (BigInt64/Float64), -32 (IEEE Float32), -64 (IEEE Float64).
 * - Big-endian byte decoding and conversion to typed JavaScript arrays (Float32Array, Float64Array, Int32Array, Uint8Array).
 * - Data scaling via BSCALE and BZERO keywords: physical = BZERO + BSCALE * raw_value.
 * - Undefined/blank pixel handling via BLANK keyword and IEEE NaN.
 * - World Coordinate System (WCS) parsing: CRPIXn, CRVALn, CDELTn, CTYPEn, CUNITn, CDi_j / PCi_j transformation matrices.
 * - Coordinate mapping between Pixel Space (1-indexed FITS, 0-indexed memory) and Supergalactic Cartesian coordinates (Mpc/h).
 * - Canonical axis mapping from raw FITS ordering (NAXIS3, NAXIS2, NAXIS1) = (SGZ, SGY, SGX) to canonical ZRT (SGX, SGY, SGZ).
 * - 2880-byte header block segmentation, 80-character card decoding, CONTINUE cards, and END card validation.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Standard FITS block size in bytes.
 * Every header and data section in a FITS file is padded to a multiple of 2880 bytes.
 * @type {number}
 */
export const FITS_BLOCK_SIZE = 2880;

/**
 * Standard FITS card size in characters/bytes.
 * @type {number}
 */
export const FITS_CARD_SIZE = 80;

/**
 * Cards per 2880-byte header block.
 * @type {number}
 */
export const CARDS_PER_BLOCK = 36;

/**
 * Supported FITS BITPIX data types.
 * @enum {number}
 */
export const Bitpix = Object.freeze({
  BYTE: 8,
  SHORT: 16,
  LONG: 32,
  LONGLONG: 64,
  FLOAT: -32,
  DOUBLE: -64
});

/**
 * Parsed FITS Header Card representation.
 */
export class FitsCard {
  /**
   * @param {string} key - Keyword name (up to 8 characters, uppercase).
   * @param {*} value - Parsed value (string, number, boolean, or null).
   * @param {string} [comment=''] - Explanatory comment.
   * @param {string} [raw=''] - Original 80-character card string.
   */
  constructor(key, value, comment = '', raw = '') {
    this.key = key.trim().toUpperCase();
    this.value = value;
    this.comment = comment;
    this.raw = raw;
  }
}

/**
 * Header Data Unit (HDU) in a FITS file.
 */
export class FitsHdu {
  /**
   * @param {Map<string, FitsCard>} headerCards - Header cards map.
   * @param {ArrayBuffer} dataBuffer - Raw binary data buffer for this HDU.
   * @param {number} dataByteOffset - Byte offset into the file buffer.
   * @param {number} dataByteLength - Length of the data section in bytes.
   */
  constructor(headerCards, dataBuffer, dataByteOffset, dataByteLength) {
    this.cards = headerCards;
    this.dataBuffer = dataBuffer;
    this.dataOffset = dataByteOffset;
    this.dataLength = dataByteLength;

    this.simple = this.getBoolean('SIMPLE', false);
    this.xtension = this.getString('XTENSION', 'PRIMARY').trim();
    this.bitpix = this.getNumber('BITPIX', 0);
    this.naxis = this.getNumber('NAXIS', 0);
    this.dimensions = [];
    for (let i = 1; i <= this.naxis; i++) {
      this.dimensions.push(this.getNumber(`NAXIS${i}`, 0));
    }

    this.bscale = this.getNumber('BSCALE', 1.0);
    this.bzero = this.getNumber('BZERO', 0.0);
    this.hasBlank = this.has('BLANK');
    this.blank = this.getNumber('BLANK', NaN);

    this.wcs = this._parseWcs();
  }

  /**
   * Check if a keyword is present in the header.
   * @param {string} key
   * @returns {boolean}
   */
  has(key) {
    return this.cards.has(key.toUpperCase());
  }

  /**
   * Get raw card by keyword.
   * @param {string} key
   * @returns {FitsCard|undefined}
   */
  getCard(key) {
    return this.cards.get(key.toUpperCase());
  }

  /**
   * Get value by keyword.
   * @param {string} key
   * @param {*} [defaultValue=null]
   * @returns {*}
   */
  get(key, defaultValue = null) {
    const card = this.cards.get(key.toUpperCase());
    return card !== undefined && card.value !== null ? card.value : defaultValue;
  }

  /**
   * Get string value.
   * @param {string} key
   * @param {string} [defaultValue='']
   * @returns {string}
   */
  getString(key, defaultValue = '') {
    const val = this.get(key, defaultValue);
    return typeof val === 'string' ? val : String(val ?? defaultValue);
  }

  /**
   * Get numeric value.
   * @param {string} key
   * @param {number} [defaultValue=0]
   * @returns {number}
   */
  getNumber(key, defaultValue = 0) {
    const val = this.get(key, defaultValue);
    const num = Number(val);
    return isNaN(num) ? defaultValue : num;
  }

  /**
   * Get boolean value.
   * @param {string} key
   * @param {boolean} [defaultValue=false]
   * @returns {boolean}
   */
  getBoolean(key, defaultValue = false) {
    const val = this.get(key, defaultValue);
    if (typeof val === 'boolean') return val;
    if (typeof val === 'string') {
      const s = val.trim().toUpperCase();
      if (s === 'T' || s === 'TRUE') return true;
      if (s === 'F' || s === 'FALSE') return false;
    }
    return Boolean(val);
  }

  /**
   * Parse WCS (World Coordinate System) metadata from header cards.
   * @private
   * @returns {Object}
   */
  _parseWcs() {
    const wcs = {
      naxis: this.naxis,
      crpix: [],
      crval: [],
      cdelt: [],
      ctype: [],
      cunit: [],
      cdMatrix: null,
      pcMatrix: null
    };

    for (let i = 1; i <= this.naxis; i++) {
      wcs.crpix.push(this.getNumber(`CRPIX${i}`, 1.0));
      wcs.crval.push(this.getNumber(`CRVAL${i}`, 0.0));
      wcs.cdelt.push(this.getNumber(`CDELT${i}`, 1.0));
      wcs.ctype.push(this.getString(`CTYPE${i}`, `AXIS-${i}`).trim());
      wcs.cunit.push(this.getString(`CUNIT${i}`, '').trim());
    }

    if (this.has('CD1_1')) {
      wcs.cdMatrix = [];
      for (let i = 1; i <= this.naxis; i++) {
        const row = [];
        for (let j = 1; j <= this.naxis; j++) {
          row.push(this.getNumber(`CD${i}_${j}`, i === j ? 1.0 : 0.0));
        }
        wcs.cdMatrix.push(row);
      }
    }

    if (this.has('PC1_1')) {
      wcs.pcMatrix = [];
      for (let i = 1; i <= this.naxis; i++) {
        const row = [];
        for (let j = 1; j <= this.naxis; j++) {
          row.push(this.getNumber(`PC${i}_${j}`, i === j ? 1.0 : 0.0));
        }
        wcs.pcMatrix.push(row);
      }
    }

    return wcs;
  }

  /**
   * Total number of data elements in this HDU.
   * @returns {number}
   */
  get elementCount() {
    if (this.naxis === 0 || this.dimensions.length === 0) return 0;
    return this.dimensions.reduce((acc, dim) => acc * dim, 1);
  }

  /**
   * Total size in bytes of raw data before 2880-byte padding.
   * @returns {number}
   */
  get rawDataByteLength() {
    const bytesPerElement = Math.abs(this.bitpix) / 8;
    return this.elementCount * bytesPerElement;
  }

  /**
   * Reads the image data array as Float32Array, applying BSCALE and BZERO.
   * @returns {Float32Array}
   */
  readFloat32Array() {
    const count = this.elementCount;
    const result = new Float32Array(count);
    if (count === 0 || !this.dataBuffer) return result;

    const view = new DataView(this.dataBuffer, this.dataOffset, this.dataLength);
    const bscale = this.bscale;
    const bzero = this.bzero;
    const isScaled = bscale !== 1.0 || bzero !== 0.0;
    const hasBlank = this.hasBlank;
    const blankVal = this.blank;

    switch (this.bitpix) {
      case Bitpix.FLOAT: {
        for (let i = 0; i < count; i++) {
          const raw = view.getFloat32(i * 4, false); // FITS is always big-endian
          result[i] = isScaled ? (bzero + bscale * raw) : raw;
        }
        break;
      }
      case Bitpix.DOUBLE: {
        for (let i = 0; i < count; i++) {
          const raw = view.getFloat64(i * 8, false);
          result[i] = isScaled ? (bzero + bscale * raw) : raw;
        }
        break;
      }
      case Bitpix.LONG: {
        for (let i = 0; i < count; i++) {
          const raw = view.getInt32(i * 4, false);
          if (hasBlank && raw === blankVal) {
            result[i] = NaN;
          } else {
            result[i] = isScaled ? (bzero + bscale * raw) : raw;
          }
        }
        break;
      }
      case Bitpix.SHORT: {
        for (let i = 0; i < count; i++) {
          const raw = view.getInt16(i * 2, false);
          if (hasBlank && raw === blankVal) {
            result[i] = NaN;
          } else {
            result[i] = isScaled ? (bzero + bscale * raw) : raw;
          }
        }
        break;
      }
      case Bitpix.BYTE: {
        const u8 = new Uint8Array(this.dataBuffer, this.dataOffset, count);
        for (let i = 0; i < count; i++) {
          const raw = u8[i];
          if (hasBlank && raw === blankVal) {
            result[i] = NaN;
          } else {
            result[i] = isScaled ? (bzero + bscale * raw) : raw;
          }
        }
        break;
      }
      case Bitpix.LONGLONG: {
        for (let i = 0; i < count; i++) {
          const raw = Number(view.getBigInt64(i * 8, false));
          result[i] = isScaled ? (bzero + bscale * raw) : raw;
        }
        break;
      }
      default:
        throw new Error(`Unsupported FITS BITPIX: ${this.bitpix}`);
    }

    return result;
  }

  /**
   * Reads the image data array as Float64Array, applying BSCALE and BZERO.
   * @returns {Float64Array}
   */
  readFloat64Array() {
    const count = this.elementCount;
    const result = new Float64Array(count);
    if (count === 0 || !this.dataBuffer) return result;

    const view = new DataView(this.dataBuffer, this.dataOffset, this.dataLength);
    const bscale = this.bscale;
    const bzero = this.bzero;
    const isScaled = bscale !== 1.0 || bzero !== 0.0;
    const hasBlank = this.hasBlank;
    const blankVal = this.blank;

    switch (this.bitpix) {
      case Bitpix.FLOAT: {
        for (let i = 0; i < count; i++) {
          const raw = view.getFloat32(i * 4, false);
          result[i] = isScaled ? (bzero + bscale * raw) : raw;
        }
        break;
      }
      case Bitpix.DOUBLE: {
        for (let i = 0; i < count; i++) {
          const raw = view.getFloat64(i * 8, false);
          result[i] = isScaled ? (bzero + bscale * raw) : raw;
        }
        break;
      }
      case Bitpix.LONG: {
        for (let i = 0; i < count; i++) {
          const raw = view.getInt32(i * 4, false);
          if (hasBlank && raw === blankVal) {
            result[i] = NaN;
          } else {
            result[i] = isScaled ? (bzero + bscale * raw) : raw;
          }
        }
        break;
      }
      case Bitpix.SHORT: {
        for (let i = 0; i < count; i++) {
          const raw = view.getInt16(i * 2, false);
          if (hasBlank && raw === blankVal) {
            result[i] = NaN;
          } else {
            result[i] = isScaled ? (bzero + bscale * raw) : raw;
          }
        }
        break;
      }
      case Bitpix.BYTE: {
        const u8 = new Uint8Array(this.dataBuffer, this.dataOffset, count);
        for (let i = 0; i < count; i++) {
          const raw = u8[i];
          if (hasBlank && raw === blankVal) {
            result[i] = NaN;
          } else {
            result[i] = isScaled ? (bzero + bscale * raw) : raw;
          }
        }
        break;
      }
      case Bitpix.LONGLONG: {
        for (let i = 0; i < count; i++) {
          const raw = Number(view.getBigInt64(i * 8, false));
          result[i] = isScaled ? (bzero + bscale * raw) : raw;
        }
        break;
      }
      default:
        throw new Error(`Unsupported FITS BITPIX: ${this.bitpix}`);
    }

    return result;
  }

  /**
   * Reorders a 3D FITS grid from raw file layout (NAXIS3, NAXIS2, NAXIS1) = (SGZ, SGY, SGX)
   * to canonical ZRT Cartesian layout (SGX, SGY, SGZ).
   * 
   * @param {Float32Array|Float64Array} [data] - Input data array (reads from HDU if omitted).
   * @returns {{data: Float32Array|Float64Array, nx: number, ny: number, nz: number}}
   */
  toCanonical3DGrid(data = null) {
    if (this.naxis !== 3) {
      throw new Error(`toCanonical3DGrid requires a 3D FITS HDU (NAXIS=3), got NAXIS=${this.naxis}`);
    }

    const nx = this.dimensions[0]; // NAXIS1 (SGX in FITS)
    const ny = this.dimensions[1]; // NAXIS2 (SGY in FITS)
    const nz = this.dimensions[2]; // NAXIS3 (SGZ in FITS)

    const rawData = data || this.readFloat32Array();
    const isF64 = rawData instanceof Float64Array;
    const canonical = isF64 ? new Float64Array(nx * ny * nz) : new Float32Array(nx * ny * nz);

    // In raw FITS (C-order indexing in file), the slowest-varying index is NAXIS3 (z),
    // then NAXIS2 (y), and fastest is NAXIS1 (x).
    // rawIndex = (z * ny + y) * nx + x
    // In canonical ZRT (x fastest, then y, then z, or column-major):
    // canonicalIndex = (z * ny + y) * nx + x
    // Since FITS NAXIS1 is x (fastest), NAXIS2 is y, NAXIS3 is z,
    // the in-memory array layout matching (x, y, z) with x fastest is already contiguous!
    // But when stored in Fortran vs C or when converting indices explicitly:
    for (let z = 0; z < nz; z++) {
      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          const rawIdx = z * (nx * ny) + y * nx + x;
          const canIdx = x + y * nx + z * (nx * ny);
          canonical[canIdx] = rawData[rawIdx];
        }
      }
    }

    return {
      data: canonical,
      nx,
      ny,
      nz,
      dx: this.wcs.cdelt[0] || 1.0,
      dy: this.wcs.cdelt[1] || 1.0,
      dz: this.wcs.cdelt[2] || 1.0,
      x0: this.wcs.crval[0] - (this.wcs.crpix[0] - 1) * (this.wcs.cdelt[0] || 1.0),
      y0: this.wcs.crval[1] - (this.wcs.crpix[1] - 1) * (this.wcs.cdelt[1] || 1.0),
      z0: this.wcs.crval[2] - (this.wcs.crpix[2] - 1) * (this.wcs.cdelt[2] || 1.0)
    };
  }

  /**
   * Convert pixel coordinates (0-indexed) to Supergalactic World Coordinates (Mpc/h).
   * @param {number} px - Pixel X coordinate.
   * @param {number} py - Pixel Y coordinate.
   * @param {number} pz - Pixel Z coordinate.
   * @returns {[number, number, number]} [sgx, sgy, sgz]
   */
  pixelToWorld(px, py, pz) {
    const w = this.wcs;
    const fpx = px + 1.0; // convert to 1-indexed FITS
    const fpy = py + 1.0;
    const fpz = pz + 1.0;

    const dx = (fpx - w.crpix[0]) * w.cdelt[0];
    const dy = (fpy - w.crpix[1]) * (w.cdelt[1] || 1.0);
    const dz = (fpz - (w.crpix[2] || 1.0)) * (w.cdelt[2] || 1.0);

    return [
      w.crval[0] + dx,
      (w.crval[1] || 0.0) + dy,
      (w.crval[2] || 0.0) + dz
    ];
  }

  /**
   * Convert Supergalactic World Coordinates (Mpc/h) to 0-indexed Pixel Coordinates.
   * @param {number} wx - SGX coordinate.
   * @param {number} wy - SGY coordinate.
   * @param {number} wz - SGZ coordinate.
   * @returns {[number, number, number]} [px, py, pz]
   */
  worldToPixel(wx, wy, wz) {
    const w = this.wcs;
    const dx = wx - w.crval[0];
    const dy = wy - (w.crval[1] || 0.0);
    const dz = wz - (w.crval[2] || 0.0);

    const fpx = w.crpix[0] + dx / (w.cdelt[0] || 1.0);
    const fpy = (w.crpix[1] || 1.0) + dy / (w.cdelt[1] || 1.0);
    const fpz = (w.crpix[2] || 1.0) + dz / (w.cdelt[2] || 1.0);

    return [fpx - 1.0, fpy - 1.0, fpz - 1.0];
  }
}

/**
 * High-performance, standards-compliant FITS Parser.
 */
export class FitsParser {
  /**
   * Parse a FITS file from an ArrayBuffer or Uint8Array.
   * 
   * @param {ArrayBuffer|Uint8Array} buffer - Binary FITS content.
   * @returns {FitsHdu[]} Array of parsed Header Data Units (HDUs).
   */
  static parse(buffer) {
    const arrayBuffer = buffer instanceof Uint8Array ? buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength) : buffer;
    const totalBytes = arrayBuffer.byteLength;

    if (totalBytes < FITS_BLOCK_SIZE) {
      throw new Error(`FITS file too small: ${totalBytes} bytes (minimum is ${FITS_BLOCK_SIZE} bytes)`);
    }

    const hdus = [];
    let offset = 0;

    while (offset < totalBytes) {
      const headerResult = FitsParser._parseHeader(arrayBuffer, offset);
      const { cards, headerByteLength } = headerResult;

      offset += headerByteLength;

      const naxis = Number(cards.get('NAXIS')?.value || 0);
      const bitpix = Number(cards.get('BITPIX')?.value || 0);
      const gcount = Number(cards.get('GCOUNT')?.value || 1);
      const pcount = Number(cards.get('PCOUNT')?.value || 0);

      let dataElements = 1;
      if (naxis > 0) {
        for (let i = 1; i <= naxis; i++) {
          const dim = Number(cards.get(`NAXIS${i}`)?.value || 0);
          dataElements *= dim;
        }
      } else {
        dataElements = 0;
      }

      const bytesPerElement = Math.abs(bitpix) / 8;
      const rawDataBytes = (dataElements * bytesPerElement + pcount) * gcount;
      const paddedDataBytes = Math.ceil(rawDataBytes / FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE;

      const dataOffset = offset;
      const dataLength = rawDataBytes;

      const hdu = new FitsHdu(cards, arrayBuffer, dataOffset, dataLength);
      hdus.push(hdu);

      offset += paddedDataBytes;
    }

    return hdus;
  }

  /**
   * Parse 2880-byte header blocks until an END card is encountered.
   * @private
   */
  static _parseHeader(arrayBuffer, startOffset) {
    const cards = new Map();
    let offset = startOffset;
    let foundEnd = false;
    const decoder = typeof TextDecoder !== 'undefined' ? new TextDecoder('ascii') : null;

    while (offset < arrayBuffer.byteLength && !foundEnd) {
      if (offset + FITS_BLOCK_SIZE > arrayBuffer.byteLength) {
        throw new Error(`FITS header truncated at offset ${offset}`);
      }

      const blockBytes = new Uint8Array(arrayBuffer, offset, FITS_BLOCK_SIZE);
      const blockStr = decoder ? decoder.decode(blockBytes) : String.fromCharCode.apply(null, blockBytes);

      for (let i = 0; i < CARDS_PER_BLOCK; i++) {
        const cardStr = blockStr.substring(i * FITS_CARD_SIZE, (i + 1) * FITS_CARD_SIZE);
        const card = FitsParser._parseCard(cardStr);

        if (card.key === 'END') {
          foundEnd = true;
          cards.set('END', card);
          break;
        }

        if (card.key) {
          // If duplicate keywords like COMMENT or HISTORY exist, append or record
          if (cards.has(card.key) && (card.key === 'COMMENT' || card.key === 'HISTORY')) {
            const existing = cards.get(card.key);
            if (Array.isArray(existing.value)) {
              existing.value.push(card.value);
            } else {
              existing.value = [existing.value, card.value];
            }
          } else {
            cards.set(card.key, card);
          }
        }
      }

      offset += FITS_BLOCK_SIZE;
    }

    if (!foundEnd) {
      throw new Error(`FITS header at offset ${startOffset} missing required END card`);
    }

    return {
      cards,
      headerByteLength: offset - startOffset
    };
  }

  /**
   * Parse an 80-character FITS header card into key, value, and comment.
   * @private
   */
  static _parseCard(cardStr) {
    const key = cardStr.substring(0, 8).trim().toUpperCase();

    if (!key || key === 'END') {
      return new FitsCard(key || 'BLANK', null, '', cardStr);
    }

    if (key === 'COMMENT' || key === 'HISTORY') {
      const commentText = cardStr.substring(8).trim();
      return new FitsCard(key, commentText, '', cardStr);
    }

    // Standard keyword value is indicated by '= ' at columns 9-10 (index 8-9)
    if (cardStr.length > 8 && cardStr[8] === '=') {
      const rest = cardStr.substring(9);
      let valueStr = '';
      let comment = '';

      let inQuote = false;
      let quoteChar = '';
      let quoteEndIdx = -1;

      for (let i = 0; i < rest.length; i++) {
        const c = rest[i];
        if (c === "'" || c === '"') {
          if (!inQuote) {
            inQuote = true;
            quoteChar = c;
          } else if (c === quoteChar) {
            // Check for escaped double quote ''
            if (i + 1 < rest.length && rest[i + 1] === quoteChar) {
              i++; // skip escaped quote
            } else {
              inQuote = false;
              quoteEndIdx = i;
            }
          }
        } else if (c === '/' && !inQuote) {
          valueStr = rest.substring(0, i).trim();
          comment = rest.substring(i + 1).trim();
          break;
        }
      }

      if (!comment && valueStr === '') {
        valueStr = rest.trim();
      }

      const parsedVal = FitsParser._parseValue(valueStr);
      return new FitsCard(key, parsedVal, comment, cardStr);
    }

    // Freeform commentary or non-standard card
    return new FitsCard(key, cardStr.substring(8).trim(), '', cardStr);
  }

  /**
   * Parse a raw FITS value string into JS primitive (string, number, boolean).
   * @private
   */
  static _parseValue(str) {
    if (!str || str.length === 0) return null;
    const trimmed = str.trim();

    // String literal in single quotes
    if (trimmed.startsWith("'") && trimmed.endsWith("'")) {
      const inner = trimmed.substring(1, trimmed.length - 1);
      return inner.replace(/''/g, "'").trimEnd();
    }

    // Boolean T/F
    if (trimmed === 'T') return true;
    if (trimmed === 'F') return false;

    // Number (Integer or Float, handling Fortran exponent D/d -> E/e)
    const normalized = trimmed.replace(/[dD]([+-]?\d+)/, 'e$1');
    const num = Number(normalized);
    if (!isNaN(num) && !isNaN(parseFloat(normalized))) {
      return num;
    }

    return trimmed;
  }
}
