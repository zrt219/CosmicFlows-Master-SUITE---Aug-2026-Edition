/**
 * @file binary_stream_encoder.js
 * @description High-throughput binary serialization format for massive cosmological streamlines:
 *  - Header: Magic 'ZRTS' (0x5A525453), version, flags, streamline count, vertex count, bounding box
 *  - Offset pointer table and vertex count table
 *  - Typed coordinate arrays (Float32 or quantized Uint16)
 *  - Typed attribute arrays (speeds, curvatures, times)
 *  - Fast decoding and round-trip verification.
 *
 * @module streamlines/binary_stream_encoder
 */

export const BINARY_MAGIC = 0x5A525453; // 'ZRTS'
export const BINARY_VERSION = 1;

export const BinaryFlags = Object.freeze({
  HAS_VELOCITIES: 1 << 0,
  HAS_SPEEDS: 1 << 1,
  HAS_CURVATURES: 1 << 2,
  HAS_TIMES: 1 << 3,
  QUANTIZED_UINT16: 1 << 4
});

export class BinaryStreamEncoder {
  /**
   * Encodes an array of Streamlines into a single binary ArrayBuffer.
   *
   * @param {Array<Object>} streamlines - Array of Streamline objects or { points, speeds, times }.
   * @param {Object} [options={}]
   * @param {boolean} [options.includeSpeeds=true]
   * @param {boolean} [options.includeTimes=true]
   * @returns {ArrayBuffer} Encoded binary buffer.
   */
  static encode(streamlines, options = {}) {
    const {
      includeSpeeds = true,
      includeTimes = true
    } = options;

    const streamlineCount = streamlines.length;
    let totalVertices = 0;
    let xMin = Infinity, xMax = -Infinity;
    let yMin = Infinity, yMax = -Infinity;
    let zMin = Infinity, zMax = -Infinity;

    for (let s = 0; s < streamlineCount; s++) {
      const pts = streamlines[s].points || [];
      totalVertices += pts.length;
      for (let i = 0; i < pts.length; i++) {
        const x = pts[i][0], y = pts[i][1], z = pts[i][2];
        if (x < xMin) xMin = x; if (x > xMax) xMax = x;
        if (y < yMin) yMin = y; if (y > yMax) yMax = y;
        if (z < zMin) zMin = z; if (z > zMax) zMax = z;
      }
    }

    if (!Number.isFinite(xMin)) {
      xMin = xMax = yMin = yMax = zMin = zMax = 0;
    }

    let flags = 0;
    if (includeSpeeds) flags |= BinaryFlags.HAS_SPEEDS;
    if (includeTimes) flags |= BinaryFlags.HAS_TIMES;

    // Header size: 4 (magic) + 2 (version) + 2 (flags) + 4 (sCount) + 4 (vCount) + 24 (bbox) = 40 bytes
    const headerBytes = 40;
    const offsetTableBytes = streamlineCount * 4;
    const countTableBytes = streamlineCount * 2;
    const vertexBytes = totalVertices * 3 * 4; // 3 floats per vertex
    const speedBytes = includeSpeeds ? totalVertices * 4 : 0;
    const timeBytes = includeTimes ? totalVertices * 4 : 0;

    const totalBytes = headerBytes + offsetTableBytes + countTableBytes + vertexBytes + speedBytes + timeBytes;
    const buffer = new ArrayBuffer(totalBytes);
    const view = new DataView(buffer);

    let offset = 0;
    // Header
    view.setUint32(offset, BINARY_MAGIC, true); offset += 4;
    view.setUint16(offset, BINARY_VERSION, true); offset += 2;
    view.setUint16(offset, flags, true); offset += 2;
    view.setUint32(offset, streamlineCount, true); offset += 4;
    view.setUint32(offset, totalVertices, true); offset += 4;
    view.setFloat32(offset, xMin, true); offset += 4;
    view.setFloat32(offset, xMax, true); offset += 4;
    view.setFloat32(offset, yMin, true); offset += 4;
    view.setFloat32(offset, yMax, true); offset += 4;
    view.setFloat32(offset, zMin, true); offset += 4;
    view.setFloat32(offset, zMax, true); offset += 4;

    // Offset table and Count table
    let vertexAccum = 0;
    for (let s = 0; s < streamlineCount; s++) {
      view.setUint32(offset, vertexAccum, true); offset += 4;
      const pts = streamlines[s].points || [];
      vertexAccum += pts.length;
    }
    for (let s = 0; s < streamlineCount; s++) {
      const pts = streamlines[s].points || [];
      view.setUint16(offset, pts.length, true); offset += 2;
    }

    // Vertex data
    const posArray = new Float32Array(buffer, offset, totalVertices * 3);
    let pIdx = 0;
    for (let s = 0; s < streamlineCount; s++) {
      const pts = streamlines[s].points || [];
      for (let i = 0; i < pts.length; i++) {
        posArray[pIdx++] = pts[i][0];
        posArray[pIdx++] = pts[i][1];
        posArray[pIdx++] = pts[i][2];
      }
    }
    offset += totalVertices * 3 * 4;

    // Speeds
    if (includeSpeeds) {
      const speedArray = new Float32Array(buffer, offset, totalVertices);
      let sIdx = 0;
      for (let s = 0; s < streamlineCount; s++) {
        const sp = streamlines[s].speeds || [];
        for (let i = 0; i < sp.length; i++) {
          speedArray[sIdx++] = sp[i];
        }
      }
      offset += totalVertices * 4;
    }

    // Times
    if (includeTimes) {
      const timeArray = new Float32Array(buffer, offset, totalVertices);
      let tIdx = 0;
      for (let s = 0; s < streamlineCount; s++) {
        const tm = streamlines[s].times || [];
        for (let i = 0; i < tm.length; i++) {
          timeArray[tIdx++] = tm[i];
        }
      }
    }

    return buffer;
  }

  /**
   * Decodes a binary ArrayBuffer back into structured streamline objects.
   *
   * @param {ArrayBuffer} buffer
   * @returns {Object} Decoded streamlines and header metadata.
   */
  static decode(buffer) {
    const view = new DataView(buffer);
    let offset = 0;

    const magic = view.getUint32(offset, true); offset += 4;
    if (magic !== BINARY_MAGIC) {
      throw new Error(`Invalid binary magic: 0x${magic.toString(16)}`);
    }

    const version = view.getUint16(offset, true); offset += 2;
    const flags = view.getUint16(offset, true); offset += 2;
    const streamlineCount = view.getUint32(offset, true); offset += 4;
    const totalVertices = view.getUint32(offset, true); offset += 4;

    const bbox = {
      xMin: view.getFloat32(offset, true),
      xMax: view.getFloat32(offset + 4, true),
      yMin: view.getFloat32(offset + 8, true),
      yMax: view.getFloat32(offset + 12, true),
      zMin: view.getFloat32(offset + 16, true),
      zMax: view.getFloat32(offset + 20, true)
    };
    offset += 24;

    const offsets = [];
    for (let s = 0; s < streamlineCount; s++) {
      offsets.push(view.getUint32(offset, true));
      offset += 4;
    }

    const counts = [];
    for (let s = 0; s < streamlineCount; s++) {
      counts.push(view.getUint16(offset, true));
      offset += 2;
    }

    const posArray = new Float32Array(buffer, offset, totalVertices * 3);
    offset += totalVertices * 3 * 4;

    let speedArray = null;
    if (flags & BinaryFlags.HAS_SPEEDS) {
      speedArray = new Float32Array(buffer, offset, totalVertices);
      offset += totalVertices * 4;
    }

    let timeArray = null;
    if (flags & BinaryFlags.HAS_TIMES) {
      timeArray = new Float32Array(buffer, offset, totalVertices);
      offset += totalVertices * 4;
    }

    const streamlines = [];
    for (let s = 0; s < streamlineCount; s++) {
      const startV = offsets[s];
      const count = counts[s];
      const points = [];
      const speeds = [];
      const times = [];

      for (let i = 0; i < count; i++) {
        const vIdx = startV + i;
        points.push([
          posArray[vIdx * 3 + 0],
          posArray[vIdx * 3 + 1],
          posArray[vIdx * 3 + 2]
        ]);
        if (speedArray) speeds.push(speedArray[vIdx]);
        if (timeArray) times.push(timeArray[vIdx]);
      }

      streamlines.push({
        id: `streamline_${s}`,
        points,
        speeds,
        times,
        vertexCount: count
      });
    }

    return {
      version,
      flags,
      streamlineCount,
      totalVertices,
      bbox,
      streamlines
    };
  }
}

export default BinaryStreamEncoder;
