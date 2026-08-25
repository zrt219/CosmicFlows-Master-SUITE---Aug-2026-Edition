/**
 * @file data_exporters.js
 * @module export/data_exporters
 * @description Research-grade format converters and serializers for ZRT Cosmicflows Workbench.
 * Supports:
 * - Basin Catalogs (CSV, ECSV with YAML headers)
 * - Streamlines (GeoJSON LineString/MultiLineString, ZRT High-Performance Compact Binary, Structured JSON)
 * - Topology Graphs (Standard XML GraphML, Node-Link JSON)
 * - Cosmological Bulk Flows (Multi-Shell CSV, JSON, VOTable XML)
 *
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Basin Catalog Exporter (CSV / ECSV / VO-Table compatible).
 */
export class BasinCatalogExporter {
  /**
   * Export basin catalog to RFC 4180 CSV with optional Astropy ECSV YAML metadata header.
   * @param {Array<Object>} basins Array of basin descriptors
   * @param {Object} [options]
   * @param {boolean} [options.includeEcsvHeader=false] If true, prepends Astropy ECSV YAML header
   * @param {string} [options.delimiter=','] Delimiter character
   * @param {number} [options.precision=6] Decimal digits for floating point values
   * @returns {string} Formatted CSV text
   */
  static exportToCsv(basins, options = {}) {
    const delimiter = options.delimiter || ',';
    const prec = options.precision !== undefined ? options.precision : 6;
    const includeEcsv = options.includeEcsvHeader || false;

    const headers = [
      'basin_id',
      'attractor_name',
      'x_mpc_h',
      'y_mpc_h',
      'z_mpc_h',
      'vx_kms',
      'vy_kms',
      'vz_kms',
      'volume_mpc3_h3',
      'mass_proxy_1e14msun',
      'galaxy_count',
      'percolation_threshold',
      'stability_index',
      'bbox_xmin',
      'bbox_xmax',
      'bbox_ymin',
      'bbox_ymax',
      'bbox_zmin',
      'bbox_zmax'
    ];

    const lines = [];

    if (includeEcsv) {
      lines.push('# %ECSV 1.0');
      lines.push('# ---');
      lines.push('# datatype:');
      lines.push('# - {name: basin_id, datatype: int32, description: Unique catchment basin identifier}');
      lines.push('# - {name: attractor_name, datatype: string, description: Identified cosmological attractor}');
      lines.push('# - {name: x_mpc_h, datatype: float64, unit: Mpc/h, description: Attractor SGX coordinate}');
      lines.push('# - {name: y_mpc_h, datatype: float64, unit: Mpc/h, description: Attractor SGY coordinate}');
      lines.push('# - {name: z_mpc_h, datatype: float64, unit: Mpc/h, description: Attractor SGZ coordinate}');
      lines.push('# - {name: vx_kms, datatype: float64, unit: km/s, description: Attractor SGX velocity}');
      lines.push('# - {name: vy_kms, datatype: float64, unit: km/s, description: Attractor SGY velocity}');
      lines.push('# - {name: vz_kms, datatype: float64, unit: km/s, description: Attractor SGZ velocity}');
      lines.push('# - {name: volume_mpc3_h3, datatype: float64, unit: (Mpc/h)^3, description: Basin watershed volume}');
      lines.push('# - {name: mass_proxy_1e14msun, datatype: float64, unit: 10^14 Msun/h, description: Integrated enclosed mass proxy}');
      lines.push('# - {name: galaxy_count, datatype: int32, description: Number of member galaxies}');
      lines.push('# - {name: percolation_threshold, datatype: float64, description: Density/potential percolation threshold}');
      lines.push('# - {name: stability_index, datatype: float64, description: Monte Carlo topological persistence stability}');
      lines.push('# - {name: bbox_xmin, datatype: float64, unit: Mpc/h}');
      lines.push('# - {name: bbox_xmax, datatype: float64, unit: Mpc/h}');
      lines.push('# - {name: bbox_ymin, datatype: float64, unit: Mpc/h}');
      lines.push('# - {name: bbox_ymax, datatype: float64, unit: Mpc/h}');
      lines.push('# - {name: bbox_zmin, datatype: float64, unit: Mpc/h}');
      lines.push('# - {name: bbox_zmax, datatype: float64, unit: Mpc/h}');
      lines.push('# meta: {software: ZRT Cosmicflows Workbench, date: ' + new Date().toISOString() + '}');
    }

    lines.push(headers.join(delimiter));

    for (const b of basins || []) {
      const bbox = b.boundingBox || b.bbox || [-100, 100, -100, 100, -100, 100];
      const pos = b.position || [b.x || 0, b.y || 0, b.z || 0];
      const vel = b.velocity || [b.vx || 0, b.vy || 0, b.vz || 0];

      const row = [
        b.id !== undefined ? b.id : (b.basin_id !== undefined ? b.basin_id : 0),
        `"${(b.name || b.attractor_name || 'Attractor').replace(/"/g, '""')}"`,
        pos[0].toFixed(prec),
        pos[1].toFixed(prec),
        pos[2].toFixed(prec),
        vel[0].toFixed(prec),
        vel[1].toFixed(prec),
        vel[2].toFixed(prec),
        (b.volume || b.volume_mpc3 || 0).toFixed(prec),
        (b.massProxy || b.mass_proxy || 0).toFixed(prec),
        b.galaxyCount !== undefined ? b.galaxyCount : (b.galaxy_count || 0),
        (b.percolationThreshold !== undefined ? b.percolationThreshold : (b.percolation_threshold || -1.0)).toFixed(prec),
        (b.stabilityIndex !== undefined ? b.stabilityIndex : (b.stability_index || 1.0)).toFixed(prec),
        bbox[0].toFixed(prec),
        bbox[1].toFixed(prec),
        bbox[2].toFixed(prec),
        bbox[3].toFixed(prec),
        bbox[4].toFixed(prec),
        bbox[5].toFixed(prec)
      ];
      lines.push(row.join(delimiter));
    }

    return lines.join('\n');
  }

  /**
   * Parses CSV string back into basin objects (supports both plain CSV and ECSV).
   * @param {string} csvText
   * @param {string} [delimiter=',']
   * @returns {Array<Object>}
   */
  static parseCsv(csvText, delimiter = ',') {
    if (!csvText || typeof csvText !== 'string') return [];
    const lines = csvText.split(/\r?\n/).filter(line => line.trim().length > 0 && !line.startsWith('#'));
    if (lines.length < 2) return [];

    const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^"|"$/g, ''));
    const results = [];

    for (let i = 1; i < lines.length; i++) {
      const rawTokens = lines[i].split(delimiter);
      if (rawTokens.length < headers.length) continue;

      const row = {};
      for (let j = 0; j < headers.length; j++) {
        let val = rawTokens[j].trim().replace(/^"|"$/g, '');
        if (!isNaN(Number(val)) && val !== '') {
          row[headers[j]] = Number(val);
        } else {
          row[headers[j]] = val;
        }
      }

      results.push({
        id: row.basin_id,
        name: row.attractor_name,
        position: [row.x_mpc_h, row.y_mpc_h, row.z_mpc_h],
        velocity: [row.vx_kms, row.vy_kms, row.vz_kms],
        volume: row.volume_mpc3_h3,
        massProxy: row.mass_proxy_1e14msun,
        galaxyCount: row.galaxy_count,
        percolationThreshold: row.percolation_threshold,
        stabilityIndex: row.stability_index,
        boundingBox: [row.bbox_xmin, row.bbox_xmax, row.bbox_ymin, row.bbox_ymax, row.bbox_zmin, row.bbox_zmax]
      });
    }

    return results;
  }
}

/**
 * Streamline Exporter & Binary Serializer.
 */
export class StreamlineExporter {
  /**
   * Exports streamlines as standard RFC 7946 GeoJSON FeatureCollection.
   * @param {Array<Object>} streamlines Array of streamline objects
   * @param {Object} [options]
   * @param {string} [options.coordinateSystem='cartesian_supergalactic'] 'cartesian_supergalactic' or 'spherical_equatorial'
   * @returns {Object} Valid GeoJSON FeatureCollection object
   */
  static exportToGeoJson(streamlines, options = {}) {
    const coordSys = options.coordinateSystem || 'cartesian_supergalactic';
    const features = [];

    for (let idx = 0; idx < (streamlines || []).length; idx++) {
      const st = streamlines[idx];
      const points = st.points || st.vertices || [];
      if (points.length < 2) continue;

      const coords = [];
      let totalArcLength = 0;
      let sumV = 0;
      let maxV = 0;

      for (let p = 0; p < points.length; p++) {
        const pt = points[p];
        const x = pt[0] !== undefined ? pt[0] : pt.x;
        const y = pt[1] !== undefined ? pt[1] : pt.y;
        const z = pt[2] !== undefined ? pt[2] : pt.z;
        const vx = pt[3] !== undefined ? pt[3] : (pt.vx || 0);
        const vy = pt[4] !== undefined ? pt[4] : (pt.vy || 0);
        const vz = pt[5] !== undefined ? pt[5] : (pt.vz || 0);

        const vMag = Math.sqrt(vx * vx + vy * vy + vz * vz);
        sumV += vMag;
        if (vMag > maxV) maxV = vMag;

        if (p > 0) {
          const prev = points[p - 1];
          const px = prev[0] !== undefined ? prev[0] : prev.x;
          const py = prev[1] !== undefined ? prev[1] : prev.y;
          const pz = prev[2] !== undefined ? prev[2] : prev.z;
          const dx = x - px, dy = y - py, dz = z - pz;
          totalArcLength += Math.sqrt(dx * dx + dy * dy + dz * dz);
        }

        if (coordSys === 'spherical_equatorial') {
          const dist = Math.sqrt(x * x + y * y + z * z);
          const ra = (Math.atan2(y, x) * 180 / Math.PI + 360) % 360;
          const dec = dist > 1e-9 ? Math.asin(Math.max(-1, Math.min(1, z / dist))) * 180 / Math.PI : 0;
          coords.push([Number(ra.toFixed(6)), Number(dec.toFixed(6)), Number(dist.toFixed(4))]);
        } else {
          coords.push([Number(x.toFixed(4)), Number(y.toFixed(4)), Number(z.toFixed(4))]);
        }
      }

      const meanV = points.length > 0 ? sumV / points.length : 0;

      features.push({
        type: 'Feature',
        id: st.id !== undefined ? st.id : `streamline_${idx}`,
        geometry: {
          type: 'LineString',
          coordinates: coords
        },
        properties: {
          line_id: st.id !== undefined ? st.id : idx,
          basin_id: st.basinId !== undefined ? st.basinId : (st.basin_id !== undefined ? st.basin_id : -1),
          attractor_name: st.attractorName || st.attractor_name || 'Unknown',
          step_count: points.length,
          arc_length_mpc_h: Number((st.arcLength || totalArcLength).toFixed(4)),
          mean_velocity_kms: Number(meanV.toFixed(2)),
          max_velocity_kms: Number(maxV.toFixed(2)),
          seed_point: st.seedPoint || [coords[0][0], coords[0][1], coords[0][2]],
          termination_reason: st.terminationReason || 'reached_attractor'
        }
      });
    }

    return {
      type: 'FeatureCollection',
      name: 'Cosmicflows_Topological_Streamlines',
      crs: {
        type: 'name',
        properties: {
          name: coordSys === 'spherical_equatorial' ? 'urn:ogc:def:crs:OGC:1.3:CRS84' : 'urn:zrt:crs:SupergalacticCartesian'
        }
      },
      features
    };
  }

  /**
   * Exports streamlines into a high-performance compact binary buffer (ZRT-BIN v1 format).
   * Format layout:
   * [12 bytes]: Magic Header "ZRT_STRM_v1\0\0"
   * [4 bytes]: uint32 Version (1)
   * [4 bytes]: uint32 NumStreamlines (N)
   * [4 bytes]: uint32 TotalPoints (P)
   * [24 bytes]: float32 Bounds [xmin, xmax, ymin, ymax, zmin, zmax]
   * [N * 16 bytes]: Index Table [offset (uint32), count (uint32), basinId (int32), arcLength (float32)]
   * [P * 28 bytes]: Vertex Payload Float32 interleaved [x, y, z, vx, vy, vz, t]
   *
   * @param {Array<Object>} streamlines
   * @returns {ArrayBuffer}
   */
  static exportToCompactBinary(streamlines = []) {
    let totalPoints = 0;
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;

    for (const st of streamlines) {
      const pts = st.points || st.vertices || [];
      totalPoints += pts.length;
      for (const pt of pts) {
        const x = pt[0] !== undefined ? pt[0] : pt.x;
        const y = pt[1] !== undefined ? pt[1] : pt.y;
        const z = pt[2] !== undefined ? pt[2] : pt.z;
        if (x < minX) minX = x; if (x > maxX) maxX = x;
        if (y < minY) minY = y; if (y > maxY) maxY = y;
        if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
      }
    }

    if (minX === Infinity) {
      minX = maxX = minY = maxY = minZ = maxZ = 0;
    }

    const numLines = streamlines.length;
    const headerSize = 12 + 4 + 4 + 4 + 24; // 48 bytes
    const indexTableSize = numLines * 16;
    const payloadSize = totalPoints * 28; // 7 floats * 4 bytes
    const totalByteLength = headerSize + indexTableSize + payloadSize;

    const buffer = new ArrayBuffer(totalByteLength);
    const view = new DataView(buffer);
    const uint8 = new Uint8Array(buffer);

    // 1. Magic Header "ZRT_STRM_v1\0\0"
    const magic = 'ZRT_STRM_v1\0\0';
    for (let i = 0; i < 12; i++) {
      uint8[i] = magic.charCodeAt(i) || 0;
    }

    // 2. Header Metadata
    view.setUint32(12, 1, true); // version
    view.setUint32(16, numLines, true);
    view.setUint32(20, totalPoints, true);
    view.setFloat32(24, minX, true);
    view.setFloat32(28, maxX, true);
    view.setFloat32(32, minY, true);
    view.setFloat32(36, maxY, true);
    view.setFloat32(40, minZ, true);
    view.setFloat32(44, maxZ, true);

    // 3. Populate Index Table and Vertex Payload
    let pointOffset = 0;
    let indexByteOffset = headerSize;
    let payloadByteOffset = headerSize + indexTableSize;

    for (let i = 0; i < numLines; i++) {
      const st = streamlines[i];
      const pts = st.points || st.vertices || [];
      const count = pts.length;
      const basinId = st.basinId !== undefined ? st.basinId : (st.basin_id !== undefined ? st.basin_id : -1);
      const arcLength = st.arcLength || 0.0;

      view.setUint32(indexByteOffset, pointOffset, true);
      view.setUint32(indexByteOffset + 4, count, true);
      view.setInt32(indexByteOffset + 8, basinId, true);
      view.setFloat32(indexByteOffset + 12, arcLength, true);
      indexByteOffset += 16;

      for (let p = 0; p < count; p++) {
        const pt = pts[p];
        const x = pt[0] !== undefined ? pt[0] : pt.x;
        const y = pt[1] !== undefined ? pt[1] : pt.y;
        const z = pt[2] !== undefined ? pt[2] : pt.z;
        const vx = pt[3] !== undefined ? pt[3] : (pt.vx || 0);
        const vy = pt[4] !== undefined ? pt[4] : (pt.vy || 0);
        const vz = pt[5] !== undefined ? pt[5] : (pt.vz || 0);
        const t = pt[6] !== undefined ? pt[6] : (pt.t || (p * 0.25));

        view.setFloat32(payloadByteOffset, x, true);
        view.setFloat32(payloadByteOffset + 4, y, true);
        view.setFloat32(payloadByteOffset + 8, z, true);
        view.setFloat32(payloadByteOffset + 12, vx, true);
        view.setFloat32(payloadByteOffset + 16, vy, true);
        view.setFloat32(payloadByteOffset + 20, vz, true);
        view.setFloat32(payloadByteOffset + 24, t, true);
        payloadByteOffset += 28;
      }

      pointOffset += count;
    }

    return buffer;
  }

  /**
   * Deserializes ZRT-BIN format buffer back into structured streamline objects.
   * @param {ArrayBuffer} buffer
   * @returns {Object} { version, bounds, streamlines: Array<Object> }
   */
  static importFromCompactBinary(buffer) {
    if (!buffer || buffer.byteLength < 48) {
      throw new Error('Invalid binary buffer: buffer is null or under minimum header length (48 bytes).');
    }

    const uint8 = new Uint8Array(buffer);
    const magicStr = String.fromCharCode(...uint8.slice(0, 10));
    if (!magicStr.startsWith('ZRT_STRM')) {
      throw new Error(`Magic header mismatch: expected ZRT_STRM, got ${magicStr}`);
    }

    const view = new DataView(buffer);
    const version = view.getUint32(12, true);
    const numLines = view.getUint32(16, true);
    const totalPoints = view.getUint32(20, true);
    const bounds = [
      view.getFloat32(24, true),
      view.getFloat32(28, true),
      view.getFloat32(32, true),
      view.getFloat32(36, true),
      view.getFloat32(40, true),
      view.getFloat32(44, true)
    ];

    const headerSize = 48;
    const indexTableSize = numLines * 16;
    const streamlines = [];

    for (let i = 0; i < numLines; i++) {
      const idxOffset = headerSize + i * 16;
      const pointOffset = view.getUint32(idxOffset, true);
      const pointCount = view.getUint32(idxOffset + 4, true);
      const basinId = view.getInt32(idxOffset + 8, true);
      const arcLength = view.getFloat32(idxOffset + 12, true);

      const points = [];
      let payloadOffset = headerSize + indexTableSize + pointOffset * 28;

      for (let p = 0; p < pointCount; p++) {
        const x = view.getFloat32(payloadOffset, true);
        const y = view.getFloat32(payloadOffset + 4, true);
        const z = view.getFloat32(payloadOffset + 8, true);
        const vx = view.getFloat32(payloadOffset + 12, true);
        const vy = view.getFloat32(payloadOffset + 16, true);
        const vz = view.getFloat32(payloadOffset + 20, true);
        const t = view.getFloat32(payloadOffset + 24, true);
        payloadOffset += 28;

        points.push([x, y, z, vx, vy, vz, t]);
      }

      streamlines.push({
        id: i,
        basinId,
        arcLength,
        points
      });
    }

    return {
      version,
      totalPoints,
      bounds,
      streamlines
    };
  }
}

/**
 * Topology Graph Exporter (GraphML XML & Node-Link JSON).
 */
export class TopologyGraphExporter {
  /**
   * Exports critical points and separatrix manifolds to standard GraphML XML format.
   * @param {Object} graphData
   * @param {Array<Object>} graphData.criticalPoints Morse critical points (minima, 1-saddles, 2-saddles, maxima)
   * @param {Array<Object>} graphData.separatrices Topological manifolds / edges connecting critical points
   * @returns {string} Standard GraphML XML string
   */
  static exportToGraphML(graphData = {}) {
    const nodes = graphData.criticalPoints || graphData.nodes || [];
    const edges = graphData.separatrices || graphData.edges || [];

    const lines = [
      '<?xml version="1.0" encoding="UTF-8"?>',
      '<graphml xmlns="http://graphml.graphdrawing.org/xmlns"',
      '         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"',
      '         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">',
      '  <!-- Node Key Definitions -->',
      '  <key id="d_name" for="node" attr.name="name" attr.type="string"/>',
      '  <key id="d_type" for="node" attr.name="type" attr.type="string"/>',
      '  <key id="d_morse_index" for="node" attr.name="morse_index" attr.type="int"/>',
      '  <key id="d_x" for="node" attr.name="x" attr.type="double"/>',
      '  <key id="d_y" for="node" attr.name="y" attr.type="double"/>',
      '  <key id="d_z" for="node" attr.name="z" attr.type="double"/>',
      '  <key id="d_potential" for="node" attr.name="potential" attr.type="double"/>',
      '  <key id="d_vmag" for="node" attr.name="velocity_magnitude" attr.type="double"/>',
      '  <key id="d_lambda1" for="node" attr.name="eigenvalue_1" attr.type="double"/>',
      '  <key id="d_lambda2" for="node" attr.name="eigenvalue_2" attr.type="double"/>',
      '  <key id="d_lambda3" for="node" attr.name="eigenvalue_3" attr.type="double"/>',
      '  <!-- Edge Key Definitions -->',
      '  <key id="e_type" for="edge" attr.name="type" attr.type="string"/>',
      '  <key id="e_weight" for="edge" attr.name="weight" attr.type="double"/>',
      '  <key id="e_length" for="edge" attr.name="length" attr.type="double"/>',
      '  <key id="e_flux" for="edge" attr.name="flux" attr.type="double"/>',
      '  <graph id="CosmologicalTopologyGraph" edgedefault="directed">'
    ];

    for (let i = 0; i < nodes.length; i++) {
      const n = nodes[i];
      const id = n.id !== undefined ? String(n.id) : `node_${i}`;
      const name = (n.name || id).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      const type = n.type || (n.morseIndex === 0 ? 'attractor' : (n.morseIndex === 3 ? 'repeller' : 'saddle'));
      const morse = n.morseIndex !== undefined ? n.morseIndex : (n.morse_index !== undefined ? n.morse_index : 0);
      const pos = n.position || [n.x || 0, n.y || 0, n.z || 0];
      const pot = n.potential !== undefined ? n.potential : 0.0;
      const vmag = n.velocityMagnitude !== undefined ? n.velocityMagnitude : 0.0;
      const eigs = n.eigenvalues || [0.0, 0.0, 0.0];

      lines.push(`    <node id="${id}">`);
      lines.push(`      <data key="d_name">${name}</data>`);
      lines.push(`      <data key="d_type">${type}</data>`);
      lines.push(`      <data key="d_morse_index">${morse}</data>`);
      lines.push(`      <data key="d_x">${pos[0]}</data>`);
      lines.push(`      <data key="d_y">${pos[1]}</data>`);
      lines.push(`      <data key="d_z">${pos[2]}</data>`);
      lines.push(`      <data key="d_potential">${pot}</data>`);
      lines.push(`      <data key="d_vmag">${vmag}</data>`);
      lines.push(`      <data key="d_lambda1">${eigs[0]}</data>`);
      lines.push(`      <data key="d_lambda2">${eigs[1]}</data>`);
      lines.push(`      <data key="d_lambda3">${eigs[2]}</data>`);
      lines.push('    </node>');
    }

    for (let j = 0; j < edges.length; j++) {
      const e = edges[j];
      const id = e.id !== undefined ? String(e.id) : `edge_${j}`;
      const source = String(e.source);
      const target = String(e.target);
      const type = e.type || 'separatrix';
      const weight = e.weight !== undefined ? e.weight : 1.0;
      const length = e.length !== undefined ? e.length : 0.0;
      const flux = e.flux !== undefined ? e.flux : 0.0;

      lines.push(`    <edge id="${id}" source="${source}" target="${target}">`);
      lines.push(`      <data key="e_type">${type}</data>`);
      lines.push(`      <data key="e_weight">${weight}</data>`);
      lines.push(`      <data key="e_length">${length}</data>`);
      lines.push(`      <data key="e_flux">${flux}</data>`);
      lines.push('    </edge>');
    }

    lines.push('  </graph>');
    lines.push('</graphml>');

    return lines.join('\n');
  }

  /**
   * Exports topology graph as Node-Link JSON compatible with D3 and Cytoscape.
   * @param {Object} graphData
   * @returns {Object} { nodes: Array, links: Array }
   */
  static exportToNodeLink(graphData = {}) {
    const nodes = (graphData.criticalPoints || graphData.nodes || []).map((n, i) => ({
      id: n.id !== undefined ? String(n.id) : `node_${i}`,
      name: n.name || `Node_${i}`,
      type: n.type || (n.morseIndex === 0 ? 'attractor' : (n.morseIndex === 3 ? 'repeller' : 'saddle')),
      morseIndex: n.morseIndex !== undefined ? n.morseIndex : 0,
      position: n.position || [n.x || 0, n.y || 0, n.z || 0],
      potential: n.potential || 0.0,
      eigenvalues: n.eigenvalues || [0.0, 0.0, 0.0]
    }));

    const links = (graphData.separatrices || graphData.edges || []).map((e, j) => ({
      id: e.id !== undefined ? String(e.id) : `edge_${j}`,
      source: String(e.source),
      target: String(e.target),
      type: e.type || 'separatrix',
      weight: e.weight !== undefined ? e.weight : 1.0,
      length: e.length || 0.0,
      flux: e.flux || 0.0
    }));

    return { nodes, links };
  }
}

/**
 * Bulk Flow Multi-Aperture Exporter.
 */
export class BulkFlowExporter {
  /**
   * Exports multi-shell bulk flow calculations into RFC 4180 CSV with covariance matrix terms.
   * @param {Array<Object>} shells
   * @returns {string} CSV text
   */
  static exportToCsv(shells = []) {
    const headers = [
      'shell_id',
      'r_min_mpc_h',
      'r_max_mpc_h',
      'r_eff_mpc_h',
      'vx_bulk_kms',
      'vy_bulk_kms',
      'vz_bulk_kms',
      'v_bulk_mag_kms',
      'v_bulk_err_kms',
      'galactic_l_deg',
      'galactic_b_deg',
      'cov_xx',
      'cov_xy',
      'cov_xz',
      'cov_yy',
      'cov_yz',
      'cov_zz',
      'galaxy_sample_size'
    ];

    const lines = [headers.join(',')];

    for (let i = 0; i < shells.length; i++) {
      const s = shells[i];
      const v = s.vBulk || [s.vx || 0, s.vy || 0, s.vz || 0];
      const cov = s.covarianceMatrix || s.cov || [
        [100, 0, 0],
        [0, 100, 0],
        [0, 0, 100]
      ];
      const vMag = s.magnitude !== undefined ? s.magnitude : Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      const vErr = s.error !== undefined ? s.error : Math.sqrt((cov[0][0] + cov[1][1] + cov[2][2]) / 3);

      const row = [
        s.id !== undefined ? s.id : i,
        (s.rMin !== undefined ? s.rMin : 0).toFixed(2),
        (s.rMax !== undefined ? s.rMax : 50).toFixed(2),
        (s.rEff !== undefined ? s.rEff : 25).toFixed(2),
        v[0].toFixed(2),
        v[1].toFixed(2),
        v[2].toFixed(2),
        vMag.toFixed(2),
        vErr.toFixed(2),
        (s.galacticL !== undefined ? s.galacticL : 0).toFixed(3),
        (s.galacticB !== undefined ? s.galacticB : 0).toFixed(3),
        cov[0][0].toFixed(3),
        cov[0][1].toFixed(3),
        cov[0][2].toFixed(3),
        cov[1][1].toFixed(3),
        cov[1][2].toFixed(3),
        cov[2][2].toFixed(3),
        s.sampleSize !== undefined ? s.sampleSize : 0
      ];
      lines.push(row.join(','));
    }

    return lines.join('\n');
  }

  /**
   * Exports bulk flow multi-shell data to IVOA VOTable 1.4 XML format.
   * @param {Array<Object>} shells
   * @returns {string} VOTable XML string
   */
  static exportToVOTable(shells = []) {
    const lines = [
      '<?xml version="1.0" encoding="UTF-8"?>',
      '<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.4">',
      '  <RESOURCE name="BulkFlowCatalog">',
      '    <TABLE name="MultiShellBulkFlows">',
      '      <DESCRIPTION>Cosmicflows Multipole Bulk Flow Decompositions</DESCRIPTION>',
      '      <FIELD name="shell_id" datatype="int" ucd="meta.id"/>',
      '      <FIELD name="r_eff" datatype="float" unit="Mpc/h" ucd="pos.distance"/>',
      '      <FIELD name="v_bulk" datatype="float" unit="km/s" ucd="phys.veloc"/>',
      '      <FIELD name="v_err" datatype="float" unit="km/s" ucd="stat.error;phys.veloc"/>',
      '      <FIELD name="galactic_l" datatype="float" unit="deg" ucd="pos.galactic.lon"/>',
      '      <FIELD name="galactic_b" datatype="float" unit="deg" ucd="pos.galactic.lat"/>',
      '      <DATA>',
      '        <TABLEDATA>'
    ];

    for (let i = 0; i < shells.length; i++) {
      const s = shells[i];
      const v = s.vBulk || [s.vx || 0, s.vy || 0, s.vz || 0];
      const vMag = s.magnitude !== undefined ? s.magnitude : Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      const vErr = s.error || 15.0;

      lines.push('          <TR>');
      lines.push(`            <TD>${s.id !== undefined ? s.id : i}</TD>`);
      lines.push(`            <TD>${(s.rEff || 0).toFixed(2)}</TD>`);
      lines.push(`            <TD>${vMag.toFixed(2)}</TD>`);
      lines.push(`            <TD>${vErr.toFixed(2)}</TD>`);
      lines.push(`            <TD>${(s.galacticL || 0).toFixed(3)}</TD>`);
      lines.push(`            <TD>${(s.galacticB || 0).toFixed(3)}</TD>`);
      lines.push('          </TR>');
    }

    lines.push('        </TABLEDATA>');
    lines.push('      </DATA>');
    lines.push('    </TABLE>');
    lines.push('  </RESOURCE>');
    lines.push('</VOTABLE>');

    return lines.join('\n');
  }
}
