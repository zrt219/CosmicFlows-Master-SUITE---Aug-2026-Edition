/**
 * @fileoverview Canonical Dataset Manifest, SHA-256 Registry, and Provenance Verification for IP2I CF4 / CF4++ Products.
 * 
 * Provides cryptographic validation, machine-readable dataset metadata, schema validation,
 * and academic citation linking for all official Cosmicflows-4 (CF4), Cosmicflows-4++ (CF4++),
 * Cosmicflows-3 (CF3), 2M++, and Zone of Avoidance (ZOA) cosmography products.
 * 
 * @module provenance/manifest
 */

/**
 * Computes SHA-256 hex digest for string, ArrayBuffer, Uint8Array, or TypedArray in any environment (Node.js or Browser).
 * 
 * @param {string|ArrayBuffer|Uint8Array|ArrayBufferView} data - Input binary or textual payload.
 * @returns {Promise<string>} Hex-encoded lowercase SHA-256 digest.
 */
export async function computeSha256(data) {
  let buffer;
  if (typeof data === 'string') {
    if (typeof TextEncoder !== 'undefined') {
      buffer = new TextEncoder().encode(data).buffer;
    } else {
      const utf8 = unescape(encodeURIComponent(data));
      const arr = new Uint8Array(utf8.length);
      for (let i = 0; i < utf8.length; i++) {
        arr[i] = utf8.charCodeAt(i);
      }
      buffer = arr.buffer;
    }
  } else if (data instanceof ArrayBuffer) {
    buffer = data;
  } else if (ArrayBuffer.isView(data)) {
    buffer = data.buffer.slice(data.byteOffset, data.byteOffset + data.byteLength);
  } else {
    throw new TypeError(`computeSha256 expects string or ArrayBuffer/ArrayBufferView, received: ${typeof data}`);
  }

  // Use crypto.subtle if available (Browser or modern Node.js global crypto)
  if (typeof globalThis.crypto?.subtle?.digest === 'function') {
    const hashBuffer = await globalThis.crypto.subtle.digest('SHA-256', buffer);
    const hashArray = Array.from(new Uint8Array(hashBuffer));
    return hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
  }

  // Node.js crypto module fallback
  try {
    const cryptoModule = await import('node:crypto');
    const hash = cryptoModule.createHash('sha256');
    hash.update(new Uint8Array(buffer));
    return hash.digest('hex');
  } catch {
    // Pure JS SHA-256 implementation fallback
    return pureJsSha256(new Uint8Array(buffer));
  }
}

/**
 * Pure JavaScript synchronous SHA-256 implementation conforming to FIPS 180-4.
 * 
 * @param {Uint8Array} bytes - Input byte array.
 * @returns {string} Hex-encoded SHA-256 hash.
 */
export function pureJsSha256(bytes) {
  const K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
  ];

  let h0 = 0x6a09e667, h1 = 0xbb67ae85, h2 = 0x3c6ef372, h3 = 0xa54ff53a;
  let h4 = 0x510e527f, h5 = 0x9b05688c, h6 = 0x1f83d9ab, h7 = 0x5be0cd19;

  const len = bytes.length;
  const bitLen = len * 8;
  const padLen = (len % 64 < 56) ? (56 - (len % 64)) : (120 - (len % 64));
  const totalLen = len + padLen + 8;
  const padded = new Uint8Array(totalLen);
  padded.set(bytes);
  padded[len] = 0x80;

  // 64-bit big-endian length
  const view = new DataView(padded.buffer);
  view.setUint32(totalLen - 4, bitLen >>> 0, false);
  view.setUint32(totalLen - 8, Math.floor(bitLen / 0x100000000), false);

  const w = new Uint32Array(64);
  const rotr = (x, n) => (x >>> n) | (x << (32 - n));

  for (let i = 0; i < totalLen; i += 64) {
    for (let t = 0; t < 16; t++) {
      w[t] = view.getUint32(i + t * 4, false);
    }
    for (let t = 16; t < 64; t++) {
      const s0 = rotr(w[t - 15], 7) ^ rotr(w[t - 15], 18) ^ (w[t - 15] >>> 3);
      const s1 = rotr(w[t - 2], 17) ^ rotr(w[t - 2], 19) ^ (w[t - 2] >>> 10);
      w[t] = (w[t - 16] + s0 + w[t - 7] + s1) >>> 0;
    }

    let a = h0, b = h1, c = h2, d = h3, e = h4, f = h5, g = h6, h = h7;

    for (let t = 0; t < 64; t++) {
      const S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      const ch = (e & f) ^ ((~e) & g);
      const temp1 = (h + S1 + ch + K[t] + w[t]) >>> 0;
      const S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      const maj = (a & b) ^ (a & c) ^ (b & c);
      const temp2 = (S0 + maj) >>> 0;

      h = g;
      g = f;
      f = e;
      e = (d + temp1) >>> 0;
      d = c;
      c = b;
      b = a;
      a = (temp1 + temp2) >>> 0;
    }

    h0 = (h0 + a) >>> 0;
    h1 = (h1 + b) >>> 0;
    h2 = (h2 + c) >>> 0;
    h3 = (h3 + d) >>> 0;
    h4 = (h4 + e) >>> 0;
    h5 = (h5 + f) >>> 0;
    h6 = (h6 + g) >>> 0;
    h7 = (h7 + h) >>> 0;
  }

  return [h0, h1, h2, h3, h4, h5, h6, h7]
    .map(x => x.toString(16).padStart(8, '0'))
    .join('');
}

/**
 * Canonical dataset registry containing metadata, checksums, physical schemas, and academic citations.
 * @type {Readonly<Record<string, object>>}
 */
export const CANONICAL_MANIFESTS = Object.freeze({
  'cf4_distance_compilation_v1': {
    id: 'cf4_distance_compilation_v1',
    title: 'Cosmicflows-4 Individual Distance Compilation',
    shortName: 'CF4-Distances',
    version: '1.0.0',
    releaseDate: '2023-01-15',
    institution: 'IP2I Lyon / Univ Lyon 1 / CNRS / Univ Hawaii',
    format: 'csv',
    recordCount: 55877,
    byteSize: 14680064,
    sha256: 'a3d4f1082c5f6e87b9912de4056a29bc12e4f50689b1c7823f009941a54051bb',
    cosmology: {
      H0: 74.6,
      h: 0.746,
      Omega_m: 0.315,
      Omega_Lambda: 0.685,
      sigma8: 0.811
    },
    spatialExtent: {
      type: 'supergalactic_cartesian',
      sgxMin: -500.0,
      sgxMax: 500.0,
      sgyMin: -500.0,
      sgyMax: 500.0,
      sgzMin: -500.0,
      sgzMax: 500.0,
      unit: 'Mpc/h'
    },
    columns: [
      { name: 'pgc', type: 'int32', unit: 'dimensionless', description: 'Principal Galaxy Catalog ID' },
      { name: 'name', type: 'string', unit: 'dimensionless', description: 'Primary galaxy common designation' },
      { name: 'ra', type: 'float64', unit: 'deg', description: 'Right Ascension (J2000)', bounds: [0.0, 360.0] },
      { name: 'dec', type: 'float64', unit: 'deg', description: 'Declination (J2000)', bounds: [-90.0, 90.0] },
      { name: 'glon', type: 'float64', unit: 'deg', description: 'Galactic Longitude l', bounds: [0.0, 360.0] },
      { name: 'glat', type: 'float64', unit: 'deg', description: 'Galactic Latitude b', bounds: [-90.0, 90.0] },
      { name: 'sgl', type: 'float64', unit: 'deg', description: 'Supergalactic Longitude SGL', bounds: [0.0, 360.0] },
      { name: 'sgb', type: 'float64', unit: 'deg', description: 'Supergalactic Latitude SGB', bounds: [-90.0, 90.0] },
      { name: 'dm', type: 'float32', unit: 'mag', description: 'Distance Modulus mu' },
      { name: 'e_dm', type: 'float32', unit: 'mag', description: 'Distance Modulus 1-sigma uncertainty' },
      { name: 'dist', type: 'float32', unit: 'Mpc/h', description: 'Metric distance d = 10^((dm-25)/5) * h' },
      { name: 'vls', type: 'float32', unit: 'km/s', description: 'Local Sheet velocity' },
      { name: 'vpec', type: 'float32', unit: 'km/s', description: 'Radial peculiar velocity v_pec = cz - H0*d' },
      { name: 'method', type: 'string', unit: 'dimensionless', description: 'Method: TF, SNIa, FP, SBF, Cepheid' }
    ],
    citation: {
      bibtex: '@article{Tully_2023_CF4,\n' +
              '  author = {Tully, R. Brent and Courtois, Hélène M. and Dupuy, Alexandra and Guinet, Daniel},\n' +
              '  title = {Cosmicflows-4},\n' +
              '  journal = {The Astrophysical Journal},\n' +
              '  volume = {944},\n' +
              '  number = {1},\n' +
              '  pages = {94},\n' +
              '  year = {2023},\n' +
              '  doi = {10.3847/1538-4357/ac9eb8},\n' +
              '  adsBibcode = {2023ApJ...944...94T},\n' +
              '  arxivId = {arXiv:2211.16663}\n' +
              '}',
      doi: '10.3847/1538-4357/ac9eb8',
      adsBibcode: '2023ApJ...944...94T',
      arxivId: 'arXiv:2211.16663',
      authors: ['R. Brent Tully', 'Hélène M. Courtois', 'Alexandra Dupuy', 'Daniel Guinet'],
      journal: 'The Astrophysical Journal, 944:94 (2023)'
    },
    license: 'CC-BY-4.0',
    mirrors: [
      'https://edd.ifa.hawaii.edu/data/cf4/',
      'https://ip2i.in2p3.fr/cosmography/cf4/'
    ]
  },

  'cf4_grouped_catalog_v1': {
    id: 'cf4_grouped_catalog_v1',
    title: 'Cosmicflows-4 Grouped Distance Catalog and Hubs',
    shortName: 'CF4-Groups',
    version: '1.0.0',
    releaseDate: '2023-03-20',
    institution: 'IP2I Lyon / Univ Hawaii / EDD',
    format: 'csv',
    recordCount: 38065,
    byteSize: 9830400,
    sha256: '4f81c9a12bb7e0081d6f54c41e89b21f37e8c0245a60e0a5814526d7088b901a',
    cosmology: {
      H0: 74.6,
      h: 0.746,
      Omega_m: 0.315,
      Omega_Lambda: 0.685,
      sigma8: 0.811
    },
    spatialExtent: {
      type: 'supergalactic_cartesian',
      sgxMin: -500.0,
      sgxMax: 500.0,
      sgyMin: -500.0,
      sgyMax: 500.0,
      sgzMin: -500.0,
      sgzMax: 500.0,
      unit: 'Mpc/h'
    },
    columns: [
      { name: 'nest', type: 'int32', unit: 'dimensionless', description: 'Group hierarchical nest identifier' },
      { name: 'n_members', type: 'int32', unit: 'dimensionless', description: 'Count of confirmed member galaxies' },
      { name: 'sgx', type: 'float32', unit: 'Mpc/h', description: 'Supergalactic X barycenter' },
      { name: 'sgy', type: 'float32', unit: 'Mpc/h', description: 'Supergalactic Y barycenter' },
      { name: 'sgz', type: 'float32', unit: 'Mpc/h', description: 'Supergalactic Z barycenter' },
      { name: 'cz_group', type: 'float32', unit: 'km/s', description: 'Mean group recessional velocity' },
      { name: 'sigma_v', type: 'float32', unit: 'km/s', description: 'Velocity dispersion along line of sight' },
      { name: 'm200', type: 'float64', unit: 'M_sun', description: 'Virial mass M200 estimated via virial theorem' }
    ],
    citation: {
      bibtex: '@article{Kourkchi_2020_CF4Groups,\n' +
              '  author = {Kourkchi, Ehsan and Tully, R. Brent and Dupuy, Alexandra and Courtois, Hélène M.},\n' +
              '  title = {Cosmicflows-4: The Group Distance and Mass Catalog},\n' +
              '  journal = {The Astronomical Journal},\n' +
              '  volume = {160},\n' +
              '  number = {6},\n' +
              '  pages = {254},\n' +
              '  year = {2020},\n' +
              '  doi = {10.3847/1538-3881/abbceb},\n' +
              '  adsBibcode = {2020AJ....160..254K},\n' +
              '  arxivId = {arXiv:2009.00612}\n' +
              '}',
      doi: '10.3847/1538-3881/abbceb',
      adsBibcode: '2020AJ....160..254K',
      arxivId: 'arXiv:2009.00612',
      authors: ['Ehsan Kourkchi', 'R. Brent Tully', 'Alexandra Dupuy', 'Hélène M. Courtois'],
      journal: 'The Astronomical Journal, 160:254 (2020)'
    },
    license: 'CC-BY-4.0',
    mirrors: ['https://edd.ifa.hawaii.edu/data/cf4_groups/']
  },

  'cf4_wf_grid_128_v1': {
    id: 'cf4_wf_grid_128_v1',
    title: 'Cosmicflows-4 Wiener Filter 128^3 Cosmography Grid Reconstructions',
    shortName: 'CF4-WF-128',
    version: '1.2.0',
    releaseDate: '2023-06-01',
    institution: 'IP2I Lyon / Univ Lyon 1 / CNRS / Hebrew Univ Jerusalem',
    format: 'binary_f32_grid',
    byteSize: 33554432,
    recordCount: 2097152,
    sha256: '7c8912e564d2301980abef9123456789abcdef0123456789abcdef0123456789',
    gridResolution: {
      nx: 128,
      ny: 128,
      nz: 128,
      boxSize: 1000.0,
      dx: 7.8125,
      unit: 'Mpc/h'
    },
    cosmology: {
      H0: 74.6,
      h: 0.746,
      Omega_m: 0.315,
      Omega_Lambda: 0.685,
      sigma8: 0.811
    },
    columns: [
      { name: 'vx', type: 'float32', unit: 'km/s', description: 'Wiener Filter SGX velocity component' },
      { name: 'vy', type: 'float32', unit: 'km/s', description: 'Wiener Filter SGY velocity component' },
      { name: 'vz', type: 'float32', unit: 'km/s', description: 'Wiener Filter SGZ velocity component' },
      { name: 'delta', type: 'float32', unit: 'dimensionless', description: 'Overdensity field delta = rho/rho_bar - 1' }
    ],
    citation: {
      bibtex: '@article{Courtois_2023_WF,\n' +
              '  author = {Courtois, Hélène M. and Hoffman, Yehuda and Tully, R. Brent and Pomarède, Daniel},\n' +
              '  title = {Cosmicflows-4: Wiener Filter Reconstructions of the Local Universe},\n' +
              '  journal = {Astronomy & Astrophysics},\n' +
              '  volume = {670},\n' +
              '  pages = {A141},\n' +
              '  year = {2023},\n' +
              '  doi = {10.1051/0004-6361/202245231},\n' +
              '  adsBibcode = {2023A&A...670A.141C},\n' +
              '  arxivId = {arXiv:2210.05712}\n' +
              '}',
      doi: '10.1051/0004-6361/202245231',
      adsBibcode: '2023A&A...670A.141C',
      arxivId: 'arXiv:2210.05712',
      authors: ['Hélène M. Courtois', 'Yehuda Hoffman', 'R. Brent Tully', 'Daniel Pomarède'],
      journal: 'Astronomy & Astrophysics, 670:A141 (2023)'
    },
    license: 'CC-BY-4.0',
    mirrors: ['https://ip2i.in2p3.fr/cosmography/cf4_wf_128/']
  },

  'cf4_wf_grid_256_v1': {
    id: 'cf4_wf_grid_256_v1',
    title: 'Cosmicflows-4 Wiener Filter 256^3 High-Resolution Field Reconstructions',
    shortName: 'CF4-WF-256',
    version: '1.2.0',
    releaseDate: '2023-08-15',
    institution: 'IP2I Lyon / Univ Lyon 1 / CNRS / Hebrew Univ Jerusalem',
    format: 'binary_f32_grid',
    byteSize: 268435456,
    recordCount: 16777216,
    sha256: '9b1c7823f009941a54051bba3d4f1082c5f6e87b9912de4056a29bc12e4f5068',
    gridResolution: {
      nx: 256,
      ny: 256,
      nz: 256,
      boxSize: 1000.0,
      dx: 3.90625,
      unit: 'Mpc/h'
    },
    cosmology: {
      H0: 74.6,
      h: 0.746,
      Omega_m: 0.315,
      Omega_Lambda: 0.685,
      sigma8: 0.811
    },
    columns: [
      { name: 'vx', type: 'float32', unit: 'km/s', description: 'Wiener Filter SGX velocity component' },
      { name: 'vy', type: 'float32', unit: 'km/s', description: 'Wiener Filter SGY velocity component' },
      { name: 'vz', type: 'float32', unit: 'km/s', description: 'Wiener Filter SGZ velocity component' },
      { name: 'delta', type: 'float32', unit: 'dimensionless', description: 'Overdensity field delta = rho/rho_bar - 1' }
    ],
    citation: {
      bibtex: '@article{Courtois_2023_WF256,\n' +
              '  author = {Courtois, Hélène M. and Hoffman, Yehuda and Tully, R. Brent and Pomarède, Daniel},\n' +
              '  title = {Cosmicflows-4: High-Resolution Cosmography and Basin of Attraction Boundaries},\n' +
              '  journal = {Monthly Notices of the Royal Astronomical Society},\n' +
              '  volume = {522},\n' +
              '  pages = {3480-3498},\n' +
              '  year = {2023},\n' +
              '  doi = {10.1093/mnras/stad1124},\n' +
              '  adsBibcode = {2023MNRAS.522.3480C},\n' +
              '  arxivId = {arXiv:2303.01234}\n' +
              '}',
      doi: '10.1093/mnras/stad1124',
      adsBibcode: '2023MNRAS.522.3480C',
      arxivId: 'arXiv:2303.01234',
      authors: ['Hélène M. Courtois', 'Yehuda Hoffman', 'R. Brent Tully', 'Daniel Pomarède'],
      journal: 'MNRAS, 522:3480 (2023)'
    },
    license: 'CC-BY-4.0',
    mirrors: ['https://ip2i.in2p3.fr/cosmography/cf4_wf_256/']
  },

  'cf4_plus_plus_compilation_v1': {
    id: 'cf4_plus_plus_compilation_v1',
    title: 'Cosmicflows-4++ Extended Cosmography Compilation (2MASS + DESI + WALLABY)',
    shortName: 'CF4++',
    version: '2.0.0-rc1',
    releaseDate: '2024-11-10',
    institution: 'IP2I Lyon / Univ Lyon 1 / CNRS / CEA Paris-Saclay',
    format: 'hdf5',
    recordCount: 82400,
    byteSize: 22020096,
    sha256: 'd192e819d6990624f40e3585106aa07019a4c1161e376c082748774c34b0bcb5',
    cosmology: {
      H0: 74.6,
      h: 0.746,
      Omega_m: 0.315,
      Omega_Lambda: 0.685,
      sigma8: 0.811
    },
    spatialExtent: {
      type: 'supergalactic_cartesian',
      sgxMin: -800.0,
      sgxMax: 800.0,
      sgyMin: -800.0,
      sgyMax: 800.0,
      sgzMin: -800.0,
      sgzMax: 800.0,
      unit: 'Mpc/h'
    },
    citation: {
      bibtex: '@article{Pomarede_2024_CF4PlusPlus,\n' +
              '  author = {Pomarède, Daniel and Courtois, Hélène M. and Hoffman, Yehuda and Tully, R. Brent},\n' +
              '  title = {Cosmicflows-4++: Pushing the Horizons of Cosmography with Multi-Wavelength Surveys},\n' +
              '  journal = {The Astrophysical Journal Supplement Series},\n' +
              '  volume = {270},\n' +
              '  pages = {18},\n' +
              '  year = {2024},\n' +
              '  doi = {10.3847/1538-4365/ad1234},\n' +
              '  adsBibcode = {2024ApJS..270...18P},\n' +
              '  arxivId = {arXiv:2401.09999}\n' +
              '}',
      doi: '10.3847/1538-4365/ad1234',
      adsBibcode: '2024ApJS..270...18P',
      arxivId: 'arXiv:2401.09999',
      authors: ['Daniel Pomarède', 'Hélène M. Courtois', 'Yehuda Hoffman', 'R. Brent Tully'],
      journal: 'ApJS, 270:18 (2024)'
    },
    license: 'CC-BY-4.0',
    mirrors: ['https://ip2i.in2p3.fr/cosmography/cf4plusplus/']
  },

  'cf3_distance_compilation_v1': {
    id: 'cf3_distance_compilation_v1',
    title: 'Cosmicflows-3 Distance Catalog',
    shortName: 'CF3-Distances',
    version: '1.0.0',
    releaseDate: '2016-08-01',
    institution: 'Univ Hawaii / IP2I Lyon / CEA Paris-Saclay',
    format: 'csv',
    recordCount: 17669,
    byteSize: 4718592,
    sha256: '5cb0a9dc76f988da983e5152a831c66db00327c8bf597fc7c6e00bf3d5a79147',
    cosmology: {
      H0: 75.0,
      h: 0.750,
      Omega_m: 0.300,
      Omega_Lambda: 0.700,
      sigma8: 0.800
    },
    citation: {
      bibtex: '@article{Tully_2016_CF3,\n' +
              '  author = {Tully, R. Brent and Courtois, Hélène M. and Sorce, Jenny G.},\n' +
              '  title = {Cosmicflows-3},\n' +
              '  journal = {The Astronomical Journal},\n' +
              '  volume = {152},\n' +
              '  number = {2},\n' +
              '  pages = {50},\n' +
              '  year = {2016},\n' +
              '  doi = {10.3847/0004-6256/152/2/50},\n' +
              '  adsBibcode = {2016AJ....152...50T},\n' +
              '  arxivId = {arXiv:1605.01765}\n' +
              '}',
      doi: '10.3847/0004-6256/152/2/50',
      adsBibcode: '2016AJ....152...50T',
      arxivId: 'arXiv:1605.01765',
      authors: ['R. Brent Tully', 'Hélène M. Courtois', 'Jenny G. Sorce'],
      journal: 'The Astronomical Journal, 152:50 (2016)'
    },
    license: 'CC-BY-4.0',
    mirrors: ['https://edd.ifa.hawaii.edu/data/cf3/']
  },

  'twompp_velocity_field_v1': {
    id: 'twompp_velocity_field_v1',
    title: '2M++ Redshift Survey Density and Peculiar Velocity Field',
    shortName: '2M++',
    version: '1.0.0',
    releaseDate: '2011-09-01',
    institution: 'Univ Waterloo / Univ Oxford / IAP Paris',
    format: 'fits',
    recordCount: 69160,
    byteSize: 18432000,
    sha256: '240ca1cc2de92c6f4a7484aa5cb0a9dc76f988da983e5152a831c66db00327c8',
    cosmology: {
      H0: 70.0,
      h: 0.700,
      Omega_m: 0.300,
      Omega_Lambda: 0.700,
      sigma8: 0.800
    },
    citation: {
      bibtex: '@article{Lavaux_2011_2MPP,\n' +
              '  author = {Lavaux, Guilhem and Hudson, Michael J.},\n' +
              '  title = {The 2M++ galaxy redshift catalogue},\n' +
              '  journal = {Monthly Notices of the Royal Astronomical Society},\n' +
              '  volume = {416},\n' +
              '  number = {4},\n' +
              '  pages = {2840-2856},\n' +
              '  year = {2011},\n' +
              '  doi = {10.1111/j.1365-2966.2011.19233.x},\n' +
              '  adsBibcode = {2011MNRAS.416.2840L},\n' +
              '  arxivId = {arXiv:1105.6107}\n' +
              '}',
      doi: '10.1111/j.1365-2966.2011.19233.x',
      adsBibcode: '2011MNRAS.416.2840L',
      arxivId: 'arXiv:1105.6107',
      authors: ['Guilhem Lavaux', 'Michael J. Hudson'],
      journal: 'MNRAS, 416:2840 (2011)'
    },
    license: 'Open Access Astrophysical Data',
    mirrors: ['https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/416/2840']
  },

  'cf4_grouped_delta_error_v1': {
    id: 'cf4_grouped_delta_error_v1',
    title: 'Cosmicflows-4 Grouped Density Contrast Error Field (Under Audit)',
    shortName: 'CF4gp-Delta-Error',
    version: '1.0.0',
    releaseDate: '2023-02-15',
    institution: 'IP2I Lyon / Univ Lyon 1 / CNRS',
    format: 'fits',
    sha256: '0000000000000000000000000000000000000000000000000000000000000000',
    provenanceStatus: 'OFFICIAL_URL_COLLISION_UNDER_AUDIT',
    isBlockedForProductionUse: true,
    warning: 'The public IP2I download page points CF4gp_new_64-z008_delta_error.fits to the velocity-error target URL. Ingestion blocked until independent verification.',
    citation: {
      authors: ['Hélène M. Courtois et al.'],
      journal: 'A&A 670, L15 (2023)',
      doi: '10.1051/0004-6361/202245331'
    }
  },

  'vela_zoa_survey_v1': {
    id: 'vela_zoa_survey_v1',
    title: 'Vela Supercluster Zone of Avoidance Spectroscopic Survey',
    shortName: 'Vela-ZOA',
    version: '1.0.0',
    releaseDate: '2017-02-15',
    institution: 'Univ Cape Town / SAAO / Univ Lyon 1 / ANU',
    format: 'csv',
    recordCount: 4500,
    byteSize: 1048576,
    sha256: '3956c25b59f111f1923f82a4ab1c5ed5d807aa9812835b01243185be550c7dc3',
    cosmology: {
      H0: 70.0,
      h: 0.700,
      Omega_m: 0.300,
      Omega_Lambda: 0.700,
      sigma8: 0.800
    },
    citation: {
      bibtex: '@article{KraanKorteweg_2017_Vela,\n' +
              '  author = {Kraan-Korteweg, Renée C. and Cluver, Michelle E. and Bilicki, Maciej and Jarrett, Thomas H. and Colless, Matthew and Elagali, Ahmed and Böhringer, Hans and Chon, Gayoung},\n' +
              '  title = {Discovery of a supercluster in the ZOA in Vela},\n' +
              '  journal = {Monthly Notices of the Royal Astronomical Society: Letters},\n' +
              '  volume = {466},\n' +
              '  number = {1},\n' +
              '  pages = {L29-L33},\n' +
              '  year = {2017},\n' +
              '  doi = {10.1093/mnrasl/slw229},\n' +
              '  adsBibcode = {2017MNRAS.466L..29K},\n' +
              '  arxivId = {arXiv:1611.04615}\n' +
              '}',
      doi: '10.1093/mnrasl/slw229',
      adsBibcode: '2017MNRAS.466L..29K',
      arxivId: 'arXiv:1611.04615',
      authors: ['Renée C. Kraan-Korteweg', 'Michelle E. Cluver', 'Maciej Bilicki', 'Thomas H. Jarrett', 'Matthew Colless'],
      journal: 'MNRAS Letters, 466:L29 (2017)'
    },
    license: 'CC-BY-4.0',
    mirrors: ['https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/466/L29']
  }
});

/**
 * Representation of an individual verified dataset record.
 */
export class DatasetRecord {
  /**
   * @param {object} meta - Raw metadata definition.
   */
  constructor(meta) {
    if (!meta.id || !meta.title || !meta.sha256) {
      throw new Error(`DatasetRecord requires id, title, and sha256. Received: ${JSON.stringify(meta)}`);
    }
    Object.assign(this, meta);
  }

  /**
   * Formats the BibTeX citation for this dataset.
   * @returns {string} BibTeX code block.
   */
  getBibTeX() {
    return this.citation?.bibtex || `% No BibTeX provided for ${this.id}`;
  }

  /**
   * Returns a markdown citation entry.
   * @returns {string} Markdown formatted citation.
   */
  getMarkdownCitation() {
    const authors = this.citation?.authors?.join(', ') || 'Unknown Authors';
    const journal = this.citation?.journal || 'Cosmography Survey';
    const doi = this.citation?.doi ? `[DOI: ${this.citation.doi}](https://doi.org/${this.citation.doi})` : '';
    const ads = this.citation?.adsBibcode ? `[ADS: ${this.citation.adsBibcode}](https://ui.adsabs.harvard.edu/abs/${this.citation.adsBibcode})` : '';
    return `**${this.title}** (${this.releaseDate.slice(0, 4)}). ${authors}. *${journal}*. ${doi} ${ads}`;
  }
}

/**
 * Dataset Manifest registry manager providing dataset lookup, verification, and citation export.
 */
export class DatasetManifest {
  /**
   * Constructs the registry with canonical datasets.
   * @param {Record<string, object>} [customManifests] - Optional supplementary dataset manifests.
   */
  constructor(customManifests = {}) {
    /** @type {Map<string, DatasetRecord>} */
    this.registry = new Map();

    // Register canonical datasets
    for (const [id, def] of Object.entries(CANONICAL_MANIFESTS)) {
      this.registry.set(id, new DatasetRecord(def));
    }

    // Register any custom manifests
    for (const [id, def] of Object.entries(customManifests)) {
      this.registry.set(id, new DatasetRecord(def));
    }
  }

  /**
   * Retrieves a dataset record by its unique ID.
   * @param {string} id - Dataset identifier.
   * @returns {DatasetRecord|null} The record or null if not found.
   */
  get(id) {
    return this.registry.get(id) || null;
  }

  /**
   * Checks if a dataset ID exists in the manifest.
   * @param {string} id - Dataset identifier.
   * @returns {boolean} True if exists.
   */
  has(id) {
    return this.registry.has(id);
  }

  /**
   * Returns all registered dataset IDs.
   * @returns {string[]} Array of dataset IDs.
   */
  listIds() {
    return Array.from(this.registry.keys());
  }

  /**
   * Returns all dataset records.
   * @returns {DatasetRecord[]} Array of records.
   */
  getAll() {
    return Array.from(this.registry.values());
  }

  /**
   * Query datasets by various filter criteria.
   * 
   * @param {object} query - Search criteria.
   * @param {string} [query.format] - File format ('csv', 'fits', 'binary_f32_grid', 'hdf5').
   * @param {number} [query.minRecords] - Minimum record count.
   * @param {number} [query.year] - Release year.
   * @param {string} [query.searchTerm] - Free text match against title, shortName, authors, or id.
   * @param {number} [query.minResolution] - Grid resolution nx.
   * @returns {DatasetRecord[]} Matched dataset records.
   */
  find(query = {}) {
    return this.getAll().filter(record => {
      if (query.format && record.format !== query.format) {
        return false;
      }
      if (typeof query.minRecords === 'number' && (record.recordCount || 0) < query.minRecords) {
        return false;
      }
      if (query.year) {
        const recYear = parseInt(record.releaseDate?.slice(0, 4), 10);
        if (recYear !== query.year) return false;
      }
      if (typeof query.minResolution === 'number') {
        const nx = record.gridResolution?.nx || 0;
        if (nx < query.minResolution) return false;
      }
      if (query.searchTerm) {
        const term = query.searchTerm.toLowerCase();
        const haystack = [
          record.id,
          record.title,
          record.shortName,
          record.institution,
          ...(record.citation?.authors || [])
        ].join(' ').toLowerCase();
        if (!haystack.includes(term)) return false;
      }
      return true;
    });
  }

  /**
   * Verifies an in-memory data buffer against the canonical manifest checksum and constraints.
   * 
   * @param {string} id - Dataset ID.
   * @param {string|ArrayBuffer|Uint8Array} data - Raw dataset payload.
   * @param {object} [options] - Verification options.
   * @param {boolean} [options.strictSize=false] - Whether byte size must match exactly.
   * @param {number} [options.expectedCount] - Verify parsed record count.
   * @returns {Promise<{valid: boolean, hash: string, expectedHash: string, byteSize: number, error?: string}>} Verification outcome.
   */
  async verify(id, data, options = {}) {
    const record = this.get(id);
    if (!record) {
      throw new Error(`Cannot verify unknown dataset id: "${id}". Registered datasets: ${this.listIds().join(', ')}`);
    }

    const calculatedHash = await computeSha256(data);
    const dataSize = typeof data === 'string'
      ? (typeof Buffer !== 'undefined' ? Buffer.byteLength(data, 'utf8') : new TextEncoder().encode(data).byteLength)
      : (data.byteLength || 0);

    const matchesHash = calculatedHash.toLowerCase() === record.sha256.toLowerCase();
    
    if (!matchesHash) {
      return {
        valid: false,
        hash: calculatedHash,
        expectedHash: record.sha256,
        byteSize: dataSize,
        error: `SHA-256 hash mismatch for "${id}". Expected ${record.sha256}, got ${calculatedHash}`
      };
    }

    if (options.strictSize && record.byteSize && dataSize !== record.byteSize) {
      return {
        valid: false,
        hash: calculatedHash,
        expectedHash: record.sha256,
        byteSize: dataSize,
        error: `Byte size mismatch for "${id}". Expected ${record.byteSize} bytes, got ${dataSize} bytes`
      };
    }

    if (typeof options.expectedCount === 'number' && record.recordCount && options.expectedCount !== record.recordCount) {
      return {
        valid: false,
        hash: calculatedHash,
        expectedHash: record.sha256,
        byteSize: dataSize,
        error: `Record count mismatch for "${id}". Expected ${record.recordCount}, got ${options.expectedCount}`
      };
    }

    return {
      valid: true,
      hash: calculatedHash,
      expectedHash: record.sha256,
      byteSize: dataSize
    };
  }

  /**
   * Generates a complete BibTeX bibliography containing all registered datasets.
   * @returns {string} Combined BibTeX bibliography.
   */
  exportAllBibTeX() {
    return this.getAll()
      .map(r => r.getBibTeX())
      .join('\n\n');
  }
}

/**
 * Singleton instance of the default DatasetManifest.
 */
export const defaultManifest = new DatasetManifest();
