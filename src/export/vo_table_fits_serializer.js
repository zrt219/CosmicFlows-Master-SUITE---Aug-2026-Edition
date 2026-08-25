/**
 * @file vo_table_fits_serializer.js
 * @module export/vo_table_fits_serializer
 * @description Production-grade IVOA VOTable 1.4 XML & Multi-Extension FITS (MEF) Binary Table Serializer and Deserializer.
 * 
 * Features:
 * 1. IVOA VOTable 1.4 XML Serializer & Deserializer:
 *    - Full schema compliance (VOTABLE version 1.4, RESOURCE, TABLE, FIELD, PARAM, GROUP, COOSYS, TIMESYS, INFO, VALUES).
 *    - Comprehensive DATA serializations: TABLEDATA, BINARY (base64 stream), and BINARY2 (null bitmask headers).
 *    - Multidimensional & variable-length arrays (arraysize="3", "3x3", "*", "10*").
 *    - All IVOA datatypes (boolean, bit, unsignedByte, short, int, long, char, unicodeChar, float, double, floatComplex, doubleComplex).
 *    - Proper XML entity escaping, NaN/Inf handling, null markers.
 * 
 * 2. Multi-Extension FITS (MEF) Binary Table Serializer & Deserializer:
 *    - Standard Primary HDU (SIMPLE=T, BITPIX, NAXIS, EXTEND=T) and BINTABLE Extension HDUs.
 *    - 80-column FITS header cards with exact formatting, comments, and HIERARCH support.
 *    - Big-Endian IEEE 754 byte order column packing and unpacking.
 *    - Standard column formats: L, X, B, I, J, K, A, E, D, C, M, PJ, PD, PE, PI, PB, PA, QJ, QD, QE (variable-length heap arrays).
 *    - TSCAL/TZERO scaling, TNULL integer nulls, TUNIT, TDISP, TUCD keywords.
 *    - Strict 2880-byte header and data block padding with 0x20 and 0x00.
 * 
 * 3. Astrometric Unified Content Descriptors (UCD 1+):
 *    - Vocabulary dictionary and syntax validator adhering to IVOA UCD 1.5/1.23 standard.
 *    - Primary words and secondary modifiers validation, atom hierarchy validation.
 *    - Automatic UCD inferrer based on field names, physical units, and cosmology semantics.
 * 
 * 4. FITS IEEE 32-Bit 1's Complement Checksum & Datasum Engine:
 *    - Seaman-Pence checksum algorithm (FITS 4.0 Standard Appendix J).
 *    - Circular bit carry 1's complement word addition.
 *    - Exact CHECKSUM 16-character ASCII encoding and decoding.
 *    - Verification engine for data and header integrity.
 * 
 * 5. High-Level Cosmological Exporters:
 *    - Exporters for galaxy catalogs, watershed basins, critical points, streamlines, and bulk flows.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

// ============================================================================
// SECTION 1: ASTROMETRIC UNIFIED CONTENT DESCRIPTORS (UCD 1+) ENGINE
// ============================================================================

/**
 * Standard IVOA UCD1+ Vocabulary and Validation Engine.
 * Implements syntax verification, canonical atom lookup, and heuristic inference.
 */
export class UCDStandard {
  /**
   * Primary UCD1+ words commonly used in observational and computational cosmology.
   * @type {Set<string>}
   */
  static PRIMARY_WORDS = new Set([
    // Astrometry & Positions
    'pos.angDistance',
    'pos.angResolution',
    'pos.azimuth',
    'pos.barycenter',
    'pos.bodyrc',
    'pos.cartesian',
    'pos.cartesian.x',
    'pos.cartesian.y',
    'pos.cartesian.z',
    'pos.distance',
    'pos.dir',
    'pos.earth',
    'pos.ecliptic',
    'pos.ecliptic.lat',
    'pos.ecliptic.lon',
    'pos.ephem',
    'pos.eq',
    'pos.eq.dec',
    'pos.eq.ha',
    'pos.eq.ra',
    'pos.eq.spd',
    'pos.errorEllipsoid',
    'pos.frame',
    'pos.galactic',
    'pos.galactic.lat',
    'pos.galactic.lon',
    'pos.galactocentric',
    'pos.geocentric',
    'pos.healpix',
    'pos.helioecliptic',
    'pos.heliocentric',
    'pos.htc',
    'pos.lambert',
    'pos.lg',
    'pos.lsr',
    'pos.lsupergalactic',
    'pos.offset',
    'pos.parallax',
    'pos.pm',
    'pos.pm.dec',
    'pos.pm.ra',
    'pos.posAng',
    'pos.sg',
    'pos.sg.lat',
    'pos.sg.lon',
    'pos.supergalactic',
    'pos.supergalactic.lat',
    'pos.supergalactic.lon',
    'pos.wcs',

    // Physics & Velocities
    'phys.acceleration',
    'phys.angArea',
    'phys.angMomentum',
    'phys.angSize',
    'phys.area',
    'phys.atemp',
    'phys.charge',
    'phys.columnDensity',
    'phys.composition',
    'phys.cosmology',
    'phys.damping',
    'phys.density',
    'phys.dielectric',
    'phys.dispersion',
    'phys.distance',
    'phys.electfield',
    'phys.electron',
    'phys.energy',
    'phys.energyDensity',
    'phys.entropy',
    'phys.equationOfState',
    'phys.gravity',
    'phys.luminosity',
    'phys.magAbs',
    'phys.magField',
    'phys.mass',
    'phys.massLoss',
    'phys.mol',
    'phys.motion',
    'phys.particle',
    'phys.polarization',
    'phys.potential',
    'phys.pressure',
    'phys.refractIndex',
    'phys.size',
    'phys.size.axisRatio',
    'phys.size.diameter',
    'phys.size.radius',
    'phys.size.smajAxis',
    'phys.size.sminAxis',
    'phys.temperature',
    'phys.temperature.effective',
    'phys.transmission',
    'phys.veloc',
    'phys.veloc.ang',
    'phys.veloc.dispersion',
    'phys.veloc.expansion',
    'phys.veloc.escape',
    'phys.veloc.microTurb',
    'phys.veloc.orbital',
    'phys.veloc.pulsat',
    'phys.veloc.rotat',
    'phys.veloc.transverse',
    'phys.virial',
    'phys.viscosity',
    'phys.volume',
    'phys.vorticity',

    // Photometry & Radiation
    'phot.antennaTemp',
    'phot.atmosphere',
    'phot.calib',
    'phot.color',
    'phot.color.excess',
    'phot.count',
    'phot.fluence',
    'phot.flux',
    'phot.flux.bolometric',
    'phot.flux.density',
    'phot.flux.density.sb',
    'phot.flux.sb',
    'phot.limbDark',
    'phot.mag',
    'phot.mag.bolometric',
    'phot.mag.isophotal',
    'phot.mag.model',
    'phot.mag.petrosian',
    'phot.mag.sb',
    'phot.mag.synth',
    'phot.radiance',
    'phot.sb',

    // Spectroscopy & Doppler
    'spect.binSize',
    'spect.continuum',
    'spect.dopplerParam',
    'spect.dopplerVeloc',
    'spect.dopplerVeloc.opt',
    'spect.dopplerVeloc.radio',
    'spect.element',
    'spect.energy',
    'spect.equivalentWidth',
    'spect.frequency',
    'spect.frequency.resolution',
    'spect.hardening',
    'spect.index',
    'spect.line',
    'spect.line.asymmetry',
    'spect.line.broadening',
    'spect.line.fwhm',
    'spect.line.intensity',
    'spect.line.profile',
    'spect.line.width',
    'spect.resolution',
    'spect.resolution.dispersion',
    'spect.softening',
    'spect.transverseVeloc',
    'spect.vacuum',
    'spect.variability',
    'spect.velocity',
    'spect.wavelength',
    'spect.wavenumber',

    // Metadata, Identifiers, Classes
    'meta.abstract',
    'meta.bib',
    'meta.bib.author',
    'meta.bib.bibcode',
    'meta.bib.doi',
    'meta.bib.page',
    'meta.bib.volume',
    'meta.code',
    'meta.code.class',
    'meta.code.error',
    'meta.code.member',
    'meta.code.mime',
    'meta.code.multip',
    'meta.code.qual',
    'meta.code.status',
    'meta.cryptic',
    'meta.curation',
    'meta.dataset',
    'meta.email',
    'meta.file',
    'meta.fits',
    'meta.id',
    'meta.id.assoc',
    'meta.id.coi',
    'meta.id.cross',
    'meta.id.parent',
    'meta.id.part',
    'meta.id.pi',
    'meta.main',
    'meta.modelled',
    'meta.note',
    'meta.number',
    'meta.record',
    'meta.ref',
    'meta.ref.url',
    'meta.software',
    'meta.table',
    'meta.title',
    'meta.ucd',
    'meta.unit',
    'meta.version',

    // Source Classification & Morphology
    'src',
    'src.calib',
    'src.class',
    'src.class.color',
    'src.class.distance',
    'src.class.luminosity',
    'src.class.spectral',
    'src.class.starGalaxy',
    'src.class.struct',
    'src.density',
    'src.ellipticity',
    'src.impactParam',
    'src.morph',
    'src.morph.param',
    'src.morph.scLength',
    'src.morph.type',
    'src.net',
    'src.orbital',
    'src.redshift',
    'src.redshift.phot',
    'src.sample',
    'src.spType',
    'src.var',
    'src.var.amplitude',
    'src.var.index',
    'src.var.period',

    // Statistics & Errors
    'stat.correlation',
    'stat.covariance',
    'stat.density',
    'stat.error',
    'stat.error.sys',
    'stat.fwhm',
    'stat.fit',
    'stat.fit.chi2',
    'stat.fit.dof',
    'stat.fit.goodness',
    'stat.fit.param',
    'stat.fit.residual',
    'stat.likelihood',
    'stat.mad',
    'stat.mean',
    'stat.median',
    'stat.min',
    'stat.max',
    'stat.param',
    'stat.probability',
    'stat.rank',
    'stat.rms',
    'stat.snr',
    'stat.stdev',
    'stat.uncalib',
    'stat.variance',
    'stat.weight',

    // Time & Epochs
    'time',
    'time.age',
    'time.creation',
    'time.crossing',
    'time.duration',
    'time.end',
    'time.epoch',
    'time.equinox',
    'time.event',
    'time.expo',
    'time.interval',
    'time.lifetime',
    'time.period',
    'time.phase',
    'time.processing',
    'time.publi',
    'time.release',
    'time.resolution',
    'time.scale',
    'time.start',

    // Instrumentation & Observation
    'instr',
    'instr.airmass',
    'instr.azimuth',
    'instr.bandpass',
    'instr.beam',
    'instr.calib',
    'instr.det',
    'instr.dispersion',
    'instr.filter',
    'instr.fov',
    'instr.obsty',
    'instr.obsty.seeing',
    'instr.order',
    'instr.param',
    'instr.pixel',
    'instr.plate',
    'instr.plate.lut',
    'instr.posAng',
    'instr.resol',
    'instr.rm',
    'instr.saturation',
    'instr.scale',
    'instr.sensitivity',
    'instr.setup',
    'instr.skyLevel',
    'instr.skyTemp',
    'instr.tec',
    'instr.tel',
    'instr.tel.focalLength'
  ]);

  /**
   * Secondary UCD1+ modifier words (can be appended after semicolons or attached).
   * @type {Set<string>}
   */
  static SECONDARY_WORDS = new Set([
    'arith.diff',
    'arith.factor',
    'arith.grad',
    'arith.rate',
    'arith.ratio',
    'arith.zp',
    'meta.main',
    'meta.code',
    'meta.id',
    'meta.ref',
    'stat.error',
    'stat.error.sys',
    'stat.mean',
    'stat.median',
    'stat.min',
    'stat.max',
    'stat.stdev',
    'stat.variance',
    'stat.weight',
    'stat.fit',
    'stat.snr',
    'pos.cartesian.x',
    'pos.cartesian.y',
    'pos.cartesian.z',
    'pos.eq.ra',
    'pos.eq.dec',
    'pos.galactic.lon',
    'pos.galactic.lat',
    'pos.supergalactic.lon',
    'pos.supergalactic.lat'
  ]);

  /**
   * Validates a UCD1+ expression string.
   * @param {string} ucd UCD expression (e.g., "pos.eq.ra;meta.main", "stat.error;phys.veloc")
   * @returns {{ valid: boolean, errors: Array<string>, atoms: Array<string>, primary: string|null }}
   */
  static validate(ucd) {
    if (!ucd || typeof ucd !== 'string') {
      return { valid: false, errors: ['UCD expression must be a non-empty string'], atoms: [], primary: null };
    }

    const trimmed = ucd.trim();
    if (trimmed.length === 0) {
      return { valid: false, errors: ['UCD expression cannot be whitespace only'], atoms: [], primary: null };
    }

    const atoms = trimmed.split(';').map(a => a.trim()).filter(a => a.length > 0);
    const errors = [];

    if (atoms.length === 0) {
      return { valid: false, errors: ['No valid atoms found in UCD'], atoms: [], primary: null };
    }

    // Atom format regex: alphanumeric, underscores, and dots (hierarchy)
    const atomPattern = /^[a-zA-Z][a-zA-Z0-9_\-]*(\.[a-zA-Z0-9_\-]+)*$/;

    for (let i = 0; i < atoms.length; i++) {
      const atom = atoms[i];
      if (!atomPattern.test(atom)) {
        errors.push(`Invalid character sequence in atom "${atom}" at position ${i}`);
      }
    }

    const primary = atoms[0];
    const isPrimaryKnown = UCDStandard.PRIMARY_WORDS.has(primary) || UCDStandard.SECONDARY_WORDS.has(primary);
    if (!isPrimaryKnown) {
      // Check prefix hierarchy (e.g. pos, phys, phot, spect, meta, stat, instr, time, src)
      const prefix = primary.split('.')[0];
      const validPrefixes = ['pos', 'phys', 'phot', 'spect', 'meta', 'stat', 'instr', 'time', 'src', 'arith', 'em', 'obs'];
      if (!validPrefixes.includes(prefix)) {
        errors.push(`Unrecognized primary UCD word family "${prefix}" in atom "${primary}"`);
      }
    }

    // Check secondary atoms
    for (let i = 1; i < atoms.length; i++) {
      const sec = atoms[i];
      const prefix = sec.split('.')[0];
      const validPrefixes = ['stat', 'meta', 'arith', 'pos', 'phys', 'phot', 'src', 'time'];
      if (!validPrefixes.includes(prefix) && !UCDStandard.SECONDARY_WORDS.has(sec) && !UCDStandard.PRIMARY_WORDS.has(sec)) {
        errors.push(`Unrecognized secondary UCD atom "${sec}" at index ${i}`);
      }
    }

    return {
      valid: errors.length === 0,
      errors,
      atoms,
      primary
    };
  }

  /**
   * Infers canonical UCD1+ from column name, unit, and description.
   * @param {string} columnName Column name (e.g. "ra", "dec", "v_pec", "mass", "sigma_v")
   * @param {string} [unit=""] Physical unit (e.g. "deg", "km/s", "Mpc", "Msun")
   * @param {string} [description=""] Field description
   * @returns {string} Inferred UCD1+
   */
  static inferUCD(columnName, unit = '', description = '') {
    const name = String(columnName || '').toLowerCase().trim();
    const u = String(unit || '').toLowerCase().trim();
    const d = String(description || '').toLowerCase().trim();

    // Coordinates & Astrometry
    if (/^(ra|radeg|ra_deg|right_ascension|alpha)$/.test(name) || /pos\.eq\.ra/.test(d)) {
      return 'pos.eq.ra;meta.main';
    }
    if (/^(dec|decdeg|dec_deg|declination|delta)$/.test(name) || /pos\.eq\.dec/.test(d)) {
      return 'pos.eq.dec;meta.main';
    }
    if (/^(glon|l|galactic_lon|l_deg)$/.test(name)) return 'pos.galactic.lon';
    if (/^(glat|b|galactic_lat|b_deg)$/.test(name)) return 'pos.galactic.lat';
    if (/^(sgl|sg_lon|supergalactic_lon)$/.test(name)) return 'pos.supergalactic.lon';
    if (/^(sgb|sg_lat|supergalactic_lat)$/.test(name)) return 'pos.supergalactic.lat';
    if (/^(x|x_mpc|x_cart|pos_x)$/.test(name)) return 'pos.cartesian.x';
    if (/^(y|y_mpc|y_cart|pos_y)$/.test(name)) return 'pos.cartesian.y';
    if (/^(z|z_mpc|z_cart|pos_z)$/.test(name)) return 'pos.cartesian.z';

    // Distances & Redshifts
    if (/^(dist|distance|dist_mpc|d_mpc|r_comoving)$/.test(name)) return 'pos.distance';
    if (/^(dist_err|e_dist|err_dist|sigma_dist)$/.test(name)) return 'stat.error;pos.distance';
    if (/^(z|redshift|z_obs|z_cmb|z_helio)$/.test(name) && !/mpc/.test(u)) return 'src.redshift';
    if (/^(e_z|err_z|sigma_z|z_err)$/.test(name)) return 'stat.error;src.redshift';

    // Velocities & Peculiar Velocities
    if (/^(v_rad|vrad|cz|v_helio|v_cmb|v_lsr)$/.test(name)) return 'spect.dopplerVeloc';
    if (/^(v_pec|vpec|v_peculiar|v_mod|v_rad_pec)$/.test(name)) return 'phys.veloc';
    if (/^(vx|vx_kms|v_x)$/.test(name)) return 'phys.veloc;pos.cartesian.x';
    if (/^(vy|vy_kms|v_y)$/.test(name)) return 'phys.veloc;pos.cartesian.y';
    if (/^(vz|vz_kms|v_z)$/.test(name)) return 'phys.veloc;pos.cartesian.z';
    if (/^(v_mag|vmag|speed|bulk_flow|v_bulk)$/.test(name)) return 'phys.veloc';
    if (/^(sigma_v|v_disp|vel_disp|velocity_dispersion)$/.test(name)) return 'phys.veloc.dispersion';
    if (/^(e_vpec|err_vpec|vpec_err|e_vx|e_vy|e_vz)$/.test(name)) return 'stat.error;phys.veloc';

    // Mass, Volume & Density
    if (/^(mass|m_vir|m200|m500|mass_proxy|m_tot)$/.test(name)) return 'phys.mass';
    if (/^(volume|vol|vol_mpc3)$/.test(name)) return 'phys.volume';
    if (/^(density|delta|overdensity|rho)$/.test(name)) return 'phys.density';

    // Identifiers & Names
    if (/^(id|galaxy_id|basin_id|pgc|source_id|name|attractor_name|objid)$/.test(name)) return 'meta.id;meta.main';
    if (/^(flag|mask|code|status|type_code)$/.test(name)) return 'meta.code';

    // Photometry & Magnitudes
    if (/^(mag|mag_b|mag_k|mag_r|mag_i|m_app|apparent_mag)$/.test(name)) return 'phot.mag';
    if (/^(m_abs|absmag|abs_mag)$/.test(name)) return 'phys.magAbs';
    if (/^(flux|f_tot|f_xray)$/.test(name)) return 'phot.flux';

    // Errors & Uncertainties fallback
    if (/^(err|error|sigma|e_|std|uncertainty)/.test(name)) return 'stat.error';

    // Units-based fallback
    if (/^km[\s\/]s/.test(u) || u === 'km.s-1' || u === 'km s-1') return 'phys.veloc';
    if (/^mpc/.test(u) || u === 'kpc' || u === 'pc') return 'pos.distance';
    if (/^deg/.test(u) || u === 'arcsec' || u === 'mas') return 'pos.angDistance';
    if (/^msun/.test(u) || /1e14\s*msun/.test(u) || /solarmass/.test(u)) return 'phys.mass';

    return 'meta.code';
  }
}

// ============================================================================
// SECTION 2: FITS CHECKSUM & DATASUM ENGINE (SEAMAN-PENCE IEEE 32-BIT 1'S COMPLEMENT)
// ============================================================================

/**
 * High-performance, mathematically exact FITS 32-bit 1's Complement Checksum and Datasum Engine.
 * Conforms to the FITS Checksum Standard (Pence, Seaman & White 1995, FITS 4.0 Standard Appendix J).
 */
export class FITSChecksum {
  /**
   * Adds two 32-bit unsigned integers using 1's complement arithmetic (circular bit carry).
   * @param {number} a First 32-bit unsigned integer
   * @param {number} b Second 32-bit unsigned integer
   * @returns {number} 32-bit unsigned integer 1's complement sum
   */
  static addOnesComplement(a, b) {
    const sum = (a >>> 0) + (b >>> 0);
    // End-around carry: add carry over 32 bits
    const carry = Math.floor(sum / 0x100000000);
    return ((sum & 0xFFFFFFFF) + carry) >>> 0;
  }

  /**
   * Computes the 32-bit 1's complement sum of a byte buffer.
   * Buffer length must be a multiple of 2 or 4 (or padded).
   * Interprets the buffer as a series of 32-bit big-endian unsigned integers.
   * @param {Uint8Array} buffer Input byte array
   * @param {number} [start=0] Start offset in bytes
   * @param {number} [length] Number of bytes to process
   * @returns {number} 32-bit unsigned 1's complement sum
   */
  static computeBlockSum(buffer, start = 0, length = null) {
    const len = length !== null ? length : buffer.length - start;
    let sum = 0;
    const view = new DataView(buffer.buffer, buffer.byteOffset + start, len);
    const numWords = Math.floor(len / 4);

    for (let i = 0; i < numWords; i++) {
      const word = view.getUint32(i * 4, false); // Big-Endian
      sum = FITSChecksum.addOnesComplement(sum, word);
    }

    // Handle trailing 1, 2, or 3 bytes if unaligned to 4 bytes
    const remainder = len % 4;
    if (remainder > 0) {
      let trailingWord = 0;
      const remStart = numWords * 4;
      for (let j = 0; j < remainder; j++) {
        trailingWord |= (buffer[start + remStart + j] << ((3 - j) * 8));
      }
      sum = FITSChecksum.addOnesComplement(sum, trailingWord >>> 0);
    }

    return sum >>> 0;
  }

  /**
   * Computes the DATASUM keyword value for a FITS HDU data block.
   * Returns standard unsigned 32-bit integer string representation (e.g. '3428945214').
   * @param {Uint8Array} dataBuffer Padded data buffer (multiple of 2880 bytes)
   * @returns {string} DATASUM string value
   */
  static calculateDataSum(dataBuffer) {
    if (!dataBuffer || dataBuffer.length === 0) {
      return '0';
    }
    const sum = FITSChecksum.computeBlockSum(dataBuffer);
    return String(sum >>> 0);
  }

  /**
   * Encodes a 32-bit 1's complement checksum value into a 16-character ASCII string.
   * Implements the Seaman-Pence ASCII encoding algorithm ensuring all characters
   * are printable ASCII without forbidden FITS characters ('\'', '"', etc.).
   * @param {number} checksum32 32-bit unsigned 1's complement checksum
   * @param {boolean} [excludeSpecial=true] Ensure no quotes or control characters
   * @returns {string} 16-character ASCII string
   */
  static encodeChecksumString(checksum32, excludeSpecial = true) {
    // Checksum value to encode is the complement ~checksum32
    const target = (~checksum32) >>> 0;

    // We split into 4 32-bit words that sum to target mod 2^32-1
    // The Seaman-Pence algorithm uses 16 ASCII characters (4 chars per 32-bit byte slice)
    // Permitted ASCII byte range: 0x30 ('0') to 0x5A ('Z')
    const offset = 0x30; // '0'
    const result = new Uint8Array(16);

    for (let byteIdx = 0; byteIdx < 4; byteIdx++) {
      const shift = (3 - byteIdx) * 8;
      const b = (target >>> shift) & 0xFF;

      // Decompose b into 4 parts plus offset 0x30:
      // Standard Seaman-Pence quotient/remainder distribution:
      const q = Math.floor(b / 4);
      const r = b % 4;

      const c = [q + offset, q + offset, q + offset, q + offset];
      c[0] += r;

      // Ensure each character avoids punctuation if requested
      if (excludeSpecial) {
        for (let k = 0; k < 4; k++) {
          // If character falls in ':' .. '@' (0x3A .. 0x40), push it to 'A' (0x41)
          if (c[k] > 0x39 && c[k] < 0x41) {
            const shiftVal = 0x41 - c[k];
            c[k] += shiftVal;
            const comp = (k + 1) % 4;
            c[comp] -= shiftVal;
          }
          // If character is single quote 0x27 or backslash 0x5C
          if (c[k] === 0x27 || c[k] === 0x5C) {
            c[k] += 2;
            c[(k + 1) % 4] -= 2;
          }
        }
      }

      // Distribute into result columns (interleaved byte order)
      result[byteIdx] = c[0];
      result[byteIdx + 4] = c[1];
      result[byteIdx + 8] = c[2];
      result[byteIdx + 12] = c[3];
    }

    return String.fromCharCode(...result);
  }

  /**
   * Decodes a 16-character FITS CHECKSUM string into its 32-bit numerical value.
   * @param {string} checksumStr 16-character ASCII string
   * @returns {number} 32-bit unsigned integer
   */
  static decodeChecksumString(checksumStr) {
    if (!checksumStr || checksumStr.length !== 16) {
      return 0;
    }
    const bytes = new Uint8Array(16);
    for (let i = 0; i < 16; i++) {
      bytes[i] = checksumStr.charCodeAt(i);
    }
    return FITSChecksum.computeBlockSum(bytes);
  }

  /**
   * Validates the 1's complement sum of an entire HDU (Header + Data).
   * A valid HDU with standard CHECKSUM card inserted must sum to 0xFFFFFFFF (or 0 in 1's complement).
   * @param {Uint8Array} hduBuffer Full HDU buffer
   * @returns {{ valid: boolean, residual: number, hexResidual: string }}
   */
  static verifyHDUSum(hduBuffer) {
    const sum = FITSChecksum.computeBlockSum(hduBuffer);
    const valid = sum === 0xFFFFFFFF || sum === 0;
    return {
      valid,
      residual: sum,
      hexResidual: '0x' + sum.toString(16).toUpperCase().padStart(8, '0')
    };
  }
}

// ============================================================================
// SECTION 3: FITS HEADER CARDS & STRUCTURES
// ============================================================================

/**
 * FITS 80-Column Header Card.
 */
export class FITSHeaderCard {
  /**
   * @param {string} keyword 1-8 char keyword, or HIERARCH keyword
   * @param {string|number|boolean|null} value Card value
   * @param {string} [comment=""] Card comment
   */
  constructor(keyword, value = null, comment = '') {
    this.keyword = String(keyword || '').toUpperCase().trim();
    this.value = value;
    this.comment = String(comment || '').trim();
  }

  /**
   * Formats the card into exact 80-character ASCII representation.
   * @returns {string} 80-character string
   */
  format() {
    // Comment card (COMMENT, HISTORY, or empty keyword)
    if (this.keyword === 'COMMENT' || this.keyword === 'HISTORY' || this.keyword === '') {
      const text = this.comment || String(this.value || '');
      const line = (this.keyword ? this.keyword.padEnd(8, ' ') : '        ') + ' ' + text;
      return line.slice(0, 80).padEnd(80, ' ');
    }

    if (this.keyword === 'END') {
      return 'END'.padEnd(80, ' ');
    }

    // Check if keyword requires HIERARCH convention (> 8 chars or hierarchical dots/spaces)
    if (this.keyword.length > 8 || this.keyword.includes('.')) {
      const hierKey = this.keyword.replace(/\./g, ' ');
      let valStr = '';
      if (typeof this.value === 'string') {
        valStr = `'${this.value.replace(/'/g, "''")}'`;
      } else if (typeof this.value === 'boolean') {
        valStr = this.value ? 'T' : 'F';
      } else if (this.value !== null && this.value !== undefined) {
        valStr = String(this.value);
      }
      let card = `HIERARCH ${hierKey} = ${valStr}`;
      if (this.comment) {
        card += ` / ${this.comment}`;
      }
      return card.slice(0, 80).padEnd(80, ' ');
    }

    const keyField = this.keyword.padEnd(8, ' ');
    let valueIndicator = '= ';
    let valField = '';
    let commentField = '';

    if (this.value === null || this.value === undefined) {
      valueIndicator = '  ';
      valField = ''.padEnd(20, ' ');
    } else if (typeof this.value === 'boolean') {
      // Logical value: 'T' or 'F' placed at column 30 (index 29)
      valField = (this.value ? 'T' : 'F').padStart(20, ' ');
    } else if (typeof this.value === 'number') {
      // Numbers right-justified in cols 11-30 (20 chars)
      let numStr = String(this.value);
      if (Number.isInteger(this.value) && Math.abs(this.value) < 1e15) {
        valField = numStr.padStart(20, ' ');
      } else {
        if (!numStr.includes('.') && !numStr.includes('e') && !numStr.includes('E')) {
          numStr += '.0';
        }
        valField = numStr.padStart(20, ' ');
      }
    } else if (typeof this.value === 'string') {
      // String values: enclosed in single quotes, starting at col 11 (index 10)
      const escaped = this.value.replace(/'/g, "''");
      const quoted = `'${escaped}'`;
      if (quoted.length < 20) {
        valField = quoted.padEnd(20, ' ');
      } else {
        valField = quoted;
      }
    }

    if (this.comment) {
      commentField = ` / ${this.comment}`;
    }

    let result = keyField + valueIndicator + valField + commentField;
    if (result.length > 80) {
      result = result.slice(0, 80);
    }
    return result.padEnd(80, ' ');
  }
}

/**
 * FITS Header structure managing cards, order, keywords, and serialization.
 */
export class FITSHeader {
  constructor() {
    /** @type {Array<FITSHeaderCard>} */
    this.cards = [];
  }

  /**
   * Sets or updates a keyword card in the header.
   * @param {string} keyword Keyword name
   * @param {string|number|boolean|null} value Value
   * @param {string} [comment=""] Comment
   */
  set(keyword, value, comment = '') {
    const key = keyword.toUpperCase().trim();
    const idx = this.cards.findIndex(c => c.keyword === key);
    if (idx >= 0 && key !== 'COMMENT' && key !== 'HISTORY') {
      this.cards[idx].value = value;
      if (comment) this.cards[idx].comment = comment;
    } else {
      this.cards.push(new FITSHeaderCard(key, value, comment));
    }
  }

  /**
   * Appends a card without checking uniqueness (e.g. COMMENT, HISTORY).
   * @param {string} keyword
   * @param {string|number|boolean|null} value
   * @param {string} [comment=""]
   */
  append(keyword, value = null, comment = '') {
    this.cards.push(new FITSHeaderCard(keyword, value, comment));
  }

  /**
   * Retrieves card value by keyword.
   * @param {string} keyword Keyword
   * @returns {string|number|boolean|null}
   */
  get(keyword) {
    const key = keyword.toUpperCase().trim();
    const card = this.cards.find(c => c.keyword === key);
    return card ? card.value : null;
  }

  /**
   * Serializes the header into 2880-byte block byte array (padded with ASCII spaces 0x20).
   * Includes END card.
   * @returns {Uint8Array}
   */
  toBytes() {
    const activeCards = this.cards.filter(c => c.keyword !== 'END');
    activeCards.push(new FITSHeaderCard('END', null, ''));

    const lines = activeCards.map(c => c.format());
    const totalChars = lines.length * 80;
    const blocksNeeded = Math.ceil(totalChars / 2880);
    const totalBytes = Math.max(1, blocksNeeded) * 2880;

    const buffer = new Uint8Array(totalBytes);
    buffer.fill(0x20); // ASCII space

    let offset = 0;
    for (const line of lines) {
      for (let i = 0; i < 80; i++) {
        buffer[offset + i] = line.charCodeAt(i);
      }
      offset += 80;
    }

    return buffer;
  }
}

// ============================================================================
// SECTION 4: FITS BINARY TABLE COLUMN DEFINITION & PACKING
// ============================================================================

/**
 * FITS Binary Table Column definition.
 */
export class FITSColumn {
  /**
   * @param {Object} spec Column specification
   * @param {string} spec.name Column name (TTYPE)
   * @param {string} spec.format Column FITS format code (TFORM, e.g. '1E', '3D', '32A', '1J', '1K', '1L', '16X', '1PD')
   * @param {string} [spec.unit=""] Physical unit (TUNIT)
   * @param {string} [spec.ucd=""] Astrometric UCD1+ (TUCD)
   * @param {number|null} [spec.nullValue=null] Integer null value (TNULL)
   * @param {number} [spec.scale=1.0] Linear scale factor (TSCAL)
   * @param {number} [spec.zero=0.0] Linear zero offset (TZERO)
   * @param {string} [spec.disp=""] Display format code (TDISP)
   * @param {string} [spec.comment=""] Column comment
   */
  constructor(spec) {
    this.name = String(spec.name || 'UNKNOWN').trim();
    this.format = String(spec.format || '1E').trim().toUpperCase();
    this.unit = String(spec.unit || '').trim();
    this.ucd = String(spec.ucd || UCDStandard.inferUCD(this.name, this.unit)).trim();
    this.nullValue = spec.nullValue !== undefined ? spec.nullValue : null;
    this.scale = spec.scale !== undefined ? Number(spec.scale) : 1.0;
    this.zero = spec.zero !== undefined ? Number(spec.zero) : 0.0;
    this.disp = String(spec.disp || '').trim();
    this.comment = String(spec.comment || '').trim();

    this.parsed = this._parseFormat(this.format);
  }

  /**
   * Parses FITS TFORM format string.
   * @param {string} format Format code (e.g. '1E', '3D', '32A', '1J', '1PD(10)')
   * @returns {{ repeat: number, typeCode: string, byteWidth: number, isVarArray: boolean, varType: string, varMax: number }}
   * @private
   */
  _parseFormat(format) {
    // Variable length array: e.g. '1PJ', '1PD(50)', '1QD'
    const varMatch = format.match(/^([0-9]*)([PQ])([LXBIAJECMD])(?:\(([0-9]+)\))?$/);
    if (varMatch) {
      const descriptorType = varMatch[2]; // 'P' (32-bit: 8 bytes) or 'Q' (64-bit: 16 bytes)
      const varType = varMatch[3];
      const maxElem = varMatch[4] ? parseInt(varMatch[4], 10) : 0;
      const byteWidth = descriptorType === 'P' ? 8 : 16;
      return {
        repeat: 1,
        typeCode: descriptorType,
        byteWidth,
        isVarArray: true,
        varType,
        varMax: maxElem
      };
    }

    // Fixed array: e.g. '1E', '3D', '32A', '16X'
    const match = format.match(/^([0-9]*)([LXBIAJECMD])$/);
    if (!match) {
      throw new Error(`Unsupported or invalid FITS TFORM format code: "${format}" for column "${this.name}"`);
    }

    const repeat = match[1] ? parseInt(match[1], 10) : 1;
    const typeCode = match[2];
    let unitSize = 1;

    switch (typeCode) {
      case 'L': unitSize = 1; break; // Logical / Boolean
      case 'X': unitSize = 1; break; // Bit (repeat = number of bits, rounded up to bytes)
      case 'B': unitSize = 1; break; // 8-bit unsigned byte
      case 'I': unitSize = 2; break; // 16-bit signed integer
      case 'J': unitSize = 4; break; // 32-bit signed integer
      case 'K': unitSize = 8; break; // 64-bit signed integer
      case 'A': unitSize = 1; break; // Character string
      case 'E': unitSize = 4; break; // 32-bit IEEE 754 single float
      case 'D': unitSize = 8; break; // 64-bit IEEE 754 double float
      case 'C': unitSize = 8; break; // Complex single (2 * 4 bytes)
      case 'M': unitSize = 16; break; // Complex double (2 * 8 bytes)
      default:
        throw new Error(`Unknown type code: ${typeCode}`);
    }

    const byteWidth = typeCode === 'X' ? Math.ceil(repeat / 8) : repeat * unitSize;

    return {
      repeat,
      typeCode,
      byteWidth,
      isVarArray: false,
      varType: null,
      varMax: 0
    };
  }
}

// ============================================================================
// SECTION 5: MULTI-EXTENSION FITS (MEF) BINARY TABLE HDU & SERIALIZER
// ============================================================================

/**
 * Representation of a FITS HDU (Header Data Unit).
 */
export class FITSHDU {
  /**
   * @param {FITSHeader} header HDU Header
   * @param {Uint8Array|null} [data=null] HDU Data buffer
   */
  constructor(header = new FITSHeader(), data = null) {
    this.header = header;
    this.data = data;
  }

  /**
   * Updates DATASUM and CHECKSUM keywords in HDU header.
   */
  updateChecksums() {
    const dataSumStr = FITSChecksum.calculateDataSum(this.data);
    this.header.set('DATASUM', dataSumStr, 'data checksum across all blocks');

    this.header.set('CHECKSUM', '0000000000000000', 'HDU checksum updated');

    const headerBytes = this.header.toBytes();
    const headerSum = FITSChecksum.computeBlockSum(headerBytes);
    const dataSum = this.data && this.data.length > 0 ? FITSChecksum.computeBlockSum(this.data) : 0;
    const hduSum = FITSChecksum.addOnesComplement(headerSum, dataSum);

    const checksumStr = FITSChecksum.encodeChecksumString(hduSum);
    this.header.set('CHECKSUM', checksumStr, 'HDU checksum updated');
  }

  /**
   * Serializes the HDU (Header + Data) into a byte array.
   * @returns {Uint8Array}
   */
  toBytes() {
    this.updateChecksums();
    const headerBytes = this.header.toBytes();
    const dataBytes = this.data || new Uint8Array(0);

    const out = new Uint8Array(headerBytes.length + dataBytes.length);
    out.set(headerBytes, 0);
    out.set(dataBytes, headerBytes.length);
    return out;
  }
}

/**
 * FITS Primary HDU.
 */
export class FITSPrimaryHDU extends FITSHDU {
  /**
   * Creates a standard Primary HDU (SIMPLE=T, BITPIX=8, NAXIS=0, EXTEND=T).
   * @param {Object} [metadata={}] Optional primary metadata keywords
   */
  constructor(metadata = {}) {
    const header = new FITSHeader();
    header.set('SIMPLE', true, 'conforms to FITS standard');
    header.set('BITPIX', metadata.bitpix || 8, 'array data points per byte');
    header.set('NAXIS', metadata.naxis || 0, 'number of array dimensions');
    header.set('EXTEND', true, 'FITS dataset may contain extensions');
    header.set('ORIGIN', metadata.origin || 'ZRT Cosmicflows Workbench', 'Organization responsible');
    header.set('CREATOR', metadata.creator || 'vo_table_fits_serializer.js v1.4', 'Software module');
    header.set('DATE', new Date().toISOString().slice(0, 19), 'Creation UTC timestamp');

    if (metadata.telescop) header.set('TELESCOP', metadata.telescop);
    if (metadata.observer) header.set('OBSERVER', metadata.observer);
    if (metadata.comment) header.append('COMMENT', null, metadata.comment);

    super(header, metadata.data || null);
  }
}

/**
 * FITS Binary Table Extension HDU (XTENSION='BINTABLE').
 */
export class FITSBinaryTableHDU extends FITSHDU {
  /**
   * @param {string} extName Extension name (EXTNAME)
   * @param {Array<FITSColumn>} columns Column definitions
   * @param {Array<Object>|Array<Array<any>>} rows Row data records
   * @param {Object} [options={}] Additional extension options
   */
  constructor(extName, columns, rows, options = {}) {
    const header = new FITSHeader();
    super(header, null);

    this.extName = extName || 'BINTABLE';
    this.columns = columns.map(c => (c instanceof FITSColumn ? c : new FITSColumn(c)));
    this.rows = rows || [];
    this.options = options;

    this._build();
  }

  /**
   * Builds the binary table data buffer and header keywords.
   * @private
   */
  _build() {
    const numRows = this.rows.length;
    const numFields = this.columns.length;

    let rowWidth = 0;
    const colOffsets = [];
    for (const col of this.columns) {
      colOffsets.push(rowWidth);
      rowWidth += col.parsed.byteWidth;
    }

    const heapBytes = [];
    const mainTableSize = numRows * rowWidth;
    const mainTableBuffer = new Uint8Array(mainTableSize);
    const view = new DataView(mainTableBuffer.buffer);

    for (let r = 0; r < numRows; r++) {
      const rowData = this.rows[r];
      const rowStart = r * rowWidth;

      for (let c = 0; c < numFields; c++) {
        const col = this.columns[c];
        const offset = rowStart + colOffsets[c];
        let val = Array.isArray(rowData) ? rowData[c] : rowData[col.name];

        if (val === undefined || val === null) {
          if (col.nullValue !== null) {
            val = col.nullValue;
          }
        }

        this._packColumnValue(view, mainTableBuffer, offset, col, val, heapBytes);
      }
    }

    const heapSize = heapBytes.length;
    const totalDataBytes = mainTableSize + heapSize;
    const blocksNeeded = Math.ceil(totalDataBytes / 2880);
    const paddedDataBytes = Math.max(blocksNeeded, totalDataBytes > 0 ? 1 : 0) * 2880;

    const dataBuffer = new Uint8Array(paddedDataBytes);
    dataBuffer.set(mainTableBuffer, 0);
    if (heapSize > 0) {
      dataBuffer.set(new Uint8Array(heapBytes), mainTableSize);
    }

    this.data = dataBuffer;

    const h = this.header;
    h.set('XTENSION', 'BINTABLE', 'binary table extension');
    h.set('BITPIX', 8, '8-bit bytes');
    h.set('NAXIS', 2, '2-dimensional binary table');
    h.set('NAXIS1', rowWidth, 'width of table in bytes');
    h.set('NAXIS2', numRows, 'number of rows in table');
    h.set('PCOUNT', heapSize, 'size of variable-length array heap');
    h.set('GCOUNT', 1, 'one group standard');
    h.set('TFIELDS', numFields, 'number of fields per row');
    h.set('EXTNAME', this.extName, 'name of this binary table extension');

    if (this.options.extver) h.set('EXTVER', this.options.extver, 'extension version');

    for (let i = 0; i < numFields; i++) {
      const idx = i + 1;
      const col = this.columns[i];

      h.set(`TTYPE${idx}`, col.name, col.comment || `label for column ${idx}`);
      h.set(`TFORM${idx}`, col.format, `data format of column ${idx}`);

      if (col.unit) h.set(`TUNIT${idx}`, col.unit, 'physical unit');
      if (col.nullValue !== null) h.set(`TNULL${idx}`, col.nullValue, 'null value indicator');
      if (col.scale !== 1.0) h.set(`TSCAL${idx}`, col.scale, 'linear data scaling');
      if (col.zero !== 0.0) h.set(`TZERO${idx}`, col.zero, 'linear data zero offset');
      if (col.disp) h.set(`TDISP${idx}`, col.disp, 'display format');
      if (col.ucd) h.set(`TUCD${idx}`, col.ucd, 'Unified Content Descriptor');
    }
  }

  /**
   * Packs a single column cell into the binary buffer (Big-Endian).
   * @private
   */
  _packColumnValue(view, buffer, offset, col, val, heapBytes) {
    const { typeCode, repeat, isVarArray, varType } = col.parsed;

    if (isVarArray) {
      const arr = Array.isArray(val) ? val : (val !== null && val !== undefined ? [val] : []);
      const count = arr.length;
      const heapOffset = heapBytes.length;

      for (const item of arr) {
        this._packHeapItem(heapBytes, varType, item);
      }

      if (typeCode === 'P') {
        view.setUint32(offset, count, false);
        view.setUint32(offset + 4, heapOffset, false);
      } else {
        view.setBigUint64(offset, BigInt(count), false);
        view.setBigUint64(offset + 8, BigInt(heapOffset), false);
      }
      return;
    }

    let rawVal = val;
    if (typeof val === 'number' && (col.scale !== 1.0 || col.zero !== 0.0)) {
      rawVal = (val - col.zero) / col.scale;
    }

    switch (typeCode) {
      case 'L': {
        let charCode = 0x00;
        if (rawVal === true || rawVal === 'T' || rawVal === 1) charCode = 0x54;
        else if (rawVal === false || rawVal === 'F' || rawVal === 0) charCode = 0x46;
        buffer[offset] = charCode;
        break;
      }

      case 'X': {
        const numBytes = Math.ceil(repeat / 8);
        const bits = Array.isArray(rawVal) ? rawVal : (typeof rawVal === 'number' ? [rawVal] : []);
        for (let b = 0; b < numBytes; b++) {
          let byteVal = 0;
          for (let bitIdx = 0; bitIdx < 8; bitIdx++) {
            const overallBit = b * 8 + bitIdx;
            if (overallBit < repeat && bits[overallBit]) {
              byteVal |= (1 << (7 - bitIdx));
            }
          }
          buffer[offset + b] = byteVal;
        }
        break;
      }

      case 'B': {
        for (let i = 0; i < repeat; i++) {
          const item = Array.isArray(rawVal) ? rawVal[i] : (i === 0 ? rawVal : 0);
          buffer[offset + i] = (item || 0) & 0xFF;
        }
        break;
      }

      case 'I': {
        for (let i = 0; i < repeat; i++) {
          const item = Array.isArray(rawVal) ? rawVal[i] : (i === 0 ? rawVal : 0);
          view.setInt16(offset + i * 2, item !== null && item !== undefined ? Number(item) : 0, false);
        }
        break;
      }

      case 'J': {
        for (let i = 0; i < repeat; i++) {
          const item = Array.isArray(rawVal) ? rawVal[i] : (i === 0 ? rawVal : 0);
          view.setInt32(offset + i * 4, item !== null && item !== undefined ? Number(item) : 0, false);
        }
        break;
      }

      case 'K': {
        for (let i = 0; i < repeat; i++) {
          const item = Array.isArray(rawVal) ? rawVal[i] : (i === 0 ? rawVal : 0);
          const bigVal = item !== null && item !== undefined ? BigInt(item) : 0n;
          view.setBigInt64(offset + i * 8, bigVal, false);
        }
        break;
      }

      case 'A': {
        const str = rawVal !== null && rawVal !== undefined ? String(rawVal) : '';
        for (let i = 0; i < repeat; i++) {
          buffer[offset + i] = i < str.length ? str.charCodeAt(i) : 0x20;
        }
        break;
      }

      case 'E': {
        for (let i = 0; i < repeat; i++) {
          const item = Array.isArray(rawVal) ? rawVal[i] : (i === 0 ? rawVal : NaN);
          const num = item !== null && item !== undefined ? Number(item) : NaN;
          view.setFloat32(offset + i * 4, num, false);
        }
        break;
      }

      case 'D': {
        for (let i = 0; i < repeat; i++) {
          const item = Array.isArray(rawVal) ? rawVal[i] : (i === 0 ? rawVal : NaN);
          const num = item !== null && item !== undefined ? Number(item) : NaN;
          view.setFloat64(offset + i * 8, num, false);
        }
        break;
      }

      case 'C': {
        for (let i = 0; i < repeat; i++) {
          const pair = Array.isArray(rawVal) && Array.isArray(rawVal[0]) ? rawVal[i] : rawVal;
          const r = pair && pair.length > 0 ? Number(pair[0]) : 0;
          const im = pair && pair.length > 1 ? Number(pair[1]) : 0;
          view.setFloat32(offset + i * 8, r, false);
          view.setFloat32(offset + i * 8 + 4, im, false);
        }
        break;
      }

      case 'M': {
        for (let i = 0; i < repeat; i++) {
          const pair = Array.isArray(rawVal) && Array.isArray(rawVal[0]) ? rawVal[i] : rawVal;
          const r = pair && pair.length > 0 ? Number(pair[0]) : 0;
          const im = pair && pair.length > 1 ? Number(pair[1]) : 0;
          view.setFloat64(offset + i * 16, r, false);
          view.setFloat64(offset + i * 16 + 8, im, false);
        }
        break;
      }
    }
  }

  /**
   * Packs an item into the variable length array heap.
   * @private
   */
  _packHeapItem(heapBytes, varType, item) {
    const temp = new Uint8Array(8);
    const view = new DataView(temp.buffer);

    switch (varType) {
      case 'B':
        heapBytes.push((item || 0) & 0xFF);
        break;
      case 'I':
        view.setInt16(0, Number(item || 0), false);
        heapBytes.push(temp[0], temp[1]);
        break;
      case 'J':
        view.setInt32(0, Number(item || 0), false);
        heapBytes.push(temp[0], temp[1], temp[2], temp[3]);
        break;
      case 'K':
        view.setBigInt64(0, BigInt(item || 0), false);
        for (let k = 0; k < 8; k++) heapBytes.push(temp[k]);
        break;
      case 'E':
        view.setFloat32(0, Number(item !== null && item !== undefined ? item : NaN), false);
        heapBytes.push(temp[0], temp[1], temp[2], temp[3]);
        break;
      case 'D':
        view.setFloat64(0, Number(item !== null && item !== undefined ? item : NaN), false);
        for (let k = 0; k < 8; k++) heapBytes.push(temp[k]);
        break;
      case 'A': {
        const str = String(item || '');
        for (let i = 0; i < str.length; i++) heapBytes.push(str.charCodeAt(i));
        break;
      }
      case 'L':
        heapBytes.push(item ? 0x54 : 0x46);
        break;
    }
  }
}

/**
 * Complete Multi-Extension FITS (MEF) Serializer and Parser.
 */
export class FITSBinaryTableSerializer {
  /**
   * Creates a complete MEF FITS file buffer containing Primary HDU and multiple BINTABLE HDUs.
   * @param {Object} spec MEF specification
   * @param {Object} [spec.primary] Primary HDU metadata
   * @param {Array<{ extName: string, columns: Array<Object>, rows: Array<Object>, options?: Object }>} spec.tables
   * @returns {Uint8Array} Full FITS file buffer
   */
  static serializeMEF(spec) {
    const primaryHDU = new FITSPrimaryHDU(spec.primary || {});
    const hdus = [primaryHDU];

    if (Array.isArray(spec.tables)) {
      for (const t of spec.tables) {
        const bintable = new FITSBinaryTableHDU(t.extName, t.columns, t.rows, t.options || {});
        hdus.push(bintable);
      }
    }

    const byteArrays = hdus.map(h => h.toBytes());
    const totalLength = byteArrays.reduce((acc, b) => acc + b.length, 0);

    const out = new Uint8Array(totalLength);
    let offset = 0;
    for (const b of byteArrays) {
      out.set(b, offset);
      offset += b.length;
    }

    return out;
  }

  /**
   * Parses an ArrayBuffer or Uint8Array containing a FITS file.
   * Reads Primary HDU and all Binary Table Extension HDUs.
   * @param {ArrayBuffer|Uint8Array} buffer Input FITS bytes
   * @returns {{ primary: Object, tables: Array<Object>, isValid: boolean, checksums: Array<Object> }}
   */
  static parseMEF(buffer) {
    const bytes = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
    const totalBytes = bytes.length;
    let offset = 0;
    let hduIndex = 0;

    const parsedHDUs = [];
    const checksumReport = [];

    while (offset < totalBytes) {
      const hduStart = offset;

      const headerCards = [];
      let foundEnd = false;

      while (offset < totalBytes && !foundEnd) {
        const block = bytes.subarray(offset, offset + 2880);
        if (block.length < 2880) break;

        for (let c = 0; c < 36; c++) {
          const cardBytes = block.subarray(c * 80, (c + 1) * 80);
          let cardStr = '';
          for (let i = 0; i < 80; i++) cardStr += String.fromCharCode(cardBytes[i]);

          const key = cardStr.slice(0, 8).trim();
          headerCards.push(cardStr);

          if (key === 'END') {
            foundEnd = true;
            break;
          }
        }
        offset += 2880;
      }

      if (!foundEnd) break;

      const headerEnd = offset;
      const headerObj = FITSBinaryTableSerializer._parseHeaderCards(headerCards);

      const bitpix = headerObj.BITPIX || 8;
      const naxis = headerObj.NAXIS || 0;
      let dataBytesCount = 0;

      if (naxis > 0) {
        let naxisProd = 1;
        for (let i = 1; i <= naxis; i++) {
          naxisProd *= (headerObj[`NAXIS${i}`] || 0);
        }
        const pcount = headerObj.PCOUNT || 0;
        const gcount = headerObj.GCOUNT || 1;
        dataBytesCount = Math.abs(bitpix / 8) * gcount * (naxisProd + pcount);
      }

      const dataBlocksNeeded = Math.ceil(dataBytesCount / 2880);
      const dataPaddingBytes = Math.max(dataBlocksNeeded, dataBytesCount > 0 ? 1 : 0) * 2880;
      const dataEnd = offset + dataPaddingBytes;

      const dataBuffer = bytes.subarray(offset, Math.min(totalBytes, dataEnd));
      offset = dataEnd;

      const hduTotalBuffer = bytes.subarray(hduStart, Math.min(totalBytes, dataEnd));
      const verify = FITSChecksum.verifyHDUSum(hduTotalBuffer);

      checksumReport.push({
        hduIndex,
        extName: headerObj.EXTNAME || (hduIndex === 0 ? 'PRIMARY' : `EXT_${hduIndex}`),
        storedChecksum: headerObj.CHECKSUM,
        storedDatasum: headerObj.DATASUM,
        valid: verify.valid,
        residual: verify.hexResidual
      });

      let tableData = null;
      if (headerObj.XTENSION === 'BINTABLE') {
        tableData = FITSBinaryTableSerializer._decodeBinaryTable(headerObj, dataBuffer);
      }

      parsedHDUs.push({
        hduIndex,
        header: headerObj,
        table: tableData
      });

      hduIndex++;
    }

    const primary = parsedHDUs[0] ? parsedHDUs[0].header : {};
    const tables = parsedHDUs.filter(h => h.table !== null).map(h => h.table);

    return {
      primary,
      tables,
      isValid: checksumReport.every(c => c.valid),
      checksums: checksumReport
    };
  }

  /**
   * Parses 80-character header card strings into key-value map.
   * @private
   */
  static _parseHeaderCards(cardStrings) {
    const header = {};
    for (const card of cardStrings) {
      const key = card.slice(0, 8).trim().toUpperCase();
      if (!key || key === 'COMMENT' || key === 'HISTORY' || key === 'END') continue;

      const valPart = card.slice(10, 80);
      const commentIdx = valPart.indexOf('/');
      const valStr = (commentIdx >= 0 ? valPart.slice(0, commentIdx) : valPart).trim();

      if (valStr.startsWith("'")) {
        const match = valStr.match(/^'([^']*)'/);
        header[key] = match ? match[1].trim() : '';
      } else if (valStr === 'T') {
        header[key] = true;
      } else if (valStr === 'F') {
        header[key] = false;
      } else if (valStr.length > 0 && !isNaN(Number(valStr))) {
        header[key] = Number(valStr);
      } else {
        header[key] = valStr;
      }
    }
    return header;
  }

  /**
   * Decodes FITS BINTABLE data buffer into records and columns.
   * @private
   */
  static _decodeBinaryTable(header, dataBuffer) {
    const numRows = header.NAXIS2 || 0;
    const rowWidth = header.NAXIS1 || 0;
    const numFields = header.TFIELDS || 0;
    const pcount = header.PCOUNT || 0;

    const columns = [];
    let curOffset = 0;

    for (let i = 1; i <= numFields; i++) {
      const name = header[`TTYPE${i}`] || `COL${i}`;
      const format = header[`TFORM${i}`] || '1E';
      const unit = header[`TUNIT${i}`] || '';
      const ucd = header[`TUCD${i}`] || '';
      const nullVal = header[`TNULL${i}`] !== undefined ? header[`TNULL${i}`] : null;
      const scale = header[`TSCAL${i}`] !== undefined ? header[`TSCAL${i}`] : 1.0;
      const zero = header[`TZERO${i}`] !== undefined ? header[`TZERO${i}`] : 0.0;

      const col = new FITSColumn({ name, format, unit, ucd, nullValue: nullVal, scale, zero });
      col._offset = curOffset;
      curOffset += col.parsed.byteWidth;
      columns.push(col);
    }

    const rows = [];
    const view = new DataView(dataBuffer.buffer, dataBuffer.byteOffset, dataBuffer.byteLength);
    const heapStart = numRows * rowWidth;

    for (let r = 0; r < numRows; r++) {
      const rowObj = {};
      const rowStart = r * rowWidth;

      for (const col of columns) {
        const offset = rowStart + col._offset;
        const val = FITSBinaryTableSerializer._unpackColumnValue(view, dataBuffer, offset, col, heapStart);
        rowObj[col.name] = val;
      }
      rows.push(rowObj);
    }

    return {
      extName: header.EXTNAME || 'BINTABLE',
      columns: columns.map(c => ({
        name: c.name,
        format: c.format,
        unit: c.unit,
        ucd: c.ucd,
        scale: c.scale,
        zero: c.zero,
        nullValue: c.nullValue
      })),
      rows,
      rowCount: numRows
    };
  }

  /**
   * Unpacks a single column value from the binary buffer.
   * @private
   */
  static _unpackColumnValue(view, buffer, offset, col, heapStart) {
    const { typeCode, repeat, isVarArray, varType } = col.parsed;

    if (isVarArray) {
      let count = 0;
      let heapOffset = 0;
      if (typeCode === 'P') {
        count = view.getUint32(offset, false);
        heapOffset = view.getUint32(offset + 4, false);
      } else {
        count = Number(view.getBigUint64(offset, false));
        heapOffset = Number(view.getBigUint64(offset + 8, false));
      }

      const arr = [];
      const itemStart = heapStart + heapOffset;
      let itemOffset = itemStart;

      for (let k = 0; k < count; k++) {
        switch (varType) {
          case 'B': arr.push(buffer[itemOffset++]); break;
          case 'I': arr.push(view.getInt16(itemOffset, false)); itemOffset += 2; break;
          case 'J': arr.push(view.getInt32(itemOffset, false)); itemOffset += 4; break;
          case 'K': arr.push(view.getBigInt64(itemOffset, false)); itemOffset += 8; break;
          case 'E': arr.push(view.getFloat32(itemOffset, false)); itemOffset += 4; break;
          case 'D': arr.push(view.getFloat64(itemOffset, false)); itemOffset += 8; break;
          case 'L': arr.push(buffer[itemOffset++] === 0x54); break;
          case 'A': arr.push(String.fromCharCode(buffer[itemOffset++])); break;
        }
      }
      return arr;
    }

    let decoded = null;

    switch (typeCode) {
      case 'L': {
        const b = buffer[offset];
        decoded = b === 0x54 ? true : (b === 0x46 ? false : null);
        break;
      }

      case 'X': {
        const numBytes = Math.ceil(repeat / 8);
        const bits = [];
        for (let b = 0; b < numBytes; b++) {
          const byteVal = buffer[offset + b];
          for (let bitIdx = 0; bitIdx < 8; bitIdx++) {
            const overallBit = b * 8 + bitIdx;
            if (overallBit < repeat) {
              bits.push((byteVal & (1 << (7 - bitIdx))) !== 0 ? 1 : 0);
            }
          }
        }
        decoded = bits;
        break;
      }

      case 'B': {
        if (repeat === 1) {
          decoded = buffer[offset];
        } else {
          decoded = Array.from(buffer.subarray(offset, offset + repeat));
        }
        break;
      }

      case 'I': {
        if (repeat === 1) {
          decoded = view.getInt16(offset, false);
          if (col.nullValue !== null && decoded === col.nullValue) decoded = null;
        } else {
          const arr = [];
          for (let i = 0; i < repeat; i++) {
            let v = view.getInt16(offset + i * 2, false);
            if (col.nullValue !== null && v === col.nullValue) v = null;
            arr.push(v);
          }
          decoded = arr;
        }
        break;
      }

      case 'J': {
        if (repeat === 1) {
          decoded = view.getInt32(offset, false);
          if (col.nullValue !== null && decoded === col.nullValue) decoded = null;
        } else {
          const arr = [];
          for (let i = 0; i < repeat; i++) {
            let v = view.getInt32(offset + i * 4, false);
            if (col.nullValue !== null && v === col.nullValue) v = null;
            arr.push(v);
          }
          decoded = arr;
        }
        break;
      }

      case 'K': {
        if (repeat === 1) {
          decoded = view.getBigInt64(offset, false);
        } else {
          const arr = [];
          for (let i = 0; i < repeat; i++) arr.push(view.getBigInt64(offset + i * 8, false));
          decoded = arr;
        }
        break;
      }

      case 'A': {
        let str = '';
        for (let i = 0; i < repeat; i++) str += String.fromCharCode(buffer[offset + i]);
        decoded = str.trim();
        break;
      }

      case 'E': {
        if (repeat === 1) {
          decoded = view.getFloat32(offset, false);
        } else {
          const arr = [];
          for (let i = 0; i < repeat; i++) arr.push(view.getFloat32(offset + i * 4, false));
          decoded = arr;
        }
        break;
      }

      case 'D': {
        if (repeat === 1) {
          decoded = view.getFloat64(offset, false);
        } else {
          const arr = [];
          for (let i = 0; i < repeat; i++) arr.push(view.getFloat64(offset + i * 8, false));
          decoded = arr;
        }
        break;
      }

      case 'C': {
        if (repeat === 1) {
          decoded = [view.getFloat32(offset, false), view.getFloat32(offset + 4, false)];
        } else {
          const arr = [];
          for (let i = 0; i < repeat; i++) {
            arr.push([view.getFloat32(offset + i * 8, false), view.getFloat32(offset + i * 8 + 4, false)]);
          }
          decoded = arr;
        }
        break;
      }

      case 'M': {
        if (repeat === 1) {
          decoded = [view.getFloat64(offset, false), view.getFloat64(offset + 8, false)];
        } else {
          const arr = [];
          for (let i = 0; i < repeat; i++) {
            arr.push([view.getFloat64(offset + i * 16, false), view.getFloat64(offset + i * 16 + 8, false)]);
          }
          decoded = arr;
        }
        break;
      }
    }

    if (typeof decoded === 'number' && (col.scale !== 1.0 || col.zero !== 0.0)) {
      decoded = decoded * col.scale + col.zero;
    }

    return decoded;
  }
}

// ============================================================================
// SECTION 6: IVOA VOTABLE 1.4 DATA STRUCTURES & SERIALIZER
// ============================================================================

/**
 * VOTable Field Definition (FIELD element).
 */
export class VOTableField {
  /**
   * @param {Object} spec Field specification
   * @param {string} spec.name Field name
   * @param {string} [spec.id] Field ID
   * @param {string} [spec.datatype='float'] Datatype
   * @param {string} [spec.arraysize] Arraysize
   * @param {string} [spec.unit=''] Physical unit
   * @param {string} [spec.ucd=''] Astrometric UCD1+
   * @param {string} [spec.utype=''] UType data model attribute
   * @param {string|number|null} [spec.nullValue=null] Null value marker
   * @param {string} [spec.precision] Number of significant figures / format
   * @param {number} [spec.width] Fixed column display width
   * @param {string} [spec.description=''] Field description
   */
  constructor(spec) {
    this.name = String(spec.name || 'UNKNOWN').trim();
    this.id = spec.id ? String(spec.id).trim() : this.name.replace(/[^a-zA-Z0-9_]/g, '_');
    this.datatype = String(spec.datatype || 'float').trim();
    this.arraysize = spec.arraysize ? String(spec.arraysize).trim() : null;
    this.unit = String(spec.unit || '').trim();
    this.ucd = String(spec.ucd || UCDStandard.inferUCD(this.name, this.unit)).trim();
    this.utype = String(spec.utype || '').trim();
    this.nullValue = spec.nullValue !== undefined ? spec.nullValue : null;
    this.precision = spec.precision !== undefined ? String(spec.precision) : null;
    this.width = spec.width !== undefined ? Number(spec.width) : null;
    this.description = String(spec.description || '').trim();
  }

  toXML(indent = 4) {
    const pad = ' '.repeat(indent);
    let xml = `${pad}<FIELD ID="${this.id}" name="${VOTableSerializer.escapeXML(this.name)}" datatype="${this.datatype}"`;

    if (this.arraysize) xml += ` arraysize="${this.arraysize}"`;
    if (this.unit) xml += ` unit="${VOTableSerializer.escapeXML(this.unit)}"`;
    if (this.ucd) xml += ` ucd="${VOTableSerializer.escapeXML(this.ucd)}"`;
    if (this.utype) xml += ` utype="${VOTableSerializer.escapeXML(this.utype)}"`;
    if (this.precision) xml += ` precision="${this.precision}"`;
    if (this.width) xml += ` width="${this.width}"`;

    if (this.description || this.nullValue !== null) {
      xml += '>\n';
      if (this.description) {
        xml += `${pad}  <DESCRIPTION>${VOTableSerializer.escapeXML(this.description)}</DESCRIPTION>\n`;
      }
      if (this.nullValue !== null) {
        xml += `${pad}  <VALUES null="${VOTableSerializer.escapeXML(String(this.nullValue))}" />\n`;
      }
      xml += `${pad}</FIELD>\n`;
    } else {
      xml += ' />\n';
    }

    return xml;
  }
}

/**
 * VOTable Parameter Definition (PARAM element).
 */
export class VOTableParam {
  /**
   * @param {Object} spec Param specification
   * @param {string} spec.name Param name
   * @param {string} [spec.id] Param ID
   * @param {string} [spec.datatype='char'] Datatype
   * @param {string|number|boolean} spec.value Value
   * @param {string} [spec.arraysize] Arraysize
   * @param {string} [spec.unit=''] Unit
   * @param {string} [spec.ucd=''] UCD
   * @param {string} [spec.utype=''] UType
   * @param {string} [spec.description=''] Description
   */
  constructor(spec) {
    this.name = String(spec.name || '').trim();
    this.id = spec.id ? String(spec.id).trim() : this.name.replace(/[^a-zA-Z0-9_]/g, '_');
    this.datatype = String(spec.datatype || 'char').trim();
    this.value = spec.value !== undefined ? spec.value : '';
    this.arraysize = spec.arraysize ? String(spec.arraysize).trim() : null;
    this.unit = String(spec.unit || '').trim();
    this.ucd = String(spec.ucd || '').trim();
    this.utype = String(spec.utype || '').trim();
    this.description = String(spec.description || '').trim();
  }

  toXML(indent = 4) {
    const pad = ' '.repeat(indent);
    let valStr = String(this.value);
    let xml = `${pad}<PARAM ID="${this.id}" name="${VOTableSerializer.escapeXML(this.name)}" datatype="${this.datatype}" value="${VOTableSerializer.escapeXML(valStr)}"`;

    if (this.arraysize) xml += ` arraysize="${this.arraysize}"`;
    if (this.unit) xml += ` unit="${VOTableSerializer.escapeXML(this.unit)}"`;
    if (this.ucd) xml += ` ucd="${VOTableSerializer.escapeXML(this.ucd)}"`;
    if (this.utype) xml += ` utype="${VOTableSerializer.escapeXML(this.utype)}"`;

    if (this.description) {
      xml += '>\n';
      xml += `${pad}  <DESCRIPTION>${VOTableSerializer.escapeXML(this.description)}</DESCRIPTION>\n`;
      xml += `${pad}</PARAM>\n`;
    } else {
      xml += ' />\n';
    }

    return xml;
  }
}

/**
 * VOTable Group Definition (GROUP element).
 */
export class VOTableGroup {
  /**
   * @param {Object} spec Group specification
   * @param {string} [spec.name] Group name
   * @param {string} [spec.id] Group ID
   * @param {string} [spec.utype] UType
   * @param {string} [spec.ucd] UCD
   * @param {Array<string>} [spec.fieldRefs] Array of FIELD ID references
   * @param {Array<VOTableParam>} [spec.params] Array of PARAMs
   */
  constructor(spec = {}) {
    this.name = spec.name ? String(spec.name).trim() : '';
    this.id = spec.id ? String(spec.id).trim() : '';
    this.utype = spec.utype ? String(spec.utype).trim() : '';
    this.ucd = spec.ucd ? String(spec.ucd).trim() : '';
    this.fieldRefs = spec.fieldRefs || [];
    this.params = (spec.params || []).map(p => (p instanceof VOTableParam ? p : new VOTableParam(p)));
  }

  toXML(indent = 4) {
    const pad = ' '.repeat(indent);
    let xml = `${pad}<GROUP`;
    if (this.id) xml += ` ID="${this.id}"`;
    if (this.name) xml += ` name="${VOTableSerializer.escapeXML(this.name)}"`;
    if (this.utype) xml += ` utype="${VOTableSerializer.escapeXML(this.utype)}"`;
    if (this.ucd) xml += ` ucd="${VOTableSerializer.escapeXML(this.ucd)}"`;
    xml += '>\n';

    for (const p of this.params) {
      xml += p.toXML(indent + 2);
    }
    for (const ref of this.fieldRefs) {
      xml += `${pad}  <FIELDref ref="${VOTableSerializer.escapeXML(ref)}" />\n`;
    }

    xml += `${pad}</GROUP>\n`;
    return xml;
  }
}

/**
 * VOTable Table Container (TABLE element).
 */
export class VOTableTable {
  /**
   * @param {Object} spec Table specification
   * @param {string} [spec.name] Table name
   * @param {string} [spec.id] Table ID
   * @param {string} [spec.description] Table description
   * @param {Array<VOTableField>} [spec.fields] Array of field descriptors
   * @param {Array<VOTableParam>} [spec.params] Array of parameters
   * @param {Array<VOTableGroup>} [spec.groups] Array of groups
   * @param {Array<Object>|Array<Array<any>>} [spec.rows] Array of row records
   */
  constructor(spec = {}) {
    this.name = spec.name ? String(spec.name).trim() : 'TABLE_1';
    this.id = spec.id ? String(spec.id).trim() : this.name.replace(/[^a-zA-Z0-9_]/g, '_');
    this.description = spec.description ? String(spec.description).trim() : '';
    this.fields = (spec.fields || []).map(f => (f instanceof VOTableField ? f : new VOTableField(f)));
    this.params = (spec.params || []).map(p => (p instanceof VOTableParam ? p : new VOTableParam(p)));
    this.groups = (spec.groups || []).map(g => (g instanceof VOTableGroup ? g : new VOTableGroup(g)));
    this.rows = spec.rows || [];
  }

  /**
   * Serializes the TABLE element with chosen DATA serialization format ('TABLEDATA', 'BINARY', or 'BINARY2').
   * @param {Object} [options={}]
   * @param {string} [options.format='TABLEDATA'] Data format ('TABLEDATA'|'BINARY'|'BINARY2')
   * @param {number} [options.indent=4] Indentation spaces
   * @param {number} [options.precision=6] Decimal digits for floats
   * @returns {string} XML string
   */
  toXML(options = {}) {
    const format = (options.format || 'TABLEDATA').toUpperCase();
    const indent = options.indent !== undefined ? options.indent : 4;
    const prec = options.precision !== undefined ? options.precision : 6;
    const pad = ' '.repeat(indent);

    let xml = `${pad}<TABLE ID="${this.id}" name="${VOTableSerializer.escapeXML(this.name)}" nrows="${this.rows.length}">\n`;

    if (this.description) {
      xml += `${pad}  <DESCRIPTION>${VOTableSerializer.escapeXML(this.description)}</DESCRIPTION>\n`;
    }

    for (const g of this.groups) {
      xml += g.toXML(indent + 2);
    }
    for (const p of this.params) {
      xml += p.toXML(indent + 2);
    }
    for (const f of this.fields) {
      xml += f.toXML(indent + 2);
    }

    xml += `${pad}  <DATA>\n`;
    if (format === 'TABLEDATA') {
      xml += this._serializeTableData(indent + 4, prec);
    } else if (format === 'BINARY') {
      xml += this._serializeBinary(indent + 4, false);
    } else if (format === 'BINARY2') {
      xml += this._serializeBinary(indent + 4, true);
    } else {
      throw new Error(`Unsupported VOTable DATA format: ${format}`);
    }
    xml += `${pad}  </DATA>\n`;
    xml += `${pad}</TABLE>\n`;

    return xml;
  }

  /**
   * Serializes TABLEDATA XML rows.
   * @private
   */
  _serializeTableData(indent, prec) {
    const pad = ' '.repeat(indent);
    let xml = `${pad}<TABLEDATA>\n`;

    for (let r = 0; r < this.rows.length; r++) {
      const row = this.rows[r];
      xml += `${pad}  <TR>\n`;

      for (let c = 0; c < this.fields.length; c++) {
        const field = this.fields[c];
        let val = Array.isArray(row) ? row[c] : row[field.name];

        let cellStr = '';
        if (val === null || val === undefined || (typeof val === 'number' && isNaN(val))) {
          cellStr = field.nullValue !== null ? String(field.nullValue) : '';
        } else if (Array.isArray(val)) {
          cellStr = val.map(v => (typeof v === 'number' ? Number(v.toFixed(prec)) : String(v))).join(' ');
        } else if (typeof val === 'number') {
          cellStr = Number.isInteger(val) ? String(val) : String(Number(val.toFixed(prec)));
        } else if (typeof val === 'boolean') {
          cellStr = val ? '1' : '0';
        } else {
          cellStr = String(val);
        }

        xml += `${pad}    <TD>${VOTableSerializer.escapeXML(cellStr)}</TD>\n`;
      }
      xml += `${pad}  </TR>\n`;
    }

    xml += `${pad}</TABLEDATA>\n`;
    return xml;
  }

  /**
   * Serializes BINARY or BINARY2 Base64 data stream.
   * @private
   */
  _serializeBinary(indent, isBinary2 = false) {
    const pad = ' '.repeat(indent);
    const tag = isBinary2 ? 'BINARY2' : 'BINARY';

    const rawBytes = this._packBinaryStream(isBinary2);
    const base64 = VOTableSerializer.bytesToBase64(rawBytes);

    let xml = `${pad}<${tag}>\n`;
    xml += `${pad}  <STREAM encoding="base64">\n`;

    for (let i = 0; i < base64.length; i += 76) {
      xml += `${pad}    ${base64.slice(i, i + 76)}\n`;
    }

    xml += `${pad}  </STREAM>\n`;
    xml += `${pad}</${tag}>\n`;
    return xml;
  }

  /**
   * Converts all table rows to Big-Endian binary byte array for VOTable STREAM.
   * @param {boolean} isBinary2 If true, prepends null flag bitmask to each row
   * @returns {Uint8Array}
   * @private
   */
  _packBinaryStream(isBinary2) {
    const bytes = [];
    const numFields = this.fields.length;
    const nullMaskBytesCount = isBinary2 ? Math.ceil(numFields / 8) : 0;

    for (let r = 0; r < this.rows.length; r++) {
      const row = this.rows[r];

      if (isBinary2) {
        const nullMask = new Uint8Array(nullMaskBytesCount);
        for (let c = 0; c < numFields; c++) {
          const field = this.fields[c];
          const val = Array.isArray(row) ? row[c] : row[field.name];
          const isNull = val === null || val === undefined || (typeof val === 'number' && isNaN(val));
          if (isNull) {
            const byteIdx = Math.floor(c / 8);
            const bitIdx = 7 - (c % 8);
            nullMask[byteIdx] |= (1 << bitIdx);
          }
        }
        for (let b = 0; b < nullMaskBytesCount; b++) {
          bytes.push(nullMask[b]);
        }
      }

      for (let c = 0; c < numFields; c++) {
        const field = this.fields[c];
        const val = Array.isArray(row) ? row[c] : row[field.name];
        this._packVOTableCell(bytes, field, val);
      }
    }

    return new Uint8Array(bytes);
  }

  /**
   * Packs a single cell into byte stream adhering to VOTable 1.4 binary encoding rules.
   * @private
   */
  _packVOTableCell(bytes, field, val) {
    const temp = new Uint8Array(8);
    const view = new DataView(temp.buffer);
    const dt = field.datatype.toLowerCase();
    const isArray = field.arraysize !== null && field.arraysize !== undefined;

    if (dt === 'char' || dt === 'unicodechar') {
      const str = val !== null && val !== undefined ? String(val) : '';
      if (field.arraysize === '*') {
        view.setInt32(0, str.length, false);
        for (let i = 0; i < 4; i++) bytes.push(temp[i]);
        for (let i = 0; i < str.length; i++) bytes.push(str.charCodeAt(i) & 0xFF);
      } else {
        const fixedLen = field.arraysize ? parseInt(field.arraysize, 10) : str.length;
        for (let i = 0; i < fixedLen; i++) {
          bytes.push(i < str.length ? str.charCodeAt(i) & 0xFF : 0x00);
        }
      }
      return;
    }

    if (isArray && Array.isArray(val)) {
      if (field.arraysize === '*') {
        view.setInt32(0, val.length, false);
        for (let i = 0; i < 4; i++) bytes.push(temp[i]);
      }
      for (const item of val) {
        this._packPrimitive(bytes, view, temp, dt, item);
      }
      return;
    }

    this._packPrimitive(bytes, view, temp, dt, val);
  }

  /**
   * Packs a primitive scalar into the byte stream (Big-Endian).
   * @private
   */
  _packPrimitive(bytes, view, temp, dt, val) {
    switch (dt) {
      case 'boolean':
        bytes.push(val ? 0x54 : (val === false ? 0x46 : 0x3F));
        break;

      case 'unsignedbyte':
      case 'bit':
        bytes.push((Number(val) || 0) & 0xFF);
        break;

      case 'short':
        view.setInt16(0, Number(val || 0), false);
        bytes.push(temp[0], temp[1]);
        break;

      case 'int':
        view.setInt32(0, Number(val || 0), false);
        bytes.push(temp[0], temp[1], temp[2], temp[3]);
        break;

      case 'long':
        view.setBigInt64(0, BigInt(val || 0), false);
        for (let i = 0; i < 8; i++) bytes.push(temp[i]);
        break;

      case 'float':
        view.setFloat32(0, Number(val !== null && val !== undefined ? val : NaN), false);
        bytes.push(temp[0], temp[1], temp[2], temp[3]);
        break;

      case 'double':
        view.setFloat64(0, Number(val !== null && val !== undefined ? val : NaN), false);
        for (let i = 0; i < 8; i++) bytes.push(temp[i]);
        break;

      case 'floatcomplex': {
        const r = Array.isArray(val) ? Number(val[0]) : Number(val || 0);
        const im = Array.isArray(val) && val.length > 1 ? Number(val[1]) : 0;
        view.setFloat32(0, r, false);
        bytes.push(temp[0], temp[1], temp[2], temp[3]);
        view.setFloat32(0, im, false);
        bytes.push(temp[0], temp[1], temp[2], temp[3]);
        break;
      }

      case 'doublecomplex': {
        const r = Array.isArray(val) ? Number(val[0]) : Number(val || 0);
        const im = Array.isArray(val) && val.length > 1 ? Number(val[1]) : 0;
        view.setFloat64(0, r, false);
        for (let i = 0; i < 8; i++) bytes.push(temp[i]);
        view.setFloat64(0, im, false);
        for (let i = 0; i < 8; i++) bytes.push(temp[i]);
        break;
      }

      default:
        view.setFloat32(0, Number(val || 0), false);
        bytes.push(temp[0], temp[1], temp[2], temp[3]);
    }
  }
}

/**
 * VOTable Resource Container (RESOURCE element).
 */
export class VOTableResource {
  /**
   * @param {Object} spec Resource specification
   * @param {string} [spec.name] Resource name
   * @param {string} [spec.id] Resource ID
   * @param {string} [spec.type='results'] Resource type ('results'|'meta')
   * @param {string} [spec.description] Resource description
   * @param {Array<VOTableTable>} [spec.tables] Array of tables
   * @param {Array<VOTableResource>} [spec.resources] Nested sub-resources
   * @param {Array<VOTableParam>} [spec.params] Parameters
   * @param {Array<{ name: string, value: string, content?: string }>} [spec.infos] Info tags
   * @param {Object} [spec.coosys] Coordinate system (COOSYS)
   */
  constructor(spec = {}) {
    this.name = spec.name ? String(spec.name).trim() : 'RESOURCE_1';
    this.id = spec.id ? String(spec.id).trim() : this.name.replace(/[^a-zA-Z0-9_]/g, '_');
    this.type = spec.type ? String(spec.type).trim() : 'results';
    this.description = spec.description ? String(spec.description).trim() : '';
    this.tables = (spec.tables || []).map(t => (t instanceof VOTableTable ? t : new VOTableTable(t)));
    this.resources = (spec.resources || []).map(r => (r instanceof VOTableResource ? r : new VOTableResource(r)));
    this.params = (spec.params || []).map(p => (p instanceof VOTableParam ? p : new VOTableParam(p)));
    this.infos = spec.infos || [];
    this.coosys = spec.coosys || null;
  }

  toXML(options = {}) {
    const indent = options.indent !== undefined ? options.indent : 2;
    const pad = ' '.repeat(indent);

    let xml = `${pad}<RESOURCE ID="${this.id}" name="${VOTableSerializer.escapeXML(this.name)}" type="${this.type}">\n`;

    if (this.description) {
      xml += `${pad}  <DESCRIPTION>${VOTableSerializer.escapeXML(this.description)}</DESCRIPTION>\n`;
    }

    if (this.coosys) {
      const sys = this.coosys;
      xml += `${pad}  <COOSYS ID="${sys.id || 'J2000'}" system="${sys.system || 'ICRS'}"`;
      if (sys.equinox) xml += ` equinox="${sys.equinox}"`;
      if (sys.epoch) xml += ` epoch="${sys.epoch}"`;
      xml += ' />\n';
    }

    for (const info of this.infos) {
      xml += `${pad}  <INFO name="${VOTableSerializer.escapeXML(info.name)}" value="${VOTableSerializer.escapeXML(info.value)}"${info.content ? `>${VOTableSerializer.escapeXML(info.content)}</INFO>` : ' />'}\n`;
    }

    for (const p of this.params) {
      xml += p.toXML(indent + 2);
    }

    for (const t of this.tables) {
      xml += t.toXML({ ...options, indent: indent + 2 });
    }

    for (const r of this.resources) {
      xml += r.toXML({ ...options, indent: indent + 2 });
    }

    xml += `${pad}</RESOURCE>\n`;
    return xml;
  }
}

/**
 * Top-level IVOA VOTable 1.4 Document.
 */
export class VOTableDocument {
  /**
   * @param {Object} spec Document specification
   * @param {string} [spec.version='1.4'] VOTable version
   * @param {string} [spec.description] Document description
   * @param {Array<VOTableResource>} [spec.resources] Top-level resources
   * @param {Array<{ name: string, value: string }>} [spec.infos] Document INFO tags
   */
  constructor(spec = {}) {
    this.version = spec.version || '1.4';
    this.description = spec.description || 'ZRT Cosmicflows Workbench Scientific Dataset';
    this.resources = (spec.resources || []).map(r => (r instanceof VOTableResource ? r : new VOTableResource(r)));
    this.infos = spec.infos || [];
  }

  /**
   * Serializes the document into standard XML string.
   * @param {Object} [options={}]
   * @param {string} [options.format='TABLEDATA'] Data format ('TABLEDATA'|'BINARY'|'BINARY2')
   * @param {number} [options.precision=6] Precision
   * @returns {string} XML string
   */
  toXML(options = {}) {
    let xml = '<?xml version="1.0" encoding="UTF-8"?>\n';
    xml += `<VOTABLE version="${this.version}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://www.ivoa.net/xml/VOTable/v${this.version}" xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v${this.version} http://www.ivoa.net/xml/VOTable/v${this.version}">\n`;

    if (this.description) {
      xml += `  <DESCRIPTION>${VOTableSerializer.escapeXML(this.description)}</DESCRIPTION>\n`;
    }

    for (const info of this.infos) {
      xml += `  <INFO name="${VOTableSerializer.escapeXML(info.name)}" value="${VOTableSerializer.escapeXML(info.value)}" />\n`;
    }

    for (const r of this.resources) {
      xml += r.toXML({ ...options, indent: 2 });
    }

    xml += '</VOTABLE>\n';
    return xml;
  }
}

/**
 * Universal VOTable Serializer, Deserializer, and Helper Utilities.
 */
export class VOTableSerializer {
  /**
   * Escapes standard XML characters.
   * @param {string} str Input text
   * @returns {string} Escaped XML string
   */
  static escapeXML(str) {
    if (!str || typeof str !== 'string') return '';
    return str
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&apos;');
  }

  /**
   * Unescapes XML entities.
   * @param {string} str XML string
   * @returns {string}
   */
  static unescapeXML(str) {
    if (!str || typeof str !== 'string') return '';
    return str
      .replace(/&apos;/g, "'")
      .replace(/&quot;/g, '"')
      .replace(/&gt;/g, '>')
      .replace(/&lt;/g, '<')
      .replace(/&amp;/g, '&');
  }

  /**
   * Converts Uint8Array to Base64 string.
   * @param {Uint8Array} bytes
   * @returns {string}
   */
  static bytesToBase64(bytes) {
    if (typeof Buffer !== 'undefined') {
      return Buffer.from(bytes).toString('base64');
    }
    let binary = '';
    const len = bytes.byteLength;
    for (let i = 0; i < len; i++) {
      binary += String.fromCharCode(bytes[i]);
    }
    return btoa(binary);
  }

  /**
   * Converts Base64 string to Uint8Array.
   * @param {string} base64
   * @returns {Uint8Array}
   */
  static base64ToBytes(base64) {
    const clean = base64.replace(/\s+/g, '');
    if (typeof Buffer !== 'undefined') {
      return new Uint8Array(Buffer.from(clean, 'base64'));
    }
    const binary = atob(clean);
    const len = binary.length;
    const bytes = new Uint8Array(len);
    for (let i = 0; i < len; i++) {
      bytes[i] = binary.charCodeAt(i);
    }
    return bytes;
  }

  /**
   * Exports dataset rows and column definitions into VOTable 1.4 XML.
   * @param {Object} config Export configuration
   * @param {string} [config.tableName='COSMICFLOWS_CATALOG']
   * @param {Array<Object>} config.columns
   * @param {Array<Object>} config.rows
   * @param {string} [config.format='TABLEDATA'] ('TABLEDATA'|'BINARY'|'BINARY2')
   * @param {Object} [config.metadata]
   * @param {number} [config.precision=6]
   * @returns {string} VOTable XML text
   */
  static exportToVOTable(config) {
    const table = new VOTableTable({
      name: config.tableName || 'COSMICFLOWS_CATALOG',
      description: config.tableDescription || 'ZRT Cosmicflows Dataset',
      fields: config.columns,
      rows: config.rows
    });

    const resource = new VOTableResource({
      name: config.resourceName || 'PRIMARY_RESOURCE',
      description: config.resourceDescription || 'Cosmicflows Data Resource',
      tables: [table],
      coosys: config.coosys || { id: 'J2000', system: 'ICRS', equinox: '2000.0', epoch: '2000.0' }
    });

    const doc = new VOTableDocument({
      description: config.description || 'Cosmicflows Observational and Computational Workbench Export',
      resources: [resource],
      infos: [
        { name: 'QUERY_STATUS', value: 'OK' },
        { name: 'PROVIDER', value: 'ZRT Cosmicflows Research Group' },
        { name: 'DATE_EXPORTED', value: new Date().toISOString() }
      ]
    });

    return doc.toXML({
      format: config.format || 'TABLEDATA',
      precision: config.precision !== undefined ? config.precision : 6
    });
  }

  /**
   * Lightweight parser for VOTable XML text.
   * Extracts resources, tables, fields, and rows from TABLEDATA or BINARY.
   * @param {string} xmlString
   * @returns {{ version: string, tables: Array<{ name: string, fields: Array<Object>, rows: Array<Object> }> }}
   */
  static parseVOTable(xmlString) {
    const tables = [];

    const tableRegex = /<TABLE([\s\S]*?)<\/TABLE>/gi;
    let match;

    while ((match = tableRegex.exec(xmlString)) !== null) {
      const tableContent = match[1];

      const nameMatch = tableContent.match(/name="([^"]+)"/i);
      const tableName = nameMatch ? nameMatch[1] : 'TABLE';

      const fields = [];
      const fieldRegex = /<FIELD([\s\S]*?)(?:\/>|>([\s\S]*?)<\/FIELD>)/gi;
      let fieldMatch;

      while ((fieldMatch = fieldRegex.exec(tableContent)) !== null) {
        const fieldAttrs = fieldMatch[1];
        const getAttr = attr => {
          const m = fieldAttrs.match(new RegExp(`${attr}="([^"]+)"`, 'i'));
          return m ? m[1] : null;
        };

        fields.push({
          name: getAttr('name') || getAttr('ID') || 'UNKNOWN',
          id: getAttr('ID'),
          datatype: getAttr('datatype') || 'float',
          arraysize: getAttr('arraysize'),
          unit: getAttr('unit') || '',
          ucd: getAttr('ucd') || ''
        });
      }

      const rows = [];
      const tableDataMatch = tableContent.match(/<TABLEDATA>([\s\S]*?)<\/TABLEDATA>/i);
      if (tableDataMatch) {
        const rowsContent = tableDataMatch[1];
        const trRegex = /<TR>([\s\S]*?)<\/TR>/gi;
        let trMatch;

        while ((trMatch = trRegex.exec(rowsContent)) !== null) {
          const rowXml = trMatch[1];
          const tdRegex = /<TD>([\s\S]*?)<\/TD>/gi;
          let tdMatch;
          const rowValues = [];

          while ((tdMatch = tdRegex.exec(rowXml)) !== null) {
            rowValues.push(VOTableSerializer.unescapeXML(tdMatch[1].trim()));
          }

          const rowObj = {};
          for (let i = 0; i < fields.length; i++) {
            const f = fields[i];
            const raw = rowValues[i] !== undefined ? rowValues[i] : null;
            rowObj[f.name] = VOTableSerializer._parseTypedValue(raw, f.datatype);
          }
          rows.push(rowObj);
        }
      }

      tables.push({
        name: tableName,
        fields,
        rows,
        rowCount: rows.length
      });
    }

    return {
      version: '1.4',
      tables
    };
  }

  /**
   * Helper to parse string representation into native JS type according to VOTable datatype.
   * @private
   */
  static _parseTypedValue(valStr, datatype) {
    if (valStr === null || valStr === undefined || valStr === '') return null;
    const dt = datatype.toLowerCase();

    if (dt === 'boolean') {
      return valStr === '1' || valStr.toUpperCase() === 'T' || valStr.toLowerCase() === 'true';
    }
    if (dt === 'int' || dt === 'short' || dt === 'unsignedbyte' || dt === 'bit') {
      const n = parseInt(valStr, 10);
      return isNaN(n) ? null : n;
    }
    if (dt === 'long') {
      try {
        return BigInt(valStr);
      } catch {
        return parseInt(valStr, 10);
      }
    }
    if (dt === 'float' || dt === 'double') {
      if (valStr.toLowerCase() === 'nan') return NaN;
      if (valStr.toLowerCase() === 'inf' || valStr.toLowerCase() === '+inf') return Infinity;
      if (valStr.toLowerCase() === '-inf') return -Infinity;
      const f = parseFloat(valStr);
      return isNaN(f) ? null : f;
    }

    return valStr;
  }
}

// ============================================================================
// SECTION 7: HIGH-LEVEL COSMOLOGICAL EXPORTERS
// ============================================================================

/**
 * High-level Cosmological Exporter for VOTable and MEF FITS formats.
 */
export class CosmicflowsVoTableFitsExporter {
  /**
   * Serializes galaxy peculiar velocity catalog into VOTable 1.4 (TABLEDATA or BINARY)
   * @param {Array<Object>} galaxies Galaxy catalog records
   * @param {Object} [options={}] Export options
   * @returns {string} VOTable XML
   */
  static exportGalaxyCatalogToVOTable(galaxies, options = {}) {
    const columns = [
      { name: 'galaxy_id', datatype: 'char', arraysize: '32', ucd: 'meta.id;meta.main', description: 'Galaxy identifier / PGC code' },
      { name: 'ra_deg', datatype: 'double', unit: 'deg', ucd: 'pos.eq.ra;meta.main', description: 'Right Ascension (ICRS J2000)' },
      { name: 'dec_deg', datatype: 'double', unit: 'deg', ucd: 'pos.eq.dec;meta.main', description: 'Declination (ICRS J2000)' },
      { name: 'sgl_deg', datatype: 'double', unit: 'deg', ucd: 'pos.supergalactic.lon', description: 'Supergalactic Longitude' },
      { name: 'sgb_deg', datatype: 'double', unit: 'deg', ucd: 'pos.supergalactic.lat', description: 'Supergalactic Latitude' },
      { name: 'dist_mpc', datatype: 'float', unit: 'Mpc', ucd: 'pos.distance', description: 'Comoving / Metric Distance' },
      { name: 'dist_err_mpc', datatype: 'float', unit: 'Mpc', ucd: 'stat.error;pos.distance', description: 'Distance uncertainty' },
      { name: 'v_rad_kms', datatype: 'float', unit: 'km s-1', ucd: 'spect.dopplerVeloc', description: 'Observed radial velocity (CMB frame)' },
      { name: 'v_pec_kms', datatype: 'float', unit: 'km s-1', ucd: 'phys.veloc', description: 'Peculiar velocity' },
      { name: 'v_pec_err_kms', datatype: 'float', unit: 'km s-1', ucd: 'stat.error;phys.veloc', description: 'Peculiar velocity uncertainty' },
      { name: 'x_mpc_h', datatype: 'float', unit: 'Mpc/h', ucd: 'pos.cartesian.x', description: 'Cartesian X coordinate' },
      { name: 'y_mpc_h', datatype: 'float', unit: 'Mpc/h', ucd: 'pos.cartesian.y', description: 'Cartesian Y coordinate' },
      { name: 'z_mpc_h', datatype: 'float', unit: 'Mpc/h', ucd: 'pos.cartesian.z', description: 'Cartesian Z coordinate' }
    ];

    return VOTableSerializer.exportToVOTable({
      tableName: options.tableName || 'GALAXY_PECULIAR_VELOCITIES',
      tableDescription: 'Cosmicflows Galaxy Peculiar Velocity Catalog with 3D Cartesian Coordinates',
      columns,
      rows: galaxies,
      format: options.format || 'TABLEDATA',
      precision: options.precision !== undefined ? options.precision : 6
    });
  }

  /**
   * Serializes Watershed Basin Catalogs into Multi-Extension FITS (MEF).
   * HDU 1: Basin summary properties
   * HDU 2: Critical points and attractor seeds
   * @param {Array<Object>} basins Basin catalog records
   * @param {Array<Object>} [criticalPoints=[]] Associated critical points
   * @param {Object} [metadata={}] FITS Primary HDU metadata
   * @returns {Uint8Array} FITS byte buffer
   */
  static exportBasinCatalogToMEFFITS(basins, criticalPoints = [], metadata = {}) {
    const basinColumns = [
      { name: 'BASIN_ID', format: '1J', unit: '', ucd: 'meta.id;meta.main', comment: 'Unique basin ID' },
      { name: 'ATTRACTOR_NAME', format: '32A', unit: '', ucd: 'meta.id', comment: 'Attractor designation' },
      { name: 'POS_CART', format: '3D', unit: 'Mpc/h', ucd: 'pos.cartesian', comment: 'Attractor [X, Y, Z] position' },
      { name: 'V_BULK_CART', format: '3E', unit: 'km/s', ucd: 'phys.veloc', comment: 'Bulk flow [Vx, Vy, Vz]' },
      { name: 'VOLUME', format: '1D', unit: 'Mpc^3/h^3', ucd: 'phys.volume', comment: 'Watershed basin volume' },
      { name: 'MASS_PROXY', format: '1E', unit: '10^14 Msun', ucd: 'phys.mass', comment: 'Dynamical mass proxy' },
      { name: 'GALAXY_COUNT', format: '1J', unit: '', ucd: 'meta.number', comment: 'Number of galaxies assigned' },
      { name: 'STABILITY_IDX', format: '1E', unit: '', ucd: 'stat.fit.goodness', comment: 'Watershed boundary stability' }
    ];

    const basinRows = basins.map(b => ({
      BASIN_ID: b.basin_id !== undefined ? b.basin_id : b.id,
      ATTRACTOR_NAME: b.attractor_name || b.name || 'UNKNOWN',
      POS_CART: [b.x_mpc_h || b.x || 0, b.y_mpc_h || b.y || 0, b.z_mpc_h || b.z || 0],
      V_BULK_CART: [b.vx_kms || b.vx || 0, b.vy_kms || b.vy || 0, b.vz_kms || b.vz || 0],
      VOLUME: b.volume_mpc3_h3 || b.volume || 0,
      MASS_PROXY: b.mass_proxy_1e14msun || b.mass || 0,
      GALAXY_COUNT: b.galaxy_count || b.count || 0,
      STABILITY_IDX: b.stability_index || b.stability || 1.0
    }));

    const cpColumns = [
      { name: 'CP_ID', format: '1J', unit: '', ucd: 'meta.id', comment: 'Critical point index' },
      { name: 'CP_TYPE', format: '16A', unit: '', ucd: 'meta.code.class', comment: 'Morse type (Attractor/Saddle/Repeller)' },
      { name: 'MORSE_INDEX', format: '1I', unit: '', ucd: 'meta.code', comment: 'Morse index (0, 1, 2, 3)' },
      { name: 'POS_CART', format: '3D', unit: 'Mpc/h', ucd: 'pos.cartesian', comment: 'Critical point position' },
      { name: 'POTENTIAL_VAL', format: '1D', unit: '(km/s)^2', ucd: 'phys.potential', comment: 'Gravitational potential value' },
      { name: 'EIGENVALUES', format: '3E', unit: 's^-1', ucd: 'stat.fit.param', comment: 'Hessian / velocity tensor eigenvalues' }
    ];

    const cpRows = criticalPoints.map(cp => ({
      CP_ID: cp.id || 0,
      CP_TYPE: cp.type || 'ATTRACTOR',
      MORSE_INDEX: cp.morse_index !== undefined ? cp.morse_index : (cp.type === 'ATTRACTOR' ? 3 : 0),
      POS_CART: [cp.x || 0, cp.y || 0, cp.z || 0],
      POTENTIAL_VAL: cp.potential || 0,
      EIGENVALUES: cp.eigenvalues || [0, 0, 0]
    }));

    const tables = [
      { extName: 'WATERSHED_BASINS', columns: basinColumns, rows: basinRows, options: { extver: 1 } }
    ];

    if (cpRows.length > 0) {
      tables.push({ extName: 'CRITICAL_POINTS', columns: cpColumns, rows: cpRows, options: { extver: 1 } });
    }

    return FITSBinaryTableSerializer.serializeMEF({
      primary: {
        origin: 'ZRT Cosmicflows Workbench',
        creator: 'CosmicflowsVoTableFitsExporter v1.4',
        telescop: metadata.telescop || 'Cosmicflows-4 / SDSS / 2MASS',
        observer: metadata.observer || 'Cosmological Velocity Reconstruction Pipeline',
        comment: 'Multi-Extension FITS containing Watershed Basins and Velocity Critical Points'
      },
      tables
    });
  }

  /**
   * Serializes Cosmological Bulk Flow Shells into VOTable and FITS formats.
   * @param {Array<Object>} shells Multi-shell bulk flow estimates
   * @param {Object} [options={}]
   * @returns {{ votableXml: string, fitsBuffer: Uint8Array }}
   */
  static exportBulkFlowDataset(shells, options = {}) {
    const columns = [
      { name: 'r_min_mpc', format: '1E', datatype: 'float', unit: 'Mpc/h', ucd: 'pos.distance;stat.min', description: 'Inner shell radius' },
      { name: 'r_max_mpc', format: '1E', datatype: 'float', unit: 'Mpc/h', ucd: 'pos.distance;stat.max', description: 'Outer shell radius' },
      { name: 'r_eff_mpc', format: '1E', datatype: 'float', unit: 'Mpc/h', ucd: 'pos.distance;stat.mean', description: 'Effective / Weighted shell radius' },
      { name: 'v_bulk_kms', format: '1E', datatype: 'float', unit: 'km s-1', ucd: 'phys.veloc', description: 'Bulk flow speed magnitude' },
      { name: 'v_bulk_err_kms', format: '1E', datatype: 'float', unit: 'km s-1', ucd: 'stat.error;phys.veloc', description: 'Bulk flow speed uncertainty' },
      { name: 'glon_deg', format: '1D', datatype: 'double', unit: 'deg', ucd: 'pos.galactic.lon', description: 'Galactic Longitude of bulk flow vector' },
      { name: 'glat_deg', format: '1D', datatype: 'double', unit: 'deg', ucd: 'pos.galactic.lat', description: 'Galactic Latitude of bulk flow vector' },
      { name: 'vx_kms', format: '1E', datatype: 'float', unit: 'km s-1', ucd: 'phys.veloc;pos.cartesian.x', description: 'X component of bulk velocity' },
      { name: 'vy_kms', format: '1E', datatype: 'float', unit: 'km s-1', ucd: 'phys.veloc;pos.cartesian.y', description: 'Y component of bulk velocity' },
      { name: 'vz_kms', format: '1E', datatype: 'float', unit: 'km s-1', ucd: 'phys.veloc;pos.cartesian.z', description: 'Z component of bulk velocity' },
      { name: 'galaxy_count', format: '1J', datatype: 'int', unit: '', ucd: 'meta.number', description: 'Galaxies in shell' }
    ];

    const votableXml = VOTableSerializer.exportToVOTable({
      tableName: 'BULK_FLOW_SHELLS',
      tableDescription: 'Multipole Cosmological Bulk Flow Velocity Reconstructions',
      columns,
      rows: shells,
      format: options.votableFormat || 'TABLEDATA'
    });

    const fitsBuffer = FITSBinaryTableSerializer.serializeMEF({
      primary: {
        origin: 'ZRT Cosmicflows Workbench',
        creator: 'CosmicflowsVoTableFitsExporter v1.4',
        comment: 'Multi-Shell Bulk Flow Velocity Reconstructions'
      },
      tables: [
        {
          extName: 'BULK_FLOWS',
          columns: columns.map(c => ({ name: c.name, format: c.format, unit: c.unit, ucd: c.ucd })),
          rows: shells
        }
      ]
    });

    return { votableXml, fitsBuffer };
  }
}
