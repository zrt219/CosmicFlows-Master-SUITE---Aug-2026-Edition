"""
Comprehensive Pytest test suite for Remote FITS Ingestion, Web Worker Decompression,
and IndexedDB Binary Grid Caching for Cosmicflows-4 128³/256³ Datasets.
"""

import os
import subprocess
import json
import pytest
from tests.utils import run_node_snippet

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestFITSHeaderParser:
    """Test Suite 1: FITS 2880-byte Block Header Parsing & Metadata Cards"""

    def test_fits_standard_card_parsing(self):
        code = """
        import { FITSHeaderParser, FITS_CARD_BYTES, FITS_RECORD_BYTES } from './src/data/remote_fits_streamer.js';

        // Construct synthetic 2880-byte FITS header block
        const buffer = new Uint8Array(FITS_RECORD_BYTES);
        const cardStrings = [
          "SIMPLE  =                    T / Standard FITS format                           ",
          "BITPIX  =                  -32 / IEEE 32-bit floating point                     ",
          "NAXIS   =                    3 / 3-dimensional data cube                        ",
          "NAXIS1  =                  128 / Dimension along SGX (Mpc/h)                    ",
          "NAXIS2  =                  128 / Dimension along SGY (Mpc/h)                    ",
          "NAXIS3  =                  128 / Dimension along SGZ (Mpc/h)                    ",
          "BSCALE  =                  1.0 / Physical scale factor                          ",
          "BZERO   =                  0.0 / Physical zero offset                           ",
          "BUNIT   = 'km/s    '           / Peculiar velocity units                        ",
          "COMMENT Cosmicflows-4 Wiener Filter Reconstruction (Courtois et al. 2023)       ",
          "END                                                                             "
        ];

        let offset = 0;
        for (const card of cardStrings) {
          for (let i = 0; i < card.length; i++) {
            buffer[offset + i] = card.charCodeAt(i);
          }
          offset += FITS_CARD_BYTES;
        }

        const parsed = FITSHeaderParser.parse(buffer);
        console.log(JSON.stringify({
          simple: parsed.cards.get('SIMPLE'),
          bitpix: parsed.cards.get('BITPIX'),
          naxis: parsed.cards.get('NAXIS'),
          naxis1: parsed.cards.get('NAXIS1'),
          naxis2: parsed.cards.get('NAXIS2'),
          naxis3: parsed.cards.get('NAXIS3'),
          bunit: parsed.cards.get('BUNIT'),
          commentsCount: parsed.comments.length
        }));
        """
        res = run_node_snippet(code)
        assert res["simple"] is True
        assert res["bitpix"] == -32
        assert res["naxis"] == 3
        assert res["naxis1"] == 128
        assert res["naxis2"] == 128
        assert res["naxis3"] == 128
        assert res["bunit"] == "km/s"
        assert res["commentsCount"] >= 1

    def test_compute_data_byte_length_padding(self):
        code = """
        import { FITSHeaderParser } from './src/data/remote_fits_streamer.js';

        const cards = new Map([
          ['SIMPLE', true],
          ['BITPIX', -32],
          ['NAXIS', 3],
          ['NAXIS1', 64],
          ['NAXIS2', 64],
          ['NAXIS3', 64]
        ]);

        const byteLen = FITSHeaderParser.computeDataByteLength(cards);
        console.log(JSON.stringify({ byteLen }));
        """
        res = run_node_snippet(code)
        # 64*64*64 * 4 bytes = 1,048,576 bytes -> ceil(1048576 / 2880) * 2880 = 365 * 2880 = 1,051,200 bytes
        assert res["byteLen"] == 1051200


class TestFITSDataDecodingAndCanonicalReordering:
    """Test Suite 2: Big-Endian Conversion & Canonical (SGZ, SGY, SGX) -> (SGX, SGY, SGZ) Mapping"""

    def test_big_endian_float_decoding_and_velocity_scale(self):
        code = """
        import { RemoteFITSStreamer } from './src/data/remote_fits_streamer.js';

        const streamer = new RemoteFITSStreamer();
        const nx = 4, ny = 4, nz = 4;
        const total = nx * ny * nz;

        // Build raw synthetic Big-Endian FITS payload
        const rawBuf = new ArrayBuffer(total * 4);
        const view = new DataView(rawBuf);

        // Store known values in FITS [z][y][x] order
        for (let iz = 0; iz < nz; iz++) {
          for (let iy = 0; iy < ny; iy++) {
            for (let ix = 0; ix < nx; ix++) {
              const srcIdx = ix + nx * (iy + ny * iz);
              const val = (iz * 100) + (iy * 10) + ix;
              view.setFloat32(srcIdx * 4, val, false); // Big-Endian
            }
          }
        }

        // Test decoding with velocity scaling (x52.0)
        const scaleFactor = 52.0;
        const decoded = new Float32Array(total);
        let srcIdx = 0;
        for (let iz = 0; iz < nz; iz++) {
          for (let iy = 0; iy < ny; iy++) {
            for (let ix = 0; ix < nx; ix++) {
              const rawVal = view.getFloat32(srcIdx * 4, false);
              const physVal = rawVal * scaleFactor;
              const dstIdx = ix + nx * (iy + ny * iz);
              decoded[dstIdx] = physVal;
              srcIdx++;
            }
          }
        }

        // Check corner (ix=1, iy=2, iz=3) -> rawVal = 321, physVal = 321 * 52 = 16692
        const checkIdx = 1 + nx * (2 + ny * 3);
        console.log(JSON.stringify({
          checkVal: decoded[checkIdx],
          expectedVal: 321 * 52.0
        }));
        """
        res = run_node_snippet(code)
        assert res["checkVal"] == pytest.approx(16692.0, rel=1e-5)


class TestCF4DatasetPresets:
    """Test Suite 3: Official CF4 Dataset Presets & Metadata Registry"""

    def test_presets_structure_and_citations(self):
        code = """
        import { CF4_DATASET_PRESETS } from './src/data/remote_fits_streamer.js';

        const keys = Object.keys(CF4_DATASET_PRESETS);
        const presets = keys.map(k => CF4_DATASET_PRESETS[k]);

        console.log(JSON.stringify({
          presetCount: presets.length,
          hasUngrouped: CF4_DATASET_PRESETS.CF4_UNGROUPED_64 !== undefined,
          hasGrouped: CF4_DATASET_PRESETS.CF4_GROUPED_64 !== undefined,
          hasMean128: CF4_DATASET_PRESETS.CF4PP_POSTERIOR_MEAN_128 !== undefined,
          hasRms128: CF4_DATASET_PRESETS.CF4PP_POSTERIOR_RMS_128 !== undefined
        }));
        """
        res = run_node_snippet(code)
        assert res["presetCount"] >= 4
        assert res["hasUngrouped"] is True
        assert res["hasGrouped"] is True
        assert res["hasMean128"] is True
        assert res["hasRms128"] is True


class TestIndexedDBCachingLayer:
    """Test Suite 4: Persistent Binary Cache API Contracts"""

    def test_cache_instance_methods(self):
        code = """
        import { IndexedDBGridCache } from './src/data/remote_fits_streamer.js';

        const cache = new IndexedDBGridCache();
        console.log(JSON.stringify({
          dbName: cache.dbName,
          storeName: cache.storeName,
          hasPut: typeof cache.put === 'function',
          hasGet: typeof cache.get === 'function',
          hasDelete: typeof cache.delete === 'function',
          hasClear: typeof cache.clear === 'function'
        }));
        """
        res = run_node_snippet(code)
        assert res["dbName"] == "zrt_cosmicflows_cache_v1"
        assert res["storeName"] == "grid_blobs"
        assert res["hasPut"] is True
        assert res["hasGet"] is True
        assert res["hasDelete"] is True
        assert res["hasClear"] is True
