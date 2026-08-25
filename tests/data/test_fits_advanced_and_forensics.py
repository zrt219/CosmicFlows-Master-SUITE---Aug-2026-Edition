import os
import pytest
import numpy as np
from tests.utils import run_node_snippet

class TestFitsAdvancedAndForensics:
    def test_fits_advanced_parser_synthetic_cube(self):
        code = """
        import { FitsAdvancedParser } from './src/data/fits_advanced_parser.js';
        
        // Construct a synthetic 2880-byte FITS buffer
        const buffer = new ArrayBuffer(5760);
        const u8 = new Uint8Array(buffer);
        const headerText = "SIMPLE  =                    T / file does conform to FITS standard             " +
                           "BITPIX  =                  -32 / number of bits per data pixel                  " +
                           "NAXIS   =                    3 / number of data axes                            " +
                           "NAXIS1  =                    2 / length of data axis 1                          " +
                           "NAXIS2  =                    2 / length of data axis 2                          " +
                           "NAXIS3  =                    2 / length of data axis 3                          " +
                           "END                                                                             ";
        for (let i = 0; i < headerText.length; i++) {
          u8[i] = headerText.charCodeAt(i);
        }
        
        // Write 8 float32 values in big-endian
        const view = new DataView(buffer, 2880);
        for (let i = 0; i < 8; i++) {
          view.setFloat32(i * 4, (i + 1) * 10.0, false);
        }
        
        const hdus = FitsAdvancedParser.parse(buffer);
        const data = hdus[0].readFloat32Array();
        
        console.log(JSON.stringify({
          naxis: hdus[0].naxis,
          dimensions: hdus[0].dimensions,
          data: Array.from(data)
        }));
        """
        res = run_node_snippet(code)
        assert res["naxis"] == 3
        assert res["dimensions"] == [2, 2, 2]
        assert len(res["data"]) == 8
        assert res["data"][0] == 10.0
        assert res["data"][7] == 80.0

    def test_dataset_forensics_statistical_moments(self):
        code = """
        import { DatasetForensics } from './src/data/dataset_forensics.js';
        
        const buf = new Float32Array([10, 20, 30, 40, 50, NaN, Infinity, -Infinity]);
        const report = DatasetForensics.analyze(buf, [2, 2, 2]);
        
        console.log(JSON.stringify(report));
        """
        res = run_node_snippet(code)
        assert res["validCount"] == 5
        assert res["nanCount"] == 1
        assert res["posInfCount"] == 1
        assert res["negInfCount"] == 1
        assert res["mean"] == 30.0
        assert res["min"] == 10.0
        assert res["max"] == 50.0
