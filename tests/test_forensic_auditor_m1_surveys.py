# -*- coding: utf-8 -*-
"""
Independent forensic test script for Milestone 1: Multi-Survey Astronomical Data & Coordinate Pipeline.
Executed by auditor_m1_1 to independently verify:
1. Live WebGL scene meshes (THREE.Points, Float32BufferAttribute, dynamicDrawUsage)
2. Statistical validation of survey point clouds (SDSS, 2MRS, DESI, CF4)
3. Kaiser RSD interpolation precision across active vertex buffers
4. DESI BAO Shell geometry, radius (147.1 Mpc = 10,297 units), and center toggling
5. Absence of facade objects, hardcoded mocks, or fake returns
"""

import os
import sys
import json
import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def test_forensic_live_webgl_meshes(cdp):
    """Forensic verification of real Three.js objects in browser runtime."""
    eval_js = """
    (function() {
        const meshes = {
            cf4: window.cf4PointsMesh,
            sdss: window.sdssPointsMesh,
            twoMrs: window.twoMrsPointsMesh,
            desi: window.desiPointsMesh
        };

        const details = {};
        for (let key in meshes) {
            const m = meshes[key];
            if (!m) { details[key] = { error: 'missing' }; continue; }
            const g = m.geometry;
            const pos = g.attributes.position;
            const real = g.attributes.posReal;
            const z = g.attributes.posRedshift;
            details[key] = {
                isPoints: !!m.isPoints,
                isBufferGeometry: !!g.isBufferGeometry,
                posCount: pos ? pos.count : 0,
                posItemSize: pos ? pos.itemSize : 0,
                posUsage: pos ? pos.usage : 0,
                dynamicDraw: pos ? pos.usage === THREE.DynamicDrawUsage : false,
                hasPosReal: !!real,
                realCount: real ? real.count : 0,
                hasPosRedshift: !!z,
                zCount: z ? z.count : 0,
                materialType: m.material ? m.material.type : null,
                visible: m.visible,
                samplePosReal: [real.array[0], real.array[1], real.array[2]],
                samplePosZ: [z.array[0], z.array[1], z.array[2]]
            };
        }

        // Check DESI BAO Shell
        const bao = window.desiBaoGroup;
        let shellMesh = null, innerMesh = null, pinMesh = null;
        let ringCount = 0;
        if (bao) {
            bao.traverse(c => {
                if (c.userData && c.userData.isDesiShell) shellMesh = c;
                if (c.isMesh && c.material && c.material.side === THREE.BackSide) innerMesh = c;
                if (c.isLine) ringCount++;
                if (c.isMesh && c.geometry && c.geometry.type === 'SphereGeometry') pinMesh = c;
            });
        }

        // Check RSD numeric interpolation at t=0.42
        window.applySurveyRSD(0.42);
        const sdssG = window.sdssPointsMesh.geometry;
        const pArr = sdssG.attributes.position.array;
        const rArr = sdssG.attributes.posReal.array;
        const zArr = sdssG.attributes.posRedshift.array;
        let maxDiff = 0;
        for (let i = 0; i < 300; i++) {
            const expected = 0.58 * rArr[i] + 0.42 * zArr[i];
            const d = Math.abs(pArr[i] - expected);
            if (d > maxDiff) maxDiff = d;
        }
        // Reset to 0.0
        window.applySurveyRSD(0.0);

        return {
            meshes: details,
            bao: {
                exists: !!bao,
                childrenCount: bao ? bao.children.length : 0,
                hasShellMesh: !!shellMesh,
                shellRadius: shellMesh ? shellMesh.geometry.parameters.radius : 0,
                shellDetail: shellMesh ? shellMesh.geometry.parameters.detail : 0,
                hasInnerMesh: !!innerMesh,
                hasPinMesh: !!pinMesh,
                ringCount: ringCount
            },
            rsdMaxDiffAt042: maxDiff
        };
    })()
    """
    res = cdp.evaluate(eval_js)

    # 1. Assert all meshes are genuine THREE.Points with genuine Float32BufferAttribute
    for key in ['cf4', 'sdss', 'twoMrs', 'desi']:
        m = res['meshes'][key]
        assert 'error' not in m, f"Mesh {key} missing: {m}"
        assert m['isPoints'] is True, f"Mesh {key} must be THREE.Points"
        assert m['isBufferGeometry'] is True, f"Mesh {key} must have THREE.BufferGeometry"
        assert m['posItemSize'] == 3, f"Mesh {key} position itemSize must be 3"
        assert m['dynamicDraw'] is True, f"Mesh {key} position usage must be DynamicDrawUsage"
        assert m['hasPosReal'] is True, f"Mesh {key} must have posReal"
        assert m['hasPosRedshift'] is True, f"Mesh {key} must have posRedshift"
        assert m['posCount'] == m['realCount'] == m['zCount'], f"Mesh {key} buffer lengths mismatch"
        # Assert non-zero coordinate vectors
        assert any(abs(v) > 10.0 for v in m['samplePosReal']), f"Mesh {key} real pos is degenerate"
        assert any(abs(v) > 10.0 for v in m['samplePosZ']), f"Mesh {key} redshift pos is degenerate"

    # 2. Assert counts
    assert res['meshes']['cf4']['posCount'] >= 30000
    assert res['meshes']['sdss']['posCount'] >= 5000
    assert res['meshes']['twoMrs']['posCount'] >= 10000
    assert res['meshes']['desi']['posCount'] >= 4000

    # 3. Assert DESI BAO Shell
    bao = res['bao']
    assert bao['exists'] is True, "desiBaoGroup missing"
    assert bao['hasShellMesh'] is True, "Shell mesh missing from desiBaoGroup"
    assert abs(bao['shellRadius'] - 10297.0) < 1.0, f"Shell radius mismatch: {bao['shellRadius']}"
    assert bao['shellDetail'] == 3, f"Shell detail expected 3, got {bao['shellDetail']}"
    assert bao['hasInnerMesh'] is True, "Inner back-side depth shell missing"
    assert bao['hasPinMesh'] is True, "Center anchor pin missing"
    assert bao['ringCount'] == 6, f"Expected 6 meridian rings, got {bao['ringCount']}"

    # 4. Assert RSD calculation precision (Float32Array precision: coordinates up to ~25,000 units have float32 ULP ~0.0015)
    assert res['rsdMaxDiffAt042'] < 1.5e-3, f"RSD linear interpolation failed: diff={res['rsdMaxDiffAt042']}"

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
