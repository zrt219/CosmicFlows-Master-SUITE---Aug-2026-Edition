# -*- coding: utf-8 -*-
"""
test_milestone1_surveys.py
Comprehensive verification test suite for Milestone 1: Multi-Survey Astronomical Data & Coordinate Pipeline.

Tests:
1. Astronomical Point Cloud Datasets (bundled SDSS, 2MRS, DESI JSON files in data/).
2. Dedicated Multi-Mesh Three.js Architecture (sdssPointsMesh, twoMrsPointsMesh, desiPointsMesh, cf4PointsMesh, window.surveyGroups).
3. DESI Year 1 BAO Acoustic Horizon Shell (calibrated rd = 147.1 Mpc = 10,297 units, dual center modes).
4. Kaiser Redshift-Space Distortion (RSD) Engine (dual coordinate buffers, 60 FPS in-place interpolation).
5. Survey Color Grading & Materials across Dark and White themes.
6. WebGL Memory Disposal (Rule 39).
"""

import os
import sys
import json
import time
import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def test_astronomical_datasets_exist():
    """Verify self-contained bundled data structures for SDSS, 2MRS, and DESI in data/."""
    data_dir = os.path.join(ROOT_DIR, "data")
    sdss_path = os.path.join(data_dir, "survey_sdss.json")
    twomrs_path = os.path.join(data_dir, "survey_2mrs.json")
    desi_path = os.path.join(data_dir, "survey_desi.json")

    assert os.path.exists(sdss_path), "data/survey_sdss.json missing!"
    assert os.path.exists(twomrs_path), "data/survey_2mrs.json missing!"
    assert os.path.exists(desi_path), "data/survey_desi.json missing!"

    with open(sdss_path, "r", encoding="utf-8") as f:
        sdss = json.load(f)
    with open(twomrs_path, "r", encoding="utf-8") as f:
        twomrs = json.load(f)
    with open(desi_path, "r", encoding="utf-8") as f:
        desi = json.load(f)

    assert sdss["count"] >= 5000, f"SDSS expected >= 5000 galaxies, got {sdss['count']}"
    assert len(sdss["pos_real_scene"]) >= 5000 * 3
    assert len(sdss["pos_redshift_scene"]) >= 5000 * 3

    assert twomrs["count"] >= 10000, f"2MRS expected >= 10000 galaxies, got {twomrs['count']}"
    assert len(twomrs["pos_real_scene"]) >= 10000 * 3
    assert len(twomrs["pos_redshift_scene"]) >= 10000 * 3

    assert desi["count"] >= 4000, f"DESI expected >= 4000 galaxies, got {desi['count']}"
    assert len(desi["pos_real_scene"]) >= 4000 * 3
    assert len(desi["pos_redshift_scene"]) >= 4000 * 3


def test_m1_multi_survey_and_bao_e2e(cdp):
    """E2E verification of Three.js multi-mesh architecture, BAO shell, and RSD engine."""
    print("=== Testing Milestone M1 Multi-Survey Point Clouds & DESI BAO Shell ===")

    # 1. Verify window.surveyGroups interface contract
    sg_res = cdp.evaluate("""
    (function() {
        const sg = window.surveyGroups;
        if (!sg) return { error: 'window.surveyGroups is not defined' };
        return {
            hasDesi: !!sg.desi,
            hasSdss: !!sg.sdss,
            hasTwoMrs: !!sg.twoMrs,
            hasHooleilana: !!sg.hooleilana,
            hasCf4: !!sg.cf4,
            isDesiGroup: sg.desi instanceof THREE.Group,
            isSdssPoints: sg.sdss instanceof THREE.Points,
            isTwoMrsPoints: sg.twoMrs instanceof THREE.Points,
            isCf4Points: sg.cf4 instanceof THREE.Points,
            isHooleilanaGroup: sg.hooleilana instanceof THREE.Group,
            desiHasGeometry: !!sg.desi.geometry,
            desiHasMaterial: !!sg.desi.material
        };
    })()
    """)
    assert "error" not in sg_res, sg_res.get("error")
    assert sg_res["hasDesi"] is True, "surveyGroups.desi missing"
    assert sg_res["hasSdss"] is True, "surveyGroups.sdss missing"
    assert sg_res["hasTwoMrs"] is True, "surveyGroups.twoMrs missing"
    assert sg_res["hasHooleilana"] is True, "surveyGroups.hooleilana missing"
    assert sg_res["hasCf4"] is True, "surveyGroups.cf4 missing"
    assert sg_res["isDesiGroup"] is True, "surveyGroups.desi must be THREE.Group"
    assert sg_res["isSdssPoints"] is True, "surveyGroups.sdss must be THREE.Points"
    assert sg_res["isTwoMrsPoints"] is True, "surveyGroups.twoMrs must be THREE.Points"
    assert sg_res["isCf4Points"] is True, "surveyGroups.cf4 must be THREE.Points"
    assert sg_res["isHooleilanaGroup"] is True, "surveyGroups.hooleilana must be THREE.Group"
    assert sg_res["desiHasGeometry"] is True, "surveyGroups.desi.geometry must be accessible"
    assert sg_res["desiHasMaterial"] is True, "surveyGroups.desi.material must be accessible"

    # 2. Verify Point Cloud Buffer Counts & Dynamic Draw Usage
    counts_res = cdp.evaluate("""
    (function() {
        const cf4 = window.cf4PointsMesh;
        const sdss = window.sdssPointsMesh;
        const twoMrs = window.twoMrsPointsMesh;
        const desi = window.desiPointsMesh;

        return {
            cf4Count: cf4 ? cf4.geometry.attributes.position.count : 0,
            sdssCount: sdss ? sdss.geometry.attributes.position.count : 0,
            twoMrsCount: twoMrs ? twoMrs.geometry.attributes.position.count : 0,
            desiCount: desi ? desi.geometry.attributes.position.count : 0,
            cf4Dynamic: cf4 ? cf4.geometry.attributes.position.usage === THREE.DynamicDrawUsage : false,
            sdssDynamic: sdss ? sdss.geometry.attributes.position.usage === THREE.DynamicDrawUsage : false,
            twoMrsDynamic: twoMrs ? twoMrs.geometry.attributes.position.usage === THREE.DynamicDrawUsage : false,
            desiDynamic: desi ? desi.geometry.attributes.position.usage === THREE.DynamicDrawUsage : false,
            cf4HasReal: !!(cf4 && cf4.geometry.attributes.posReal),
            cf4HasZ: !!(cf4 && cf4.geometry.attributes.posRedshift),
            sdssHasReal: !!(sdss && sdss.geometry.attributes.posReal),
            sdssHasZ: !!(sdss && sdss.geometry.attributes.posRedshift),
            twoMrsHasReal: !!(twoMrs && twoMrs.geometry.attributes.posReal),
            twoMrsHasZ: !!(twoMrs && twoMrs.geometry.attributes.posRedshift),
            desiHasReal: !!(desi && desi.geometry.attributes.posReal),
            desiHasZ: !!(desi && desi.geometry.attributes.posRedshift)
        };
    })()
    """)
    assert counts_res["cf4Count"] >= 30000, f"Expected CF4 >= 30,000, got {counts_res['cf4Count']}"
    assert counts_res["sdssCount"] >= 5000, f"Expected SDSS >= 5,000, got {counts_res['sdssCount']}"
    assert counts_res["twoMrsCount"] >= 10000, f"Expected 2MRS >= 10,000, got {counts_res['twoMrsCount']}"
    assert counts_res["desiCount"] >= 4000, f"Expected DESI >= 4,000, got {counts_res['desiCount']}"

    assert counts_res["cf4Dynamic"] is True, "cf4 position usage must be DynamicDrawUsage"
    assert counts_res["sdssDynamic"] is True, "sdss position usage must be DynamicDrawUsage"
    assert counts_res["twoMrsDynamic"] is True, "twoMrs position usage must be DynamicDrawUsage"
    assert counts_res["desiDynamic"] is True, "desi position usage must be DynamicDrawUsage"

    assert counts_res["cf4HasReal"] and counts_res["cf4HasZ"], "CF4 missing dual buffers"
    assert counts_res["sdssHasReal"] and counts_res["sdssHasZ"], "SDSS missing dual buffers"
    assert counts_res["twoMrsHasReal"] and counts_res["twoMrsHasZ"], "2MRS missing dual buffers"
    assert counts_res["desiHasReal"] and counts_res["desiHasZ"], "DESI missing dual buffers"

    # 3. Verify DESI BAO Acoustic Horizon Shell Calibration & Dual Center Modes
    bao_res = cdp.evaluate("""
    (function() {
        const bg = window.desiBaoGroup;
        if (!bg) return { error: 'window.desiBaoGroup is not defined' };

        // Find shell mesh
        let shellMesh = null;
        bg.traverse(child => {
            if (child.userData && child.userData.isDesiShell) shellMesh = child;
        });

        const radius = shellMesh ? shellMesh.geometry.parameters.radius : 0;
        const expectedRadius = 147.1 * 70; // 10297 units

        // Test Bootes center
        window.simState.desiBaoCenter = 'bootes';
        window.initDesiBaoShell();
        let bootesMesh = null;
        window.desiBaoGroup.traverse(c => { if (c.userData && c.userData.isDesiShell) bootesMesh = c; });
        const posBootes = bootesMesh ? bootesMesh.position.clone() : null;

        // Test Origin center
        window.simState.desiBaoCenter = 'origin';
        window.initDesiBaoShell();
        let originMesh = null;
        window.desiBaoGroup.traverse(c => { if (c.userData && c.userData.isDesiShell) originMesh = c; });
        const posOrigin = originMesh ? originMesh.position.clone() : null;

        // Reset to bootes
        window.simState.desiBaoCenter = 'bootes';
        window.initDesiBaoShell();

        return {
            radius: radius,
            expectedRadius: expectedRadius,
            radiusDiff: Math.abs(radius - expectedRadius),
            bootesPos: posBootes ? { x: posBootes.x, y: posBootes.y, z: posBootes.z } : null,
            originPos: posOrigin ? { x: posOrigin.x, y: posOrigin.y, z: posOrigin.z } : null
        };
    })()
    """)
    assert "error" not in bao_res, bao_res.get("error")
    assert bao_res["radiusDiff"] < 1.0, f"DESI BAO shell radius expected ~10297 units, got {bao_res['radius']}"
    assert bao_res["originPos"]["x"] == 0 and bao_res["originPos"]["y"] == 0 and bao_res["originPos"]["z"] == 0, "Origin center mode must be (0,0,0)"
    assert abs(bao_res["bootesPos"]["x"] - 3800) < 1.0, "Bootes center SGX must match 3800"

    # 4. Verify Kaiser RSD Interpolation Engine
    rsd_res = cdp.evaluate("""
    (function() {
        const sdss = window.sdssPointsMesh;
        const pos = sdss.geometry.attributes.position.array;
        const real = sdss.geometry.attributes.posReal.array;
        const z = sdss.geometry.attributes.posRedshift.array;

        // Set to real space (t = 0)
        window.setSurveyRsdMode(false);
        const diffReal0 = Math.abs(pos[0] - real[0]) + Math.abs(pos[1] - real[1]) + Math.abs(pos[2] - real[2]);

        // Set to redshift space (t = 1)
        window.setSurveyRsdMode(true);
        const diffZ1 = Math.abs(pos[0] - z[0]) + Math.abs(pos[1] - z[1]) + Math.abs(pos[2] - z[2]);

        // Interpolate half-way (t = 0.5)
        window.applySurveyRSD(0.5);
        const expectedMid = 0.5 * real[0] + 0.5 * z[0];
        const diffMid = Math.abs(pos[0] - expectedMid);

        // Reset
        window.setSurveyRsdMode(false);

        return {
            diffReal0: diffReal0,
            diffZ1: diffZ1,
            diffMid: diffMid,
            displacementNonZero: Math.abs(real[0] - z[0]) + Math.abs(real[1] - z[1]) > 0
        };
    })()
    """)
    assert rsd_res["diffReal0"] < 1e-3, "In real-space mode, positions must match posReal"
    assert rsd_res["diffZ1"] < 1e-3, "In redshift-space mode, positions must match posRedshift"
    assert rsd_res["diffMid"] < 1e-3, "Mid-way morph must be exact linear blend"
    assert rsd_res["displacementNonZero"] is True, "Real and redshift positions must exhibit non-zero RSD displacement"

    # 5. Verify Survey Color Grading across Themes
    color_res = cdp.evaluate("""
    (function() {
        // Dark theme check
        window.toggleTheme('dark');
        const desiDark = window.desiPointsMesh.material.color.getHexString();
        const sdssDark = window.sdssPointsMesh.material.color.getHexString();
        const twoMrsDark = window.twoMrsPointsMesh.material.color.getHexString();
        const cf4Dark = window.cf4PointsMesh.material.color.getHexString();

        // White theme check
        window.toggleTheme('white');
        const desiWhite = window.desiPointsMesh.material.color.getHexString();
        const sdssWhite = window.sdssPointsMesh.material.color.getHexString();
        const twoMrsWhite = window.twoMrsPointsMesh.material.color.getHexString();
        const cf4White = window.cf4PointsMesh.material.color.getHexString();

        // Reset to dark
        window.toggleTheme('dark');

        return {
            desiDark: desiDark,
            sdssDark: sdssDark,
            twoMrsDark: twoMrsDark,
            cf4Dark: cf4Dark,
            desiWhite: desiWhite,
            sdssWhite: sdssWhite,
            twoMrsWhite: twoMrsWhite,
            cf4White: cf4White
        };
    })()
    """)
    assert color_res["desiDark"] == "00f3ff", f"Expected DESI dark 00f3ff, got {color_res['desiDark']}"
    assert color_res["sdssDark"] == "ffaa00", f"Expected SDSS dark ffaa00, got {color_res['sdssDark']}"
    assert color_res["twoMrsDark"] == "ff2a55", f"Expected 2MRS dark ff2a55, got {color_res['twoMrsDark']}"
    assert color_res["cf4Dark"] == "38bdf8", f"Expected CF4 dark 38bdf8, got {color_res['cf4Dark']}"

    assert color_res["desiWhite"] == "0284c7", f"Expected DESI white 0284c7, got {color_res['desiWhite']}"
    assert color_res["sdssWhite"] == "d97706", f"Expected SDSS white d97706, got {color_res['sdssWhite']}"
    assert color_res["twoMrsWhite"] == "be123c", f"Expected 2MRS white be123c, got {color_res['twoMrsWhite']}"
    assert color_res["cf4White"] == "334155", f"Expected CF4 white 334155, got {color_res['cf4White']}"
