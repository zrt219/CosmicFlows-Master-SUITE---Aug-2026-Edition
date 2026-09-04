#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scripts/verify_m1_adversarial.py
Empirical Adversarial Verification Suite for Milestone 1:
Multi-Survey Cosmological Data Expansion & Acoustic Shell Geometry.

Challenger: challenger_m1_2
Archetype: EMPIRICAL CHALLENGER (critic, specialist)
"""

import os
import sys
import json
import math
import numpy as np

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi

# Standard de Vaucouleurs Supergalactic to Galactic rotation matrix (RC3 / IAU)
def build_sg_to_gal_matrix():
    l_sgp = 47.37 * DEG2RAD
    b_sgp = 6.32 * DEG2RAD
    l_sg0 = 137.37 * DEG2RAD

    z_x = math.cos(b_sgp) * math.cos(l_sgp)
    z_y = math.cos(b_sgp) * math.sin(l_sgp)
    z_z = math.sin(b_sgp)

    x_x = math.cos(l_sg0)
    x_y = math.sin(l_sg0)
    x_z = 0.0

    y_x = z_y * x_z - z_z * x_y
    y_y = z_z * x_x - z_x * x_z
    y_z = z_x * x_y - z_y * x_x
    y_norm = math.hypot(y_x, y_y, y_z)

    R = np.array([
        [x_x, y_x / y_norm, z_x],
        [x_y, y_y / y_norm, z_y],
        [x_z, y_z / y_norm, z_z]
    ])
    return R

ROT_SG_TO_GAL = build_sg_to_gal_matrix()

def sg_to_gal(sgx, sgy, sgz):
    v_sg = np.array([sgx, sgy, sgz])
    v_gal = ROT_SG_TO_GAL @ v_sg
    r = np.linalg.norm(v_gal)
    if r < 1e-9:
        return 0.0, 0.0, 0.0
    b = math.asin(np.clip(v_gal[2] / r, -1.0, 1.0)) * RAD2DEG
    l = math.atan2(v_gal[1], v_gal[0]) * RAD2DEG
    if l < 0:
        l += 360.0
    return l, b, r


class EmpiricalM1Challenge:
    def __init__(self):
        self.results = {}
        self.failures = []

    def log_result(self, category, check_name, passed, detail=""):
        if category not in self.results:
            self.results[category] = []
        self.results[category].append({
            "check": check_name,
            "passed": passed,
            "detail": detail
        })
        status = "[PASS]" if passed else "[FAIL]"
        print(f"  {status} {check_name}: {detail}")
        if not passed:
            self.failures.append(f"{category} -> {check_name}: {detail}")

    # =========================================================================
    # 1. SURVEY JSON DATASET AUDIT
    # =========================================================================
    def audit_survey_datasets(self):
        print("\n=== [PART 1] EMPIRICAL ASTROMETRIC DATASET AUDIT ===")
        data_dir = os.path.join(ROOT_DIR, "data")
        files = {
            "SDSS": (os.path.join(data_dir, "survey_sdss.json"), 4000),
            "2MRS": (os.path.join(data_dir, "survey_2mrs.json"), 8000),
            "DESI": (os.path.join(data_dir, "survey_desi.json"), 3000)
        }

        loaded_data = {}
        for name, (path, min_count) in files.items():
            if not os.path.exists(path):
                self.log_result("Dataset Existence", f"{name} file present", False, f"{path} not found")
                continue
            
            size_kb = os.path.getsize(path) / 1024.0
            try:
                with open(path, "r", encoding="utf-8") as f:
                    d = json.load(f)
                loaded_data[name] = d
                self.log_result("Dataset Parsing", f"{name} JSON parse", True, f"{size_kb:.1f} KB, parsed successfully")
            except Exception as e:
                self.log_result("Dataset Parsing", f"{name} JSON parse", False, str(e))
                continue

            # Vertex count assertions
            count = d.get("count", 0)
            passed_count = count > min_count
            self.log_result(
                "Vertex Count",
                f"{name} vertex count requirement (>{min_count})",
                passed_count,
                f"count={count} (required > {min_count})"
            )

            # Array length and schema alignment
            pos_real_mpch = d.get("pos_real_mpch", [])
            pos_real_scene = d.get("pos_real_scene", [])
            pos_redshift_scene = d.get("pos_redshift_scene", [])
            v_los = d.get("v_los", [])
            v_thermal = d.get("v_thermal", [])
            target_class = d.get("target_class", [])
            clusters = d.get("cluster_association", [])

            align_passed = (
                len(pos_real_mpch) == count * 3 and
                len(pos_real_scene) == count * 3 and
                len(pos_redshift_scene) == count * 3 and
                len(v_los) == count and
                len(v_thermal) == count and
                len(target_class) == count and
                len(clusters) == count
            )
            self.log_result(
                "Array Alignment",
                f"{name} 1:1 buffer alignment across all 7 channels",
                align_passed,
                f"real_mpch={len(pos_real_mpch)}, real_scene={len(pos_real_scene)}, z_scene={len(pos_redshift_scene)}, v_los={len(v_los)}"
            )

            # Check for NaN / Inf
            all_floats = pos_real_mpch + pos_real_scene + pos_redshift_scene + v_los + v_thermal
            has_nan_inf = any(math.isnan(x) or math.isinf(x) for x in all_floats)
            self.log_result("Numerical Integrity", f"{name} zero NaN/Infinity", not has_nan_inf, f"Checked {len(all_floats)} floats")

            # Physical scale factor: 1 Mpc/h = 74.6 scene units (H0 = 74.6)
            h0 = d.get("h0", 74.6)
            scale_matches = True
            for i in range(min(count, 200)):
                mpch_x, mpch_y, mpch_z = pos_real_mpch[i*3 : i*3+3]
                # Three.js canonical mapping: sg_to_three(sgx, sgy, sgz) = [sgx, sgz, -sgy]
                expected_three_x = round(mpch_x * 74.6, 1)
                expected_three_y = round(mpch_z * 74.6, 1)
                expected_three_z = round(-mpch_y * 74.6, 1)
                actual_x, actual_y, actual_z = pos_real_scene[i*3 : i*3+3]
                if abs(actual_x - expected_three_x) > 0.5 or abs(actual_y - expected_three_y) > 0.5 or abs(actual_z - expected_three_z) > 0.5:
                    scale_matches = False
                    break
            self.log_result(
                "Physical Transformation",
                f"{name} comoving Mpc/h -> Three.js (SGX, SGZ, -SGY) * 74.6 scale",
                scale_matches,
                f"Checked first 200 galaxies against canonical axis mapping"
            )

        # =====================================================================
        # 2MRS ASTROMETRY: 8 OCTANTS & ZONE OF AVOIDANCE (|b| < 15 deg)
        # =====================================================================
        if "2MRS" in loaded_data:
            d2 = loaded_data["2MRS"]
            mpch = d2["pos_real_mpch"]
            scene = d2["pos_real_scene"]
            n = d2["count"]

            # Octants in Supergalactic coordinates
            octants_sg = set()
            octant_counts_sg = {}
            # Octants in Scene coordinates
            octants_scene = set()

            zoa_galaxies_15 = 0
            zoa_galaxies_10 = 0
            zoa_galaxies_5 = 0
            b_list = []

            for i in range(n):
                gx = mpch[i*3]
                gy = mpch[i*3+1]
                gz = mpch[i*3+2]
                code_sg = f"{'+' if gx>=0 else '-'}{'+' if gy>=0 else '-'}{'+' if gz>=0 else '-'}"
                octants_sg.add(code_sg)
                octant_counts_sg[code_sg] = octant_counts_sg.get(code_sg, 0) + 1

                sx = scene[i*3]
                sy = scene[i*3+1]
                sz = scene[i*3+2]
                code_scene = f"{'+' if sx>=0 else '-'}{'+' if sy>=0 else '-'}{'+' if sz>=0 else '-'}"
                octants_scene.add(code_scene)

                l, b, r = sg_to_gal(gx, gy, gz)
                b_list.append(b)
                if abs(b) < 15.0:
                    zoa_galaxies_15 += 1
                if abs(b) < 10.0:
                    zoa_galaxies_10 += 1
                if abs(b) < 5.0:
                    zoa_galaxies_5 += 1

            self.log_result(
                "2MRS All-Sky Coverage",
                "2MRS penetrates all 8 Supergalactic octants",
                len(octants_sg) == 8,
                f"Octants populated: {len(octants_sg)}/8. Counts: {octant_counts_sg}"
            )
            self.log_result(
                "2MRS All-Sky Coverage",
                "2MRS penetrates all 8 Three.js scene octants",
                len(octants_scene) == 8,
                f"Scene octants populated: {len(octants_scene)}/8"
            )
            self.log_result(
                "2MRS Zone of Avoidance",
                "2MRS penetrates ZoA (|b| < 15 deg)",
                zoa_galaxies_15 > 500,
                f"{zoa_galaxies_15} galaxies ({zoa_galaxies_15/n*100:.1f}%) within |b| < 15 deg"
            )
            self.log_result(
                "2MRS Zone of Avoidance",
                "2MRS deep core penetration (|b| < 10 deg and |b| < 5 deg)",
                zoa_galaxies_10 > 200 and zoa_galaxies_5 > 50,
                f"|b| < 10 deg: {zoa_galaxies_10} galaxies, |b| < 5 deg: {zoa_galaxies_5} galaxies (min |b| = {min(abs(x) for x in b_list):.3f} deg)"
            )

            # Check ZoA cluster presence
            clusters_set = set(d2.get("cluster_association", []))
            has_puppi = "Puppis Cluster (ZoA)" in clusters_set
            has_ophi = "Ophiuchus Cluster (ZoA)" in clusters_set
            self.log_result(
                "2MRS ZoA Clusters",
                "ZoA hidden benchmark clusters (Puppis & Ophiuchus) present",
                has_puppi and has_ophi,
                f"Puppis present: {has_puppi}, Ophiuchus present: {has_ophi}"
            )

        # =====================================================================
        # SDSS ASTROMETRY: NORTHERN GALACTIC CAP & SLOAN GREAT WALL SPINE
        # =====================================================================
        if "SDSS" in loaded_data:
            ds = loaded_data["SDSS"]
            mpch = ds["pos_real_mpch"]
            clusters = ds["cluster_association"]
            n = ds["count"]

            sgw_spine_count = sum(1 for c in clusters if c == "SGW_Spine")
            ngc_filament_count = sum(1 for c in clusters if c == "NGC_Filament")
            cluster_gal_count = n - sgw_spine_count - ngc_filament_count

            # Galactic latitude distribution: SDSS NGC is primarily b > 20 deg
            b_list = []
            for i in range(n):
                gx, gy, gz = mpch[i*3 : i*3+3]
                l, b, r = sg_to_gal(gx, gy, gz)
                b_list.append(b)

            b_arr = np.array(b_list)
            ngc_fraction = np.sum(b_arr > 0.0) / float(n)

            # Spine centroid
            spine_pts = [mpch[i*3:i*3+3] for i in range(n) if clusters[i] == "SGW_Spine"]
            spine_mean_mpch = np.mean(spine_pts, axis=0) if spine_pts else [0,0,0]
            spine_mean_kms = spine_mean_mpch * 74.6

            self.log_result(
                "SDSS Astrometry",
                "SDSS points reside predominantly in Northern Galactic Cap (b > 0 deg)",
                ngc_fraction > 0.85,
                f"NGC fraction: {ngc_fraction*100:.1f}% (mean b = {np.mean(b_arr):.1f} deg, min b = {np.min(b_arr):.1f} deg)"
            )
            self.log_result(
                "SDSS Astrometry",
                "Sloan Great Wall spine centroid matches published astrometry",
                abs(spine_mean_kms[0] - 12250) < 1500 and abs(spine_mean_kms[1] - 14000) < 1500 and abs(spine_mean_kms[2] - 11000) < 1500,
                f"Spine mean: ({spine_mean_kms[0]:.0f}, {spine_mean_kms[1]:.0f}, {spine_mean_kms[2]:.0f}) km/s"
            )
            self.log_result(
                "SDSS Morphology",
                "SDSS spine and rich clusters partition (Abell 2142, Abell 2199, Hercules, Coma)",
                sgw_spine_count >= 2500 and cluster_gal_count >= 800,
                f"Spine galaxies: {sgw_spine_count}, NGC sheets: {ngc_filament_count}, Clusters: {cluster_gal_count}"
            )

        # =====================================================================
        # DESI ASTROMETRY: BGS/LRG FIBERS & BAO SOUND HORIZON SHELL
        # =====================================================================
        if "DESI" in loaded_data:
            dd = loaded_data["DESI"]
            mpch = dd["pos_real_mpch"]
            t_classes = dd["target_class"]
            clusters = dd["cluster_association"]
            n = dd["count"]

            has_lrg = "LRG" in t_classes
            has_bgs = "BGS_BRIGHT" in t_classes
            lrg_count = sum(1 for c in t_classes if c == "LRG")
            bgs_count = sum(1 for c in t_classes if c == "BGS_BRIGHT")
            elg_count = sum(1 for c in t_classes if c == "ELG")

            # BAO shell radial calibration from Boötes: (3800, 16200, 8900) km/s
            bootes_kms = np.array([3800.0, 16200.0, 8900.0])
            shell_radii_kms = []
            for i in range(n):
                if clusters[i] == "DESI_BAO_Shell":
                    p_kms = np.array(mpch[i*3:i*3+3]) * 74.6
                    r_b = np.linalg.norm(p_kms - bootes_kms)
                    shell_radii_kms.append(r_b)

            shell_radii_kms = np.array(shell_radii_kms)
            mean_rd_kms = np.mean(shell_radii_kms) if len(shell_radii_kms) > 0 else 0
            mean_rd_mpc = mean_rd_kms / 74.6
            std_rd_kms = np.std(shell_radii_kms) if len(shell_radii_kms) > 0 else 0

            self.log_result(
                "DESI Target Classes",
                "DESI fiber targets include LRG and BGS_BRIGHT spectroscopic classes",
                has_lrg and has_bgs,
                f"LRG={lrg_count}, BGS_BRIGHT={bgs_count}, ELG={elg_count} (total={n})"
            )
            self.log_result(
                "DESI BAO Shell Astrometry",
                "DESI BAO shell points calibrated to r_d = 147.1 Mpc from Boötes core",
                abs(mean_rd_mpc - 147.1) < 5.0 and len(shell_radii_kms) >= 2000,
                f"Shell points: {len(shell_radii_kms)}, mean radius: {mean_rd_mpc:.2f} Mpc ({mean_rd_kms:.1f} km/s), std: {std_rd_kms:.1f} km/s"
            )

    # =========================================================================
    # 2. IN-BROWSER HEADLESS CDP VERIFICATION (THREE.JS RUNTIME)
    # =========================================================================
    def run_cdp_runtime_checks(self):
        print("\n=== [PART 2] IN-BROWSER HEADLESS CDP RUNTIME CHECKS ===")
        from tests.cdp_client import ChromeCDPClient
        from tests.conftest import ensure_http_server

        port = 9222
        try:
            ensure_http_server(8000)
            cdp = ChromeCDPClient(port=port, url="http://localhost:8000", spawn_headless=True)
            cdp.start()
            cdp.wait_for_condition("document.readyState === 'complete' && !!document.querySelector('#scene canvas')", timeout=15.0)
            print(f"Connected to Chrome CDP on port {port}")
        except Exception as e:
            self.log_result("CDP Connection", "Connect to Chrome CDP", False, str(e))
            return

        # 1. DESI BAO Shell Geometry, Dual Centers & Meridian Wave Rings
        bao_eval = cdp.evaluate("""
        (function() {
            try {
                const bg = window.desiBaoGroup;
                if (!bg) return { error: "window.desiBaoGroup is undefined" };

                // 1. Traverse and identify components
                let shellMesh = null;
                let innerMesh = null;
                let rings = [];
                let pin = null;

                bg.traverse(c => {
                    if (c.userData && c.userData.isDesiShell) shellMesh = c;
                    else if (c.isMesh && c.geometry && c.geometry.type === 'IcosahedronGeometry') innerMesh = c;
                    else if (c.isLine) rings.push(c);
                    else if (c.isMesh && c.geometry && c.geometry.type === 'SphereGeometry') pin = c;
                });

                if (!shellMesh) return { error: "DESI outer shell mesh not found" };

                const outerRadius = shellMesh.geometry.parameters ? shellMesh.geometry.parameters.radius : null;
                const innerRadius = (innerMesh && innerMesh.geometry.parameters) ? innerMesh.geometry.parameters.radius : null;

                // Meridian ring verification
                const ringCount = rings.length;
                let ringPointsValid = true;
                let ringRadii = [];
                rings.forEach(r => {
                    const pos = r.geometry.attributes.position;
                    if (!pos || pos.count < 60) ringPointsValid = false;
                    // Compute radius of ring from first vertex
                    const x = pos.getX(0), y = pos.getY(0), z = pos.getZ(0);
                    ringRadii.push(Math.sqrt(x*x + y*y + z*z));
                });

                // 2. Dual Center Modes Verification: 'bootes' vs 'origin'
                // Test Boötes
                window.simState.desiBaoCenter = 'bootes';
                window.initDesiBaoShell();
                let bShell = null;
                window.desiBaoGroup.traverse(c => { if (c.userData && c.userData.isDesiShell) bShell = c; });
                const bPos = bShell ? { x: bShell.position.x, y: bShell.position.y, z: bShell.position.z } : null;

                // Test Origin
                window.simState.desiBaoCenter = 'origin';
                window.initDesiBaoShell();
                let oShell = null;
                window.desiBaoGroup.traverse(c => { if (c.userData && c.userData.isDesiShell) oShell = c; });
                const oPos = oShell ? { x: oShell.position.x, y: oShell.position.y, z: oShell.position.z } : null;

                // Reset to bootes
                window.simState.desiBaoCenter = 'bootes';
                window.initDesiBaoShell();

                return {
                    outerRadius: outerRadius,
                    innerRadius: innerRadius,
                    ringCount: ringCount,
                    ringPointsValid: ringPointsValid,
                    ringMeanRadius: ringRadii.length > 0 ? (ringRadii.reduce((a,b)=>a+b,0)/ringRadii.length) : 0,
                    hasPin: !!pin,
                    bootesPos: bPos,
                    originPos: oPos
                };
            } catch(e) {
                return { error: e.message };
            }
        })()
        """)

        if "error" in bao_eval:
            self.log_result("DESI BAO Shell", "Evaluate BAO shell in Three.js scene", False, bao_eval["error"])
        else:
            outer_r = bao_eval.get("outerRadius")
            expected_r = 147.1 * 70.0  # 10,297.0 units
            r_diff = abs(outer_r - expected_r) if outer_r is not None else 999
            self.log_result(
                "DESI BAO Shell Radius",
                "DESI BAO shell radius exactly 10,297 ± 5 units",
                r_diff <= 5.0,
                f"outerRadius={outer_r} (expected {expected_r}, diff={r_diff:.2f} units)"
            )

            inner_r = bao_eval.get("innerRadius")
            expected_inner = expected_r * 0.985
            inner_diff = abs(inner_r - expected_inner) if inner_r is not None else 999
            self.log_result(
                "DESI BAO Inner Shell",
                "DESI inner depth shell radius matches 0.985 * r_d",
                inner_diff <= 5.0,
                f"innerRadius={inner_r} (expected {expected_inner:.1f})"
            )

            ring_count = bao_eval.get("ringCount", 0)
            ring_valid = bao_eval.get("ringPointsValid", False)
            ring_r = bao_eval.get("ringMeanRadius", 0)
            self.log_result(
                "DESI BAO Meridian Rings",
                "6 longitudinal meridian wave rings with 64-segment resolution at r = 10,297",
                ring_count == 6 and ring_valid and abs(ring_r - expected_r) <= 5.0,
                f"ringCount={ring_count}, ringValid={ring_valid}, meanRadius={ring_r:.1f}"
            )

            b_pos = bao_eval.get("bootesPos") or {}
            # In Three.js coordinates: sgToThree(3800, 16200, 8900) = (3800, 8900, -16200)
            bootes_correct = (
                abs(b_pos.get("x", 0) - 3800) < 1.0 and
                abs(b_pos.get("y", 0) - 8900) < 1.0 and
                abs(b_pos.get("z", 0) - (-16200)) < 1.0
            )
            self.log_result(
                "DESI BAO Dual Centers",
                "Boötes center mode maps to sgToThree(3800, 16200, 8900)",
                bootes_correct,
                f"bootesPos={b_pos}"
            )

            o_pos = bao_eval.get("originPos") or {}
            origin_correct = (
                abs(o_pos.get("x", 0)) < 1e-4 and
                abs(o_pos.get("y", 0)) < 1e-4 and
                abs(o_pos.get("z", 0)) < 1e-4
            )
            self.log_result(
                "DESI BAO Dual Centers",
                "Origin center mode maps to (0, 0, 0)",
                origin_correct,
                f"originPos={o_pos}"
            )

        # 2. WebGL Resource Disposal (`disposeHierarchy`)
        disposal_eval = cdp.evaluate("""
        (function() {
            try {
                if (typeof window.disposeHierarchy !== 'function') {
                    return { error: "window.disposeHierarchy is not exposed" };
                }

                // Create synthetic nested hierarchy with Mesh, Points, Lines, Groups, Textures, Materials
                const testGroup = new THREE.Group();
                const geo1 = new THREE.BufferGeometry();
                geo1.setAttribute('position', new THREE.Float32BufferAttribute([0,0,0, 1,1,1, 2,2,2], 3));
                const mat1 = new THREE.PointsMaterial({ size: 10, color: 0xff0000 });
                const pts = new THREE.Points(geo1, mat1);
                testGroup.add(pts);

                const geo2 = new THREE.IcosahedronGeometry(100, 2);
                const mat2 = new THREE.MeshBasicMaterial({ color: 0x00ff00 });
                const mesh = new THREE.Mesh(geo2, mat2);
                testGroup.add(mesh);

                let geo1Disposed = false;
                let mat1Disposed = false;
                let geo2Disposed = false;
                let mat2Disposed = false;

                geo1.addEventListener('dispose', () => { geo1Disposed = true; });
                mat1.addEventListener('dispose', () => { mat1Disposed = true; });
                geo2.addEventListener('dispose', () => { geo2Disposed = true; });
                mat2.addEventListener('dispose', () => { mat2Disposed = true; });

                // Execute disposal
                window.disposeHierarchy(testGroup);

                const childrenCleared = testGroup.children.length === 0;

                // Also test disposing desiBaoGroup safely and re-initializing it
                let desiBaoDisposalSuccess = false;
                try {
                    const bg = window.desiBaoGroup;
                    window.disposeHierarchy(bg);
                    desiBaoDisposalSuccess = true;
                    window.initDesiBaoShell();
                } catch(e) {
                    desiBaoDisposalSuccess = false;
                }

                return {
                    geo1Disposed,
                    mat1Disposed,
                    geo2Disposed,
                    mat2Disposed,
                    childrenCleared,
                    desiBaoDisposalSuccess
                };
            } catch(e) {
                return { error: e.message };
            }
        })()
        """)

        if "error" in disposal_eval:
            self.log_result("WebGL Disposal", "disposeHierarchy execution", False, disposal_eval["error"])
        else:
            all_disposed = (
                disposal_eval.get("geo1Disposed") and
                disposal_eval.get("mat1Disposed") and
                disposal_eval.get("geo2Disposed") and
                disposal_eval.get("mat2Disposed") and
                disposal_eval.get("childrenCleared") and
                disposal_eval.get("desiBaoDisposalSuccess")
            )
            self.log_result(
                "WebGL Resource Disposal",
                "disposeHierarchy() recursively cleans geometries, materials, and children with 0 exceptions",
                all_disposed,
                f"geoDisposed={disposal_eval.get('geo1Disposed')}, matDisposed={disposal_eval.get('mat1Disposed')}, childrenCleared={disposal_eval.get('childrenCleared')}, liveShellDisposal={disposal_eval.get('desiBaoDisposalSuccess')}"
            )

        # 3. Dynamic Dual Buffer RSD & Performance Check
        rsd_eval = cdp.evaluate("""
        (function() {
            try {
                const meshes = [
                    { name: 'SDSS', mesh: window.sdssPointsMesh },
                    { name: '2MRS', mesh: window.twoMrsPointsMesh },
                    { name: 'DESI', mesh: window.desiPointsMesh },
                    { name: 'CF4', mesh: window.cf4PointsMesh }
                ];

                let allBuffersPresent = true;
                let details = {};

                meshes.forEach(m => {
                    if (!m.mesh || !m.mesh.geometry) {
                        allBuffersPresent = false;
                        details[m.name] = "missing";
                        return;
                    }
                    const attrs = m.mesh.geometry.attributes;
                    const hasPos = !!attrs.position;
                    const hasReal = !!attrs.posReal;
                    const hasZ = !!attrs.posRedshift;
                    const isDynamic = hasPos && attrs.position.usage === THREE.DynamicDrawUsage;
                    if (!hasPos || !hasReal || !hasZ || !isDynamic) {
                        allBuffersPresent = false;
                    }
                    details[m.name] = {
                        count: attrs.position ? attrs.position.count : 0,
                        hasReal, hasZ, isDynamic
                    };
                });

                // Benchmark 100 RSD morph iterations for 60 FPS performance (< 16.6ms / frame)
                const t0 = performance.now();
                for (let i = 0; i < 60; i++) {
                    window.applySurveyRSD(i / 60.0);
                }
                const t1 = performance.now();
                const avgFrameTimeMs = (t1 - t0) / 60.0;

                // Reset to real-space
                window.setSurveyRsdMode(false);

                return {
                    allBuffersPresent,
                    details,
                    avgFrameTimeMs
                };
            } catch(e) {
                return { error: e.message };
            }
        })()
        """)

        if "error" in rsd_eval:
            self.log_result("Kaiser RSD Engine", "Evaluate RSD dual buffers and morphing", False, rsd_eval["error"])
        else:
            self.log_result(
                "Kaiser RSD Dual Buffers",
                "All 4 survey meshes possess posReal, posRedshift, and DynamicDrawUsage",
                rsd_eval.get("allBuffersPresent", False),
                f"Buffers: {rsd_eval.get('details')}"
            )
            avg_t = rsd_eval.get("avgFrameTimeMs", 999)
            self.log_result(
                "Kaiser RSD 60 FPS Budget",
                "applySurveyRSD in-place morphing executes well within 16.6 ms frame budget",
                avg_t < 5.0,
                f"avgFrameTime = {avg_t:.3f} ms / frame (budget < 16.6 ms)"
            )

        try:
            cdp.close()
        except Exception:
            pass

    def summary(self):
        print("\n=================================================================")
        print("                 CHALLENGER VERIFICATION SUMMARY                ")
        print("=================================================================")
        total_checks = sum(len(v) for v in self.results.values())
        total_passed = sum(sum(1 for c in v if c["passed"]) for v in self.results.values())
        total_failed = len(self.failures)

        print(f"Total Checks: {total_checks}")
        print(f"Passed:       {total_passed}")
        print(f"Failed:       {total_failed}")

        if total_failed == 0:
            print("\nVERDICT: >>> APPROVE <<<")
            return "APPROVE"
        else:
            print(f"\nVERDICT: >>> FAIL <<< ({total_failed} failures)")
            for f in self.failures:
                print(f"  - {f}")
            return "FAIL"


if __name__ == "__main__":
    challenger = EmpiricalM1Challenge()
    challenger.audit_survey_datasets()
    challenger.run_cdp_runtime_checks()
    verdict = challenger.summary()
    sys.exit(0 if verdict == "APPROVE" else 1)
