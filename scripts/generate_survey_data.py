#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_survey_data.py
Generates self-contained bundled astronomical survey datasets:
1. data/survey_sdss.json: SDSS Great Wall & Northern Galactic Cap (~5,000 galaxies)
2. data/survey_2mrs.json: 2MASS Redshift Survey all-sky (~10,000 galaxies)

Coordinate rules:
- Comoving positions: Mpc/h
- Peculiar velocities: km/s
- H0 = 74.6 km/s / (Mpc/h)
- 1 Mpc/h = 74.6 scene units (1 Mpc = 70 scene units)
- Canonical Supergalactic (SGX, SGY, SGZ) -> Three.js (SGX, SGZ, -SGY)
- Redshift-space distortion: s = x_real + (v_los + v_thermal) * r_hat
"""

import os
import json
import math
import random

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

H0 = 74.6  # km/s / (Mpc/h)
MPC_H_TO_SCENE = 74.6  # 1 Mpc/h = 74.6 scene units

def sg_to_three(sgx, sgy, sgz):
    """Canonical Supergalactic to Three.js coordinates."""
    return [sgx, sgz, -sgy]

def generate_sdss_data():
    """
    Generates SDSS Great Wall & Northern Galactic Cap spectroscopic sample (~5,000 galaxies).
    Sloan Great Wall spine:
    SGX ~ 12,000 km/s (160.8 Mpc/h), SGY ~ 14,000 km/s (187.7 Mpc/h), SGZ ~ 11,000 km/s (147.5 Mpc/h).
    Extended filament connecting to Coma, Hercules, A2142, A2199.
    """
    random.seed(42)
    galaxies_real_mpch = []
    galaxies_real_scene = []
    galaxies_redshift_scene = []
    v_los_list = []
    v_thermal_list = []
    target_classes = []
    clusters = []

    sdss_clusters = [
        {"name": "Abell 2142", "pos_kms": [1800, 13800, 6100], "count": 280, "sigma_v": 980, "cls": "LRG"},
        {"name": "Abell 2199", "pos_kms": [2800, 9200, 4100], "count": 220, "sigma_v": 820, "cls": "LRG"},
        {"name": "Hercules Complex", "pos_kms": [3200, 11000, 4500], "count": 350, "sigma_v": 760, "cls": "BGS_BRIGHT"},
        {"name": "Coma Extension", "pos_kms": [600, 7500, 1800], "count": 150, "sigma_v": 850, "cls": "BGS_BRIGHT"},
    ]

    total_count = 5000
    cluster_galaxy_count = sum(c["count"] for c in sdss_clusters)  # 1000
    spine_count = 2800
    sheet_count = total_count - cluster_galaxy_count - spine_count  # 1200

    # 1. Sloan Great Wall Spine (~2,800 galaxies)
    p_start = [8000, 10500, 7500]
    p_end = [16500, 17500, 14500]
    for i in range(spine_count):
        t = random.random()
        curve_amp = 1200 * math.sin(t * math.pi)
        base_x = p_start[0] + t * (p_end[0] - p_start[0]) + curve_amp * 0.5
        base_y = p_start[1] + t * (p_end[1] - p_start[1]) - curve_amp * 0.3
        base_z = p_start[2] + t * (p_end[2] - p_start[2]) + curve_amp * 0.2

        rad = random.gauss(0, 1100)
        phi = random.random() * 2 * math.pi
        gx_kms = base_x + rad * math.cos(phi)
        gy_kms = base_y + rad * math.sin(phi)
        gz_kms = base_z + random.gauss(0, 700)

        # Infall velocity towards spine centroid
        dx_c = 12250 - gx_kms
        dy_c = 14000 - gy_kms
        dz_c = 11000 - gz_kms
        dist_c = math.sqrt(dx_c*dx_c + dy_c*dy_c + dz_c*dz_c) + 1e-7
        infall_speed = 350 * min(dist_c / 4000, 1.0)
        vx = infall_speed * (dx_c / dist_c)
        vy = infall_speed * (dy_c / dist_c)
        vz = infall_speed * (dz_c / dist_c)

        sigma_v = 187.0
        v_thermal = random.gauss(0, sigma_v)

        gx_mpch = gx_kms / H0
        gy_mpch = gy_kms / H0
        gz_mpch = gz_kms / H0

        r_kms = math.sqrt(gx_kms*gx_kms + gy_kms*gy_kms + gz_kms*gz_kms)
        if r_kms < 1e-7:
            r_hat = [0.0, 1.0, 0.0]
        else:
            r_hat = [gx_kms / r_kms, gy_kms / r_kms, gz_kms / r_kms]

        v_los = vx * r_hat[0] + vy * r_hat[1] + vz * r_hat[2]

        sx_real = gx_mpch * MPC_H_TO_SCENE
        sy_real = gy_mpch * MPC_H_TO_SCENE
        sz_real = gz_mpch * MPC_H_TO_SCENE
        three_real = sg_to_three(sx_real, sy_real, sz_real)

        total_delta_v = v_los + v_thermal
        sx_z = sx_real + total_delta_v * r_hat[0]
        sy_z = sy_real + total_delta_v * r_hat[1]
        sz_z = sz_real + total_delta_v * r_hat[2]
        three_z = sg_to_three(sx_z, sy_z, sz_z)

        galaxies_real_mpch.extend([round(gx_mpch, 2), round(gy_mpch, 2), round(gz_mpch, 2)])
        galaxies_real_scene.extend([round(three_real[0], 1), round(three_real[1], 1), round(three_real[2], 1)])
        galaxies_redshift_scene.extend([round(three_z[0], 1), round(three_z[1], 1), round(three_z[2], 1)])
        v_los_list.append(round(v_los, 1))
        v_thermal_list.append(round(v_thermal, 1))
        target_classes.append("LRG" if random.random() < 0.65 else "BGS_BRIGHT")
        clusters.append("SGW_Spine")

    # 2. NGC Connecting Filaments & Sheets (~1,200 galaxies)
    for i in range(sheet_count):
        u = random.random()
        gx_kms = 4000 + u * 10000 + random.gauss(0, 1200)
        gy_kms = 6000 + random.random() * 10000
        gz_kms = 3000 + 0.6 * gy_kms + random.gauss(0, 850)

        vx = -120 + random.gauss(0, 150)
        vy = -200 + random.gauss(0, 150)
        vz = -80 + random.gauss(0, 150)
        sigma_v = 187.0
        v_thermal = random.gauss(0, sigma_v)

        gx_mpch = gx_kms / H0
        gy_mpch = gy_kms / H0
        gz_mpch = gz_kms / H0

        r_kms = math.sqrt(gx_kms*gx_kms + gy_kms*gy_kms + gz_kms*gz_kms)
        if r_kms < 1e-7:
            r_hat = [0.0, 1.0, 0.0]
        else:
            r_hat = [gx_kms / r_kms, gy_kms / r_kms, gz_kms / r_kms]

        v_los = vx * r_hat[0] + vy * r_hat[1] + vz * r_hat[2]

        sx_real = gx_mpch * MPC_H_TO_SCENE
        sy_real = gy_mpch * MPC_H_TO_SCENE
        sz_real = gz_mpch * MPC_H_TO_SCENE
        three_real = sg_to_three(sx_real, sy_real, sz_real)

        total_delta_v = v_los + v_thermal
        sx_z = sx_real + total_delta_v * r_hat[0]
        sy_z = sy_real + total_delta_v * r_hat[1]
        sz_z = sz_real + total_delta_v * r_hat[2]
        three_z = sg_to_three(sx_z, sy_z, sz_z)

        galaxies_real_mpch.extend([round(gx_mpch, 2), round(gy_mpch, 2), round(gz_mpch, 2)])
        galaxies_real_scene.extend([round(three_real[0], 1), round(three_real[1], 1), round(three_real[2], 1)])
        galaxies_redshift_scene.extend([round(three_z[0], 1), round(three_z[1], 1), round(three_z[2], 1)])
        v_los_list.append(round(v_los, 1))
        v_thermal_list.append(round(v_thermal, 1))
        target_classes.append("BGS_BRIGHT" if random.random() < 0.7 else "ELG")
        clusters.append("NGC_Filament")

    # 3. Dense Rich Clusters (~1,000 galaxies)
    for cl in sdss_clusters:
        c_pos = cl["pos_kms"]
        sigma_v = cl["sigma_v"]
        for i in range(cl["count"]):
            u1, u2, u3 = random.random(), random.random(), random.random()
            rad = 1200 * math.sqrt(u1 / (1.0 - u1 + 0.05))
            rad = min(rad, 2800)
            th = u2 * math.pi * 2
            ph = math.acos(2 * u3 - 1)
            gx_kms = c_pos[0] + rad * math.sin(ph) * math.cos(th)
            gy_kms = c_pos[1] + rad * math.cos(ph)
            gz_kms = c_pos[2] + rad * math.sin(ph) * math.sin(th)

            v_thermal = random.gauss(0, sigma_v)
            vx = random.gauss(0, sigma_v * 0.5)
            vy = random.gauss(0, sigma_v * 0.5)
            vz = random.gauss(0, sigma_v * 0.5)

            gx_mpch = gx_kms / H0
            gy_mpch = gy_kms / H0
            gz_mpch = gz_kms / H0

            r_kms = math.sqrt(gx_kms*gx_kms + gy_kms*gy_kms + gz_kms*gz_kms)
            if r_kms < 1e-7:
                r_hat = [0.0, 1.0, 0.0]
            else:
                r_hat = [gx_kms / r_kms, gy_kms / r_kms, gz_kms / r_kms]

            v_los = vx * r_hat[0] + vy * r_hat[1] + vz * r_hat[2]

            sx_real = gx_mpch * MPC_H_TO_SCENE
            sy_real = gy_mpch * MPC_H_TO_SCENE
            sz_real = gz_mpch * MPC_H_TO_SCENE
            three_real = sg_to_three(sx_real, sy_real, sz_real)

            total_delta_v = v_los + v_thermal
            sx_z = sx_real + total_delta_v * r_hat[0]
            sy_z = sy_real + total_delta_v * r_hat[1]
            sz_z = sz_real + total_delta_v * r_hat[2]
            three_z = sg_to_three(sx_z, sy_z, sz_z)

            galaxies_real_mpch.extend([round(gx_mpch, 2), round(gy_mpch, 2), round(gz_mpch, 2)])
            galaxies_real_scene.extend([round(three_real[0], 1), round(three_real[1], 1), round(three_real[2], 1)])
            galaxies_redshift_scene.extend([round(three_z[0], 1), round(three_z[1], 1), round(three_z[2], 1)])
            v_los_list.append(round(v_los, 1))
            v_thermal_list.append(round(v_thermal, 1))
            target_classes.append(cl["cls"])
            clusters.append(cl["name"])

    data = {
        "survey": "SDSS",
        "name": "SDSS Great Wall & Northern Galactic Cap Filaments",
        "count": len(v_los_list),
        "h0": H0,
        "units": {
            "spatial": "Mpc/h",
            "velocity": "km/s",
            "scene": "scene_units (1 Mpc/h = 74.6 units)"
        },
        "pos_real_mpch": galaxies_real_mpch,
        "pos_real_scene": galaxies_real_scene,
        "pos_redshift_scene": galaxies_redshift_scene,
        "v_los": v_los_list,
        "v_thermal": v_thermal_list,
        "target_class": target_classes,
        "cluster_association": clusters
    }

    out_file = os.path.join(DATA_DIR, "survey_sdss.json")
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(',', ':'))
    print(f"[SDSS] Wrote {len(v_los_list)} galaxies to {out_file} ({os.path.getsize(out_file) / 1024:.1f} KB)")


def generate_2mrs_data():
    """
    Generates 2MASS Redshift Survey (2MRS) all-sky bright galaxy sample (~10,000 galaxies).
    Complete to Ks <= 11.75 across |b| > 5 deg out to cz ~ 10,000 km/s (134 Mpc/h).
    Penetrates the Zone of Avoidance.
    """
    random.seed(1234)
    galaxies_real_mpch = []
    galaxies_real_scene = []
    galaxies_redshift_scene = []
    v_los_list = []
    v_thermal_list = []
    target_classes = []
    clusters = []

    two_mrs_structures = [
        {"name": "Virgo Cluster", "pos_kms": [-280, 1300, -100], "count": 480, "sigma_v": 720, "spread": 600},
        {"name": "Centaurus (A3526)", "pos_kms": [-4200, 1200, 3100], "count": 450, "sigma_v": 850, "spread": 1000},
        {"name": "Norma / GA (A3627)", "pos_kms": [-4800, -850, 3900], "count": 520, "sigma_v": 920, "spread": 1100},
        {"name": "Hydra (A1060)", "pos_kms": [-3800, -2100, 2400], "count": 380, "sigma_v": 680, "spread": 850},
        {"name": "Perseus-Pisces Core", "pos_kms": [4800, -2500, -500], "count": 550, "sigma_v": 890, "spread": 1200},
        {"name": "Coma Cluster", "pos_kms": [500, 7000, 1500], "count": 500, "sigma_v": 1020, "spread": 1100},
        {"name": "Pavo-Indus Complex", "pos_kms": [-2900, -4500, -1200], "count": 420, "sigma_v": 740, "spread": 950},
        {"name": "Fornax Cluster", "pos_kms": [-1200, -1600, -800], "count": 320, "sigma_v": 450, "spread": 600},
        {"name": "Puppis Cluster (ZoA)", "pos_kms": [-6800, -3200, 500], "count": 360, "sigma_v": 780, "spread": 900},
        {"name": "Ophiuchus Cluster (ZoA)", "pos_kms": [-6500, 2800, 8200], "count": 340, "sigma_v": 850, "spread": 950},
        {"name": "Antlia Cluster", "pos_kms": [-2400, -900, 1800], "count": 280, "sigma_v": 580, "spread": 700},
        {"name": "Leo Cluster (A1367)", "pos_kms": [450, 6200, 2800], "count": 350, "sigma_v": 750, "spread": 900},
    ]

    total_count = 10000
    cluster_count = sum(s["count"] for s in two_mrs_structures)  # 4950
    field_count = total_count - cluster_count  # 5050

    # 1. Structure Galaxies (Virgo, GA, Coma, PP, ZoA)
    for st in two_mrs_structures:
        c_pos = st["pos_kms"]
        sigma_v = st["sigma_v"]
        spread = st["spread"]
        for i in range(st["count"]):
            u1, u2, u3 = random.random(), random.random(), random.random()
            rad = spread * math.pow(-2 * math.log(max(u1, 0.0001)), 0.35) * 0.8
            th = u2 * math.pi * 2
            ph = math.acos(2 * u3 - 1)
            gx_kms = c_pos[0] + rad * math.sin(ph) * math.cos(th)
            gy_kms = c_pos[1] + rad * math.cos(ph)
            gz_kms = c_pos[2] + rad * math.sin(ph) * math.sin(th)

            v_thermal = random.gauss(0, sigma_v)
            vx = random.gauss(0, sigma_v * 0.4)
            vy = random.gauss(0, sigma_v * 0.4)
            vz = random.gauss(0, sigma_v * 0.4)

            gx_mpch = gx_kms / H0
            gy_mpch = gy_kms / H0
            gz_mpch = gz_kms / H0

            r_kms = math.sqrt(gx_kms*gx_kms + gy_kms*gy_kms + gz_kms*gz_kms)
            if r_kms < 1e-7:
                r_hat = [0.0, 1.0, 0.0]
            else:
                r_hat = [gx_kms / r_kms, gy_kms / r_kms, gz_kms / r_kms]

            v_los = vx * r_hat[0] + vy * r_hat[1] + vz * r_hat[2]

            sx_real = gx_mpch * MPC_H_TO_SCENE
            sy_real = gy_mpch * MPC_H_TO_SCENE
            sz_real = gz_mpch * MPC_H_TO_SCENE
            three_real = sg_to_three(sx_real, sy_real, sz_real)

            total_delta_v = v_los + v_thermal
            sx_z = sx_real + total_delta_v * r_hat[0]
            sy_z = sy_real + total_delta_v * r_hat[1]
            sz_z = sz_real + total_delta_v * r_hat[2]
            three_z = sg_to_three(sx_z, sy_z, sz_z)

            galaxies_real_mpch.extend([round(gx_mpch, 2), round(gy_mpch, 2), round(gz_mpch, 2)])
            galaxies_real_scene.extend([round(three_real[0], 1), round(three_real[1], 1), round(three_real[2], 1)])
            galaxies_redshift_scene.extend([round(three_z[0], 1), round(three_z[1], 1), round(three_z[2], 1)])
            v_los_list.append(round(v_los, 1))
            v_thermal_list.append(round(v_thermal, 1))
            target_classes.append("2MRS_SPEC")
            clusters.append(st["name"])

    # 2. All-Sky Field Galaxies (penetrating ZoA)
    for i in range(field_count):
        u1, u2, u3 = random.random(), random.random(), random.random()
        r = 600 + math.pow(u1, 0.55) * 9400
        th = u2 * math.pi * 2
        ph = math.acos(2 * u3 - 1)

        gx_kms = r * math.sin(ph) * math.cos(th)
        gy_kms = r * math.cos(ph)
        gz_kms = r * math.sin(ph) * math.sin(th)

        vx = -180 + random.gauss(0, 160)
        vy = -120 + random.gauss(0, 160)
        vz = 140 + random.gauss(0, 160)
        sigma_v = 187.0
        v_thermal = random.gauss(0, sigma_v)

        gx_mpch = gx_kms / H0
        gy_mpch = gy_kms / H0
        gz_mpch = gz_kms / H0

        r_kms = math.sqrt(gx_kms*gx_kms + gy_kms*gy_kms + gz_kms*gz_kms)
        if r_kms < 1e-7:
            r_hat = [0.0, 1.0, 0.0]
        else:
            r_hat = [gx_kms / r_kms, gy_kms / r_kms, gz_kms / r_kms]

        v_los = vx * r_hat[0] + vy * r_hat[1] + vz * r_hat[2]

        sx_real = gx_mpch * MPC_H_TO_SCENE
        sy_real = gy_mpch * MPC_H_TO_SCENE
        sz_real = gz_mpch * MPC_H_TO_SCENE
        three_real = sg_to_three(sx_real, sy_real, sz_real)

        total_delta_v = v_los + v_thermal
        sx_z = sx_real + total_delta_v * r_hat[0]
        sy_z = sy_real + total_delta_v * r_hat[1]
        sz_z = sz_real + total_delta_v * r_hat[2]
        three_z = sg_to_three(sx_z, sy_z, sz_z)

        galaxies_real_mpch.extend([round(gx_mpch, 2), round(gy_mpch, 2), round(gz_mpch, 2)])
        galaxies_real_scene.extend([round(three_real[0], 1), round(three_real[1], 1), round(three_real[2], 1)])
        galaxies_redshift_scene.extend([round(three_z[0], 1), round(three_z[1], 1), round(three_z[2], 1)])
        v_los_list.append(round(v_los, 1))
        v_thermal_list.append(round(v_thermal, 1))
        target_classes.append("2MRS_SPEC")
        clusters.append("2MRS_Field")

    data = {
        "survey": "2MRS",
        "name": "2MASS Redshift Survey All-Sky Galaxy Sample",
        "count": len(v_los_list),
        "h0": H0,
        "units": {
            "spatial": "Mpc/h",
            "velocity": "km/s",
            "scene": "scene_units (1 Mpc/h = 74.6 units)"
        },
        "pos_real_mpch": galaxies_real_mpch,
        "pos_real_scene": galaxies_real_scene,
        "pos_redshift_scene": galaxies_redshift_scene,
        "v_los": v_los_list,
        "v_thermal": v_thermal_list,
        "target_class": target_classes,
        "cluster_association": clusters
    }

    out_file = os.path.join(DATA_DIR, "survey_2mrs.json")
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(',', ':'))
    print(f"[2MRS] Wrote {len(v_los_list)} galaxies to {out_file} ({os.path.getsize(out_file) / 1024:.1f} KB)")


def generate_desi_data():
    """
    Generates DESI Year 1 spectroscopic footprint (~4,000 galaxies)
    covering the Boötes core, surrounding acoustic horizon shell, and NGC/SGC footprints.
    Boötes Core: (3800, 16200, 8900) km/s
    """
    random.seed(9876)
    galaxies_real_mpch = []
    galaxies_real_scene = []
    galaxies_redshift_scene = []
    v_los_list = []
    v_thermal_list = []
    target_classes = []
    clusters = []

    bootes_kms = [3800, 16200, 8900]
    total_count = 4000
    shell_count = 2400
    field_count = total_count - shell_count  # 1600

    # 1. Galaxies delineating the DESI BAO sound horizon shell (radius rd = 147.1 Mpc * 74.6 = 10973 km/s)
    rd_kms = 147.1 * 74.6  # ~10973.66 km/s
    for i in range(shell_count):
        u1, u2 = random.random(), random.random()
        th = u1 * math.pi * 2
        ph = math.acos(2 * u2 - 1)
        # Shell thickness +/- 750 km/s (~10 Mpc)
        r_shell = rd_kms + random.gauss(0, 650)
        gx_kms = bootes_kms[0] + r_shell * math.sin(ph) * math.cos(th)
        gy_kms = bootes_kms[1] + r_shell * math.cos(ph)
        gz_kms = bootes_kms[2] + r_shell * math.sin(ph) * math.sin(th)

        # Infall / expansion velocity perturbations
        vx = random.gauss(0, 180)
        vy = random.gauss(0, 180)
        vz = random.gauss(0, 180)
        sigma_v = 187.0
        v_thermal = random.gauss(0, sigma_v)

        gx_mpch = gx_kms / H0
        gy_mpch = gy_kms / H0
        gz_mpch = gz_kms / H0

        r_kms = math.sqrt(gx_kms*gx_kms + gy_kms*gy_kms + gz_kms*gz_kms)
        if r_kms < 1e-7:
            r_hat = [0.0, 1.0, 0.0]
        else:
            r_hat = [gx_kms / r_kms, gy_kms / r_kms, gz_kms / r_kms]

        v_los = vx * r_hat[0] + vy * r_hat[1] + vz * r_hat[2]

        sx_real = gx_mpch * MPC_H_TO_SCENE
        sy_real = gy_mpch * MPC_H_TO_SCENE
        sz_real = gz_mpch * MPC_H_TO_SCENE
        three_real = sg_to_three(sx_real, sy_real, sz_real)

        total_delta_v = v_los + v_thermal
        sx_z = sx_real + total_delta_v * r_hat[0]
        sy_z = sy_real + total_delta_v * r_hat[1]
        sz_z = sz_real + total_delta_v * r_hat[2]
        three_z = sg_to_three(sx_z, sy_z, sz_z)

        galaxies_real_mpch.extend([round(gx_mpch, 2), round(gy_mpch, 2), round(gz_mpch, 2)])
        galaxies_real_scene.extend([round(three_real[0], 1), round(three_real[1], 1), round(three_real[2], 1)])
        galaxies_redshift_scene.extend([round(three_z[0], 1), round(three_z[1], 1), round(three_z[2], 1)])
        v_los_list.append(round(v_los, 1))
        v_thermal_list.append(round(v_thermal, 1))
        target_classes.append("LRG" if random.random() < 0.6 else "ELG")
        clusters.append("DESI_BAO_Shell")

    # 2. Galaxies in the Boötes core and DESI NGC survey footprint
    for i in range(field_count):
        u1, u2, u3 = random.random(), random.random(), random.random()
        r = 1500 + math.pow(u1, 0.5) * 16000
        th = u2 * math.pi * 2
        ph = 0.2 + u3 * (math.pi * 0.6)  # NGC footprint
        gx_kms = r * math.sin(ph) * math.cos(th)
        gy_kms = 4000 + r * math.cos(ph)
        gz_kms = 2000 + r * math.sin(ph) * math.sin(th)

        vx = random.gauss(0, 200)
        vy = random.gauss(0, 200)
        vz = random.gauss(0, 200)
        sigma_v = 187.0
        v_thermal = random.gauss(0, sigma_v)

        gx_mpch = gx_kms / H0
        gy_mpch = gy_kms / H0
        gz_mpch = gz_kms / H0

        r_kms = math.sqrt(gx_kms*gx_kms + gy_kms*gy_kms + gz_kms*gz_kms)
        if r_kms < 1e-7:
            r_hat = [0.0, 1.0, 0.0]
        else:
            r_hat = [gx_kms / r_kms, gy_kms / r_kms, gz_kms / r_kms]

        v_los = vx * r_hat[0] + vy * r_hat[1] + vz * r_hat[2]

        sx_real = gx_mpch * MPC_H_TO_SCENE
        sy_real = gy_mpch * MPC_H_TO_SCENE
        sz_real = gz_mpch * MPC_H_TO_SCENE
        three_real = sg_to_three(sx_real, sy_real, sz_real)

        total_delta_v = v_los + v_thermal
        sx_z = sx_real + total_delta_v * r_hat[0]
        sy_z = sy_real + total_delta_v * r_hat[1]
        sz_z = sz_real + total_delta_v * r_hat[2]
        three_z = sg_to_three(sx_z, sy_z, sz_z)

        galaxies_real_mpch.extend([round(gx_mpch, 2), round(gy_mpch, 2), round(gz_mpch, 2)])
        galaxies_real_scene.extend([round(three_real[0], 1), round(three_real[1], 1), round(three_real[2], 1)])
        galaxies_redshift_scene.extend([round(three_z[0], 1), round(three_z[1], 1), round(three_z[2], 1)])
        v_los_list.append(round(v_los, 1))
        v_thermal_list.append(round(v_thermal, 1))
        target_classes.append("BGS_BRIGHT" if random.random() < 0.5 else "LRG")
        clusters.append("DESI_NGC_Footprint")

    data = {
        "survey": "DESI",
        "name": "DESI Year 1 Spectroscopic Sample (Boötes & BAO Footprint)",
        "count": len(v_los_list),
        "h0": H0,
        "units": {
            "spatial": "Mpc/h",
            "velocity": "km/s",
            "scene": "scene_units (1 Mpc/h = 74.6 units)"
        },
        "pos_real_mpch": galaxies_real_mpch,
        "pos_real_scene": galaxies_real_scene,
        "pos_redshift_scene": galaxies_redshift_scene,
        "v_los": v_los_list,
        "v_thermal": v_thermal_list,
        "target_class": target_classes,
        "cluster_association": clusters
    }

    out_file = os.path.join(DATA_DIR, "survey_desi.json")
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(',', ':'))
    print(f"[DESI] Wrote {len(v_los_list)} galaxies to {out_file} ({os.path.getsize(out_file) / 1024:.1f} KB)")


if __name__ == "__main__":
    generate_sdss_data()
    generate_2mrs_data()
    generate_desi_data()

