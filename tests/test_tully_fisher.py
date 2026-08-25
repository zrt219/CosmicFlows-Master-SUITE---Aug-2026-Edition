"""
Automated Test Suite for Tully-Fisher Kinematic & Distance Resolver Subsystem
Validates:
- Multi-band TFR calibrations (WISE W1, 2MASS Ks, Optical I, Optical B)
- Baryonic Tully-Fisher (BTFR) mass scaling
- Inclination deprojection and turbulent motion correction
- Extinction-corrected distance modulus and metric distance inversion
- Radial peculiar velocity derivation (v_pec = cz - H0*d)
- Benchmark galaxy presets
- Interactive UI DOM reactivity
- 3D WebGL probe seeding in Three.js scene
- Publication LaTeX caption modal generation
"""

import pytest
import math

class TestTullyFisherMath:
    def test_tfr_api_exposure(self, cdp):
        res = cdp.evaluate("""
        (function() {
            return {
                hasTF: !!(window.cosmicflows && window.cosmicflows.tullyFisher),
                hasSolve: typeof window.solveTullyFisher === 'function',
                hasCalibrations: !!(window.cosmicflows && window.cosmicflows.tullyFisher && window.cosmicflows.tullyFisher.CALIBRATIONS),
                hasPresets: !!(window.cosmicflows && window.cosmicflows.tullyFisher && window.cosmicflows.tullyFisher.PRESETS),
                hasBTFR: !!(window.cosmicflows && window.cosmicflows.tullyFisher && window.cosmicflows.tullyFisher.BTFR)
            };
        })()
        """)
        assert res["hasTF"] is True
        assert res["hasSolve"] is True
        assert res["hasCalibrations"] is True
        assert res["hasPresets"] is True
        assert res["hasBTFR"] is True

    def test_tfr_wise_w1_standard_calculation(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const sol = window.solveTullyFisher({
                wObs: 450,
                incl: 90,
                mApp: 8.12,
                band: 'W1',
                aExt: 0.15,
                cz: 528,
                H0: 74.6
            });
            return sol;
        })()
        """)
        # wCorr = 450/sin(90) - 6.5 = 443.5 km/s
        # vRot = 443.5 / 2 = 221.75 km/s
        # log10(443.5) = 2.64689
        # mAbs = -9.75 * (2.64689 - 2.5) - 21.84 = -9.75 * 0.14689 - 21.84 = -1.432 - 21.84 = -23.272
        # mu = 8.12 - 0.15 - (-23.272) = 31.242
        # d = 10^((31.242-25)/5) = 10^1.2484 = 17.72 Mpc
        # vHubble = 17.72 * 74.6 = 1321.9 km/s
        # vPec = 528 - 1321.9 = -793.9 km/s (Infall)
        assert abs(res["wCorr"] - 443.5) < 0.1
        assert abs(res["vRot"] - 221.75) < 0.1
        assert abs(res["mAbs"] - (-23.272)) < 0.1
        assert res["distMpc"] > 0
        assert res["vPec"] < 0
        assert res["mBaryonic"] > 1e11

    def test_tfr_multiband_calibrations(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const bands = ['W1', 'Ks', 'I', 'B'];
            const results = {};
            for (let b of bands) {
                results[b] = window.solveTullyFisher({
                    wObs: 400,
                    incl: 80,
                    mApp: 10.0,
                    band: b,
                    cz: 1000,
                    H0: 75.0
                });
            }
            return results;
        })()
        """)
        assert "W1" in res and "Ks" in res and "I" in res and "B" in res
        # Check that NIR/IR bands give brighter absolute magnitudes for standard spirals
        assert res["W1"]["mAbs"] < res["B"]["mAbs"]
        for b in ["W1", "Ks", "I", "B"]:
            assert res[b]["distMpc"] > 0
            assert math.isfinite(res[b]["vPec"])

    def test_tfr_inclination_deprojection_and_clamping(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const sol90 = window.solveTullyFisher({ wObs: 300, incl: 90 });
            const sol45 = window.solveTullyFisher({ wObs: 300, incl: 45 });
            const solLow = window.solveTullyFisher({ wObs: 300, incl: 5 }); // Should clamp to >= 15 deg
            return {
                wCorr90: sol90.wCorr,
                wCorr45: sol45.wCorr,
                wCorrLow: solLow.wCorr,
                inclLow: solLow.inclDeg
            };
        })()
        """)
        assert res["wCorr45"] > res["wCorr90"]
        assert res["inclLow"] == 15.0
        assert res["wCorrLow"] > 0

    def test_tfr_baryonic_btfr_scaling(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const sol1 = window.solveTullyFisher({ wObs: 200, incl: 90 });
            const sol2 = window.solveTullyFisher({ wObs: 400, incl: 90 });
            return {
                mb1: sol1.mBaryonic,
                mb2: sol2.mBaryonic,
                vRot1: sol1.vRot,
                vRot2: sol2.vRot,
                ratio: sol2.mBaryonic / sol1.mBaryonic
            };
        })()
        """)
        # If rotation speed doubles (~2x), mass scales as 2^4 = 16x
        vRatio = res["vRot2"] / res["vRot1"]
        expectedRatio = math.pow(vRatio, 4.0)
        assert abs(res["ratio"] - expectedRatio) < 0.1

    def test_tfr_presets_validation(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const presets = window.cosmicflows.tullyFisher.PRESETS;
            const keys = Object.keys(presets);
            const verified = {};
            for (let k of keys) {
                const p = presets[k];
                const sol = window.solveTullyFisher({
                    wObs: p.wObs,
                    incl: p.incl,
                    mApp: p.mApp,
                    band: p.band,
                    aExt: p.aExt,
                    cz: p.cz
                });
                verified[k] = {
                    name: p.name,
                    distMpc: sol.distMpc,
                    vPec: sol.vPec,
                    vRot: sol.vRot
                };
            }
            return { count: keys.length, verified };
        })()
        """)
        assert res["count"] >= 7
        assert "ngc891" in res["verified"]
        assert "m31" in res["verified"]
        assert "ngc5128" in res["verified"]
        # M31 has negative recession velocity (blueshifted local group approaching)
        assert res["verified"]["m31"]["vPec"] < 0


class TestTullyFisherUIAnd3DScene:
    def test_tfr_tab_and_inputs_exist_in_dom(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const tabBtn = document.querySelector('.tab-btn[data-tab="tab-tfr"]');
            const pane = document.getElementById('tab-tfr');
            const selPreset = document.getElementById('sel-tfr-preset');
            const rngW = document.getElementById('rng-tfr-w');
            const rngIncl = document.getElementById('rng-tfr-incl');
            const rngMapp = document.getElementById('rng-tfr-mapp');
            const btnSeed = document.getElementById('btn-tfr-seed');
            const btnLatex = document.getElementById('btn-tfr-latex');

            return {
                hasTabBtn: !!tabBtn,
                hasPane: !!pane,
                hasPreset: !!selPreset,
                hasRngW: !!rngW,
                hasRngIncl: !!rngIncl,
                hasRngMapp: !!rngMapp,
                hasBtnSeed: !!btnSeed,
                hasBtnLatex: !!btnLatex
            };
        })()
        """)
        assert res["hasTabBtn"] is True
        assert res["hasPane"] is True
        assert res["hasPreset"] is True
        assert res["hasRngW"] is True
        assert res["hasRngIncl"] is True
        assert res["hasRngMapp"] is True
        assert res["hasBtnSeed"] is True
        assert res["hasBtnLatex"] is True

    def test_tfr_ui_preset_change_updates_outputs(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const sel = document.getElementById('sel-tfr-preset');
            sel.value = 'm31';
            sel.dispatchEvent(new Event('change'));
            
            const vRotText = document.getElementById('res-tfr-vrot').textContent;
            const distText = document.getElementById('res-tfr-dist').textContent;
            const vPecText = document.getElementById('res-tfr-vpec').textContent;

            return { vRotText, distText, vPecText };
        })()
        """)
        assert "km/s" in res["vRotText"]
        assert "Mpc" in res["distText"]
        assert "Infall" in res["vPecText"] or "Inflow" in res["vPecText"]

    def test_tfr_seed_probe_creates_3d_object_in_scene(self, cdp):
        res = cdp.evaluate("""
        (function() {
            window.cosmicflows.tullyFisher.seedProbe();
            const probeGroup = window.cosmicflows.scene.getObjectByName('tfrProbeGroup');
            return {
                hasGroup: !!probeGroup,
                childrenCount: probeGroup ? probeGroup.children.length : 0
            };
        })()
        """)
        assert res["hasGroup"] is True
        # Sphere, Ring Aura, Streamline Line, LOS Dashed Line
        assert res["childrenCount"] >= 3

    def test_tfr_latex_export_modal_generation(self, cdp):
        res = cdp.evaluate("""
        (function() {
            window.cosmicflows.tullyFisher.exportLatex();
            const modal = document.getElementById('caption-modal');
            const code = document.getElementById('latex-code').value;
            const isVisible = modal.style.display === 'block';
            // Clean up
            modal.style.display = 'none';
            return {
                isVisible,
                hasLatexCode: code.includes('\\\\label{eq:tfr_calibration}') && code.includes('v_{\\\\text{pec}}')
            };
        })()
        """)
        assert res["isVisible"] is True
        assert res["hasLatexCode"] is True
