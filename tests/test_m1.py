# -*- coding: utf-8 -*-
import pytest

def test_m1_astrophysics_math_and_engines(cdp):
    print("=== Starting Milestone M1 Astrophysics Math & Cosmographic Engines Test Suite ===")
    
    # 1. Verify window.cosmicflows Programmatic API Bridge Interface Contract
    cf_exists = cdp.evaluate("""
        (function() {
            const cf = window.cosmicflows;
            if (!cf) return false;
            return typeof cf.getVelocitySg === 'function' &&
                   typeof cf.getDivergenceSg === 'function' &&
                   typeof cf.getVorticitySg === 'function' &&
                   typeof cf.integrateRK4 === 'function' &&
                   typeof cf.switchScienceEngine === 'function' &&
                   typeof cf.toggleTheme === 'function' &&
                   typeof cf.exportPublicationFigure === 'function';
        })()
    """)
    assert cf_exists is True, "window.cosmicflows missing required API methods!"

    # 2. 3D Velocity vector evaluation
    v_vec = cdp.evaluate("""
        (function() {
            window.cosmicflows.switchScienceEngine('cf4-wf');
            const v = new THREE.Vector3();
            window.cosmicflows.getVelocitySg(new THREE.Vector3(1000, 2000, 3000), v);
            return { x: v.x, y: v.y, z: v.z, isFinite: isFinite(v.x) && isFinite(v.y) && isFinite(v.z) };
        })()
    """)
    assert v_vec["isFinite"] is True, "Velocity vector must be finite"

    # 3. Engines switch
    for engine in ['cf4-wf', 'cf4-hmc', 'cf4-zoa', 'vweb', 'twomrs', 'gadget4', 'nusser']:
        res = cdp.evaluate(f"window.cosmicflows.switchScienceEngine('{engine}')")
        assert res is not False, f"Engine {engine} switch failed"
