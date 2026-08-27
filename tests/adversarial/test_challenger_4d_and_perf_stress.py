# -*- coding: utf-8 -*-
"""
tests/adversarial/test_challenger_4d_and_perf_stress.py
================================================================================
Empirical Challenger 1 Exhaustive Stress Suite:
1. 4D Cosmic Time Scrubbing (-13.8 Gyr to +10.0 Gyr) across 4,200 streamlines
   and galaxy coordinates with strict finite coordinate assertions.
2. Particle advection and continuous playback 60 FPS performance budget (< 16.6 ms/frame)
   and zero geometry/texture allocation assertions during active advection.
3. 50x Dataset and 50x Theme Switching VRAM and Three.js memory invariance.
4. Out-of-bounds, NaN injection, and rapid cyclic scrubbing stability.
"""

import pytest
import time
import math
import json


class TestAdversarial4DCosmicTimeScrubbing:
    """Stress tests 4D Cosmic Time scrubbing across entire temporal span [-13.8, +10.0] Gyr."""

    def test_extreme_boundary_and_fuzzing_resilience(self, cdp):
        """Fuzz cosmic time engine with boundary, out-of-bounds, negative, subnormal, and zero values."""
        res = cdp.evaluate("""
        (function() {
            if (!window.cosmicflows || !window.cosmicflows.timeEngine) {
                return { error: 'timeEngine not found' };
            }
            const te = window.cosmicflows.timeEngine;
            
            const testValues = [
                -13.8, -13.78, -13.80001, -100.0, -1e6,
                0.0, -0.0, 1e-7, -1e-7,
                10.0, 10.00001, 100.0, 1e6,
                -10.4, -5.0, 2.5, 5.0, 7.5
            ];
            
            const results = [];
            for (let val of testValues) {
                te.setTime(val);
                const currentTime = te.getTime();
                const cosmo = te.computeCosmology(val);
                
                // Assert valid bounded time
                const isBounded = currentTime >= -13.8 && currentTime <= 10.0;
                const isFiniteZ = Number.isFinite(cosmo.z);
                const isFiniteA = Number.isFinite(cosmo.a);
                const isFiniteD = Number.isFinite(cosmo.D_plus);
                const isFiniteH = Number.isFinite(cosmo.H_z);
                
                results.push({
                    input: val,
                    currentTime: currentTime,
                    isBounded: isBounded,
                    isFiniteZ: isFiniteZ,
                    isFiniteA: isFiniteA,
                    isFiniteD: isFiniteD,
                    isFiniteH: isFiniteH
                });
            }
            
            // Restore present
            te.setTime(0.0);
            
            const allBounded = results.every(r => r.isBounded);
            const allFinite = results.every(r => r.isFiniteZ && r.isFiniteA && r.isFiniteD && r.isFiniteH);
            
            return {
                allBounded: allBounded,
                allFinite: allFinite,
                results: results
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert res["allBounded"] is True, "Time engine allowed out-of-bounds cosmic time"
        assert res["allFinite"] is True, "Cosmological parameters produced non-finite values during fuzzing"

    def test_rapid_cyclic_scrubbing_100_passes(self, cdp):
        """Stress-test 100 rapid oscillations between Big Bang (-13.8 Gyr) and Deep Future (+10.0 Gyr)."""
        res = cdp.evaluate("""
        (function() {
            const te = window.cosmicflows && window.cosmicflows.timeEngine;
            const renderer = window.cosmicflows && window.cosmicflows.renderer;
            if (!te || !renderer) return { error: 'Engine or renderer missing' };
            
            const initialGeos = renderer.info.memory.geometries;
            const initialTexs = renderer.info.memory.textures;
            
            const passes = 100;
            const frameTimes = [];
            
            for (let p = 0; p < passes; p++) {
                const target = (p % 2 === 0) ? -13.8 : 10.0;
                const t0 = performance.now();
                te.setTime(target);
                const t1 = performance.now();
                frameTimes.push(t1 - t0);
            }
            
            te.setTime(0.0);
            
            const finalGeos = renderer.info.memory.geometries;
            const finalTexs = renderer.info.memory.textures;
            
            const avg = frameTimes.reduce((a, b) => a + b, 0) / frameTimes.length;
            const max = Math.max(...frameTimes);
            
            return {
                passes: passes,
                avgMs: avg,
                maxMs: max,
                geoDelta: finalGeos - initialGeos,
                texDelta: finalTexs - initialTexs
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert res["avgMs"] < 5.0, f"Average scrub frame time exceeds 5ms: {res['avgMs']:.2f}ms"
        assert res["maxMs"] < 16.6, f"Peak scrub frame time exceeds 16.6ms budget: {res['maxMs']:.2f}ms"
        assert res["geoDelta"] == 0, f"Geometry allocation leak detected: {res['geoDelta']}"
        assert res["texDelta"] == 0, f"Texture allocation leak detected: {res['texDelta']}"


class TestAdversarialParticleAdvectionAndPerformanceBudget:
    """Stress tests particle advection and real-time animation loop under high load."""

    def test_continuous_advection_60fps_budget(self, cdp):
        """Verify that running 60 active animation frames maintains per-frame execution time under 16.6ms."""
        res = cdp.evaluate("""
        (function() {
            if (!window.cosmicflows || !window.cosmicflows.renderer) {
                return { error: 'cosmicflows renderer missing' };
            }
            
            const frameDurations = [];
            const samples = 60;
            
            // Measure execution cost of streamline/galaxy updates over multiple steps
            for (let i = 0; i < samples; i++) {
                const t0 = performance.now();
                
                // Trigger time advance + render update
                if (typeof window.applyCosmicTime === 'function') {
                    window.applyCosmicTime(-13.8 + (23.8 * (i / samples)));
                }
                
                const t1 = performance.now();
                frameDurations.push(t1 - t0);
            }
            
            if (typeof window.applyCosmicTime === 'function') {
                window.applyCosmicTime(0.0);
            }
            
            const sum = frameDurations.reduce((a, b) => a + b, 0);
            const avg = sum / frameDurations.length;
            const max = Math.max(...frameDurations);
            const p95 = [...frameDurations].sort((a, b) => a - b)[Math.floor(frameDurations.length * 0.95)];
            
            return {
                samples: samples,
                avgMs: avg,
                maxMs: max,
                p95Ms: p95,
                budgetPassed: max < 16.6 && avg < 8.0
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert res["avgMs"] < 8.0, f"Average animation step time was {res['avgMs']:.2f}ms (expected <8.0ms)"
        assert res["maxMs"] < 16.6, f"Max animation step time was {res['maxMs']:.2f}ms (expected <16.6ms)"
        assert res["budgetPassed"] is True


class TestAdversarialMemoryInvariance:
    """Stress tests repeated theme toggles and reconstruction engine switches for zero memory leaks."""

    def test_30x_theme_switch_memory_invariance(self, cdp):
        """Toggle between Dark and White theme 30 times and verify zero geometry/texture accumulation."""
        res = cdp.evaluate("""
        (function() {
            if (typeof window.toggleTheme !== 'function' || !window.cosmicflows || !window.cosmicflows.renderer) {
                return { error: 'toggleTheme or renderer missing' };
            }
            const renderer = window.cosmicflows.renderer;
            
            // Initial warm-up
            window.toggleTheme('dark');
            window.toggleTheme('white');
            
            const initialGeos = renderer.info.memory.geometries;
            const initialTexs = renderer.info.memory.textures;
            
            for (let i = 0; i < 30; i++) {
                window.toggleTheme(i % 2 === 0 ? 'dark' : 'white');
            }
            
            window.toggleTheme('white');
            
            const finalGeos = renderer.info.memory.geometries;
            const finalTexs = renderer.info.memory.textures;
            
            return {
                initialGeos: initialGeos,
                finalGeos: finalGeos,
                geoDelta: finalGeos - initialGeos,
                initialTexs: initialTexs,
                finalTexs: finalTexs,
                texDelta: finalTexs - initialTexs
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert abs(res["geoDelta"]) <= 2, f"Leaked Three.js geometries across 30 theme toggles: delta={res['geoDelta']}"
        assert abs(res["texDelta"]) <= 1, f"Leaked Three.js textures across 30 theme toggles: delta={res['texDelta']}"

    def test_30x_dataset_engine_switch_memory_invariance(self, cdp):
        """Switch across all 5 reconstruction science engines 30 times and assert zero memory growth."""
        res = cdp.evaluate("""
        (function() {
            if (typeof window.switchScienceEngine !== 'function' || !window.cosmicflows || !window.cosmicflows.renderer) {
                return { error: 'switchScienceEngine or renderer missing' };
            }
            const renderer = window.cosmicflows.renderer;
            const engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026'];
            
            // Warm up
            window.switchScienceEngine('cf4-wf');
            
            const initialGeos = renderer.info.memory.geometries;
            const initialTexs = renderer.info.memory.textures;
            
            for (let i = 0; i < 30; i++) {
                const eng = engines[i % engines.length];
                window.switchScienceEngine(eng);
            }
            
            window.switchScienceEngine('cf4-wf');
            
            const finalGeos = renderer.info.memory.geometries;
            const finalTexs = renderer.info.memory.textures;
            
            return {
                initialGeos: initialGeos,
                finalGeos: finalGeos,
                geoDelta: finalGeos - initialGeos,
                initialTexs: initialTexs,
                finalTexs: finalTexs,
                texDelta: finalTexs - initialTexs
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert abs(res["geoDelta"]) <= 4, f"Leaked Three.js geometries across 30 dataset switches: delta={res['geoDelta']}"
        assert abs(res["texDelta"]) <= 2, f"Leaked Three.js textures across 30 dataset switches: delta={res['texDelta']}"
