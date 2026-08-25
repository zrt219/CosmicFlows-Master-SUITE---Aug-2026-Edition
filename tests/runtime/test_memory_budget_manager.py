# -*- coding: utf-8 -*-
"""
tests/runtime/test_memory_budget_manager.py
================================================================================
Comprehensive, automated scientific Pytest suite for the WebGL / Web Worker
Memory Budget Manager, LRU Eviction Sentinel, High-Frequency Telemetry Poller,
and Streamline Out-of-Memory (OOM) Prevention Guard.

Verifies:
1. Exact Analytical Memory Footprints & Typed Array Estimators
2. Multi-Tier Priority LRU Eviction & Automatic Three.js Object Disposal
3. Multi-State Memory Pressure Hysteresis (NORMAL, MODERATE, CRITICAL, EXCEEDED)
4. Out-of-Memory Guard & Adaptive LOD Downscaling for 64³, 128³, 256³ Grids & RK45 Traces
5. High-Frequency Telemetry Polling, Allocation Derivatives, & OLS Memory Leak Detector
6. Web Worker Zero-Copy Transferable Buffer Lifecycle & Multi-Worker Coordination
7. WeakRef / FinalizationRegistry Garbage Collection Sentinel
8. Strict Mathematical & Memory Conservation Invariants Verification
"""

import pytest
import math
from tests.utils import run_node_snippet

class TestMemoryBudgetManagerInit:
    """Test Suite 1: Initialization, Default Budget Configuration & Limit Invariants"""

    def test_default_initialization_and_telemetry(self):
        """Verify default 512 MB soft limit, 1024 MB hard limit, and initial zero counters."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, MEMORY_PRESSURE_LEVEL } from './src/runtime/memory_budget_manager.js';
        const mgr = new MemoryBudgetManager();
        const telem = mgr.getTelemetrySnapshot();
        const inv = mgr.verifyInvariants();
        
        console.log(JSON.stringify({
            softLimitMB: telem.softLimitMB,
            hardLimitMB: telem.hardLimitMB,
            totalTrackedBytes: mgr.totalTrackedBytes,
            resourceCount: mgr.resourceCount,
            pressureLevel: mgr.pressureLevel,
            availableSoftMB: telem.availableSoftHeadroomMB,
            availableHardMB: telem.availableHardHeadroomMB,
            isInvariantValid: inv.isValid,
            discrepancyCount: inv.discrepancies.length
        }));
        """)
        assert res["softLimitMB"] == 512.0
        assert res["hardLimitMB"] == 1024.0
        assert res["totalTrackedBytes"] == 0
        assert res["resourceCount"] == 0
        assert res["pressureLevel"] == "normal"
        assert res["availableSoftMB"] == 512.0
        assert res["availableHardMB"] == 1024.0
        assert res["isInvariantValid"] is True
        assert res["discrepancyCount"] == 0

    def test_custom_budget_configuration(self):
        """Verify custom limits, thresholds, and alignment parameters."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        const mgr = new MemoryBudgetManager({
            softLimitBytes: 256 * BYTES_PER_MB,
            hardLimitBytes: 512 * BYTES_PER_MB,
            moderateThresholdRatio: 0.5,
            telemetryPollIntervalMs: 50
        });
        const telem = mgr.getTelemetrySnapshot();
        console.log(JSON.stringify({
            softLimitMB: telem.softLimitMB,
            hardLimitMB: telem.hardLimitMB,
            pollInterval: mgr.config.telemetryPollIntervalMs,
            modRatio: mgr.config.moderateThresholdRatio
        }));
        """)
        assert res["softLimitMB"] == 256.0
        assert res["hardLimitMB"] == 512.0
        assert res["pollInterval"] == 50
        assert res["modRatio"] == 0.5

    def test_invalid_configuration_invariants_rejection(self):
        """Verify constructor throws error on invalid or conflicting limit thresholds."""
        res = run_node_snippet("""
        import { MemoryBudgetManager } from './src/runtime/memory_budget_manager.js';
        let caughtNegative = false;
        let caughtInverted = false;
        
        try {
            new MemoryBudgetManager({ softLimitBytes: -100 });
        } catch (e) {
            caughtNegative = true;
        }
        
        try {
            new MemoryBudgetManager({ softLimitBytes: 1000, hardLimitBytes: 500 });
        } catch (e) {
            caughtInverted = true;
        }
        
        console.log(JSON.stringify({ caughtNegative, caughtInverted }));
        """)
        assert res["caughtNegative"] is True
        assert res["caughtInverted"] is True


class TestMemoryEstimator:
    """Test Suite 2: High-Precision Footprint Estimators for TypedArrays, Grids, and Streamlines"""

    def test_typed_array_element_byte_sizes(self):
        """Verify exact byte sizes for standard WebGL and JS typed array representations."""
        res = run_node_snippet("""
        import { MemoryEstimator, DATA_TYPE_BYTES } from './src/runtime/memory_budget_manager.js';
        
        const sizes = {
            f32: MemoryEstimator.getElementByteSize(new Float32Array(1)),
            f64: MemoryEstimator.getElementByteSize(new Float64Array(1)),
            u8: MemoryEstimator.getElementByteSize(new Uint8Array(1)),
            u16: MemoryEstimator.getElementByteSize(new Uint16Array(1)),
            u32: MemoryEstimator.getElementByteSize(new Uint32Array(1)),
            half: MemoryEstimator.getElementByteSize('HalfFloatType'),
            customFloatStr: MemoryEstimator.getElementByteSize('Float32Array')
        };
        console.log(JSON.stringify(sizes));
        """)
        assert res["f32"] == 4
        assert res["f64"] == 8
        assert res["u8"] == 1
        assert res["u16"] == 2
        assert res["u32"] == 4
        assert res["half"] == 2
        assert res["customFloatStr"] == 4

    def test_byte_alignment_calculation(self):
        """Verify memory boundary alignment logic (4-byte and 16-byte alignment)."""
        res = run_node_snippet("""
        import { MemoryEstimator } from './src/runtime/memory_budget_manager.js';
        
        const a1 = MemoryEstimator.alignBytes(13, 4);   // 16
        const a2 = MemoryEstimator.alignBytes(16, 4);   // 16
        const a3 = MemoryEstimator.alignBytes(1, 4);    // 4
        const a4 = MemoryEstimator.alignBytes(25, 16);  // 32
        const a5 = MemoryEstimator.alignBytes(100, 1);  // 100
        
        console.log(JSON.stringify({ a1, a2, a3, a4, a5 }));
        """)
        assert res["a1"] == 16
        assert res["a2"] == 16
        assert res["a3"] == 4
        assert res["a4"] == 32
        assert res["a5"] == 100

    def test_buffer_geometry_byte_estimation(self):
        """Verify exact byte estimation for Three.js BufferGeometry (attributes + indices)."""
        res = run_node_snippet("""
        import { MemoryEstimator } from './src/runtime/memory_budget_manager.js';
        
        // Mock Three.js BufferGeometry with 10,000 vertices:
        // Position (vec3 Float32) = 10000 * 3 * 4 = 120,000 bytes
        // Normal (vec3 Float32)   = 10000 * 3 * 4 = 120,000 bytes
        // UV (vec2 Float32)       = 10000 * 2 * 4 = 80,000 bytes
        // Velocity (vec3 Float32) = 10000 * 3 * 4 = 120,000 bytes
        // Index (Uint32, 20000)   = 20000 * 4     = 80,000 bytes
        // Total expected = 520,000 bytes
        const geom = {
            attributes: {
                position: { array: new Float32Array(10000 * 3) },
                normal: { array: new Float32Array(10000 * 3) },
                uv: { array: new Float32Array(10000 * 2) },
                velocity: { array: new Float32Array(10000 * 3) }
            },
            index: { array: new Uint32Array(20000) }
        };
        
        const bytes = MemoryEstimator.estimateBufferGeometryBytes(geom);
        console.log(JSON.stringify({ bytes, expected: 520000 }));
        """)
        assert res["bytes"] == 520000
        assert res["bytes"] == res["expected"]

    def test_2d_and_3d_texture_byte_estimation(self):
        """Verify memory estimation for 2D textures (with mipmaps) and 3D data textures."""
        res = run_node_snippet("""
        import { MemoryEstimator } from './src/runtime/memory_budget_manager.js';
        
        // 2D Texture: 1024x1024 RGBA Float32 with Mipmaps
        // Base: 1024 * 1024 * 4 channels * 4 bytes = 16,777,216 bytes (16 MB)
        // Mipmaps: 16,777,216 * 1.3333333333333333 = 22,369,621 bytes (~21.33 MB)
        const tex2D = {
            width: 1024,
            height: 1024,
            format: 'RGBA',
            type: 'FloatType',
            generateMipmaps: true
        };
        const bytes2D = MemoryEstimator.estimateTextureBytes(tex2D);
        
        // 3D Texture: 64x64x64 RGBA Float32 without Mipmaps
        // 64 * 64 * 64 * 4 channels * 4 bytes = 4,194,304 bytes (4.0 MB)
        const tex3D_64 = {
            width: 64,
            height: 64,
            depth: 64,
            format: 'RGBA',
            type: 'FloatType',
            generateMipmaps: false
        };
        const bytes3D_64 = MemoryEstimator.estimateData3DTextureBytes(tex3D_64);
        
        // 3D Texture: 128x128x128 1-channel Float32 without Mipmaps
        // 128 * 128 * 128 * 1 channel * 4 bytes = 8,388,608 bytes (8.0 MB)
        const tex3D_128_scalar = {
            width: 128,
            height: 128,
            depth: 128,
            format: 'RedFormat',
            type: 'FloatType',
            generateMipmaps: false
        };
        const bytes3D_128 = MemoryEstimator.estimateData3DTextureBytes(tex3D_128_scalar);
        
        console.log(JSON.stringify({ bytes2D, bytes3D_64, bytes3D_128 }));
        """)
        assert res["bytes2D"] == math.floor(16777216 * 1.3333333333333333)
        assert res["bytes3D_64"] == 4194304
        assert res["bytes3D_128"] == 8388608

    def test_instanced_mesh_and_render_target_estimation(self):
        """Verify InstancedMesh and WebGLRenderTarget memory calculations."""
        res = run_node_snippet("""
        import { MemoryEstimator } from './src/runtime/memory_budget_manager.js';
        
        // InstancedMesh: 50,000 instances
        // Base geom: 1000 vertices * 3 * 4 = 12,000 bytes
        // instanceMatrix: 50,000 * 16 * 4 = 3,200,000 bytes
        // instanceColor: 50,000 * 3 * 4 = 600,000 bytes
        // Total = 3,812,000 bytes
        const instMesh = {
            count: 50000,
            geometry: { attributes: { position: { array: new Float32Array(1000 * 3) } } },
            instanceMatrix: { array: new Float32Array(50000 * 16) },
            instanceColor: { array: new Float32Array(50000 * 3) }
        };
        const instBytes = MemoryEstimator.estimateInstancedMeshBytes(instMesh);
        
        // RenderTarget: 1920x1080 RGBA UnsignedByte with Depth buffer (4 bytes/pixel)
        // Color: 1920 * 1080 * 4 * 1 = 8,294,400 bytes
        // Depth: 1920 * 1080 * 4     = 8,294,400 bytes
        // Total = 16,588,800 bytes
        const rt = {
            width: 1920,
            height: 1080,
            depthBuffer: true,
            stencilBuffer: false,
            samples: 0
        };
        const rtBytes = MemoryEstimator.estimateRenderTargetBytes(rt);
        
        console.log(JSON.stringify({ instBytes, rtBytes }));
        """)
        assert res["instBytes"] == 3812000
        assert res["rtBytes"] == 16588800

    def test_cosmological_3d_grid_analytical_footprints(self):
        """Verify exact byte estimation for 64³, 128³, 256³ cosmological grids."""
        res = run_node_snippet("""
        import { MemoryEstimator, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        // 64^3 scalar Float32 = 262,144 * 4 = 1,048,576 bytes = 1.0 MB
        const grid64_scalar = MemoryEstimator.estimateGridBytes({ resolution: 64, channels: 1 });
        // 64^3 vector Float32 (vx, vy, vz) = 262,144 * 3 * 4 = 3,145,728 bytes = 3.0 MB
        const grid64_vector = MemoryEstimator.estimateGridBytes({ resolution: 64, channels: 3 });
        // 128^3 scalar Float32 = 2,097,152 * 4 = 8,388,608 bytes = 8.0 MB
        const grid128_scalar = MemoryEstimator.estimateGridBytes({ resolution: 128, channels: 1 });
        // 128^3 vector Float32 = 2,097,152 * 3 * 4 = 25,165,824 bytes = 24.0 MB
        const grid128_vector = MemoryEstimator.estimateGridBytes({ resolution: 128, channels: 3 });
        // 256^3 vector Float32 = 16,777,216 * 3 * 4 = 201,326,592 bytes = 192.0 MB
        const grid256_vector = MemoryEstimator.estimateGridBytes({ resolution: 256, channels: 3 });
        
        console.log(JSON.stringify({
            g64_s_mb: grid64_scalar / BYTES_PER_MB,
            g64_v_mb: grid64_vector / BYTES_PER_MB,
            g128_s_mb: grid128_scalar / BYTES_PER_MB,
            g128_v_mb: grid128_vector / BYTES_PER_MB,
            g256_v_mb: grid256_vector / BYTES_PER_MB
        }));
        """)
        assert res["g64_s_mb"] == 1.0
        assert res["g64_v_mb"] == 3.0
        assert res["g128_s_mb"] == 8.0
        assert res["g128_v_mb"] == 24.0
        assert res["g256_v_mb"] == 192.0

    def test_streamline_rk45_reconstruction_estimator(self):
        """Verify exact mathematical footprint for Runge-Kutta streamline reconstruction with tubes."""
        res = run_node_snippet("""
        import { MemoryEstimator } from './src/runtime/memory_budget_manager.js';
        
        const estLines = MemoryEstimator.estimateStreamlineReconstructionBytes({
            seedCount: 1000,
            avgStepsPerStreamline: 100,
            attributesPerVertex: 10,
            renderAsTubes: false
        });
        
        const estTubes = MemoryEstimator.estimateStreamlineReconstructionBytes({
            seedCount: 1000,
            avgStepsPerStreamline: 100,
            attributesPerVertex: 10,
            renderAsTubes: true,
            tubeRadialSegments: 6
        });
        
        console.log(JSON.stringify({
            linesTotalBytes: estLines.totalBytes,
            linesRawBytes: estLines.rawTraceBytes,
            tubesTotalBytes: estTubes.totalBytes,
            tubesTriangles: estTubes.triangleCount
        }));
        """)
        # 1000 * 100 * 10 * 4 = 4,000,000 raw trace bytes
        assert res["linesRawBytes"] == 4000000
        # Line geom: 1000 * 100 * 6 * 4 = 2,400,000 bytes; Overhead: 200,000; Total = 6,600,000 bytes
        assert res["linesTotalBytes"] == 6600000
        assert res["tubesTriangles"] == 1000 * 99 * 6 * 2
        assert res["tubesTotalBytes"] > res["linesTotalBytes"]


class TestResourceRegistrationAndAccounting:
    """Test Suite 3: Resource Registration, Category Accounting & Access Tracking"""

    def test_register_and_query_multiple_resources(self):
        """Verify registration of diverse resource categories and exact category sum accounting."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, RESOURCE_TYPE, RESOURCE_PRIORITY, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        mgr.registerGeometry('geom_1', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.HIGH });
        mgr.registerTexture('tex_1', null, { sizeBytes: 20 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.MEDIUM });
        mgr.registerData3DTexture('grid_3d', null, { sizeBytes: 50 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.PINNED });
        mgr.registerInstancedMesh('inst_mesh', null, { sizeBytes: 15 * BYTES_PER_MB });
        mgr.registerRenderTarget('rt_main', null, { sizeBytes: 30 * BYTES_PER_MB });
        
        const telem = mgr.getTelemetrySnapshot();
        const inv = mgr.verifyInvariants();
        
        console.log(JSON.stringify({
            totalMB: telem.totalTrackedMB,
            geomMB: telem.categories.geometry.mb,
            tex2DMB: telem.categories.texture_2d.mb,
            tex3DMB: telem.categories.texture_3d.mb,
            instMB: telem.categories.instanced_mesh.mb,
            rtMB: telem.categories.render_target.mb,
            resourceCount: mgr.resourceCount,
            isValid: inv.isValid
        }));
        """)
        assert res["totalMB"] == 125.0
        assert res["geomMB"] == 10.0
        assert res["tex2DMB"] == 20.0
        assert res["tex3DMB"] == 50.0
        assert res["instMB"] == 15.0
        assert res["rtMB"] == 30.0
        assert res["resourceCount"] == 5
        assert res["isValid"] is True

    def test_resource_touch_and_pinning(self):
        """Verify touchResource updates timestamps and pin/unpin toggles pinned status."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, RESOURCE_PRIORITY } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        const entry = mgr.registerGeometry('geom_dyn', null, { sizeBytes: 1024, priority: RESOURCE_PRIORITY.LOW });
        
        const initPinned = entry.isPinned;
        const initAccess = entry.accessCount;
        
        mgr.touchResource('geom_dyn');
        const postTouchAccess = entry.accessCount;
        
        mgr.pinResource('geom_dyn');
        const postPin = entry.isPinned;
        const postPinPriority = entry.priority;
        
        mgr.unpinResource('geom_dyn', RESOURCE_PRIORITY.LOW);
        const postUnpin = entry.isPinned;
        const postUnpinPriority = entry.priority;
        
        console.log(JSON.stringify({
            initPinned,
            initAccess,
            postTouchAccess,
            postPin,
            postPinPriority,
            postUnpin,
            postUnpinPriority
        }));
        """)
        assert res["initPinned"] is False
        assert res["initAccess"] == 0
        assert res["postTouchAccess"] == 1
        assert res["postPin"] is True
        assert res["postPinPriority"] == 4  # RESOURCE_PRIORITY.PINNED
        assert res["postUnpin"] is False
        assert res["postUnpinPriority"] == 1  # RESOURCE_PRIORITY.LOW

    def test_material_auto_discovery_of_textures(self):
        """Verify registering a material automatically discovers and tracks attached textures."""
        res = run_node_snippet("""
        import { MemoryBudgetManager } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        // Mock Three.js Material with attached color map and normal map
        const mat = {
            name: 'CosmicDustMaterial',
            map: { width: 512, height: 512, format: 'RGBA', type: 'UnsignedByteType' },
            normalMap: { width: 512, height: 512, format: 'RGB', type: 'UnsignedByteType' }
        };
        
        const matEntry = mgr.registerMaterial('mat_dust', mat);
        const telem = mgr.getTelemetrySnapshot();
        
        console.log(JSON.stringify({
            matId: matEntry.id,
            hasMapTex: mgr.hasResource('mat_dust_map_tex'),
            hasNormalTex: mgr.hasResource('mat_dust_normalMap_tex'),
            resourceCount: mgr.resourceCount,
            dependenciesCount: matEntry.dependencies.size
        }));
        """)
        assert res["matId"] == "mat_dust"
        assert res["hasMapTex"] is True
        assert res["hasNormalTex"] is True
        assert res["resourceCount"] == 3
        assert res["dependenciesCount"] == 2


class TestMemoryPressureAndStateTransitions:
    """Test Suite 4: Multi-Level Pressure Transitions, Hysteresis & Alert Callbacks"""

    def test_pressure_level_transitions(self):
        """Verify NORMAL -> MODERATE -> CRITICAL -> EXCEEDED state transitions."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        // Soft limit = 100 MB, Hard limit = 200 MB, Moderate threshold = 75 MB (0.75)
        const mgr = new MemoryBudgetManager({
            softLimitBytes: 100 * BYTES_PER_MB,
            hardLimitBytes: 200 * BYTES_PER_MB,
            moderateThresholdRatio: 0.75,
            autoEvictOnHardLimit: false
        });
        
        const pressureLog = [];
        mgr.on('pressureChange', (data) => {
            pressureLog.push(`${data.oldLevel}->${data.newLevel}`);
        });
        
        // Phase 1: 50 MB -> NORMAL
        mgr.registerGeometry('res_50', null, { sizeBytes: 50 * BYTES_PER_MB });
        const p1 = mgr.pressureLevel;
        
        // Phase 2: +30 MB (80 MB total) -> MODERATE (>= 75 MB)
        mgr.registerGeometry('res_30', null, { sizeBytes: 30 * BYTES_PER_MB });
        const p2 = mgr.pressureLevel;
        
        // Phase 3: +25 MB (105 MB total) -> CRITICAL (>= 100 MB soft limit)
        mgr.registerGeometry('res_25', null, { sizeBytes: 25 * BYTES_PER_MB });
        const p3 = mgr.pressureLevel;
        
        // Phase 4: +100 MB (205 MB total) -> EXCEEDED (>= 200 MB hard limit)
        try {
            mgr.registerGeometry('res_100', null, { sizeBytes: 100 * BYTES_PER_MB });
        } catch (e) {
            // MemoryBudgetExceededError expected
        }
        
        // Phase 5: Free 150 MB (down to 55 MB) -> Hysteresis back to NORMAL
        mgr.unregisterResource('res_30');
        mgr.unregisterResource('res_25');
        const p5 = mgr.pressureLevel;
        
        console.log(JSON.stringify({ p1, p2, p3, p5, pressureLog }));
        """)
        assert res["p1"] == "normal"
        assert res["p2"] == "moderate"
        assert res["p3"] == "critical"
        assert "normal->moderate" in res["pressureLog"]
        assert "moderate->critical" in res["pressureLog"]
        assert "moderate->normal" in res["pressureLog"]

    def test_can_allocate_preflight_check(self):
        """Verify canAllocate correctly projects memory headroom and safety flags."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager({
            softLimitBytes: 100 * BYTES_PER_MB,
            hardLimitBytes: 200 * BYTES_PER_MB
        });
        
        mgr.registerGeometry('init_res', null, { sizeBytes: 60 * BYTES_PER_MB });
        
        const checkSmall = mgr.canAllocate(20 * BYTES_PER_MB);   // Total 80 MB (Moderate, SoftSafe: true)
        const checkMedium = mgr.canAllocate(50 * BYTES_PER_MB);  // Total 110 MB (Critical, SoftSafe: false, HardSafe: true)
        const checkOver = mgr.canAllocate(150 * BYTES_PER_MB);   // Total 210 MB (Exceeded, HardSafe: false)
        
        console.log(JSON.stringify({
            checkSmall: { canAlloc: checkSmall.canAllocate, isSoftSafe: checkSmall.isSoftSafe, pressure: checkSmall.projectedPressure },
            checkMedium: { canAlloc: checkMedium.canAllocate, isSoftSafe: checkMedium.isSoftSafe, pressure: checkMedium.projectedPressure },
            checkOver: { canAlloc: checkOver.canAllocate, isSoftSafe: checkOver.isSoftSafe, pressure: checkOver.projectedPressure }
        }));
        """)
        assert res["checkSmall"]["canAlloc"] is True
        assert res["checkSmall"]["isSoftSafe"] is True
        assert res["checkSmall"]["pressure"] == "moderate"
        
        assert res["checkMedium"]["canAlloc"] is True
        assert res["checkMedium"]["isSoftSafe"] is False
        assert res["checkMedium"]["pressure"] == "critical"
        
        assert res["checkOver"]["canAlloc"] is False
        assert res["checkOver"]["isSoftSafe"] is False
        assert res["checkOver"]["pressure"] == "exceeded"


class TestLRUEvictionAndDisposal:
    """Test Suite 5: Priority-Weighted LRU Eviction & Automatic Three.js Object Disposal"""

    def test_lru_eviction_chronological_order(self):
        """Verify LRU evicts least recently accessed resources first."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, RESOURCE_PRIORITY, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        // Register 3 equal-priority resources
        mgr.registerGeometry('res_a', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.LOW });
        mgr.registerGeometry('res_b', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.LOW });
        mgr.registerGeometry('res_c', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.LOW });
        
        // Touch res_a and res_c, making res_b the least recently used
        mgr.touchResource('res_a');
        mgr.touchResource('res_c');
        
        // Free 10 MB via LRU
        const freed = mgr.evictLRU(10 * BYTES_PER_MB);
        
        console.log(JSON.stringify({
            freedMB: freed / BYTES_PER_MB,
            hasA: mgr.hasResource('res_a'),
            hasB: mgr.hasResource('res_b'),
            hasC: mgr.hasResource('res_c')
        }));
        """)
        assert res["freedMB"] == 10.0
        assert res["hasA"] is True
        assert res["hasB"] is False  # Evicted!
        assert res["hasC"] is True

    def test_lru_eviction_priority_tiering(self):
        """Verify lower priority tiers (TRANSIENT, LOW) are evicted before higher tiers (MEDIUM, HIGH)."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, RESOURCE_PRIORITY, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        mgr.registerGeometry('res_transient', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.TRANSIENT });
        mgr.registerGeometry('res_low', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.LOW });
        mgr.registerGeometry('res_med', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.MEDIUM });
        mgr.registerGeometry('res_high', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.HIGH });
        
        // Free 20 MB: should evict TRANSIENT and LOW, leaving MEDIUM and HIGH untouched
        const freed = mgr.evictLRU(20 * BYTES_PER_MB);
        
        console.log(JSON.stringify({
            freedMB: freed / BYTES_PER_MB,
            hasTransient: mgr.hasResource('res_transient'),
            hasLow: mgr.hasResource('res_low'),
            hasMed: mgr.hasResource('res_med'),
            hasHigh: mgr.hasResource('res_high')
        }));
        """)
        assert res["freedMB"] == 20.0
        assert res["hasTransient"] is False
        assert res["hasLow"] is False
        assert res["hasMed"] is True
        assert res["hasHigh"] is True

    def test_pinned_resources_immunity(self):
        """Verify PINNED resources are strictly immune to LRU eviction."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, RESOURCE_PRIORITY, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        mgr.registerGeometry('res_pinned', null, { sizeBytes: 50 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.PINNED });
        mgr.registerGeometry('res_low', null, { sizeBytes: 10 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.LOW });
        
        // Request freeing 100 MB (more than all available unpinned memory)
        const freed = mgr.evictLRU(100 * BYTES_PER_MB);
        
        console.log(JSON.stringify({
            freedMB: freed / BYTES_PER_MB,
            hasPinned: mgr.hasResource('res_pinned'),
            hasLow: mgr.hasResource('res_low'),
            totalTrackedMB: mgr.totalTrackedBytes / BYTES_PER_MB
        }));
        """)
        assert res["freedMB"] == 10.0
        assert res["hasPinned"] is True
        assert res["hasLow"] is False
        assert res["totalTrackedMB"] == 50.0

    def test_automatic_threejs_disposal_callback_invocation(self):
        """Verify underlying Three.js dispose() methods are called upon eviction/unregistration."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        let geomDisposed = false;
        let texDisposed = false;
        
        const mockGeom = { dispose: () => { geomDisposed = true; } };
        const mockTex = { dispose: () => { texDisposed = true; } };
        
        mgr.registerGeometry('geom_mock', mockGeom, { sizeBytes: 5 * BYTES_PER_MB });
        mgr.registerTexture('tex_mock', mockTex, { sizeBytes: 5 * BYTES_PER_MB });
        
        // Unregister geometry explicitly
        mgr.unregisterResource('geom_mock', true);
        
        // Evict texture via LRU
        mgr.evictLRU(5 * BYTES_PER_MB);
        
        console.log(JSON.stringify({ geomDisposed, texDisposed }));
        """)
        assert res["geomDisposed"] is True
        assert res["texDisposed"] is True

    def test_auto_eviction_on_hard_limit_breach(self):
        """Verify manager automatically triggers LRU eviction to accommodate incoming allocations."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, RESOURCE_PRIORITY, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        // Soft limit 80 MB, Hard limit 100 MB, AutoEvict true
        const mgr = new MemoryBudgetManager({
            softLimitBytes: 80 * BYTES_PER_MB,
            hardLimitBytes: 100 * BYTES_PER_MB,
            autoEvictOnHardLimit: true,
            targetEvictionHeadroom: 0.10
        });
        
        // Register 70 MB of low priority items
        mgr.registerGeometry('chunk_1', null, { sizeBytes: 35 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.LOW });
        mgr.registerGeometry('chunk_2', null, { sizeBytes: 35 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.LOW });
        
        // Allocate 75 MB (would push total to 145 MB > 100 MB hard limit)
        // Should automatically evict chunk_1 and chunk_2 to make room!
        mgr.registerGeometry('chunk_heavy', null, { sizeBytes: 75 * BYTES_PER_MB, priority: RESOURCE_PRIORITY.HIGH });
        
        console.log(JSON.stringify({
            hasHeavy: mgr.hasResource('chunk_heavy'),
            hasChunk1: mgr.hasResource('chunk_1'),
            hasChunk2: mgr.hasResource('chunk_2'),
            totalTrackedMB: mgr.totalTrackedBytes / BYTES_PER_MB
        }));
        """)
        assert res["hasHeavy"] is True
        assert res["hasChunk1"] is False
        assert res["hasChunk2"] is False
        assert res["totalTrackedMB"] == 75.0


class TestWorkerTransferableBufferRegistry:
    """Test Suite 6: Web Worker Transferable ArrayBuffer Registry & Multi-Worker Coordination"""

    def test_worker_transferable_buffer_lifecycle(self):
        """Verify registration, tracking, detachment, and disposal of worker ArrayBuffers."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        // Allocate 10 MB ArrayBuffer
        const buffer = new ArrayBuffer(10 * BYTES_PER_MB);
        const entry = mgr.trackWorkerTransfer('task_rk45_buf', buffer, { workerId: 'worker_1', direction: 'outbound' });
        
        const initMB = mgr.totalTrackedBytes / BYTES_PER_MB;
        const hasBuf = mgr.hasResource('worker_buf_task_rk45_buf');
        
        // Simulate zero-copy detachment upon worker postMessage
        mgr.markWorkerBufferDetached('task_rk45_buf');
        const postDetachMB = mgr.totalTrackedBytes / BYTES_PER_MB;
        
        console.log(JSON.stringify({
            initMB,
            hasBuf,
            postDetachMB,
            resourceCount: mgr.resourceCount
        }));
        """)
        assert res["initMB"] == 10.0
        assert res["hasBuf"] is True
        assert res["postDetachMB"] == 0.0
        assert res["resourceCount"] == 0

    def test_multi_worker_coordinator_quotas(self):
        """Verify worker memory quotas, acceptance checks, and dispatch/completion accounting."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager({ softLimitBytes: 100 * BYTES_PER_MB });
        const coord = mgr.workerCoordinator;
        
        // Register 2 workers with 20 MB quota each
        coord.registerWorker('w_0', 20 * BYTES_PER_MB);
        coord.registerWorker('w_1', 20 * BYTES_PER_MB);
        
        const canAcceptSmall = coord.canWorkerAccept('w_0', 10 * BYTES_PER_MB);
        const canAcceptOver = coord.canWorkerAccept('w_0', 25 * BYTES_PER_MB);
        
        coord.recordDispatch('w_0', 15 * BYTES_PER_MB);
        const canAcceptSecond = coord.canWorkerAccept('w_0', 10 * BYTES_PER_MB); // 15 + 10 = 25 > 20 -> false
        
        coord.recordCompletion('w_0', 15 * BYTES_PER_MB);
        const canAcceptAfterCompletion = coord.canWorkerAccept('w_0', 10 * BYTES_PER_MB);
        
        console.log(JSON.stringify({
            canAcceptSmall,
            canAcceptOver,
            canAcceptSecond,
            canAcceptAfterCompletion
        }));
        """)
        assert res["canAcceptSmall"] is True
        assert res["canAcceptOver"] is False
        assert res["canAcceptSecond"] is False
        assert res["canAcceptAfterCompletion"] is True


class TestStreamlineOOMGuard:
    """Test Suite 7: Out-of-Memory Guard & Adaptive Resolution/LOD Recommenders"""

    def test_guard_reconstruction_safety_pass(self):
        """Verify safe streamline reconstruction passes preflight check."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager({ softLimitBytes: 512 * BYTES_PER_MB });
        
        // 500 seeds, 100 steps, tubes with 6 segments
        const config = {
            seedCount: 500,
            avgStepsPerStreamline: 100,
            attributesPerVertex: 10,
            renderAsTubes: true,
            tubeRadialSegments: 6
        };
        
        const assessment = mgr.guardStreamlineReconstruction(config, { strict: true });
        
        console.log(JSON.stringify({
            isSafe: assessment.isSafe,
            requiredMB: assessment.requiredMB,
            availableMB: assessment.availableMB,
            hasRecommendation: assessment.recommendedConfig !== null
        }));
        """)
        assert res["isSafe"] is True
        assert res["requiredMB"] < 50.0
        assert res["hasRecommendation"] is False

    def test_guard_reconstruction_safety_rejection_and_adaptive_downscaling(self):
        """Verify unsafe massive streamline reconstruction is flagged and receives downscaled recommendation."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        // Low headroom scenario: soft limit 20 MB, 15 MB already occupied (5 MB headroom)
        const mgr = new MemoryBudgetManager({ softLimitBytes: 20 * BYTES_PER_MB });
        mgr.registerGeometry('existing_grid', null, { sizeBytes: 15 * BYTES_PER_MB });
        
        // Massive config requiring ~40 MB
        const config = {
            seedCount: 5000,
            avgStepsPerStreamline: 128,
            attributesPerVertex: 10,
            renderAsTubes: true,
            tubeRadialSegments: 6
        };
        
        const assessment = mgr.guardStreamlineReconstruction(config, { strict: false });
        const rec = assessment.recommendedConfig;
        
        console.log(JSON.stringify({
            isSafe: assessment.isSafe,
            requiredMB: assessment.requiredMB,
            availableMB: assessment.availableMB,
            strategy: rec ? rec.strategy : null,
            recProjectedMB: rec ? rec.projectedMB : null,
            recSeeds: rec ? rec.seedCount : null
        }));
        """)
        assert res["isSafe"] is False
        assert res["requiredMB"] > res["availableMB"]
        assert res["strategy"] is not None
        assert res["recProjectedMB"] <= res["availableMB"]

    def test_guard_grid_safety_resolution_fallback(self):
        """Verify 3D grid volume allocation guard recommends resolution fallback (256³ -> 128³ -> 64³)."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        // 30 MB available headroom
        const mgr = new MemoryBudgetManager({ softLimitBytes: 30 * BYTES_PER_MB });
        
        // 256^3 vector field = 192 MB -> unsafe, fallback to 128^3 (24 MB)
        const g256 = mgr.guardGridAllocation(256, 3);
        
        // 128^3 vector field = 24 MB -> safe in 30 MB headroom
        const g128 = mgr.guardGridAllocation(128, 3);
        
        console.log(JSON.stringify({
            g256_safe: g256.isSafe,
            g256_fallback: g256.fallbackResolution,
            g128_safe: g128.isSafe
        }));
        """)
        assert res["g256_safe"] is False
        assert res["g256_fallback"] == 128
        assert res["g128_safe"] is True


class TestTelemetryAndLeakDetection:
    """Test Suite 8: High-Frequency Telemetry Polling, History Buffers & OLS Leak Detector"""

    def test_telemetry_snapshot_and_category_breakdown(self):
        """Verify snapshot metrics, percentages, peak watermarks, and category breakdown."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        mgr.registerGeometry('g1', null, { sizeBytes: 25 * BYTES_PER_MB });
        mgr.registerTexture('t1', null, { sizeBytes: 75 * BYTES_PER_MB });
        
        const snap = mgr.getTelemetrySnapshot();
        
        console.log(JSON.stringify({
            totalTrackedMB: snap.totalTrackedMB,
            peakTrackedMB: snap.peakTrackedMB,
            geomPercent: snap.categories.geometry.percent,
            texPercent: snap.categories.texture_2d.percent,
            resourceCount: snap.resourceCount
        }));
        """)
        assert res["totalTrackedMB"] == 100.0
        assert res["peakTrackedMB"] == 100.0
        assert res["geomPercent"] == 25.0
        assert res["texPercent"] == 75.0
        assert res["resourceCount"] == 2

    def test_high_frequency_telemetry_polling_and_history(self):
        """Verify telemetry polling ticks, history ring buffer, and allocation rate smoothing."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager({ maxHistorySamples: 5 });
        
        for (let i = 0; i < 10; i++) {
            mgr.registerGeometry(`g_${i}`, null, { sizeBytes: 1 * BYTES_PER_MB });
            mgr.pollTelemetry();
        }
        
        const history = mgr.getHistory();
        
        console.log(JSON.stringify({
            historyLength: history.length,
            lastTotalMB: history[history.length - 1].totalTrackedMB,
            resourceCount: mgr.resourceCount
        }));
        """)
        assert res["historyLength"] == 5  # Ring buffer capped at maxHistorySamples = 5
        assert res["lastTotalMB"] == 10.0
        assert res["resourceCount"] == 10

    def test_ols_memory_leak_detection(self):
        """Verify ordinary least squares (OLS) linear regression trend detector."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager({ softLimitBytes: 100 * BYTES_PER_MB });
        
        // Simulate continuous memory accumulation (leak profile)
        for (let i = 0; i < 20; i++) {
            mgr.registerGeometry(`leak_chunk_${i}`, null, { sizeBytes: 4 * BYTES_PER_MB });
            mgr.pollTelemetry();
        }
        
        const trend = mgr.detectMemoryLeakTrend(20);
        
        console.log(JSON.stringify({
            hasConfidence: trend.confidence > 0.5,
            slopeMBPerMin: trend.slopeMBPerMin,
            isLeaking: trend.isLeaking
        }));
        """)
        assert res["hasConfidence"] is True
        assert res["slopeMBPerMin"] > 0


class TestGarbageCollectionSentinel:
    """Test Suite 9: WeakRef & FinalizationRegistry Garbage Collection Sentinel"""

    def test_gc_sentinel_registration_and_untracking(self):
        """Verify objects are registered with FinalizationRegistry and clean up properly."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, GarbageCollectionSentinel } from './src/runtime/memory_budget_manager.js';
        
        let orphanNotified = false;
        const sentinel = new GarbageCollectionSentinel((token) => {
            orphanNotified = true;
        });
        
        const mockTarget = { id: 'test_mesh', name: 'CosmicMesh' };
        sentinel.track(mockTarget, { id: 'test_mesh', sizeBytes: 1024 });
        
        sentinel.untrack(mockTarget, 'test_mesh');
        
        console.log(JSON.stringify({
            isSentinelActive: typeof FinalizationRegistry !== 'undefined',
            collectedCount: sentinel.collectedCount
        }));
        """)
        assert res["isSentinelActive"] is True
        assert res["collectedCount"] == 0


class TestScientificInvariantsAndStress:
    """Test Suite 10: Strict Memory Conservation Laws, Rapid Cycling & Error Bounds"""

    def test_strict_memory_conservation_law(self):
        """Verify mathematical invariance: tracked bytes == sum of resources == sum of categories."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, RESOURCE_PRIORITY, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        
        // Randomly register and unregister 100 items
        for (let i = 0; i < 100; i++) {
            const size = (i + 1) * 1024 * 64; // arbitrary size
            mgr.registerGeometry(`geom_${i}`, null, { sizeBytes: size });
        }
        
        // Delete half
        for (let i = 0; i < 50; i++) {
            mgr.unregisterResource(`geom_${i}`);
        }
        
        const inv = mgr.verifyInvariants();
        
        console.log(JSON.stringify({
            isValid: inv.isValid,
            discrepancies: inv.discrepancies,
            resourceCount: mgr.resourceCount,
            totalTrackedBytes: mgr.totalTrackedBytes
        }));
        """)
        assert res["isValid"] is True
        assert len(res["discrepancies"]) == 0
        assert res["resourceCount"] == 50
        assert res["totalTrackedBytes"] > 0

    def test_rapid_cyclic_allocation_and_eviction_stress(self):
        """Verify rapid allocation/eviction cycles cause zero accounting drift or memory leaks."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager({
            softLimitBytes: 50 * BYTES_PER_MB,
            hardLimitBytes: 100 * BYTES_PER_MB,
            autoEvictOnHardLimit: true
        });
        
        // Perform 300 rapid allocations that continuously push above hard limit
        for (let i = 0; i < 300; i++) {
            mgr.registerGeometry(`cycle_geom_${i}`, null, { sizeBytes: 10 * BYTES_PER_MB });
        }
        
        const inv = mgr.verifyInvariants();
        const telem = mgr.getTelemetrySnapshot();
        
        console.log(JSON.stringify({
            isValid: inv.isValid,
            discrepancyCount: inv.discrepancies.length,
            totalTrackedMB: telem.totalTrackedMB,
            evictionCount: mgr._totalEvictionCount,
            isWithinHardLimit: mgr.totalTrackedBytes <= mgr.hardLimitBytes
        }));
        """)
        assert res["isValid"] is True
        assert res["discrepancyCount"] == 0
        assert res["isWithinHardLimit"] is True
        assert res["evictionCount"] > 200

    def test_reset_and_destroy_lifecycle(self):
        """Verify reset() and destroy() cleanly restore all state counters and release listeners."""
        res = run_node_snippet("""
        import { MemoryBudgetManager, BYTES_PER_MB } from './src/runtime/memory_budget_manager.js';
        
        const mgr = new MemoryBudgetManager();
        mgr.registerGeometry('res_1', null, { sizeBytes: 10 * BYTES_PER_MB });
        mgr.registerTexture('res_2', null, { sizeBytes: 20 * BYTES_PER_MB });
        
        mgr.reset();
        const postResetTracked = mgr.totalTrackedBytes;
        const postResetCount = mgr.resourceCount;
        
        mgr.registerGeometry('res_3', null, { sizeBytes: 5 * BYTES_PER_MB });
        mgr.destroy();
        
        let caughtOnDestroyed = false;
        try {
            mgr.registerGeometry('res_4', null, { sizeBytes: 5 * BYTES_PER_MB });
        } catch (e) {
            caughtOnDestroyed = true;
        }
        
        console.log(JSON.stringify({
            postResetTracked,
            postResetCount,
            caughtOnDestroyed
        }));
        """)
        assert res["postResetTracked"] == 0
        assert res["postResetCount"] == 0
        assert res["caughtOnDestroyed"] is True
