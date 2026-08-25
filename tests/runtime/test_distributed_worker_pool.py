# -*- coding: utf-8 -*-
"""
Automated Pytest Suite for Distributed Worker Pool and Zero-Copy Task Scheduling.
"""

import pytest
from tests.utils import run_node_snippet

class TestDistributedWorkerPool:
    def test_worker_pool_initialization(self):
        res = run_node_snippet("""
        import { DistributedWorkerPool } from './src/runtime/distributed_worker_pool.js';
        const pool = new DistributedWorkerPool({ maxWorkers: 8 });
        pool.initialize();
        const telem = pool.getTelemetry();
        console.log(JSON.stringify({
            isInit: telem.isInitialized,
            totalWorkers: telem.totalWorkers,
            idleCount: telem.idleWorkerCount
        }));
        """)
        assert res["isInit"] is True
        assert res["totalWorkers"] == 8
        assert res["idleCount"] == 8

    def test_streamline_batch_dispatch(self):
        res = run_node_snippet("""
        import { DistributedWorkerPool } from './src/runtime/distributed_worker_pool.js';
        const pool = new DistributedWorkerPool({ maxWorkers: 4 });
        const seeds = [
            { id: 1, sgx: 10, sgy: 20, sgz: 30 },
            { id: 2, sgx: -10, sgy: 40, sgz: -20 },
            { id: 3, sgx: 50, sgy: -30, sgz: 10 }
        ];
        
        const results = await pool.dispatchStreamlineBatch(seeds, {}, { chunkSize: 2 });
        const telem = pool.getTelemetry();
        console.log(JSON.stringify({
            resultCount: results.length,
            completedTasks: telem.completedTasks,
            firstResultId: results[0].id
        }));
        """)
        assert res["resultCount"] == 3
        assert res["completedTasks"] == 3
        assert res["firstResultId"] == 1
