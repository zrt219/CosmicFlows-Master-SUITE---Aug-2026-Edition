self.importScripts('./trilinear.js');

self.onmessage = function(e) {
  if (e.data.type === 'analyze') {
    try {
      const { gridBuffer, N, halfExtent, streamlineResults } = e.data;
      const velocityGrid = new Float32Array(gridBuffer);

      // --- Helper Functions ---
      function postProgress(percent, status) {
        self.postMessage({ type: 'progress', percent, status });
      }

      function getV(ix, iy, iz) {
        // Clamp indices
        ix = Math.max(0, Math.min(N - 1, ix));
        iy = Math.max(0, Math.min(N - 1, iy));
        iz = Math.max(0, Math.min(N - 1, iz));
        const idx = (iz * N * N + iy * N + ix) * 3;
        return [velocityGrid[idx], velocityGrid[idx + 1], velocityGrid[idx + 2]];
      }

      // Linear interpolation inside cell for continuous velocity
      function sampleVelocity(x, y, z) {
        // Convert to grid coordinates [0, N-1]
        const gx = (x + halfExtent) / (2 * halfExtent) * (N - 1);
        const gy = (y + halfExtent) / (2 * halfExtent) * (N - 1);
        const gz = (z + halfExtent) / (2 * halfExtent) * (N - 1);

        const ix0 = Math.floor(gx), ix1 = ix0 + 1;
        const iy0 = Math.floor(gy), iy1 = iy0 + 1;
        const iz0 = Math.floor(gz), iz1 = iz0 + 1;

        const tx = gx - ix0;
        const ty = gy - iy0;
        const tz = gz - iz0;

        const c000 = getV(ix0, iy0, iz0);
        const c100 = getV(ix1, iy0, iz0);
        const c010 = getV(ix0, iy1, iz0);
        const c110 = getV(ix1, iy1, iz0);
        const c001 = getV(ix0, iy0, iz1);
        const c101 = getV(ix1, iy0, iz1);
        const c011 = getV(ix0, iy1, iz1);
        const c111 = getV(ix1, iy1, iz1);

        const v = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
          const c00 = c000[i] * (1 - tx) + c100[i] * tx;
          const c01 = c001[i] * (1 - tx) + c101[i] * tx;
          const c10 = c010[i] * (1 - tx) + c110[i] * tx;
          const c11 = c011[i] * (1 - tx) + c111[i] * tx;
          
          const c0 = c00 * (1 - ty) + c10 * ty;
          const c1 = c01 * (1 - ty) + c11 * ty;
          
          v[i] = c0 * (1 - tz) + c1 * tz;
        }
        return v;
      }
      
      function computeJacobian(x, y, z) {
        const dx = halfExtent / N;
        const vx1 = sampleVelocity(x + dx, y, z);
        const vx0 = sampleVelocity(x - dx, y, z);
        const vy1 = sampleVelocity(x, y + dx, z);
        const vy0 = sampleVelocity(x, y - dx, z);
        const vz1 = sampleVelocity(x, y, z + dx);
        const vz0 = sampleVelocity(x, y, z - dx);
        
        return [
          [(vx1[0]-vx0[0])/(2*dx), (vy1[0]-vy0[0])/(2*dx), (vz1[0]-vz0[0])/(2*dx)],
          [(vx1[1]-vx0[1])/(2*dx), (vy1[1]-vy0[1])/(2*dx), (vz1[1]-vz0[1])/(2*dx)],
          [(vx1[2]-vx0[2])/(2*dx), (vy1[2]-vy0[2])/(2*dx), (vz1[2]-vz0[2])/(2*dx)]
        ];
      }

      function invert3x3(m) {
        const det = m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2]) -
                    m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
                    m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
        if (Math.abs(det) < 1e-8) return null;
        const invDet = 1 / det;
        return [
          [(m[1][1]*m[2][2]-m[2][1]*m[1][2])*invDet, (m[0][2]*m[2][1]-m[0][1]*m[2][2])*invDet, (m[0][1]*m[1][2]-m[0][2]*m[1][1])*invDet],
          [(m[1][2]*m[2][0]-m[1][0]*m[2][2])*invDet, (m[0][0]*m[2][2]-m[0][2]*m[2][0])*invDet, (m[1][0]*m[0][2]-m[0][0]*m[1][2])*invDet],
          [(m[1][0]*m[2][1]-m[2][0]*m[1][1])*invDet, (m[2][0]*m[0][1]-m[0][0]*m[2][1])*invDet, (m[0][0]*m[1][1]-m[1][0]*m[0][1])*invDet]
        ];
      }
      
      // Calculate eigenvalues of a symmetric 3x3 matrix (Strain rate tensor)
      // Using an analytical approach or iterative for 3x3 symmetric
      function calculateEigen(sym) {
        // Simple Jacobi eigenvalue algorithm for 3x3
        let V = [[1,0,0], [0,1,0], [0,0,1]];
        let A = [[sym[0][0], sym[0][1], sym[0][2]],
                 [sym[1][0], sym[1][1], sym[1][2]],
                 [sym[2][0], sym[2][1], sym[2][2]]];
        
        for (let iter = 0; iter < 50; iter++) {
          let maxVal = 0, p = 0, q = 1;
          for (let i = 0; i < 3; i++) {
            for (let j = i+1; j < 3; j++) {
              if (Math.abs(A[i][j]) > maxVal) {
                maxVal = Math.abs(A[i][j]);
                p = i; q = j;
              }
            }
          }
          if (maxVal < 1e-6) break;
          
          let theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
          let t = Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta*theta + 1));
          if (theta === 0) t = 1;
          
          let c = 1 / Math.sqrt(t*t + 1);
          let s = t * c;
          
          // Apply rotation
          let App = c*c*A[p][p] - 2*s*c*A[p][q] + s*s*A[q][q];
          let Aqq = s*s*A[p][p] + 2*s*c*A[p][q] + c*c*A[q][q];
          A[p][q] = A[q][p] = 0;
          A[p][p] = App;
          A[q][q] = Aqq;
          
          for (let i = 0; i < 3; i++) {
            if (i !== p && i !== q) {
              let Aip = A[i][p], Aiq = A[i][q];
              A[i][p] = A[p][i] = c*Aip - s*Aiq;
              A[i][q] = A[q][i] = s*Aip + c*Aiq;
            }
            let Vip = V[i][p], Viq = V[i][q];
            V[i][p] = c*Vip - s*Viq;
            V[i][q] = s*Vip + c*Viq;
          }
        }
        
        return {
          eigenvalues: [A[0][0], A[1][1], A[2][2]],
          eigenvectors: V
        };
      }


      // --- Analysis 1: Watershed Basin Classification (0-50%) ---
      postProgress(0, 'Starting Watershed Basin Classification');
      const M = 64; // 64^3 grid for basin
      const basinGrid = new Int32Array(M * M * M);
      
      const numTerminals = streamlineResults.length;
      const terminalPos = new Float32Array(numTerminals * 3);
      const terminalBasin = new Int32Array(numTerminals);
      
      let maxBasinId = -1;
      
      for (let i = 0; i < numTerminals; i++) {
        let res = streamlineResults[i];
        terminalPos[i*3] = res.terminalPos[0];
        terminalPos[i*3+1] = res.terminalPos[1];
        terminalPos[i*3+2] = res.terminalPos[2];
        terminalBasin[i] = res.basinId;
        if (res.basinId > maxBasinId) maxBasinId = res.basinId;
      }
      
      // Assign nearest terminal to each 64^3 cell
      for (let z = 0; z < M; z++) {
        const pz = (z + 0.5) / M * (2 * halfExtent) - halfExtent;
        for (let y = 0; y < M; y++) {
          const py = (y + 0.5) / M * (2 * halfExtent) - halfExtent;
          for (let x = 0; x < M; x++) {
            const px = (x + 0.5) / M * (2 * halfExtent) - halfExtent;
            
            let minDist = Infinity;
            let bestBasin = -1;
            for (let i = 0; i < numTerminals; i++) {
              let dx = px - terminalPos[i*3];
              let dy = py - terminalPos[i*3+1];
              let dz = pz - terminalPos[i*3+2];
              let dist2 = dx*dx + dy*dy + dz*dz;
              if (dist2 < minDist) {
                minDist = dist2;
                bestBasin = terminalBasin[i];
              }
            }
            basinGrid[z * M * M + y * M + x] = bestBasin;
          }
        }
        if (z % 8 === 0) postProgress((z / M) * 20, 'Assigning nearest terminals...');
      }
      
      postProgress(20, 'Smoothing basin boundaries...');
      // Smooth with 3x3x3 majority filter
      const basinGridSmoothed = new Int32Array(M * M * M);
      for (let z = 0; z < M; z++) {
        for (let y = 0; y < M; y++) {
          for (let x = 0; x < M; x++) {
            let counts = {};
            for (let dz = -1; dz <= 1; dz++) {
              for (let dy = -1; dy <= 1; dy++) {
                for (let dx = -1; dx <= 1; dx++) {
                  let nz = z + dz, ny = y + dy, nx = x + dx;
                  if (nx >= 0 && nx < M && ny >= 0 && ny < M && nz >= 0 && nz < M) {
                    let b = basinGrid[nz * M * M + ny * M + nx];
                    counts[b] = (counts[b] || 0) + 1;
                  }
                }
              }
            }
            let bestB = -1, maxC = -1;
            for (let b in counts) {
              if (counts[b] > maxC) { maxC = counts[b]; bestB = Number(b); }
            }
            basinGridSmoothed[z * M * M + y * M + x] = bestB;
          }
        }
      }
      
      postProgress(30, 'Computing basin statistics...');
      const numBasins = maxBasinId + 1;
      const basinStats = Array(numBasins).fill(null).map(() => ({
        galaxyCount: 0,
        volumeCount: 0,
        sumSpeed: 0,
        sumTerminalPos: [0,0,0],
        terminalCount: 0
      }));
      
      for (let i = 0; i < streamlineResults.length; i++) {
        let res = streamlineResults[i];
        let b = res.basinId;
        if (b >= 0 && b < numBasins) {
          basinStats[b].galaxyCount++;
          let speed = res.initialSpeed || Math.sqrt(res.initialVelocity ? (res.initialVelocity[0]**2 + res.initialVelocity[1]**2 + res.initialVelocity[2]**2) : 0);
          basinStats[b].sumSpeed += speed;
          basinStats[b].sumTerminalPos[0] += res.terminalPos[0];
          basinStats[b].sumTerminalPos[1] += res.terminalPos[1];
          basinStats[b].sumTerminalPos[2] += res.terminalPos[2];
          basinStats[b].terminalCount++;
        }
      }
      
      for (let i = 0; i < M * M * M; i++) {
        let b = basinGridSmoothed[i];
        if (b >= 0 && b < numBasins) basinStats[b].volumeCount++;
      }
      
      const basinSummary = basinStats.map((stat, id) => {
        return {
          basinId: id,
          galaxyCount: stat.galaxyCount,
          volumeFraction: stat.volumeCount / (M * M * M),
          meanInflowSpeed: stat.galaxyCount > 0 ? stat.sumSpeed / stat.galaxyCount : 0,
          attractorPosition: stat.terminalCount > 0 ? [
            stat.sumTerminalPos[0] / stat.terminalCount,
            stat.sumTerminalPos[1] / stat.terminalCount,
            stat.sumTerminalPos[2] / stat.terminalCount
          ] : [0,0,0]
        };
      }).filter(b => b.galaxyCount > 0 || b.volumeFraction > 0);
      
      postProgress(40, 'Identifying basin boundaries...');
      const boundaryVoxels = [];
      const interBasinFluxMap = new Map();
      
      for (let z = 0; z < M; z++) {
        for (let y = 0; y < M; y++) {
          for (let x = 0; x < M; x++) {
            let b = basinGridSmoothed[z * M * M + y * M + x];
            let isBoundary = false;
            for (let dz = -1; dz <= 1; dz++) {
              for (let dy = -1; dy <= 1; dy++) {
                for (let dx = -1; dx <= 1; dx++) {
                  if (dx===0 && dy===0 && dz===0) continue;
                  let nz = z + dz, ny = y + dy, nx = x + dx;
                  if (nx >= 0 && nx < M && ny >= 0 && ny < M && nz >= 0 && nz < M) {
                    let nb = basinGridSmoothed[nz * M * M + ny * M + nx];
                    if (nb !== b && nb !== -1) {
                      isBoundary = true;
                    }
                  }
                }
              }
            }
            if (isBoundary) {
              const px = (x + 0.5) / M * (2 * halfExtent) - halfExtent;
              const py = (y + 0.5) / M * (2 * halfExtent) - halfExtent;
              const pz = (z + 0.5) / M * (2 * halfExtent) - halfExtent;
              boundaryVoxels.push([px, py, pz, b]);
            }
          }
        }
      }
      
      // Calculate flux - simple version based on streamlines
      for (let i = 0; i < streamlineResults.length; i++) {
        let res = streamlineResults[i];
        if (res.path) {
          // If streamline path is available, trace it through basins
          let lastBasin = -1;
          for (let p of res.path) {
            let px = p[0], py = p[1], pz = p[2];
            let ix = Math.floor((px + halfExtent) / (2 * halfExtent) * M);
            let iy = Math.floor((py + halfExtent) / (2 * halfExtent) * M);
            let iz = Math.floor((pz + halfExtent) / (2 * halfExtent) * M);
            ix = Math.max(0, Math.min(M-1, ix));
            iy = Math.max(0, Math.min(M-1, iy));
            iz = Math.max(0, Math.min(M-1, iz));
            
            let b = basinGridSmoothed[iz * M * M + iy * M + ix];
            if (lastBasin !== -1 && b !== lastBasin && b !== -1) {
              let key = lastBasin < b ? `${lastBasin}-${b}` : `${b}-${lastBasin}`;
              interBasinFluxMap.set(key, (interBasinFluxMap.get(key) || 0) + 1);
            }
            lastBasin = b;
          }
        }
      }
      
      const interBasinFlux = Array.from(interBasinFluxMap.entries()).map(([key, count]) => {
        let parts = key.split('-');
        return { basinA: parseInt(parts[0]), basinB: parseInt(parts[1]), fluxCount: count };
      });


      // --- Analysis 2: Morse Theory Critical Point Finder (50-100%) ---
      postProgress(50, 'Starting Morse Theory Critical Point Finder...');
      
      const candidates = [];
      // Search 128^3 grid for zero-crossing candidates
      for (let z = 0; z < N - 1; z++) {
        for (let y = 0; y < N - 1; y++) {
          for (let x = 0; x < N - 1; x++) {
            let minVx = Infinity, maxVx = -Infinity;
            let minVy = Infinity, maxVy = -Infinity;
            let minVz = Infinity, maxVz = -Infinity;
            
            for (let dz = 0; dz <= 1; dz++) {
              for (let dy = 0; dy <= 1; dy++) {
                for (let dx = 0; dx <= 1; dx++) {
                  let v = getV(x+dx, y+dy, z+dz);
                  minVx = Math.min(minVx, v[0]); maxVx = Math.max(maxVx, v[0]);
                  minVy = Math.min(minVy, v[1]); maxVy = Math.max(maxVy, v[1]);
                  minVz = Math.min(minVz, v[2]); maxVz = Math.max(maxVz, v[2]);
                }
              }
            }
            
            if (minVx < 0 && maxVx > 0 && minVy < 0 && maxVy > 0 && minVz < 0 && maxVz > 0) {
              const px = (x + 0.5) / (N - 1) * (2 * halfExtent) - halfExtent;
              const py = (y + 0.5) / (N - 1) * (2 * halfExtent) - halfExtent;
              const pz = (z + 0.5) / (N - 1) * (2 * halfExtent) - halfExtent;
              candidates.push([px, py, pz]);
            }
          }
        }
        if (z % 16 === 0) postProgress(50 + (z / N) * 20, 'Scanning for zero-crossings...');
      }
      
      postProgress(70, `Refining ${candidates.length} candidates via Newton's method...`);
      
      let rawCriticalPoints = [];
      let cpStats = { nAttractors: 0, nRepellers: 0, nSaddleFilament: 0, nSaddleWall: 0 };
      
      for (let i = 0; i < candidates.length; i++) {
        let pos = candidates[i];
        let refinedPos = [...pos];
        let valid = true;
        
        for (let iter = 0; iter < 5; iter++) {
          let v = sampleVelocity(refinedPos[0], refinedPos[1], refinedPos[2]);
          let J = computeJacobian(refinedPos[0], refinedPos[1], refinedPos[2]);
          let Jinv = invert3x3(J);
          
          if (!Jinv) { valid = false; break; }
          
          refinedPos[0] -= (Jinv[0][0]*v[0] + Jinv[0][1]*v[1] + Jinv[0][2]*v[2]);
          refinedPos[1] -= (Jinv[1][0]*v[0] + Jinv[1][1]*v[1] + Jinv[1][2]*v[2]);
          refinedPos[2] -= (Jinv[2][0]*v[0] + Jinv[2][1]*v[1] + Jinv[2][2]*v[2]);
          
          if (Math.abs(refinedPos[0]) > halfExtent || Math.abs(refinedPos[1]) > halfExtent || Math.abs(refinedPos[2]) > halfExtent) {
            valid = false; break;
          }
        }
        
        if (valid) {
          let vFinal = sampleVelocity(refinedPos[0], refinedPos[1], refinedPos[2]);
          let vMag = Math.sqrt(vFinal[0]**2 + vFinal[1]**2 + vFinal[2]**2);
          
          let J = computeJacobian(refinedPos[0], refinedPos[1], refinedPos[2]);
          // Strain rate tensor S_ij = (J_ij + J_ji)/2
          let S = [
            [J[0][0], (J[0][1]+J[1][0])/2, (J[0][2]+J[2][0])/2],
            [(J[1][0]+J[0][1])/2, J[1][1], (J[1][2]+J[2][1])/2],
            [(J[2][0]+J[0][2])/2, (J[2][1]+J[1][2])/2, J[2][2]]
          ];
          
          let eigen = calculateEigen(S);
          let e = eigen.eigenvalues;
          
          let numNeg = 0;
          if (e[0] < 0) numNeg++;
          if (e[1] < 0) numNeg++;
          if (e[2] < 0) numNeg++;
          
          let type = '';
          if (numNeg === 3) type = 'attractor';
          else if (numNeg === 0) type = 'repeller';
          else if (numNeg === 2) type = 'saddle-filament';
          else if (numNeg === 1) type = 'saddle-wall';
          
          rawCriticalPoints.push({
            position: refinedPos,
            type: type,
            eigenvalues: e,
            eigenvectors: eigen.eigenvectors,
            vMag: vMag
          });
        }
        
        if (i % 1000 === 0) postProgress(70 + (i / candidates.length) * 10, 'Refining critical points...');
      }
      
      postProgress(80, 'De-duplicating critical points...');
      // De-duplicate within 500 km/s (assuming spatial units relate to this, or simply 500 distance units)
      const DEDUP_DIST = 500;
      let criticalPoints = [];
      
      // Sort by velocity magnitude (ascending) to keep the one with smallest |v|
      rawCriticalPoints.sort((a, b) => a.vMag - b.vMag);
      
      for (let i = 0; i < rawCriticalPoints.length; i++) {
        let cp1 = rawCriticalPoints[i];
        let isDup = false;
        for (let j = 0; j < criticalPoints.length; j++) {
          let cp2 = criticalPoints[j];
          let dx = cp1.position[0] - cp2.position[0];
          let dy = cp1.position[1] - cp2.position[1];
          let dz = cp1.position[2] - cp2.position[2];
          let dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
          if (dist < DEDUP_DIST) {
            isDup = true;
            break;
          }
        }
        if (!isDup) {
          criticalPoints.push(cp1);
          if (cp1.type === 'attractor') cpStats.nAttractors++;
          if (cp1.type === 'repeller') cpStats.nRepellers++;
          if (cp1.type === 'saddle-filament') cpStats.nSaddleFilament++;
          if (cp1.type === 'saddle-wall') cpStats.nSaddleWall++;
        }
      }
      
      postProgress(90, 'Constructing skeleton graph...');
      let attractors = criticalPoints.filter(cp => cp.type === 'attractor');
      let saddleFilaments = criticalPoints.filter(cp => cp.type === 'saddle-filament');
      
      let skeletonEdges = [];
      
      for (let i = 0; i < saddleFilaments.length; i++) {
        let saddle = saddleFilaments[i];
        // Connect to 2 nearest attractors along the positive eigenvector direction
        // The positive eigenvalue is the one > 0
        let posIdx = saddle.eigenvalues.findIndex(val => val > 0);
        if (posIdx !== -1) {
          let dir = [saddle.eigenvectors[0][posIdx], saddle.eigenvectors[1][posIdx], saddle.eigenvectors[2][posIdx]];
          
          // Find closest attractor in +dir and -dir
          let bestPosAttractor = null, bestNegAttractor = null;
          let minDistPos = Infinity, minDistNeg = Infinity;
          
          for (let j = 0; j < attractors.length; j++) {
            let att = attractors[j];
            let dx = att.position[0] - saddle.position[0];
            let dy = att.position[1] - saddle.position[1];
            let dz = att.position[2] - saddle.position[2];
            let dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
            if (dist < 1e-5) continue;
            
            // Dot product to check direction
            let dot = (dx*dir[0] + dy*dir[1] + dz*dir[2]) / dist;
            
            if (dot > 0.5) { // roughly in positive direction
              if (dist < minDistPos) { minDistPos = dist; bestPosAttractor = att; }
            } else if (dot < -0.5) { // roughly in negative direction
              if (dist < minDistNeg) { minDistNeg = dist; bestNegAttractor = att; }
            }
          }
          
          if (bestPosAttractor) {
            skeletonEdges.push({ source: saddle.position, target: bestPosAttractor.position });
          }
          if (bestNegAttractor) {
            skeletonEdges.push({ source: saddle.position, target: bestNegAttractor.position });
          }
        }
      }
      
      postProgress(100, 'Analysis complete');
      
      self.postMessage({
        type: 'complete',
        basins: {
          summary: basinSummary,
          boundaryVoxels: boundaryVoxels,
          interBasinFlux: interBasinFlux
        },
        criticalPoints: {
          points: criticalPoints.map(cp => ({
            position: cp.position,
            type: cp.type,
            eigenvalues: cp.eigenvalues,
            eigenvectors: cp.eigenvectors
          })),
          skeletonEdges: skeletonEdges
        },
        stats: cpStats
      });
      
    } catch (error) {
      self.postMessage({ type: 'error', error: error.message, stack: error.stack });
    }
  }
};
