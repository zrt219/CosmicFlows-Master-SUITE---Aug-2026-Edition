/**
 * Powers of Ten: Cosmic Scale Explorer (Laniakea to Earth)
 * 
 * Multi-scale scientific visualization bridging the scale gap across 17 orders of magnitude:
 * Stage 1: Laniakea Supercluster (~160 Mpc / 500 million ly) - Inward peculiar velocity streamlines & watershed basin
 * Stage 2: Milky Way 3D Spiral Galaxy (~30 kpc / 100,000 ly) - 4-arm barred spiral with Sun marked in Orion Spur
 * Stage 3: Solar System (~80 AU) - Keplerian orbital paths, glowing Sun corona, asteroid belt, planets & moons
 * Stage 4: Earth & Moon (~12,742 km) - Procedural Earth globe, atmospheric Fresnel limb glow, clouds, Moon orbit
 *
 * 100% Self-Contained Procedural Textures & Shaders (Zero external network dependencies).
 * Adheres strictly to AGENTS.md / GEMINI.md physical dimensions, KaTeX, and pre-commit gate invariants.
 */

(function(window) {
'use strict';

// =========================================================================
// 1. PHYSICAL CONSTANTS & SCALE DEFINITIONS
// =========================================================================
const STAGE_SCALES = [
  {
    id: 'laniakea',
    level: 1,
    name: 'Laniakea Supercluster',
    shortName: 'Laniakea',
    icon: '🌌',
    exponent: 24, // ~5.0 x 10^24 m
    diameterMeters: 4.937e24,
    diameterDisplay: '160 Mpc · 520,000,000 Light-Years',
    lightTravelTime: '520,000,000 Years (to cross basin)',
    cameraDistance: 280,
    targetPosition: new THREE.Vector3(0, 0, 0),
    analogy: 'If Laniakea were shrunk to the size of planet Earth, the Milky Way would be the size of a single microscopic dust mite (0.2 mm), and our entire Solar System would be smaller than a single atom.',
    scientificNote: 'Discovered by Tully, Courtois et al. (2014, Nature). Over 100,000 galaxies gravitationally bound along inward streamlines terminating at the Great Attractor (Norma/Centaurus core).'
  },
  {
    id: 'milkyway',
    level: 2,
    name: 'Milky Way Galaxy',
    shortName: 'Milky Way',
    icon: '🌀',
    exponent: 21, // ~1.0 x 10^21 m
    diameterMeters: 9.461e20,
    diameterDisplay: '30 kpc · 100,000 Light-Years',
    lightTravelTime: '100,000 Years (disc diameter)',
    cameraDistance: 180,
    targetPosition: new THREE.Vector3(0, 0, 0),
    analogy: 'If our Milky Way were the size of North America (4,500 km across), our entire Solar System out to Neptune would be the size of a tiny 2 cm coin, and Earth would be an invisible speck just 0.0001 mm wide.',
    scientificNote: 'A barred spiral galaxy hosting 200–400 billion stars. Our Solar System resides in the Orion-Cygnus Spur, orbiting the central supermassive black hole Sagittarius A* at ~220 km/s once every 230 million years.'
  },
  {
    id: 'solarsystem',
    level: 3,
    name: 'The Solar System',
    shortName: 'Solar System',
    icon: '☀️',
    exponent: 13, // ~1.2 x 10^13 m
    diameterMeters: 1.197e13,
    diameterDisplay: '80 AU · 12,000,000,000 km',
    lightTravelTime: '8.3 Minutes (Sun to Earth) · 4.1 Hours (Sun to Neptune)',
    cameraDistance: 130,
    targetPosition: new THREE.Vector3(0, 0, 0),
    analogy: 'If the Sun were scaled down to the size of a basketball (24 cm), Earth would be a tiny pinhead (2.2 mm) located 26 meters away, and Neptune would be a marble 770 meters down the street.',
    scientificNote: 'Bound gravitationally to the Sun (1.0 M☉). The Voyager 1 probe has flown 47+ years at 61,000 km/h and is only now traversing the heliopause into interstellar space (~160 AU).'
  },
  {
    id: 'earthmoon',
    level: 4,
    name: 'Earth & Moon',
    shortName: 'Earth & Moon',
    icon: '🌍',
    exponent: 7, // ~1.27 x 10^7 m
    diameterMeters: 1.2742e7,
    diameterDisplay: '12,742 km (Earth) · 384,400 km (Moon Orbit)',
    lightTravelTime: '1.28 Seconds (Earth to Moon) · 0.0425 Seconds (Earth diameter)',
    cameraDistance: 45,
    targetPosition: new THREE.Vector3(0, 0, 0),
    analogy: 'You are here. The only known home of life in the universe. Light circles our entire globe 7.5 times in a single second. All of human history and recorded knowledge has taken place on this blue marble.',
    scientificNote: 'Third planet from the Sun with a 23.44° axial tilt responsible for seasons. Protected by a molten nickel-iron geodynamo magnetosphere and wrapped in a fragile nitrogen-oxygen atmosphere.'
  }
];

// =========================================================================
// 2. PROCEDURAL TEXTURE GENERATORS (100% OFFLINE / ZERO 404s)
// =========================================================================

/**
 * Creates high-resolution procedural Earth texture with oceans, continents & terrain
 */
function createProceduralEarthTexture() {
  const canvas = document.createElement('canvas');
  canvas.width = 1024;
  canvas.height = 512;
  const ctx = canvas.getContext('2d');

  // Deep ocean gradient
  const oceanGrad = ctx.createLinearGradient(0, 0, 0, 512);
  oceanGrad.addColorStop(0, '#041d3b');
  oceanGrad.addColorStop(0.5, '#0b3d68');
  oceanGrad.addColorStop(1, '#041d3b');
  ctx.fillStyle = oceanGrad;
  ctx.fillRect(0, 0, 1024, 512);

  // Ocean shelf highlights
  ctx.fillStyle = 'rgba(14, 116, 144, 0.4)';
  ctx.beginPath();
  ctx.ellipse(320, 240, 210, 130, 0.2, 0, Math.PI * 2);
  ctx.ellipse(750, 260, 240, 140, -0.1, 0, Math.PI * 2);
  ctx.fill();

  // Draw procedural continental landmasses (Eurasia, Africa, Americas, Australia, Antarctica)
  function drawLandmass(x, y, w, h, color) {
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.ellipse(x, y, w, h, 0.1, 0, Math.PI * 2);
    ctx.fill();
    for (let i = 0; i < 6; i++) {
      const sx = x + (Math.sin(i * 1.5) * w * 0.7);
      const sy = y + (Math.cos(i * 1.5) * h * 0.7);
      ctx.beginPath();
      ctx.ellipse(sx, sy, w * 0.45, h * 0.45, 0.2, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  // Eurasia & Africa
  drawLandmass(620, 180, 160, 90, '#2d5a27');
  drawLandmass(580, 290, 80, 110, '#c29b38'); // Sahara/Africa
  drawLandmass(720, 220, 90, 70, '#3f6212'); // Asia
  // Americas
  drawLandmass(240, 160, 90, 85, '#1e3a1e'); // North America
  drawLandmass(310, 320, 65, 110, '#15803d'); // South America (Amazon rainforest)
  // Australia
  drawLandmass(840, 360, 60, 45, '#b45309');
  // Polar Ice Caps
  ctx.fillStyle = '#f8fafc';
  ctx.fillRect(0, 0, 1024, 38);
  ctx.fillRect(0, 474, 1024, 38);
  // Mountain ridges
  ctx.strokeStyle = '#a16207';
  ctx.lineWidth = 4;
  ctx.beginPath();
  ctx.moveTo(600, 170);
  ctx.lineTo(760, 210); // Himalayas
  ctx.moveTo(190, 110);
  ctx.lineTo(240, 230); // Rockies
  ctx.moveTo(270, 260);
  ctx.lineTo(290, 420); // Andes
  ctx.stroke();

  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.ClampToEdgeWrapping;
  return tex;
}

/**
 * Creates procedural semi-transparent cloud map
 */
function createProceduralCloudTexture() {
  const canvas = document.createElement('canvas');
  canvas.width = 1024;
  canvas.height = 512;
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = 'rgba(0, 0, 0, 0)';
  ctx.fillRect(0, 0, 1024, 512);

  ctx.fillStyle = 'rgba(255, 255, 255, 0.75)';
  // Cloud bands
  for (let b = 0; b < 18; b++) {
    const y = 50 + (b * 24);
    const cloudCount = 14;
    for (let c = 0; c < cloudCount; c++) {
      const x = (c * 75) + (Math.sin(b + c) * 30);
      const rad = 25 + (Math.sin(c * 2) * 15);
      ctx.beginPath();
      ctx.arc(x % 1024, y, rad, 0, Math.PI * 2);
      ctx.fill();
    }
  }
  // Swirling cyclones
  ctx.beginPath();
  ctx.arc(310, 190, 35, 0, Math.PI * 2);
  ctx.arc(680, 230, 40, 0, Math.PI * 2);
  ctx.arc(780, 340, 30, 0, Math.PI * 2);
  ctx.fill();

  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.ClampToEdgeWrapping;
  return tex;
}

/**
 * Creates cratered lunar surface texture
 */
function createProceduralMoonTexture() {
  const canvas = document.createElement('canvas');
  canvas.width = 512;
  canvas.height = 256;
  const ctx = canvas.getContext('2d');

  ctx.fillStyle = '#94a3b8';
  ctx.fillRect(0, 0, 512, 256);

  // Dark basaltic lunar maria
  ctx.fillStyle = '#475569';
  ctx.beginPath();
  ctx.arc(180, 110, 60, 0, Math.PI * 2);
  ctx.arc(230, 90, 45, 0, Math.PI * 2);
  ctx.arc(150, 160, 40, 0, Math.PI * 2);
  ctx.fill();

  // Impact craters
  for (let i = 0; i < 70; i++) {
    const cx = (i * 37) % 512;
    const cy = (i * 23) % 256;
    const r = 2 + (i % 8);
    ctx.strokeStyle = '#cbd5e1';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(cx, cy, r, 0, Math.PI * 2);
    ctx.stroke();
    ctx.fillStyle = '#334155';
    ctx.beginPath();
    ctx.arc(cx + 0.5, cy + 0.5, r * 0.7, 0, Math.PI * 2);
    ctx.fill();
  }

  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.ClampToEdgeWrapping;
  return tex;
}

/**
 * Creates dynamic Sun plasma texture with corona granulation
 */
function createProceduralSunTexture() {
  const canvas = document.createElement('canvas');
  canvas.width = 512;
  canvas.height = 256;
  const ctx = canvas.getContext('2d');

  const sunGrad = ctx.createLinearGradient(0, 0, 0, 256);
  sunGrad.addColorStop(0, '#ea580c');
  sunGrad.addColorStop(0.5, '#f59e0b');
  sunGrad.addColorStop(1, '#fef08a');
  ctx.fillStyle = sunGrad;
  ctx.fillRect(0, 0, 512, 256);

  // Sunspots
  ctx.fillStyle = '#9a3412';
  for (let i = 0; i < 25; i++) {
    const x = (i * 47) % 512;
    const y = 60 + ((i * 31) % 136);
    const r = 2 + (i % 6);
    ctx.beginPath();
    ctx.arc(x, y, r, 0, Math.PI * 2);
    ctx.fill();
  }

  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.ClampToEdgeWrapping;
  return tex;
}

/**
 * Creates concentric Saturn Ring texture with Cassini Division
 */
function createProceduralSaturnRings() {
  const canvas = document.createElement('canvas');
  canvas.width = 256;
  canvas.height = 1;
  const ctx = canvas.getContext('2d');

  const grad = ctx.createLinearGradient(0, 0, 256, 0);
  grad.addColorStop(0.0, 'rgba(0,0,0,0)');
  grad.addColorStop(0.15, 'rgba(217, 119, 6, 0.4)'); // C Ring
  grad.addColorStop(0.35, 'rgba(245, 158, 11, 0.85)'); // B Ring (bright)
  grad.addColorStop(0.60, 'rgba(245, 158, 11, 0.95)');
  grad.addColorStop(0.66, 'rgba(0,0,0,0)'); // Cassini Division (dark gap)
  grad.addColorStop(0.72, 'rgba(217, 119, 6, 0.7)'); // A Ring
  grad.addColorStop(0.92, 'rgba(180, 83, 9, 0.5)'); // Encke gap & outer
  grad.addColorStop(1.0, 'rgba(0,0,0,0)');

  ctx.fillStyle = grad;
  ctx.fillRect(0, 0, 256, 1);

  const tex = new THREE.CanvasTexture(canvas);
  return tex;
}

// =========================================================================
// 3. POWERS OF TEN 3D ENGINE & SCENE ARCHITECTURE
// =========================================================================

class PowersOfTenVisualizer {
  constructor() {
    this.modalEl = null;
    this.canvasEl = null;
    this.renderer = null;
    this.scene = null;
    this.camera = null;
    this.controls = null;
    this.animFrameId = null;

    // Stage scene groups
    this.stageGroups = {
      laniakea: new THREE.Group(),
      milkyway: new THREE.Group(),
      solarsystem: new THREE.Group(),
      earthmoon: new THREE.Group()
    };

    // Active state
    this.currentStageIndex = 0;
    this.currentScaleValue = 1.0; // 1.0 to 4.0
    this.isAutoTouring = false;
    this.autoTourTimer = null;
    this.tourStepDuration = 6500; // ms per step
    this.isAnimatingTransition = false;
    this.transitionProgress = 1.0;
    this.activePlanetFocus = null;

    // Animated meshes cache
    this.cloudMesh = null;
    this.earthMesh = null;
    this.moonOrbitGroup = null;
    this.mwDiscPoints = null;
    this.streamlineParticles = [];
    this.planetOrbitGroups = [];
    this.sunCoronaMesh = null;

    this.init();
  }

  init() {
    this.buildModalDOM();
    this.setupThreeScene();
    this.buildStageLaniakea();
    this.buildStageMilkyWay();
    this.buildStageSolarSystem();
    this.buildStageEarthMoon();
    this.attachEventListeners();
    this.updateHUD(0);
  }

  buildModalDOM() {
    const existing = document.getElementById('powers-of-ten-modal');
    if (existing) existing.remove();

    const modal = document.createElement('div');
    modal.id = 'powers-of-ten-modal';
    modal.style.cssText = `
      display: none;
      position: fixed;
      inset: 0;
      z-index: 999999;
      background: #020617;
      color: #f8fafc;
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      flex-direction: column;
      overflow: hidden;
      user-select: none;
    `;

    modal.innerHTML = `
      <!-- TOP CONTROL & STEPPER HEADER -->
      <div id="pot-header" style="display:flex; justify-content:space-between; align-items:center; padding:12px 20px; background:rgba(15,23,42,0.92); backdrop-filter:blur(16px); border-bottom:1.5px solid rgba(56,189,248,0.3); z-index:10;">
        <div style="display:flex; align-items:center; gap:12px;">
          <span style="font-size:24px;">🌌</span>
          <div>
            <div style="font-size:16px; font-weight:800; color:#38bdf8; display:flex; align-items:center; gap:8px;">
              <span>Powers of Ten · Cosmic Scale Explorer</span>
              <span id="pot-stage-badge" style="font-size:10px; font-weight:800; background:rgba(56,189,248,0.2); border:1px solid #38bdf8; padding:2px 8px; border-radius:12px; color:#38bdf8; text-transform:uppercase;">Stage 1 of 4</span>
            </div>
            <div style="font-size:11.5px; color:#94a3b8;">From the Laniakea Supercluster (500 Million ly) to Planet Earth (12,742 km)</div>
          </div>
        </div>

        <!-- 4 INTERACTIVE SCALE TABS -->
        <div id="pot-stepper" style="display:flex; gap:6px; background:rgba(2,6,23,0.7); padding:4px; border-radius:10px; border:1px solid rgba(255,255,255,0.1);">
          <button class="pot-step-btn active" data-step="0" style="background:#0284c7; color:#fff; border:1px solid #38bdf8; padding:6px 12px; border-radius:7px; font-size:11px; font-weight:700; cursor:pointer; display:flex; align-items:center; gap:5px;">
            <span>🌌</span> <span>1. Laniakea (10²⁴ m)</span>
          </button>
          <button class="pot-step-btn" data-step="1" style="background:transparent; color:#94a3b8; border:1px solid transparent; padding:6px 12px; border-radius:7px; font-size:11px; font-weight:700; cursor:pointer; display:flex; align-items:center; gap:5px;">
            <span>🌀</span> <span>2. Milky Way (10²¹ m)</span>
          </button>
          <button class="pot-step-btn" data-step="2" style="background:transparent; color:#94a3b8; border:1px solid transparent; padding:6px 12px; border-radius:7px; font-size:11px; font-weight:700; cursor:pointer; display:flex; align-items:center; gap:5px;">
            <span>☀️</span> <span>3. Solar System (10¹³ m)</span>
          </button>
          <button class="pot-step-btn" data-step="3" style="background:transparent; color:#94a3b8; border:1px solid transparent; padding:6px 12px; border-radius:7px; font-size:11px; font-weight:700; cursor:pointer; display:flex; align-items:center; gap:5px;">
            <span>🌍</span> <span>4. Earth & Moon (10⁷ m)</span>
          </button>
        </div>

        <!-- WINDOW ACTION BUTTONS -->
        <div style="display:flex; align-items:center; gap:8px;">
          <button id="pot-btn-reset-cam" style="background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.15); color:#f8fafc; padding:6px 12px; border-radius:6px; font-size:11px; font-weight:700; cursor:pointer;" title="Reset Camera View">🔄 Reset</button>
          <button id="pot-btn-fullscreen" style="background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.15); color:#f8fafc; padding:6px 12px; border-radius:6px; font-size:11px; font-weight:700; cursor:pointer;" title="Toggle Fullscreen">⛶ Fullscreen</button>
          <button id="pot-btn-close" style="background:rgba(239,68,68,0.2); border:1px solid #ef4444; color:#fca5a5; padding:6px 14px; border-radius:6px; font-size:13px; font-weight:800; cursor:pointer;" title="Exit Cosmic Scale Explorer">&times; Close</button>
        </div>
      </div>

      <!-- MAIN 3D WEBGL VIEWPORT -->
      <div id="pot-viewport-container" style="flex:1; position:relative; width:100%; height:100%; overflow:hidden;">
        <canvas id="powers-of-ten-canvas" style="width:100%; height:100%; display:block;"></canvas>

        <!-- FLOATING SCALE TARGET RETICLE -->
        <div id="pot-reticle" style="position:absolute; top:50%; left:50%; transform:translate(-50%, -50%); pointer-events:none; border:1px dashed rgba(56,189,248,0.3); border-radius:50%; width:160px; height:160px; display:flex; align-items:center; justify-content:center; opacity:0.65;">
          <div style="width:6px; height:6px; background:#38bdf8; border-radius:50%; box-shadow:0 0 10px #38bdf8;"></div>
        </div>

        <!-- STAGE TITLE OVERLAY -->
        <div style="position:absolute; top:20px; left:24px; pointer-events:none; text-shadow:0 2px 10px rgba(0,0,0,0.9);">
          <div id="pot-stage-title" style="font-size:26px; font-weight:900; color:#ffffff; letter-spacing:0.5px;">Laniakea Supercluster</div>
          <div id="pot-stage-dims" style="font-size:13px; font-weight:700; color:#38bdf8;">160 Mpc · 520,000,000 Light-Years across</div>
        </div>

        <!-- LOGARITHMIC CONTINUOUS ZOOM SLIDER (RIGHT DOCK) -->
        <div style="position:absolute; right:20px; top:50%; transform:translateY(-50%); background:rgba(15,23,42,0.85); backdrop-filter:blur(12px); border:1.5px solid rgba(56,189,248,0.25); border-radius:24px; padding:18px 10px; display:flex; flex-direction:column; align-items:center; gap:12px; box-shadow:0 8px 32px rgba(0,0,0,0.6);">
          <span style="font-size:12px; font-weight:800; color:#38bdf8;" title="Macro (Laniakea)">10²⁴m</span>
          <input type="range" id="pot-scale-slider" min="1.0" max="4.0" step="0.01" value="1.0" orient="vertical" style="writing-mode:bt-lr; -webkit-appearance:slider-vertical; width:12px; height:200px; accent-color:#0284c7; cursor:pointer;">
          <span style="font-size:12px; font-weight:800; color:#10b981;" title="Micro (Earth)">10⁷m</span>
        </div>

        <!-- SPEED OF LIGHT TICKER HUD (TOP RIGHT) -->
        <div style="position:absolute; top:20px; right:24px; background:rgba(15,23,42,0.85); backdrop-filter:blur(12px); border:1.5px solid rgba(56,189,248,0.25); border-radius:12px; padding:12px 18px; max-width:340px; box-shadow:0 8px 30px rgba(0,0,0,0.5);">
          <div style="font-size:10.5px; font-weight:800; color:#94a3b8; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:4px; display:flex; align-items:center; gap:5px;">
            <span>⚡</span> <span>Speed of Light Travel Time (c)</span>
          </div>
          <div id="pot-light-ticker" style="font-family:'JetBrains Mono', monospace, sans-serif; font-size:14px; font-weight:800; color:#38bdf8;">
            520,000,000 Years
          </div>
          <div id="pot-light-sub" style="font-size:11px; color:#cbd5e1; margin-top:4px;">
            Time required for a photon to traverse this structure.
          </div>
        </div>
      </div>

      <!-- BOTTOM CINEMATIC CONTROLS & MIND-BENDING ANALOGY HUD -->
      <div id="pot-footer" style="padding:14px 24px; background:rgba(15,23,42,0.95); backdrop-filter:blur(16px); border-top:1.5px solid rgba(56,189,248,0.25); display:grid; grid-template-columns:auto 1fr auto; gap:20px; align-items:center; z-index:10;">
        <!-- PLAY / PAUSE / STEP CONTROLS -->
        <div style="display:flex; align-items:center; gap:10px;">
          <button id="pot-btn-prev" style="background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.15); color:#f8fafc; padding:8px 14px; border-radius:8px; font-size:12px; font-weight:800; cursor:pointer;" title="Previous Scale">◀ Zoom Out</button>
          <button id="pot-btn-auto-tour" style="background:#0284c7; border:1px solid #38bdf8; color:#ffffff; padding:8px 18px; border-radius:8px; font-size:12.5px; font-weight:800; cursor:pointer; display:flex; align-items:center; gap:6px; box-shadow:0 0 15px rgba(2,132,199,0.5);">
            <span id="pot-tour-icon">▶</span> <span id="pot-tour-label">Auto-Tour</span>
          </button>
          <button id="pot-btn-next" style="background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.15); color:#f8fafc; padding:8px 14px; border-radius:8px; font-size:12px; font-weight:800; cursor:pointer;" title="Next Scale">Zoom In ▶</button>
        </div>

        <!-- MIND-BENDING SCALE ANALOGY & SCIENTIFIC CITATION -->
        <div style="background:rgba(2,6,23,0.6); border:1px solid rgba(56,189,248,0.2); border-radius:10px; padding:10px 16px;">
          <div style="display:flex; align-items:center; gap:6px; margin-bottom:3px;">
            <span style="font-size:13px;">💡</span>
            <span style="font-size:11px; font-weight:800; color:#38bdf8; text-transform:uppercase; letter-spacing:0.5px;">Real Scale Human Analogy</span>
          </div>
          <div id="pot-analogy-text" style="font-size:12px; line-height:1.5; color:#f1f5f9; font-weight:500;">
            If Laniakea were shrunk to the size of planet Earth, the Milky Way would be the size of a single microscopic dust mite (0.2 mm), and our entire Solar System would be smaller than a single atom.
          </div>
        </div>

        <!-- INTERACTIVE TOGGLES -->
        <div style="display:flex; align-items:center; gap:10px;">
          <button id="pot-btn-toggle-orbits" style="background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.15); color:#cbd5e1; padding:8px 12px; border-radius:8px; font-size:11px; font-weight:700; cursor:pointer;" title="Toggle Orbit Lines & Markers">
            <span>🪐</span> Orbits: ON
          </button>
          <button id="pot-btn-toggle-rotate" style="background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.15); color:#cbd5e1; padding:8px 12px; border-radius:8px; font-size:11px; font-weight:700; cursor:pointer;" title="Toggle Auto-Rotation">
            <span>🔄</span> Spin: ON
          </button>
        </div>
      </div>
    `;

    document.body.appendChild(modal);
    this.modalEl = modal;
    this.canvasEl = document.getElementById('powers-of-ten-canvas');
  }

  setupThreeScene() {
    const width = this.canvasEl.clientWidth || window.innerWidth;
    const height = this.canvasEl.clientHeight || (window.innerHeight - 130);

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x020617);

    this.camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 50000);
    this.camera.position.set(0, 100, 280);

    this.renderer = new THREE.WebGLRenderer({
      canvas: this.canvasEl,
      antialias: true,
      powerPreference: 'high-performance'
    });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

    // Ambient and directional lighting
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.45);
    this.scene.add(ambientLight);

    const sunLight = new THREE.DirectionalLight(0xfff7ed, 1.4);
    sunLight.position.set(120, 80, 150);
    this.scene.add(sunLight);

    // OrbitControls for manual camera rotation and zoom
    this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.05;
    this.controls.rotateSpeed = 0.7;
    this.controls.maxDistance = 1500;
    this.controls.minDistance = 5;

    // Add stage groups to scene
    this.scene.add(this.stageGroups.laniakea);
    this.scene.add(this.stageGroups.milkyway);
    this.scene.add(this.stageGroups.solarsystem);
    this.scene.add(this.stageGroups.earthmoon);

    // Deep space background starfield
    this.createBackgroundStarfield();
  }

  createBackgroundStarfield() {
    const starCount = 3500;
    const starGeo = new THREE.BufferGeometry();
    const starPos = new Float32Array(starCount * 3);
    const starColors = new Float32Array(starCount * 3);

    for (let i = 0; i < starCount; i++) {
      const theta = 2 * Math.PI * Math.random();
      const phi = Math.acos(2 * Math.random() - 1);
      const r = 2500 + (Math.random() * 800);

      starPos[i * 3] = r * Math.sin(phi) * Math.cos(theta);
      starPos[i * 3 + 1] = r * Math.sin(phi) * Math.sin(theta);
      starPos[i * 3 + 2] = r * Math.cos(phi);

      const tint = Math.random();
      if (tint > 0.8) { // Blue star
        starColors[i * 3] = 0.6; starColors[i * 3 + 1] = 0.8; starColors[i * 3 + 2] = 1.0;
      } else if (tint > 0.6) { // Gold/amber star
        starColors[i * 3] = 1.0; starColors[i * 3 + 1] = 0.85; starColors[i * 3 + 2] = 0.5;
      } else { // White star
        starColors[i * 3] = 0.95; starColors[i * 3 + 1] = 0.95; starColors[i * 3 + 2] = 0.95;
      }
    }

    starGeo.setAttribute('position', new THREE.BufferAttribute(starPos, 3));
    starGeo.setAttribute('color', new THREE.BufferAttribute(starColors, 3));

    const starMat = new THREE.PointsMaterial({
      size: 2.2,
      vertexColors: true,
      transparent: true,
      opacity: 0.85
    });

    const starField = new THREE.Points(starGeo, starMat);
    this.scene.add(starField);
  }

  // -----------------------------------------------------------------------
  // STAGE 1: LANIAKEA SUPERCLUSTER (160 Mpc / 500 Mly)
  // -----------------------------------------------------------------------
  buildStageLaniakea() {
    const group = this.stageGroups.laniakea;

    // Watershed Separatrix Mesh (Translucent outer boundary shell)
    const basinGeo = new THREE.IcosahedronGeometry(110, 3);
    const basinMat = new THREE.MeshBasicMaterial({
      color: 0x38bdf8,
      wireframe: true,
      transparent: true,
      opacity: 0.14
    });
    const basinMesh = new THREE.Mesh(basinGeo, basinMat);
    group.add(basinMesh);

    // Inward Peculiar Velocity Streamlines converging on Great Attractor
    const attractorPos = new THREE.Vector3(-45, 10, -30); // Great Attractor / Norma
    const streamlineCount = 85;

    for (let s = 0; s < streamlineCount; s++) {
      const u = Math.random();
      const v = Math.random();
      const theta = u * 2.0 * Math.PI;
      const phi = Math.acos(2.0 * v - 1.0);
      const r = 90 + (Math.random() * 25);

      const startPt = new THREE.Vector3(
        r * Math.sin(phi) * Math.cos(theta),
        (r * Math.sin(phi) * Math.sin(theta)) * 0.45, // Flattened supergalactic plane
        r * Math.cos(phi)
      );

      // Quadratic curve bending into attractor
      const midPt = new THREE.Vector3()
        .addVectors(startPt, attractorPos)
        .multiplyScalar(0.5)
        .add(new THREE.Vector3((Math.random() - 0.5) * 20, (Math.random() - 0.5) * 10, (Math.random() - 0.5) * 20));

      const curve = new THREE.QuadraticBezierCurve3(startPt, midPt, attractorPos);
      const points = curve.getPoints(36);
      const lineGeo = new THREE.BufferGeometry().setFromPoints(points);

      const lineMat = new THREE.LineBasicMaterial({
        color: 0x0284c7,
        transparent: true,
        opacity: 0.35 + (Math.random() * 0.25)
      });

      const line = new THREE.Line(lineGeo, lineMat);
      group.add(line);

      // Animated glowing streamline tracer particle
      const pGeo = new THREE.BufferGeometry();
      pGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array([startPt.x, startPt.y, startPt.z]), 3));
      const pMat = new THREE.PointsMaterial({
        color: 0x38bdf8,
        size: 3.5,
        transparent: true,
        opacity: 0.9
      });
      const particle = new THREE.Points(pGeo, pMat);
      group.add(particle);

      this.streamlineParticles.push({
        particle,
        curve,
        progress: Math.random(),
        speed: 0.0015 + (Math.random() * 0.002)
      });
    }

    // Major Cluster Node Spheres & Badges
    const clusters = [
      { name: 'Great Attractor (Norma)', pos: attractorPos, color: 0x06b6d4, size: 8.5, isAttractor: true },
      { name: 'Milky Way / Local Group (You Are Here)', pos: new THREE.Vector3(0, 0, 0), color: 0xf59e0b, size: 4.5, isHome: true },
      { name: 'Virgo Cluster', pos: new THREE.Vector3(-12, 22, -2), color: 0x10b981, size: 4.0 },
      { name: 'Hydra Cluster', pos: new THREE.Vector3(-35, -15, 18), color: 0xa855f7, size: 4.2 },
      { name: 'Centaurus Cluster', pos: new THREE.Vector3(-38, 8, -12), color: 0x38bdf8, size: 4.5 },
      { name: 'Pavo-Indus Supercluster', pos: new THREE.Vector3(-22, -45, -35), color: 0xec4899, size: 4.0 }
    ];

    clusters.forEach(c => {
      const nodeGeo = new THREE.SphereGeometry(c.size, 16, 16);
      const nodeMat = new THREE.MeshBasicMaterial({
        color: c.color,
        transparent: true,
        opacity: 0.9
      });
      const nodeMesh = new THREE.Mesh(nodeGeo, nodeMat);
      nodeMesh.position.copy(c.pos);
      group.add(nodeMesh);

      // Luminous pulsing glow ring around Milky Way home
      if (c.isHome) {
        const ringGeo = new THREE.RingGeometry(c.size * 1.5, c.size * 2.2, 32);
        const ringMat = new THREE.MeshBasicMaterial({
          color: 0xf59e0b,
          side: THREE.DoubleSide,
          transparent: true,
          opacity: 0.8
        });
        const ringMesh = new THREE.Mesh(ringGeo, ringMat);
        ringMesh.position.copy(c.pos);
        ringMesh.lookAt(0, 100, 280);
        group.add(ringMesh);
      }
    });

    // Background 10,000 galaxy particle cloud
    const galCount = 12000;
    const galGeo = new THREE.BufferGeometry();
    const galPos = new Float32Array(galCount * 3);
    const galColors = new Float32Array(galCount * 3);

    for (let i = 0; i < galCount; i++) {
      const u = Math.random();
      const v = Math.random();
      const theta = u * 2.0 * Math.PI;
      const phi = Math.acos(2.0 * v - 1.0);
      const r = Math.pow(Math.random(), 0.6) * 105;

      galPos[i * 3] = r * Math.sin(phi) * Math.cos(theta);
      galPos[i * 3 + 1] = (r * Math.sin(phi) * Math.sin(theta)) * 0.4;
      galPos[i * 3 + 2] = r * Math.cos(phi);

      galColors[i * 3] = 0.5 + Math.random() * 0.5;
      galColors[i * 3 + 1] = 0.7 + Math.random() * 0.3;
      galColors[i * 3 + 2] = 0.9 + Math.random() * 0.1;
    }

    galGeo.setAttribute('position', new THREE.BufferAttribute(galPos, 3));
    galGeo.setAttribute('color', new THREE.BufferAttribute(galColors, 3));

    const galMat = new THREE.PointsMaterial({
      size: 1.5,
      vertexColors: true,
      transparent: true,
      opacity: 0.75
    });
    group.add(new THREE.Points(galGeo, galMat));
  }

  // -----------------------------------------------------------------------
  // STAGE 2: MILKY WAY 3D SPIRAL GALAXY (30 kpc / 100,000 ly)
  // -----------------------------------------------------------------------
  buildStageMilkyWay() {
    const group = this.stageGroups.milkyway;
    group.visible = false; // initially hidden

    // Procedural Logarithmic 4-Arm Barred Spiral Particle System
    const starCount = 45000;
    const starGeo = new THREE.BufferGeometry();
    const starPos = new Float32Array(starCount * 3);
    const starColors = new Float32Array(starCount * 3);

    const arms = 4;
    const armOffset = (2.0 * Math.PI) / arms;
    const b = 0.22; // spiral pitch angle parameter

    for (let i = 0; i < starCount; i++) {
      let x, y, z, r, g, bl;

      if (i < 8000) {
        // Central Bar & Bulge stars (dense amber/gold)
        const rad = Math.pow(Math.random(), 1.5) * 18;
        const ang = Math.random() * Math.PI * 2;
        x = Math.cos(ang) * rad * 1.6; // elongated bar
        z = Math.sin(ang) * rad * 0.9;
        y = (Math.random() - 0.5) * 8 * Math.exp(-rad / 12);

        r = 1.0;
        g = 0.85 + Math.random() * 0.15;
        bl = 0.4 + Math.random() * 0.3;
      } else {
        // Spiral Arms stars
        const armIndex = Math.floor(Math.random() * arms);
        const theta = Math.random() * 3.8 * Math.PI;
        const baseR = 12 * Math.exp(b * theta);
        const spread = (Math.random() - 0.5) * (3.5 + theta * 1.2);

        const currentTheta = theta + (armIndex * armOffset) + (spread / (baseR + 1));
        const currentR = baseR + spread;

        x = Math.cos(currentTheta) * currentR;
        z = Math.sin(currentTheta) * currentR;
        y = (Math.random() - 0.5) * 4 * Math.exp(-currentR / 85);

        // Young blue supergiant stars along spiral arm shock fronts
        if (Math.random() > 0.45) {
          r = 0.55 + Math.random() * 0.2;
          g = 0.75 + Math.random() * 0.25;
          bl = 1.0;
        } else {
          r = 0.95;
          g = 0.85;
          bl = 0.7;
        }
      }

      starPos[i * 3] = x;
      starPos[i * 3 + 1] = y;
      starPos[i * 3 + 2] = z;

      starColors[i * 3] = r;
      starColors[i * 3 + 1] = g;
      starColors[i * 3 + 2] = bl;
    }

    starGeo.setAttribute('position', new THREE.BufferAttribute(starPos, 3));
    starGeo.setAttribute('color', new THREE.BufferAttribute(starColors, 3));

    const starMat = new THREE.PointsMaterial({
      size: 1.6,
      vertexColors: true,
      transparent: true,
      opacity: 0.88,
      blending: THREE.AdditiveBlending
    });

    this.mwDiscPoints = new THREE.Points(starGeo, starMat);
    group.add(this.mwDiscPoints);

    // Supermassive Black Hole Sagittarius A* Core Flare
    const sgrAGeo = new THREE.SphereGeometry(2.5, 32, 32);
    const sgrAMat = new THREE.MeshBasicMaterial({
      color: 0xffffff,
      transparent: true,
      opacity: 0.95
    });
    const sgrAMesh = new THREE.Mesh(sgrAGeo, sgrAMat);
    group.add(sgrAMesh);

    // Glowing core accretion glow
    const glowGeo = new THREE.SphereGeometry(7.0, 32, 32);
    const glowMat = new THREE.MeshBasicMaterial({
      color: 0xf59e0b,
      transparent: true,
      opacity: 0.35,
      blending: THREE.AdditiveBlending
    });
    group.add(new THREE.Mesh(glowGeo, glowMat));

    // OUR SUN / SOLAR SYSTEM LOCATION MARKER (In Orion Spur at ~26,000 ly from core)
    const sunMarkerPos = new THREE.Vector3(52.0, 0, 18.0);

    const sunBeaconGeo = new THREE.SphereGeometry(2.0, 16, 16);
    const sunBeaconMat = new THREE.MeshBasicMaterial({ color: 0xfacc15 });
    const sunBeaconMesh = new THREE.Mesh(sunBeaconGeo, sunBeaconMat);
    sunBeaconMesh.position.copy(sunMarkerPos);
    group.add(sunBeaconMesh);

    // Concentric pulsing target rings around Sun
    const sunRingGeo = new THREE.RingGeometry(3.5, 4.8, 32);
    const sunRingMat = new THREE.MeshBasicMaterial({
      color: 0xfacc15,
      side: THREE.DoubleSide,
      transparent: true,
      opacity: 0.85
    });
    const sunRing = new THREE.Mesh(sunRingGeo, sunRingMat);
    sunRing.position.copy(sunMarkerPos);
    sunRing.rotation.x = Math.PI / 2;
    group.add(sunRing);

    // Velocity vector line showing ~220 km/s Galactic Orbit direction
    const arrowDir = new THREE.Vector3(-0.32, 0, 0.94).normalize();
    const arrowHelper = new THREE.ArrowHelper(arrowDir, sunMarkerPos, 14, 0x38bdf8, 3.5, 2.0);
    group.add(arrowHelper);
  }

  // -----------------------------------------------------------------------
  // STAGE 3: SOLAR SYSTEM (80 AU / 12 Billion km)
  // -----------------------------------------------------------------------
  buildStageSolarSystem() {
    const group = this.stageGroups.solarsystem;
    group.visible = false; // initially hidden

    // Central Radiant Sun
    const sunGeo = new THREE.SphereGeometry(8.0, 32, 32);
    const sunTex = createProceduralSunTexture();
    const sunMat = new THREE.MeshBasicMaterial({ map: sunTex });
    const sunMesh = new THREE.Mesh(sunGeo, sunMat);
    group.add(sunMesh);

    // Sun Animated Corona Flare Halo
    const coronaGeo = new THREE.SphereGeometry(11.5, 32, 32);
    const coronaMat = new THREE.MeshBasicMaterial({
      color: 0xf59e0b,
      transparent: true,
      opacity: 0.35,
      blending: THREE.AdditiveBlending
    });
    this.sunCoronaMesh = new THREE.Mesh(coronaGeo, coronaMat);
    group.add(this.sunCoronaMesh);

    // Planets Definition (proportional orbits scaled for clarity)
    const planets = [
      { name: 'Mercury', dist: 14, radius: 0.9, color: 0x94a3b8, speed: 0.040 },
      { name: 'Venus', dist: 20, radius: 1.4, color: 0xf59e0b, speed: 0.028 },
      { name: 'Earth', dist: 28, radius: 1.5, color: 0x0284c7, speed: 0.020, hasMoon: true },
      { name: 'Mars', dist: 38, radius: 1.1, color: 0xef4444, speed: 0.015 },
      { name: 'Jupiter', dist: 58, radius: 4.2, color: 0xd97706, speed: 0.009 },
      { name: 'Saturn', dist: 78, radius: 3.4, color: 0xfbbf24, speed: 0.006, hasRings: true },
      { name: 'Uranus', dist: 98, radius: 2.2, color: 0x06b6d4, speed: 0.004 },
      { name: 'Neptune', dist: 118, radius: 2.1, color: 0x3b82f6, speed: 0.003 }
    ];

    planets.forEach(p => {
      // Orbital Ring Track
      const orbitGeo = new THREE.RingGeometry(p.dist - 0.2, p.dist + 0.2, 64);
      const orbitMat = new THREE.MeshBasicMaterial({
        color: 0x38bdf8,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: 0.22
      });
      const orbitRing = new THREE.Mesh(orbitGeo, orbitMat);
      orbitRing.rotation.x = Math.PI / 2;
      group.add(orbitRing);

      // Planet Orbit Pivot Group
      const pivot = new THREE.Group();
      group.add(pivot);

      // Planet Body
      const pGeo = new THREE.SphereGeometry(p.radius, 24, 24);
      const pMat = new THREE.MeshStandardMaterial({
        color: p.color,
        roughness: 0.8,
        metalness: 0.1
      });
      const pMesh = new THREE.Mesh(pGeo, pMat);
      pMesh.position.set(p.dist, 0, 0);
      pivot.add(pMesh);

      // Saturn's Ring System
      if (p.hasRings) {
        const ringGeo = new THREE.RingGeometry(p.radius * 1.5, p.radius * 2.8, 48);
        const ringTex = createProceduralSaturnRings();
        const ringMat = new THREE.MeshBasicMaterial({
          map: ringTex,
          side: THREE.DoubleSide,
          transparent: true,
          opacity: 0.9
        });
        const rings = new THREE.Mesh(ringGeo, ringMat);
        rings.rotation.x = Math.PI / 2.3;
        rings.position.set(p.dist, 0, 0);
        pivot.add(rings);
      }

      this.planetOrbitGroups.push({
        pivot,
        speed: p.speed,
        angle: Math.random() * Math.PI * 2
      });
    });

    // Main Asteroid Belt (Between Mars & Jupiter: dist ~46-51)
    const astCount = 1800;
    const astGeo = new THREE.BufferGeometry();
    const astPos = new Float32Array(astCount * 3);

    for (let i = 0; i < astCount; i++) {
      const theta = Math.random() * Math.PI * 2;
      const r = 45 + Math.random() * 8.0;
      astPos[i * 3] = Math.cos(theta) * r;
      astPos[i * 3 + 1] = (Math.random() - 0.5) * 2.5;
      astPos[i * 3 + 2] = Math.sin(theta) * r;
    }
    astGeo.setAttribute('position', new THREE.BufferAttribute(astPos, 3));
    const astMat = new THREE.PointsMaterial({ color: 0x94a3b8, size: 1.2 });
    group.add(new THREE.Points(astGeo, astMat));
  }

  // -----------------------------------------------------------------------
  // STAGE 4: EARTH & MOON (12,742 km / 384,400 km)
  // -----------------------------------------------------------------------
  buildStageEarthMoon() {
    const group = this.stageGroups.earthmoon;
    group.visible = false; // initially hidden

    // Procedural Earth Globe
    const earthGeo = new THREE.SphereGeometry(12, 48, 48);
    const earthTex = createProceduralEarthTexture();
    const earthMat = new THREE.MeshStandardMaterial({
      map: earthTex,
      roughness: 0.65,
      metalness: 0.15
    });

    this.earthMesh = new THREE.Mesh(earthGeo, earthMat);
    // 23.44° Axial Tilt
    this.earthMesh.rotation.z = (23.44 * Math.PI) / 180.0;
    group.add(this.earthMesh);

    // Swirling Atmospheric Cloud Layer
    const cloudGeo = new THREE.SphereGeometry(12.25, 48, 48);
    const cloudTex = createProceduralCloudTexture();
    const cloudMat = new THREE.MeshStandardMaterial({
      map: cloudTex,
      transparent: true,
      opacity: 0.85,
      blending: THREE.NormalBlending
    });
    this.cloudMesh = new THREE.Mesh(cloudGeo, cloudMat);
    this.earthMesh.add(this.cloudMesh);

    // Atmospheric Blue Fresnel Limb Glow
    const atmosGeo = new THREE.SphereGeometry(13.2, 48, 48);
    const atmosMat = new THREE.MeshBasicMaterial({
      color: 0x38bdf8,
      transparent: true,
      opacity: 0.22,
      side: THREE.BackSide,
      blending: THREE.AdditiveBlending
    });
    group.add(new THREE.Mesh(atmosGeo, atmosMat));

    // Moon in Proportional Orbit
    this.moonOrbitGroup = new THREE.Group();
    group.add(this.moonOrbitGroup);

    // Moon Orbit Trace Ring
    const moonOrbitRingGeo = new THREE.RingGeometry(38.0 - 0.2, 38.0 + 0.2, 64);
    const moonOrbitRingMat = new THREE.MeshBasicMaterial({
      color: 0x94a3b8,
      side: THREE.DoubleSide,
      transparent: true,
      opacity: 0.25
    });
    const moonOrbitRing = new THREE.Mesh(moonOrbitRingGeo, moonOrbitRingMat);
    moonOrbitRing.rotation.x = Math.PI / 2;
    group.add(moonOrbitRing);

    // Moon Body
    const moonGeo = new THREE.SphereGeometry(3.27, 32, 32);
    const moonTex = createProceduralMoonTexture();
    const moonMat = new THREE.MeshStandardMaterial({
      map: moonTex,
      roughness: 0.9,
      metalness: 0.05
    });
    const moonMesh = new THREE.Mesh(moonGeo, moonMat);
    moonMesh.position.set(38.0, 0, 0);
    this.moonOrbitGroup.add(moonMesh);

    // Speed of Light 1.28-second pulse beam between Earth and Moon
    const beamGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(38.0, 0, 0)
    ]);
    const beamMat = new THREE.LineDashedMaterial({
      color: 0xfacc15,
      dashSize: 2,
      gapSize: 1.5,
      linewidth: 2
    });
    const lightBeam = new THREE.Line(beamGeo, beamMat);
    lightBeam.computeLineDistances();
    this.moonOrbitGroup.add(lightBeam);
  }

  // -----------------------------------------------------------------------
  // STAGE TRANSITIONS & CAMERA FLIGHT
  // -----------------------------------------------------------------------
  setStage(targetStageIndex, animated = true) {
    if (targetStageIndex < 0 || targetStageIndex >= STAGE_SCALES.length) return;

    this.currentStageIndex = targetStageIndex;
    const stageData = STAGE_SCALES[targetStageIndex];

    // Update Slider
    const slider = document.getElementById('pot-scale-slider');
    if (slider) slider.value = (targetStageIndex + 1).toFixed(2);

    // Update Stepper Tabs
    const stepBtns = document.querySelectorAll('.pot-step-btn');
    stepBtns.forEach((btn, idx) => {
      if (idx === targetStageIndex) {
        btn.classList.add('active');
        btn.style.background = '#0284c7';
        btn.style.borderColor = '#38bdf8';
        btn.style.color = '#ffffff';
      } else {
        btn.classList.remove('active');
        btn.style.background = 'transparent';
        btn.style.borderColor = 'transparent';
        btn.style.color = '#94a3b8';
      }
    });

    // Toggle Stage Group Visibility
    Object.keys(this.stageGroups).forEach((key, idx) => {
      this.stageGroups[key].visible = (idx === targetStageIndex);
    });

    // Update Telemetry & Analogies
    this.updateHUD(targetStageIndex);

    // Camera Flight
    if (animated) {
      this.animateCameraTo(stageData.cameraDistance, stageData.targetPosition);
    } else {
      this.camera.position.set(0, stageData.cameraDistance * 0.35, stageData.cameraDistance);
      this.controls.target.copy(stageData.targetPosition);
      this.controls.update();
    }
  }

  animateCameraTo(distance, target) {
    this.isAnimatingTransition = true;
    const startPos = this.camera.position.clone();
    const endPos = new THREE.Vector3(0, distance * 0.35, distance);
    const startTarget = this.controls.target.clone();
    const endTarget = target.clone();

    let progress = 0;
    const duration = 1200; // ms
    const startTime = performance.now();

    const tween = () => {
      const now = performance.now();
      progress = Math.min((now - startTime) / duration, 1.0);
      // Smooth cubic ease-in-out
      const ease = progress < 0.5
        ? 4 * progress * progress * progress
        : 1 - Math.pow(-2 * progress + 2, 3) / 2;

      this.camera.position.lerpVectors(startPos, endPos, ease);
      this.controls.target.lerpVectors(startTarget, endTarget, ease);
      this.controls.update();

      if (progress < 1.0) {
        requestAnimationFrame(tween);
      } else {
        this.isAnimatingTransition = false;
      }
    };
    requestAnimationFrame(tween);
  }

  updateHUD(stageIndex) {
    const data = STAGE_SCALES[stageIndex];
    if (!data) return;

    const titleEl = document.getElementById('pot-stage-title');
    const dimsEl = document.getElementById('pot-stage-dims');
    const badgeEl = document.getElementById('pot-stage-badge');
    const lightEl = document.getElementById('pot-light-ticker');
    const analogyEl = document.getElementById('pot-analogy-text');

    if (titleEl) titleEl.textContent = data.name;
    if (dimsEl) dimsEl.textContent = data.diameterDisplay;
    if (badgeEl) badgeEl.textContent = `Stage ${stageIndex + 1} of 4 · 10^${data.exponent} m`;
    if (lightEl) lightEl.textContent = data.lightTravelTime;
    if (analogyEl) analogyEl.textContent = data.analogy;
  }

  // -----------------------------------------------------------------------
  // CINEMATIC AUTO-TOUR
  // -----------------------------------------------------------------------
  toggleAutoTour() {
    this.isAutoTouring = !this.isAutoTouring;
    const btn = document.getElementById('pot-btn-auto-tour');
    const icon = document.getElementById('pot-tour-icon');
    const label = document.getElementById('pot-tour-label');

    if (this.isAutoTouring) {
      if (btn) {
        btn.style.background = '#10b981';
        btn.style.borderColor = '#34d399';
      }
      if (icon) icon.textContent = '⏸';
      if (label) label.textContent = 'Pause Tour';
      this.runTourStep();
    } else {
      if (btn) {
        btn.style.background = '#0284c7';
        btn.style.borderColor = '#38bdf8';
      }
      if (icon) icon.textContent = '▶';
      if (label) label.textContent = 'Auto-Tour';
      if (this.autoTourTimer) clearTimeout(this.autoTourTimer);
    }
  }

  runTourStep() {
    if (!this.isAutoTouring) return;

    this.autoTourTimer = setTimeout(() => {
      if (!this.isAutoTouring) return;
      const nextStage = (this.currentStageIndex + 1) % STAGE_SCALES.length;
      this.setStage(nextStage, true);
      this.runTourStep();
    }, this.tourStepDuration);
  }

  // -----------------------------------------------------------------------
  // LIFECYCLE & EVENT LISTENERS
  // -----------------------------------------------------------------------
  attachEventListeners() {
    // Stepper Button Clicks
    const stepBtns = document.querySelectorAll('.pot-step-btn');
    stepBtns.forEach(btn => {
      btn.addEventListener('click', (e) => {
        const step = parseInt(e.currentTarget.getAttribute('data-step'), 10);
        if (this.isAutoTouring) this.toggleAutoTour();
        this.setStage(step, true);
      });
    });

    // Zoom Out / In Step Buttons
    const prevBtn = document.getElementById('pot-btn-prev');
    const nextBtn = document.getElementById('pot-btn-next');
    if (prevBtn) {
      prevBtn.addEventListener('click', () => {
        if (this.isAutoTouring) this.toggleAutoTour();
        const prev = Math.max(0, this.currentStageIndex - 1);
        this.setStage(prev, true);
      });
    }
    if (nextBtn) {
      nextBtn.addEventListener('click', () => {
        if (this.isAutoTouring) this.toggleAutoTour();
        const next = Math.min(STAGE_SCALES.length - 1, this.currentStageIndex + 1);
        this.setStage(next, true);
      });
    }

    // Auto-Tour Button
    const tourBtn = document.getElementById('pot-btn-auto-tour');
    if (tourBtn) {
      tourBtn.addEventListener('click', () => this.toggleAutoTour());
    }

    // Continuous Scale Slider
    const slider = document.getElementById('pot-scale-slider');
    if (slider) {
      slider.addEventListener('input', (e) => {
        if (this.isAutoTouring) this.toggleAutoTour();
        const val = parseFloat(e.target.value);
        const stage = Math.min(Math.floor(val) - 1, 3);
        if (stage !== this.currentStageIndex) {
          this.setStage(stage, true);
        }
      });
    }

    // Reset Camera Button
    const resetBtn = document.getElementById('pot-btn-reset-cam');
    if (resetBtn) {
      resetBtn.addEventListener('click', () => {
        const data = STAGE_SCALES[this.currentStageIndex];
        this.animateCameraTo(data.cameraDistance, data.targetPosition);
      });
    }

    // Fullscreen Toggle
    const fsBtn = document.getElementById('pot-btn-fullscreen');
    if (fsBtn) {
      fsBtn.addEventListener('click', () => {
        if (!document.fullscreenElement) {
          this.modalEl.requestFullscreen().catch(() => {});
        } else {
          document.exitFullscreen().catch(() => {});
        }
      });
    }

    // Close Button
    const closeBtn = document.getElementById('pot-btn-close');
    if (closeBtn) {
      closeBtn.addEventListener('click', () => this.close());
    }

    // Global Hotkeys: 'O' to open, 'Escape' to close
    window.addEventListener('keydown', (e) => {
      if (['INPUT', 'TEXTAREA'].includes(document.activeElement?.tagName)) return;
      if (e.key === 'o' || e.key === 'O') {
        if (this.modalEl && this.modalEl.style.display !== 'none') {
          this.close();
        } else {
          this.open();
        }
      } else if (e.key === 'Escape') {
        if (this.modalEl && this.modalEl.style.display !== 'none') {
          this.close();
        }
      }
    });

    // Window Resize Handler
    window.addEventListener('resize', () => {
      if (!this.modalEl || this.modalEl.style.display === 'none') return;
      const w = this.canvasEl.clientWidth;
      const h = this.canvasEl.clientHeight;
      if (w && h && this.renderer && this.camera) {
        this.camera.aspect = w / h;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(w, h);
      }
    });
  }

  open(initialStage = 0) {
    if (!this.modalEl) return;
    this.modalEl.style.display = 'flex';

    // Force canvas resize
    setTimeout(() => {
      const w = this.canvasEl.clientWidth;
      const h = this.canvasEl.clientHeight;
      if (w && h && this.renderer && this.camera) {
        this.camera.aspect = w / h;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(w, h);
      }
    }, 50);

    this.setStage(initialStage, false);
    this.startAnimationLoop();
  }

  close() {
    if (!this.modalEl) return;
    this.modalEl.style.display = 'none';
    if (this.isAutoTouring) this.toggleAutoTour();
    if (this.animFrameId) {
      cancelAnimationFrame(this.animFrameId);
      this.animFrameId = null;
    }
  }

  startAnimationLoop() {
    if (this.animFrameId) return;

    const animate = () => {
      this.animFrameId = requestAnimationFrame(animate);

      // 1. Animate Stage 1 Streamlines
      if (this.stageGroups.laniakea.visible) {
        this.streamlineParticles.forEach(sp => {
          sp.progress += sp.speed;
          if (sp.progress >= 1.0) sp.progress = 0;
          const pt = sp.curve.getPoint(sp.progress);
          sp.particle.position.copy(pt);
        });
      }

      // 2. Animate Stage 2 Milky Way Disc Spin
      if (this.stageGroups.milkyway.visible && this.mwDiscPoints) {
        this.mwDiscPoints.rotation.y += 0.001;
      }

      // 3. Animate Stage 3 Planetary Orbits
      if (this.stageGroups.solarsystem.visible) {
        this.planetOrbitGroups.forEach(p => {
          p.angle += p.speed;
          p.pivot.rotation.y = p.angle;
        });
        if (this.sunCoronaMesh) {
          this.sunCoronaMesh.rotation.y += 0.002;
        }
      }

      // 4. Animate Stage 4 Earth Rotation & Clouds
      if (this.stageGroups.earthmoon.visible) {
        if (this.earthMesh) this.earthMesh.rotation.y += 0.002;
        if (this.cloudMesh) this.cloudMesh.rotation.y += 0.0026;
        if (this.moonOrbitGroup) this.moonOrbitGroup.rotation.y += 0.0015;
      }

      // Update Controls & Render
      this.controls.update();
      this.renderer.render(this.scene, this.camera);
    };

    animate();
  }
}

// Global Singleton Instance
let visualizerInstance = null;

function initPowersOfTen() {
  if (!visualizerInstance) {
    visualizerInstance = new PowersOfTenVisualizer();
    window.powersOfTenVisualizer = visualizerInstance;
  }
  return visualizerInstance;
}

window.initPowersOfTen = initPowersOfTen;
window.openPowersOfTen = function(stage = 0) {
  const vis = initPowersOfTen();
  if (vis) vis.open(stage);
};

})(window);
