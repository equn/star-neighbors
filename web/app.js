/* The Sun's Changing Neighborhood
 *
 * Every star's heliocentric offset is stored as a Chebyshev polynomial in time,
 * fitted to an orbit integrated through a Milky Way potential. The coefficients
 * live in a floating-point texture and the polynomial is evaluated in the vertex
 * shader, so the whole catalog moves on real dynamics at frame rate without
 * any per-frame work on the CPU.
 */
'use strict';

const LY = 3.261563777;
/* Model constants are read from meta.json at boot so the client cannot silently
 * disagree with the data it was handed; these are only fallbacks. The GLSL and
 * the Chebyshev scratch buffer are both built after that read. */
let TSPAN = 5.0;            // Myr, the integration and fit window
let DEG = 12;
let NC = DEG + 1;

const S = {
  n: 0, meta: null, d: null, coeffs: null,
  t: 0, play: 0,        // 0 paused, +1 forwards, -1 backwards
  lastDir: 1,           // direction Space resumes in
  sel: -1, hover: -1,
  cam: { yaw: 0.6, pitch: 0.32, dist: 46, target: [0, 0, 0] },
  // "From Earth": the camera stands at the Sun and looks outward, and stars
  // are drawn at their apparent magnitude instead of their distance.
  sky: { on: false, fov: 1.2, magLim: 6.5, restore: null },
  filt: { hideC: true, onlyFeat: false, radius: 20, maxSpread: 2.0 },
  tab: 'near', query: '',
  visible: null, pulse: -1e9,
  dirty: true,          // something changed; redraw this frame
  preFocus: null,      // view to restore when the selection is cleared
  camDirty: false,     // set once the user moves the camera themselves
};

/* ---------------------------------------------------------------- utilities */
const $ = (id) => document.getElementById(id);
const fmtT = (myr) => {
  const y = myr * 1e6;
  if (Math.abs(y) < 1) return 'today';
  const s = y < 0 ? '−' : '+';
  const a = Math.abs(y);
  if (a < 1e5) return `${s}${Math.round(a / 100) / 10} kyr`;
  if (a < 1e6) return `${s}${Math.round(a / 1000)} kyr`;
  return `${s}${(a / 1e6).toFixed(3)} Myr`;
};
const fmtD = (pc) => `${(pc * LY).toFixed(2)} ly`;

/* Chebyshev basis at u = t / TSPAN, evaluated once per frame on the CPU. */
let basis = new Float32Array(NC);
function chebBasis(t) {
  if (basis.length !== NC) basis = new Float32Array(NC);
  const u = Math.max(-1, Math.min(1, t / TSPAN));
  basis[0] = 1; basis[1] = u;
  for (let k = 2; k < NC; k++) basis[k] = 2 * u * basis[k - 1] - basis[k - 2];
  return basis;
}
/* CPU evaluation for one star (used by the chart and picking). */
function posAt(i, t) {
  const b = chebBasis(t), o = i * 3 * NC;
  let x = 0, y = 0, z = 0;
  for (let k = 0; k < NC; k++) {
    const w = b[k];
    x += S.coeffs[o + k] * w;
    y += S.coeffs[o + NC + k] * w;
    z += S.coeffs[o + 2 * NC + k] * w;
  }
  return [x, y, z];
}
function distAt(i, t) { const p = posAt(i, t); return Math.hypot(p[0], p[1], p[2]); }

/* Apparent magnitude from the Sun, m = M + 5 log10(d/10 pc). The distance is
 * whatever the trajectory says it is at time t, so this is what the star looks
 * like from Earth at that moment; 99 means the catalog has no photometry. */
function appMag(i, t) {
  const M = S.d.abs_g[i];
  if (M === null || M === undefined || M > 90) return 99;
  return M + 5 * Math.log10(Math.max(distAt(i, t), 1e-6)) - 5;
}
/* The n brightest stars right now, brightest first. Recomputed per frame while
 * time is moving, which is 12k cheap polynomial evaluations - the same work the
 * GPU is doing anyway. */
function skyBrightest(n, t = S.t) {
  const out = [];
  for (let i = 0; i < S.n; i++) {
    if (!S.visible[i]) continue;
    const m = appMag(i, t);
    if (m > S.sky.magLim) continue;
    out.push({ i, m });
  }
  out.sort((a, b) => a.m - b.m);
  return out.slice(0, n);
}

/* ------------------------------------------------------------------- webgl */
const canvas = $('gl');
/* The 3-D view needs WebGL2 (float textures + gl_VertexID). Everything else -
 * the distance-time chart, the encounter tables, search, the Monte Carlo
 * figures - is Canvas 2-D and DOM, so a machine without WebGL2 still gets a
 * working tool rather than a blank page. */
const gl = canvas.getContext('webgl2', { antialias: true, alpha: false });
const HAS_GL = !!gl;

const VERT = () => `#version 300 es
precision highp float;
precision highp int;
layout(location=0) in vec3 a_color;
layout(location=1) in float a_size;
layout(location=3) in float a_absg;   // absolute G; 99 where there is none
layout(location=4) in float a_vis;    // 1 if the filters keep this star

uniform sampler2D u_coeffs;      // NC wide, n tall; rgb = x,y,z coefficient
uniform float u_basis[${NC}];
uniform mat4 u_vp;
uniform float u_pxScale;
uniform float u_radius;
uniform int u_sel;
uniform int u_sky;               // 1 = standing on Earth looking out
uniform float u_magLim;          // faintest apparent magnitude drawn
uniform float u_skyScale;

out vec3 v_color;
out float v_alpha;
flat out int v_id;

void main() {
  int id = gl_VertexID;
  v_id = id;
  vec3 p = vec3(0.0);
  for (int k = 0; k < ${NC}; k++) {
    p += texelFetch(u_coeffs, ivec2(k, id), 0).rgb * u_basis[k];
  }

  float d = length(p);

  // The star you explicitly selected is always drawn. Culling it was how
  // "Focus" could fly the camera to a star that was not on screen: a third of
  // the shipped catalog lies beyond any reachable view radius.
  bool sel = (id == u_sel);
  // Whether the filters keep this star is decided once, on the CPU, and
  // uploaded. Restating the rules here is what let "Max spread" apply to the
  // lists and not to the view.
  bool hidden = a_vis < 0.5;

  // Apparent magnitude as seen from the Sun at this instant. The distance is
  // the one the trajectory is already reporting, so this costs nothing: the
  // same Chebyshev evaluation that places the star also brightens it.
  float m = a_absg + 5.0 * log(max(d, 1e-6)) * 0.4342944819 - 5.0;
  bool noPhot = a_absg > 90.0;

  bool cull = u_sky == 1
    ? (!sel && (noPhot || m > u_magLim || hidden))
    : (!sel && (d > u_radius || hidden));

  gl_Position = u_vp * vec4(p, 1.0);
  if (cull || gl_Position.w <= 0.0) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); }

  float s;
  if (u_sky == 1) {
    // On the sky a star is a point source: its size on the page stands for
    // brightness, not for anything physical, so it comes from the magnitude
    // and not from the distance the perspective divide would give.
    float b = u_magLim - m;                 // magnitudes above the limit
    s = (1.5 + b * 1.3) * u_skyScale;
    v_alpha = noPhot ? 0.0 : clamp(b * 0.45, 0.16, 1.0);
  } else {
    // apparent size: brighter (smaller a_size index) and nearer stars draw larger
    s = a_size * u_pxScale / max(gl_Position.w, 0.15);
    v_alpha = sel ? 1.0 : clamp(1.35 - d / u_radius * 0.45, 0.42, 1.0);
  }
  if (id == u_sel) { s *= 2.2; v_alpha = 1.0; }
  // The floor is in device pixels, so on a 2x display 2.0 was under one CSS
  // pixel - too few samples for the star's color to be anything but a guess
  // at the far end of the field. 2.8 is still a dot, with enough of it to
  // carry a hue.
  gl_PointSize = clamp(s, 2.8, 64.0);

  v_color = a_color;
}`;

const FRAG = () => `#version 300 es
precision highp float;
precision highp int;
in vec3 v_color;
in float v_alpha;
out vec4 outColor;
void main() {
  vec2 c = gl_PointCoord * 2.0 - 1.0;
  float r = dot(c, c);
  if (r > 1.0) discard;
  // bright core with a soft halo, so crowded regions bloom naturally
  float core = exp(-r * 5.5);
  float halo = exp(-r * 1.7) * 0.35;
  // core + halo peaks at 1.35, and the output is premultiplied (col * a), so
  // an unclamped a drove all three channels past 1.0 at the center of every
  // star and the write clipped them to pure white before blending - which is
  // why the middle of a point was white however carefully col was computed.
  float a = min((core + halo) * v_alpha, 1.0);

  // Lean the hot center toward white rather than summing with it. Adding a
  // constant to all three channels cannot coexist with hue here, because
  // teff_to_rgb pins red at 255 for everything cooler than 6600 K: the old
  // core * 0.55 term drove green and blue up to meet it and the whole point
  // turned white on the 2-4 px discs this view mostly draws. Ross 248 is
  // rgb(255,142,29) in the data and was reaching the screen as
  // rgb(255,255,134). Cubing the core also confines the whitening to the
  // middle ~30% of the radius, so only the very center reads as hot.
  vec3 col = mix(v_color, vec3(1.0), core * core * core * 0.28);

  // No selection tint. This used to be a 55% mix toward blue, which repainted
  // the one star the interface is asking you to look at - the tour narrates
  // Gliese 710 as "distinctly orange" over a point rendered rgb(220,243,250).
  // The reticle, the name label, a 2.2x size boost and full alpha already say
  // which star is selected, and none of them lie about its temperature.
  outColor = vec4(col * a, a);
}`;

const LVERT = `#version 300 es
precision highp float;
layout(location=0) in vec3 a_pos;
uniform mat4 u_vp;
// Screen-space nudge, in normalized device units. gl.lineWidth() is capped at
// 1 on essentially every desktop GL driver - this context reports a maximum of
// exactly 1 - so the only way to draw a guide ring thicker than one physical
// pixel is to draw it several times, shifted by a fraction of a pixel. The
// shift has to happen after projection and scale with w, or it would grow and
// shrink with distance instead of staying a constant width on screen.
uniform vec2 u_nudge;
void main() {
  vec4 p = u_vp * vec4(a_pos, 1.0);
  p.xy += u_nudge * p.w;
  gl_Position = p;
}`;

const LFRAG = `#version 300 es
precision highp float;
uniform vec4 u_col;
out vec4 outColor;
void main() { outColor = u_col; }`;

function compile(src, type) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src); gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}
function program(vs, fs) {
  const p = gl.createProgram();
  gl.attachShader(p, compile(vs, gl.VERTEX_SHADER));
  gl.attachShader(p, compile(fs, gl.FRAGMENT_SHADER));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}

let progStars, progLines, vaoStars, vaoLines, texCoeffs, uni = {}, lineCount = 0;
let oortStart = 0, oortCount = 0;
let vaoTrail, bufTrail, trailCount = 0;
let bufVis = null;

/* Push the CPU's visibility decision to the GPU.
 *
 * The shader used to re-derive it from the flag bits, duplicating the logic in
 * recomputeVisible(). The two drifted, as duplicated logic does: the CPU also
 * applies the "Max spread" filter and the shader never knew about it, so
 * dragging that slider removed 3,587 stars from the list, the chart and the
 * sky labels while changing the 3-D view by exactly zero pixels. Uploading the
 * result instead of restating the rules means a filter cannot apply to one and
 * not the other, including any filter added later. */
function uploadVis() {
  // Called from both initGL (buffer exists, filters not yet run) and
  // recomputeVisible (the reverse, before the GL context is up on a reload),
  // so both halves have to be present before there is anything to send.
  if (!HAS_GL || !bufVis || !S.visible) return;
  const v = new Float32Array(S.n);
  for (let i = 0; i < S.n; i++) v[i] = S.visible[i];
  gl.bindBuffer(gl.ARRAY_BUFFER, bufVis);
  gl.bufferSubData(gl.ARRAY_BUFFER, 0, v);
}
const TRAIL_N = 320;
// Data flag bits 1..128 come from meta.json; this one is added by the renderer
// and must stay clear of them.

function initGL() {
  if (!HAS_GL) return;
  progStars = program(VERT(), FRAG());
  progLines = program(LVERT, LFRAG);
  for (const k of ['u_coeffs', 'u_basis', 'u_vp', 'u_pxScale', 'u_radius',
                   'u_sel',
                   'u_sky', 'u_magLim', 'u_skyScale']) {
    uni[k] = gl.getUniformLocation(progStars, k);
  }
  uni.l_vp = gl.getUniformLocation(progLines, 'u_vp');
  uni.l_col = gl.getUniformLocation(progLines, 'u_col');
  uni.l_nudge = gl.getUniformLocation(progLines, 'u_nudge');

  // --- coefficient texture: NC wide, n tall, RGBA32F -------------------------
  const tex = new Float32Array(S.n * NC * 4);
  for (let i = 0; i < S.n; i++) {
    const src = i * 3 * NC, dst = i * NC * 4;
    for (let k = 0; k < NC; k++) {
      tex[dst + k * 4 + 0] = S.coeffs[src + k];
      tex[dst + k * 4 + 1] = S.coeffs[src + NC + k];
      tex[dst + k * 4 + 2] = S.coeffs[src + 2 * NC + k];
      tex[dst + k * 4 + 3] = 0;
    }
  }
  texCoeffs = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, texCoeffs);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, NC, S.n, 0, gl.RGBA, gl.FLOAT, tex);

  // --- per-star attributes ---------------------------------------------------
  const col = new Float32Array(S.n * 3);
  const size = new Float32Array(S.n);
  const absg = new Float32Array(S.n);
  for (let i = 0; i < S.n; i++) {
    col[i * 3] = S.d.color[i * 3] / 255;
    col[i * 3 + 1] = S.d.color[i * 3 + 1] / 255;
    col[i * 3 + 2] = S.d.color[i * 3 + 2] / 255;
    const M = S.d.abs_g[i];
    // absolute magnitude -> a rough physical size; giants get visibly bigger
    size[i] = (M > 90) ? 1.4 : Math.max(0.55, Math.min(7.0, 2.9 - 0.26 * M));
    absg[i] = (M === null || M === undefined) ? 99 : M;
  }
  vaoStars = gl.createVertexArray();
  gl.bindVertexArray(vaoStars);
  const mk = (data, loc, n) => {
    const b = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, b);
    gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, n, gl.FLOAT, false, 0, 0);
  };
  mk(col, 0, 3); mk(size, 1, 1); mk(absg, 3, 1);
  // Visibility is uploaded rather than re-derived in the shader; see uploadVis.
  bufVis = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, bufVis);
  gl.bufferData(gl.ARRAY_BUFFER, S.n * 4, gl.DYNAMIC_DRAW);
  gl.enableVertexAttribArray(4);
  gl.vertexAttribPointer(4, 1, gl.FLOAT, false, 0, 0);
  gl.bindVertexArray(null);
  uploadVis();

  buildGuides();

  // trail buffer, rewritten whenever the selection changes
  vaoTrail = gl.createVertexArray();
  gl.bindVertexArray(vaoTrail);
  bufTrail = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, bufTrail);
  gl.bufferData(gl.ARRAY_BUFFER, (TRAIL_N + 1) * 3 * 4, gl.DYNAMIC_DRAW);
  gl.enableVertexAttribArray(0);
  gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
  gl.bindVertexArray(null);
}

/* The path a star actually traces past the Sun. It is visibly curved: over
 * megayears the Galactic tide and differential rotation bend it away from the
 * straight line a naive model would draw. */
function buildTrail(i) {
  if (!HAS_GL || i < 0) { trailCount = 0; return; }
  const v = new Float32Array((TRAIL_N + 1) * 3);
  for (let k = 0; k <= TRAIL_N; k++) {
    const t = -TSPAN + 2 * TSPAN * k / TRAIL_N;
    const p = posAt(i, t);
    v[k * 3] = p[0]; v[k * 3 + 1] = p[1]; v[k * 3 + 2] = p[2];
  }
  gl.bindBuffer(gl.ARRAY_BUFFER, bufTrail);
  gl.bufferSubData(gl.ARRAY_BUFFER, 0, v);
  trailCount = TRAIL_N + 1;
}

/* Reference rings at 5/10/20 pc plus the Galactic axes. */
function buildGuides() {
  const v = [];
  const ring = (r, ax) => {
    const N = 128;
    for (let i = 0; i < N; i++) {
      for (const j of [i, (i + 1) % N]) {
        const a = j / N * Math.PI * 2, c = Math.cos(a) * r, s = Math.sin(a) * r;
        v.push(ax === 2 ? c : (ax === 1 ? c : 0), ax === 2 ? s : (ax === 1 ? 0 : c),
               ax === 0 ? s : (ax === 1 ? s : 0));
      }
    }
  };
  for (const r of RING_PC) ring(r, 2);
  // No ray toward the Galactic center. It used to be drawn out to 45 pc with a
  // label on it, which put a marker for something 8,000 pc away inside a map
  // whose widest ring is 40 - and at close zoom the label landed beside the
  // Oort cloud. Orientation belongs to the axis gizmo, which is plainly a
  // compass rather than a map, and says GC without implying a distance.
  lineCount = v.length / 3;
  // The Oort cloud is roughly spherical, so three orthogonal great circles
  // rather than one ring in the plane - a single circle would read as another
  // distance marker instead of a shell you can pass inside of. Kept after the
  // rest of the geometry so it can be drawn separately, in its own color and
  // only when it is worth drawing.
  oortStart = lineCount;
  for (const ax of [0, 1, 2]) ring(OORT_PC, ax);
  oortCount = v.length / 3 - oortStart;
  vaoLines = gl.createVertexArray();
  gl.bindVertexArray(vaoLines);
  const b = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, b);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(v), gl.STATIC_DRAW);
  gl.enableVertexAttribArray(0);
  gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
  gl.bindVertexArray(null);
}

/* ------------------------------------------------------------------ matrix */
function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2), nf = 1 / (near - far);
  return [f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0];
}
function lookAt(eye, ctr, up) {
  const z = norm(sub(eye, ctr)), x = norm(cross(up, z)), y = cross(z, x);
  return [x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0,
          -dot(x, eye), -dot(y, eye), -dot(z, eye), 1];
}
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a) => { const l = Math.hypot(...a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
function mul(a, b) {
  const o = new Array(16);
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
    o[i * 4 + j] = a[j] * b[i * 4] + a[4 + j] * b[i * 4 + 1] + a[8 + j] * b[i * 4 + 2] + a[12 + j] * b[i * 4 + 3];
  }
  return o;
}
function camEye() {
  const { yaw, pitch, dist, target } = S.cam;
  return [target[0] + dist * Math.cos(pitch) * Math.cos(yaw),
          target[1] + dist * Math.cos(pitch) * Math.sin(yaw),
          target[2] + dist * Math.sin(pitch)];
}
/* Unit vector the sky camera is pointing along, from the same yaw/pitch the
 * orbiting camera uses, so switching modes keeps you facing the same way. */
function skyDir() {
  const { yaw, pitch } = S.cam;
  return [Math.cos(pitch) * Math.cos(yaw),
          Math.cos(pitch) * Math.sin(yaw),
          Math.sin(pitch)];
}
function viewProj() {
  const asp = canvas.width / canvas.height;
  if (S.sky.on) {
    // Standing at the Sun. The near plane has to clear a star mid-encounter:
    // Gliese 710 comes inside 0.06 pc and would otherwise clip away at the
    // exact moment it is worth looking at.
    const d = skyDir();
    return mul(perspective(S.sky.fov, asp, 1e-4, 4000),
               lookAt([0, 0, 0], d, [0, 0, 1]));
  }
  return mul(perspective(Math.PI / 4, asp, 0.05, 4000), lookAt(camEye(), S.cam.target, [0, 0, 1]));
}

/* -------------------------------------------------------------------- draw */
let vp = null;
function render() {
  if (!HAS_GL) return;
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const w = Math.floor(canvas.clientWidth * dpr), h = Math.floor(canvas.clientHeight * dpr);
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
  gl.viewport(0, 0, w, h);
  gl.clearColor(0.02, 0.027, 0.051, 1);
  gl.clear(gl.COLOR_BUFFER_BIT);
  gl.disable(gl.DEPTH_TEST);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE);      // additive: stars glow through each other

  vp = viewProj();

  // The distance rings and the trajectory trail are scaffolding for the map.
  // From inside the scene they would be a cage around the observer.
  if (!S.sky.on) {
  gl.useProgram(progLines);
  gl.uniformMatrix4fv(uni.l_vp, false, vp);
  gl.bindVertexArray(vaoLines);

  // Draw each guide four times, offset by fractions of a physical pixel, to
  // build a line about two pixels wide out of one-pixel primitives. Blending
  // is additive, so the passes also brighten it: the rings used to be a
  // 0.13-alpha blue that vanished on any screen with a bit of glare.
  const PASSES = [[0, 0], [0.8, 0], [0, 0.8], [0.8, 0.8]];
  const thick = (col, first, count) => {
    gl.uniform4f(uni.l_col, col[0], col[1], col[2], col[3]);
    for (const [dx, dy] of PASSES) {
      gl.uniform2f(uni.l_nudge, dx * 2 / w, dy * 2 / h);
      gl.drawArrays(gl.LINES, first, count);
    }
    gl.uniform2f(uni.l_nudge, 0, 0);
  };

  thick([0.74, 0.82, 0.95, 0.16], 0, lineCount);          // distance rings

  if (oortCount && viewAcross() <= OORT_SHOW_PC) {
    thick([0.62, 0.85, 1.0, 0.13], oortStart, oortCount); // Oort shell
  }

  if (trailCount) {
    gl.uniform4f(uni.l_col, 0.44, 0.76, 1.0, 0.55);
    gl.bindVertexArray(vaoTrail);
    gl.drawArrays(gl.LINE_STRIP, 0, trailCount);
    gl.bindVertexArray(vaoLines);
  }
  }

  gl.useProgram(progStars);
  gl.activeTexture(gl.TEXTURE0);
  gl.bindTexture(gl.TEXTURE_2D, texCoeffs);
  gl.uniform1i(uni.u_coeffs, 0);
  gl.uniform1fv(uni.u_basis, chebBasis(S.t));
  gl.uniformMatrix4fv(uni.u_vp, false, vp);
  // convert a star's nominal size into pixels through the actual projection,
  // so points scale correctly with zoom and display density
  const focal = 0.5 * h / Math.tan(Math.PI / 8);
  gl.uniform1f(uni.u_pxScale, focal * 0.055);
  gl.uniform1f(uni.u_radius, S.filt.radius);
  gl.uniform1i(uni.u_sel, S.sel);
  gl.uniform1i(uni.u_sky, S.sky.on ? 1 : 0);
  gl.uniform1f(uni.u_magLim, S.sky.magLim);
  // Narrowing the field magnifies the sky but not a point source, so grow the
  // dots a little as you zoom in or the view empties out.
  gl.uniform1f(uni.u_skyScale, dpr * Math.pow(1.2 / S.sky.fov, 0.35));
  gl.bindVertexArray(vaoStars);
  gl.drawArrays(gl.POINTS, 0, S.n);
  gl.bindVertexArray(null);

  drawSunMarker();
}

/* The Sun is drawn in 2-D over the scene so it always reads as the origin. */
function drawSunMarker() {
  if (!HAS_GL) return;
  // In the sky view the Sun is behind you, not in front of you.
  if (S.sky.on) { sunEl.style.display = 'none'; return; }
  const o = project([0, 0, 0]);
  if (!o) { sunEl.style.display = 'none'; return; }
  // on narrow layouts the panels cover most of the canvas; a marker floating
  // over the star list reads as a bug rather than as the Sun
  for (const id of ['rail', 'dock', 'info', 'title', 'nav', 'about', 'tourbar']) {
    const el = $(id);
    if (!el || el.hidden || getComputedStyle(el).display === 'none') continue;
    const r = el.getBoundingClientRect();
    if (o[0] >= r.left - 6 && o[0] <= r.right + 6
        && o[1] >= r.top - 6 && o[1] <= r.bottom + 6) {
      sunEl.style.display = 'none';
      return;
    }
  }
  sunEl.style.display = 'block';
  sunEl.style.left = `${o[0]}px`;
  sunEl.style.top = `${o[1]}px`;
}
function project(p) {
  const v = vp;
  const x = v[0] * p[0] + v[4] * p[1] + v[8] * p[2] + v[12];
  const y = v[1] * p[0] + v[5] * p[1] + v[9] * p[2] + v[13];
  const w = v[3] * p[0] + v[7] * p[1] + v[11] * p[2] + v[15];
  if (w <= 0) return null;
  return [(x / w * 0.5 + 0.5) * canvas.clientWidth, (1 - (y / w * 0.5 + 0.5)) * canvas.clientHeight];
}

/* ------------------------------------------------------- orientation & zoom */
/* Once the view has been tumbled it is easy to lose track of which way the
 * Galaxy lies, so the camera's orientation is always shown explicitly against
 * the three Galactic axes, and the plane rings are labeled with their radius
 * to keep a sense of scale. */

const overlay = $('overlay');
const ox = overlay.getContext('2d');

const PULSE_SECS = 1.5;                 // selection ring animation
const DIST_MIN = 0.6, DIST_MAX = 300;
const DIST_RATIO = DIST_MAX / DIST_MIN;
const distToU = (d) => 1 - Math.log(d / DIST_MIN) / Math.log(DIST_RATIO);
const uToDist = (u) => DIST_MIN * Math.pow(DIST_RATIO, 1 - u);
const FOV = Math.PI / 4;
/* width of the view at the Sun's distance - the most intuitive scale readout */
const viewAcross = () => 2 * S.cam.dist * Math.tan(FOV / 2);

const RING_PC = [5, 10, 20, 40];
/* The Oort cloud has no edge, only a consensus about where it thins out. 1.5 ly
 * is the round number for the outer boundary (about 95,000 AU) and is what the
 * chart shades, so both read it from here rather than each carrying a literal.
 * Drawn only when the view is tight enough for it to mean anything: at the
 * default 124 ly across it would be four pixels of noise. */
const OORT_LY = 1.5;
const OORT_PC = OORT_LY / LY;
const OORT_SHOW_PC = 12;                // show once the view is this wide or less

/* Camera axes in world space. */
function camBasis() {
  const eye = camEye();
  const fwd = norm(sub(eye, S.cam.target));      // points back toward the eye
  const right = norm(cross([0, 0, 1], fwd));
  const up = cross(fwd, right);
  return { right, up, fwd };
}

/* --- smooth moves between viewpoints, so presets do not teleport ---------- */
/* The CSS media query only reaches CSS animation; the camera fly-to and the
 * selection pulse are drawn on canvas and ignored it entirely. Someone who asks
 * for reduced motion should not get a 700 ms swoop across the sky. */
const REDUCED_MOTION = matchMedia('(prefers-reduced-motion: reduce)');
let camAnim = null;
function flyTo(yaw, pitch, dist, target = [0, 0, 0], ms = 720) {
  if (REDUCED_MOTION.matches) {          // arrive, do not travel
    camAnim = null;
    S.cam.yaw = yaw; S.cam.pitch = pitch; S.cam.dist = dist;
    S.cam.target = target.slice();
    markCamDirty(); invalidate();
    return;
  }
  // take the short way round in yaw
  let dy = ((yaw - S.cam.yaw + Math.PI) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI) - Math.PI;
  camAnim = {
    t0: performance.now(), ms,
    from: { yaw: S.cam.yaw, pitch: S.cam.pitch, dist: S.cam.dist, target: S.cam.target.slice() },
    to: { yaw: S.cam.yaw + dy, pitch, dist, target },
  };
}
function stepCamAnim(now) {
  if (!camAnim) return;
  invalidate();
  const k = Math.min(1, (now - camAnim.t0) / camAnim.ms);
  const e = k < 0.5 ? 4 * k * k * k : 1 - Math.pow(-2 * k + 2, 3) / 2;   // ease in-out
  const { from, to } = camAnim;
  S.cam.yaw = from.yaw + (to.yaw - from.yaw) * e;
  S.cam.pitch = from.pitch + (to.pitch - from.pitch) * e;
  S.cam.dist = from.dist * Math.pow(to.dist / from.dist, e);             // log-lerp
  for (let a = 0; a < 3; a++) {
    S.cam.target[a] = from.target[a] + (to.target[a] - from.target[a]) * e;
  }
  if (k >= 1) camAnim = null;
}

/* Frame a star together with the Sun.
 *
 * Alpha Centauri A and B sit 7 screen pixels apart in the default view, so
 * selecting one gives almost no visual signal: they read as a single dot and
 * the info card is in the far corner. Flying the camera in and marking the star
 * where it actually is turns "which of those specks did I just click" into
 * something you can see. The look-at point becomes the midpoint of the Sun-star
 * line and the camera swings perpendicular to it, so both ends stay in frame.
 */
function snapshotCam() {
  return { yaw: S.cam.yaw, pitch: S.cam.pitch, dist: S.cam.dist,
           target: S.cam.target.slice(), radius: S.filt.radius };
}
/* Any hand-driven camera move means the stored "before" view is no longer what
 * the user wants to come back to, so stop offering to restore it. */
function markCamDirty() { S.camDirty = true; invalidate(); }

/* The slider was linear 2-60 pc, but the bundle ships stars out to ~970 pc -
 * 4,164 of them beyond 60 - so a third of the catalog was unreachable and
 * Focus silently failed on it. Logarithmic, like the zoom control, so the
 * crowded 2-20 pc range keeps its resolution while the far field is reachable. */
const RAD_MIN = 2, RAD_MAX = 1000;
const radToU = (pc) => Math.log(pc / RAD_MIN) / Math.log(RAD_MAX / RAD_MIN);
const uToRad = (u) => Math.round(RAD_MIN * Math.pow(RAD_MAX / RAD_MIN, u));
const fmtRadius = (pc) => (pc >= 1000 ? `${(pc / 1000).toFixed(2)} kpc` : `${pc} pc`);

/* Switching between the map and the view from inside it. The orbiting camera's
 * distance is meaningless once you are standing at the Sun, so it is put aside
 * and handed back on the way out; the heading is kept, so you leave facing the
 * way you were looking. */
function setSky(on) {
  if (on === S.sky.on) return;
  camAnim = null;
  if (on) {
    S.sky.restore = { dist: S.cam.dist, target: S.cam.target.slice() };
    S.cam.target = [0, 0, 0];
  } else if (S.sky.restore) {
    S.cam.dist = S.sky.restore.dist;
    S.cam.target = S.sky.restore.target;
    S.sky.restore = null;
  }
  S.sky.on = on;
  const b = $('vSky');
  b.setAttribute('aria-pressed', String(on));
  b.classList.toggle('on', on);
  document.body.classList.toggle('skymode', on);
  buildTrail(S.sky.on ? -1 : S.sel);
  invalidate();
}

const SKY_FOV_MIN = 0.20, SKY_FOV_MAX = 2.10;      // radians, ~11 to ~120 deg
const fovToU = (f) => 1 - Math.log(f / SKY_FOV_MIN)
                          / Math.log(SKY_FOV_MAX / SKY_FOV_MIN);
const uToFov = (u) => SKY_FOV_MIN * Math.pow(SKY_FOV_MAX / SKY_FOV_MIN, 1 - u);
function setFov(f) {
  if (!Number.isFinite(f)) return;          // see setT: NaN here is permanent
  S.sky.fov = Math.max(SKY_FOV_MIN, Math.min(SKY_FOV_MAX, f));
  camAnim = null;
  invalidate();
}

function setRadius(pc) {
  S.filt.radius = pc;
  $('frad').value = radToU(pc);
  $('fradv').textContent = fmtRadius(pc);
}

/* A tour shot has to be framed into the part of the canvas the caption is not
 * covering, and that part is a different size on every layout: the succession
 * beat's text is three lines at 840 px wide and eight lines on a phone, which
 * leaves 63% of the height in one case and 43% in the other. So measure it
 * rather than assume it. Returns the lift as a fraction of the view height -
 * enough to move the pair's center to the middle of what remains visible -
 * and the extra distance that keeps both ends of the pair inside it. */
const VIEW_H_PER_DIST = 2 * Math.tan(Math.PI / 8);   // 45-degree vertical field
function tourFraming() {
  const h = canvas.clientHeight;
  const bar = $('tourbar');
  if (!document.body.classList.contains('touring') || !bar || bar.hidden || !h) {
    return { lift: 0, pull: 1 };
  }
  // The clear strip runs from under the title bar to the top of the caption.
  // Both ends matter: centering on the caption alone over-lifted, and put the
  // Sun a pixel above the title bar in the widest shot on a phone.
  const y0 = canvas.getBoundingClientRect().top;
  const head = $('title');
  const hi = head ? head.getBoundingClientRect().bottom - y0 : 0;
  const lo = bar.getBoundingClientRect().top - y0;
  const f = Math.max(0.3, Math.min(1, (lo - hi) / h));   // usable fraction
  // The pair's height on screen goes as 1/dist, so the pull has to go as 1/f,
  // not fall off linearly: a linear version was fine at the 63% a wide window
  // leaves but far too weak at the 30% left by a ten-line caption on a phone,
  // where the Sun ended up above the title bar. 0.85 is set so that the widest
  // shot in the tour still fits, and reproduces the 1.35x that reads well on
  // desktop.
  return { lift: 0.5 - (hi + lo) / (2 * h), pull: Math.max(1, Math.min(2.8, 0.85 / f)) };
}

function focusStar(i, ms = 760) {
  if (i < 0 || !HAS_GL) return;
  const p0 = posAt(i, S.t);
  // From the Sun there is nowhere to fly to: focusing means turning to face it.
  if (S.sky.on) {
    if (!S.preFocus) S.preFocus = snapshotCam();
    S.camDirty = false;
    const r0 = Math.hypot(p0[0], p0[1], p0[2]) || 1;
    flyTo(Math.atan2(p0[1], p0[0]),
          Math.max(-1.5, Math.min(1.5, Math.asin(p0[2] / r0))),
          S.cam.dist, [0, 0, 0], ms);
    S.pulse = performance.now();
    return;
  }
  // remember where we came from, but only the first time in a browsing run
  if (!S.preFocus) S.preFocus = snapshotCam();
  S.camDirty = false;
  const p = posAt(i, S.t);
  const r = Math.hypot(p[0], p[1], p[2]);

  // a star outside the view radius would be culled the moment we arrived
  if (r > S.filt.radius) {
    setRadius(Math.min(RAD_MAX, Math.ceil(r * 1.15)));
    refilter();
  }

  const yaw = Math.atan2(p[1], p[0]) + Math.PI / 2;   // look across the pair
  const pitch = 0.30;

  // Framing centered on the canvas is wrong during the tour, because the
  // caption covers the bottom of it: measured over the succession beat,
  // Ross 248 sat behind the text in 48 of 84 samples and settled at y=669 on
  // a 930 px canvas whose caption starts at 587. Pull back and lower the
  // look-at point, which lifts the pair into the part of the frame the viewer
  // can actually see. The eye follows the target, so dropping the target in z
  // translates the whole rig down and the subject rises on screen; dividing by
  // cos(pitch) accounts for the camera not being level with the world z axis.
  const fr = tourFraming();
  const dist = Math.max(1.1, Math.min(300, (r * 1.9 + 0.5) * fr.pull));
  const mid = [p[0] / 2, p[1] / 2, p[2] / 2];
  mid[2] -= fr.lift * VIEW_H_PER_DIST * dist / Math.cos(pitch);
  flyTo(yaw, pitch, dist, mid, ms);
  S.pulse = performance.now();
}

/* --- axis gizmo ----------------------------------------------------------- */
const GIZMO_AXES = [
  { v: [1, 0, 0], label: 'GC', color: '#ffc46b', major: true },
  { v: [-1, 0, 0], label: 'AC', color: '#8a7a5e', major: false },
  { v: [0, 1, 0], label: 'ROT', color: '#8fe6a0', major: true },
  { v: [0, -1, 0], label: '', color: '#5c7d63', major: false },
  { v: [0, 0, 1], label: 'NGP', color: '#7fcdff', major: true },
  { v: [0, 0, -1], label: 'SGP', color: '#5b7e96', major: false },
];

function drawGizmo() {
  if (!HAS_GL) return;
  const el = $('gizmo');
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const w = el.clientWidth || 84, h = el.clientHeight || 84;
  const bw = Math.floor(w * dpr), bh = Math.floor(h * dpr);   // see drawOverlay
  if (el.width !== bw || el.height !== bh) { el.width = bw; el.height = bh; }
  const g = el.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);

  const cxp = w / 2, cyp = h / 2, R = Math.min(w, h) / 2 - 13;
  const { right, up, fwd } = camBasis();

  g.strokeStyle = 'rgba(140,170,220,.16)';
  g.beginPath(); g.arc(cxp, cyp, R + 5, 0, 7); g.stroke();

  // draw far axes first so the near ones sit on top
  const drawn = GIZMO_AXES.map((a) => ({
    ...a,
    x: cxp + dot(a.v, right) * R,
    y: cyp - dot(a.v, up) * R,
    depth: dot(a.v, fwd),
  })).sort((p, q) => p.depth - q.depth);

  g.lineCap = 'round';
  for (const a of drawn) {
    const front = a.depth > -0.05;
    g.globalAlpha = a.major ? (front ? 1 : 0.45) : (front ? 0.5 : 0.25);
    g.strokeStyle = a.color;
    g.lineWidth = a.major ? 1.8 : 1.2;
    g.beginPath(); g.moveTo(cxp, cyp); g.lineTo(a.x, a.y); g.stroke();
    g.fillStyle = a.color;
    g.beginPath(); g.arc(a.x, a.y, a.major ? 3 : 2, 0, 7); g.fill();
    if (a.label) {
      g.font = a.major ? '600 8.5px ui-monospace, monospace' : '8px ui-monospace, monospace';
      g.textAlign = 'center'; g.textBaseline = 'middle';
      // push the label outward from the center so it clears the dot
      const len = Math.hypot(a.x - cxp, a.y - cyp) || 1;
      g.fillText(a.label, cxp + (a.x - cxp) / len * (len + 8),
                          cyp + (a.y - cyp) / len * (len + 8));
    }
  }
  g.globalAlpha = 1;
  g.fillStyle = 'rgba(255,233,168,.95)';
  g.beginPath(); g.arc(cxp, cyp, 2.4, 0, 7); g.fill();
}

/* --- in-scene labels: ring radii and the Galactic-center axis -------------- */
/* Labels are placed wherever the screen is actually clear. A fixed offset angle
 * puts them behind the side panels as soon as the camera turns, so each ring is
 * sampled all the way round and the label goes to the visible point nearest the
 * open middle of the view. */
const PANEL_IDS = ['title', 'rail', 'dock', 'nav', 'info', 'skybar', 'tourbar'];

function panelRects() {
  const out = [];
  for (const id of PANEL_IDS) {
    const el = $(id);
    if (!el || getComputedStyle(el).display === 'none') continue;
    const r = el.getBoundingClientRect();
    if (r.width > 0 && r.height > 0) out.push(r);
  }
  return out;
}
function clearOf(x, y, rects, halfW = 0, pad = 8) {
  return !rects.some((r) => x + halfW >= r.left - pad && x - halfW <= r.right + pad
                         && y >= r.top - pad && y <= r.bottom + pad);
}
/* Clearance test for center-aligned text, using its measured width. */
function clearText(x, y, text, rects) {
  return clearOf(x, y, rects, ox.measureText(text).width / 2 + 2);
}
/* Center of whatever screen space the panels leave over. */
function freeAnchor(rects, w, h) {
  const rail = rects.find((r) => r.left < w * 0.5 && r.height > h * 0.3);
  const dock = rects.find((r) => r.top > h * 0.6);
  const left = rail ? rail.right : 0;
  const bottom = dock ? dock.top : h;
  return [(left + w) / 2, bottom / 2];
}

function drawOverlay() {
  if (!HAS_GL) return;                 // nothing is projected without a camera
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const w = overlay.clientWidth, h = overlay.clientHeight;
  // Both dimensions decide this, not just the width. Entering full screen from
  // an already-maximised window keeps the width and only grows the height - the
  // menu bar, tab strip and window chrome go away - so a width-only test left
  // this canvas at its old backing height while the GL canvas beside it
  // resized. The browser then stretched a 1400-tall buffer over an 1680-tall
  // box, and every label and reticle drifted downward in proportion to its own
  // y: 85 px at mid-screen, which is what put the marker below its star.
  // Floor to match the GL canvas exactly; a fractional assignment would be
  // truncated and then never compare equal, reallocating every frame.
  const bw = Math.floor(w * dpr), bh = Math.floor(h * dpr);
  if (overlay.width !== bw || overlay.height !== bh) { overlay.width = bw; overlay.height = bh; }
  ox.setTransform(dpr, 0, 0, dpr, 0, 0);
  ox.clearRect(0, 0, w, h);
  ox.font = '9.5px ui-monospace, monospace';
  ox.textAlign = 'center'; ox.textBaseline = 'middle';

  const rects = panelRects();
  const [ax, ay] = freeAnchor(rects, w, h);
  const onScreen = (p) => p && p[0] > 20 && p[0] < w - 20 && p[1] > 14 && p[1] < h - 14;

  // Ring radii. Sample the whole ring and label the clear point nearest the
  // anchor. Seen edge-on every ring projects onto the same line, so labels also
  // have to dodge each other or they stack into an unreadable pile.
  const sunPt = project([0, 0, 0]);
  const placed = [];
  const freeOfPlaced = (x, y, hw) => !placed.some(
    (q) => Math.abs(q.y - y) < 12 && Math.abs(q.x - x) < q.hw + hw + 6);

  const N = 48;
  // The Oort shell is labeled by the same rule as the distance rings, and is
  // in the list only while it is actually drawn.
  const guides = S.sky.on ? [] : RING_PC.map((r) => ({ r, label: `${r} pc` }));
  if (!S.sky.on && oortCount && viewAcross() <= OORT_SHOW_PC) {
    guides.push({ r: OORT_PC, label: 'Oort cloud', oort: true });
  }
  for (const g of guides) {
    const r = g.r;
    if (r > S.cam.dist * 3) continue;
    const label = g.label;
    const hw = ox.measureText(label).width / 2 + 2;

    // A ring's radius only reads correctly at its widest point on screen. Seen
    // edge-on the ring collapses to a line and every other point is
    // foreshortened, which would put "20 pc" inside "5 pc". So find the widest
    // extent across the whole ring first - including points that are off-screen
    // or behind a panel - and only then look for a usable spot near it. If the
    // ring's true extent is not visible, it gets no label at all.
    const pts = [];
    let maxR = 0;
    for (let k = 0; k < N; k++) {
      const a = k / N * Math.PI * 2;
      const p = project([r * Math.cos(a), r * Math.sin(a), 0]);
      const sr = (p && sunPt) ? Math.hypot(p[0] - sunPt[0], p[1] - sunPt[1]) : -1;
      if (sr > maxR) maxR = sr;
      pts.push({ p, sr });
    }
    let best = null, bestD = Infinity;
    for (const c of pts) {
      if (c.sr < maxR * 0.92) continue;                  // must read as the radius
      if (!onScreen(c.p) || !clearText(c.p[0], c.p[1], label, rects)) continue;
      if (!freeOfPlaced(c.p[0], c.p[1], hw)) continue;
      const d = Math.hypot(c.p[0] - ax, c.p[1] - ay);
      if (d < bestD) { bestD = d; best = c.p; }
    }
    if (!best) continue;
    placed.push({ x: best[0], y: best[1], hw });
    ox.fillStyle = g.oort ? 'rgba(111,195,255,.62)' : 'rgba(150,180,225,.55)';
    ox.fillText(label, best[0], best[1]);
  }

  // In the sky view the useful labels are the stars you can actually see, so
  // name the brightest handful in view rather than the geometry.
  if (S.sky.on) {
    for (const b of skyBrightest(14)) {
      if (b.i === S.sel) continue;          // the reticle already names it
      const p = project(posAt(b.i, S.t));
      if (!onScreen(p)) continue;
      const y = p[1] + 13;
      const label = S.d.names[b.i];
      const hw = ox.measureText(label).width / 2 + 2;
      if (!clearText(p[0], y, label, rects) || !freeOfPlaced(p[0], y, hw)) continue;
      placed.push({ x: p[0], y, hw });
      ox.fillStyle = b.m < 1 ? 'rgba(226,238,255,.95)' : 'rgba(174,196,232,.7)';
      ox.fillText(label, p[0], y);
    }
    drawSelectionMarker(rects);
    return;
  }

  // Nothing here labels the Galactic center. It sits 8,000 pc away and the
  // widest ring on this map is 40, so any in-scene marker for it is off by two
  // orders of magnitude no matter where it is placed. The axis gizmo carries
  // the direction instead: it reads as a compass, not as a position.
  const sp = sunPt;

  if (sp && onScreen(sp) && clearText(sp[0], sp[1] + 15, 'Sun', rects)
      && freeOfPlaced(sp[0], sp[1] + 15, 11)) {
    ox.fillStyle = 'rgba(255,233,168,.9)';
    ox.fillText('Sun', sp[0], sp[1] + 15);
  }

  drawSelectionMarker(rects);
}

/* Reticle on the selected star: a short expanding pulse to catch the eye, then
 * a steady ring so it stays identifiable in a crowded field. */
function drawSelectionMarker(rects) {
  const i = S.sel;
  if (i < 0 || !S.visible[i]) return;
  const p = posAt(i, S.t);   // no radius test: the selection is always marked
  const q = project(p);
  const w = overlay.clientWidth, h = overlay.clientHeight;
  if (!q || q[0] < -40 || q[0] > w + 40 || q[1] < -40 || q[1] > h + 40) return;

  // expanding pulse; the steady reticle below it still marks the selection
  const age = (performance.now() - (S.pulse || -1e9)) / 1000;
  if (!REDUCED_MOTION.matches && age >= 0 && age < PULSE_SECS) {
    for (let k = 0; k < 3; k++) {
      const a = age - k * 0.22;
      if (a <= 0) continue;
      const f = Math.min(1, a / PULSE_SECS);
      ox.beginPath();
      ox.arc(q[0], q[1], 9 + f * 62, 0, 7);
      ox.strokeStyle = `rgba(111,195,255,${0.55 * (1 - f) ** 2})`;
      ox.lineWidth = 2 * (1 - f) + 0.5;
      ox.stroke();
    }
  }

  // steady reticle
  const R = 13;
  ox.strokeStyle = 'rgba(111,195,255,.9)';
  ox.lineWidth = 1.4;
  ox.beginPath(); ox.arc(q[0], q[1], R, 0, 7); ox.stroke();
  ox.beginPath();
  for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
    ox.moveTo(q[0] + dx * (R + 3), q[1] + dy * (R + 3));
    ox.lineTo(q[0] + dx * (R + 7), q[1] + dy * (R + 7));
  }
  ox.stroke();

  const label = S.d.names[i];
  ox.font = '600 11px system-ui';
  const hw = ox.measureText(label).width / 2 + 4;
  let ly = q[1] - R - 9;
  if (!clearOf(q[0], ly, rects, hw)) ly = q[1] + R + 14;
  ox.fillStyle = 'rgba(9,14,24,.75)';
  ox.fillRect(q[0] - hw, ly - 8, hw * 2, 15);
  ox.fillStyle = '#8fd4ff';
  ox.textAlign = 'center'; ox.textBaseline = 'middle';
  ox.fillText(label, q[0], ly);
}

const sunEl = document.createElement('div');
sunEl.style.cssText = `position:absolute;width:14px;height:14px;margin:-7px 0 0 -7px;border-radius:50%;
  background:radial-gradient(circle,#fff 0%,#ffe9a8 45%,rgba(255,200,90,0) 72%);
  pointer-events:none;display:none;z-index:2;`;

/* ------------------------------------------------------------------ picking */
function pick(mx, my) {
  if (!HAS_GL) return -1;
  let best = -1, bestD = 18;
  const t = S.t;
  for (let i = 0; i < S.n; i++) {
    if (!S.visible[i]) continue;
    const p = posAt(i, t);
    if (Math.hypot(p[0], p[1], p[2]) > S.filt.radius) continue;
    const s = project(p);
    if (!s) continue;
    const d = Math.hypot(s[0] - mx, s[1] - my);
    if (d < bestD) { bestD = d; best = i; }
  }
  return best;
}
/* Relative width of the resampled range. Grade and this are different things:
 * a star can have flawless astrometry (grade A) and still be extrapolated so far
 * that its range is wider than the prediction itself. */
function spread(i) {
  const med = S.d.mc_dmin_med[i];
  if (!(med > 0)) return 0;
  return (S.d.mc_dmin_hi[i] - S.d.mc_dmin_lo[i]) / med;
}

function recomputeVisible() {
  S.visible = new Uint8Array(S.n);
  for (let i = 0; i < S.n; i++) {
    if (S.filt.hideC && S.d.grade[i] === 'C') continue;
    if (S.filt.onlyFeat && !(S.d.flags[i] & 16)) continue;
    if (S.filt.maxSpread < 2.0 && spread(i) > S.filt.maxSpread) continue;
    S.visible[i] = 1;
  }
  // Inside the function rather than beside its callers, so a filter added
  // later cannot reach the lists without also reaching the view.
  uploadVis();
}

/* -------------------------------------------------------------------- chart */
/* Both reference charts use linear axes over +/-80 kyr. Over +/-5 Myr that
 * collapses every famous encounter into a spike at the origin, so the time axis
 * is symmetric-log (linear within +/-10 kyr, logarithmic outside) and distance
 * is logarithmic. Gliese 710 at 0.17 ly and Ross 154 at 6.4 ly then both read
 * clearly on one plot, which no linear pair of axes can do. */
const chart = $('chart');
const cx = chart.getContext('2d');
let chartStars = [];
let labelSet = new Set();      // the subset of chartStars whose names are drawn
let chartHover = -1;
let timeDrag = false, hoverCursor = false;

/* Time <-> screen mapping. symlog() takes Myr to [-1,1]; invSymlog() comes back.
 * Everything that positions itself along the time axis - the chart cursor, the
 * transport slider, playback - works in these units, so the slider thumb and the
 * cursor line always sit at the same place. A slider linear in Myr would crawl
 * for four megayears and then cross the entire interesting region in a pixel. */
const T0 = 0.01;                       // Myr; linear region half-width
/* Recomputed once meta.json supplies the real window, like TSPAN itself. */
let symK = Math.log10(1 + TSPAN / T0);
const symlog = (t) => Math.sign(t) * Math.log10(1 + Math.abs(t) / T0) / symK;
const invSymlog = (u) => Math.sign(u) * T0 * (Math.pow(10, Math.abs(u) * symK) - 1);

const D_LO = 0.06, D_HI = 34;          // ly, log distance range
const logd = (ly) => (Math.log10(Math.max(ly, D_LO)) - Math.log10(D_LO))
                   / (Math.log10(D_HI) - Math.log10(D_LO));

const T_TICKS_FULL = [-5, -3, -1, -0.3, -0.1, -0.03, -0.01, 0,
                      0.01, 0.03, 0.1, 0.3, 1, 3, 5];
const T_TICKS_SPARSE = [-5, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 5];
const T_TICKS_TINY = [-5, -1, -0.1, 0, 0.1, 1, 5];
/* Fifteen labels need roughly 700 px of plot; below that they collide into an
 * unreadable smear, so drop to every other one, and on a phone drop the
 * ten-thousand-year pair as well - "−10k now +10k" ran together. */
const timeTicks = (pw) => (pw < 420 ? T_TICKS_TINY
                         : pw < 700 ? T_TICKS_SPARSE : T_TICKS_FULL);
const D_TICKS = [0.1, 0.3, 1, 3, 10, 30];
const tickLabel = (t) => {
  if (t === 0) return 'now';
  const a = Math.abs(t), s = t < 0 ? '−' : '+';
  return a < 1 ? `${s}${Math.round(a * 1000)}k` : `${s}${a}M`;
};

/* Survey designations run to 28 characters, which no plot label can carry
 * without swamping the curve it names. Selected and hovered stars keep the full
 * name; the rest are truncated, and the full name is a hover or a list row away. */
function chartLabel(i, full) {
  const n = S.d.names[i];
  return (full || n.length <= 20) ? n : n.slice(0, 19) + '…';
}

function pickChartStars() {
  // notable stars that actually have an encounter inside the window, plus the
  // closest approaches overall, so the plot stays legible instead of hairy
  const feat = [], other = [];
  for (let i = 0; i < S.n; i++) {
    if (!S.visible[i]) continue;
    if (isCensored(i) || Math.abs(rankT(i)) > TSPAN * 0.999) continue;
    (S.d.flags[i] & 16 ? feat : other).push(i);
  }
  feat.sort((a, b) => rankD(a) - rankD(b));
  other.sort((a, b) => rankD(a) - rankD(b));
  // a narrow plot cannot carry 36 labeled curves; show the closest few
  const narrow = chart.clientWidth < 700;
  const set = new Set([...feat.slice(0, narrow ? 6 : 18),
                       ...other.slice(0, narrow ? 3 : 10)]);
  if (S.sel >= 0) set.add(S.sel);
  chartStars = [...set];
  // Curves are cheap to read as texture; names are not. Two dozen of them
  // collide into a smear around the present epoch, so only the closest few get
  // labeled, plus whatever the reader is pointing at.
  const budget = [...set].filter((i) => !isCensored(i))
                         .sort((a, b) => rankD(a) - rankD(b))
                         .slice(0, narrow ? 3 : 6);
  labelSet = new Set(budget);   // the selected and hovered stars are added when drawn
}

function chartGeom() {
  const w = chart.clientWidth, h = chart.clientHeight;
  const L = 46, R = 14, T = 12, B = 22;   // L/R mirror --plot-l/--plot-r in the CSS
  return { w, h, L, R, T, B, pw: w - L - R, ph: h - T - B };
}
function chartX(t) { const g = chartGeom(); return g.L + (symlog(t) * 0.5 + 0.5) * g.pw; }
/* Screen x back to a time, clamped to the plotted range. */
function chartT(px) {
  const g = chartGeom();
  return Math.max(-TSPAN, Math.min(TSPAN,
    invSymlog(Math.max(-1, Math.min(1, ((px - g.L) / g.pw) * 2 - 1)))));
}
const SCRUB_H = 22;            // grab strip along the top of the plot
function onScrubBand(mx, my) {
  const g = chartGeom();
  return mx >= g.L - 8 && mx <= g.L + g.pw + 8 && my >= g.T - 10 && my <= g.T + SCRUB_H;
}
function onCursor(mx, my) {
  const g = chartGeom();
  return my >= g.T - 10 && my <= g.T + g.ph && Math.abs(mx - chartX(S.t)) <= 10;
}
function chartY(ly) { const g = chartGeom(); return g.T + g.ph - logd(ly) * g.ph; }

function drawChart() {
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const g = chartGeom();
  const cbw = Math.floor(g.w * dpr), cbh = Math.floor(g.h * dpr);   // see drawOverlay
  if (chart.width !== cbw || chart.height !== cbh) { chart.width = cbw; chart.height = cbh; }
  cx.setTransform(dpr, 0, 0, dpr, 0, 0);
  cx.clearRect(0, 0, g.w, g.h);

  const X = chartX, Y = chartY;

  // --- grid ---------------------------------------------------------------
  cx.font = '10px ui-monospace, monospace';
  cx.strokeStyle = 'rgba(140,170,220,.10)'; cx.lineWidth = 1;
  cx.fillStyle = '#5b6880'; cx.textAlign = 'right';
  for (const d of D_TICKS) {
    const y = Y(d);
    cx.beginPath(); cx.moveTo(g.L, y); cx.lineTo(g.L + g.pw, y); cx.stroke();
    cx.fillText(d < 1 ? `${d}` : `${d} ly`, g.L - 6, y + 3);
  }
  cx.textAlign = 'center';
  for (const t of timeTicks(g.pw)) {
    const x = X(t);
    cx.strokeStyle = t === 0 ? 'rgba(140,170,220,.28)' : 'rgba(140,170,220,.10)';
    cx.beginPath(); cx.moveTo(x, g.T); cx.lineTo(x, g.T + g.ph); cx.stroke();
    cx.fillStyle = t === 0 ? '#8595ad' : '#5b6880';
    cx.fillText(tickLabel(t), x, g.h - 7);
  }
  if (g.pw >= 420) {          // otherwise it lands on top of the +5M tick
    cx.fillStyle = '#5b6880'; cx.textAlign = 'left';
    cx.fillText('yr', g.L + g.pw + 2, g.h - 7);
  }

  // --- Oort cloud band ----------------------------------------------------
  const oort = OORT_LY;                    // ly, outer edge; shared with the 3-D shell
  cx.fillStyle = 'rgba(111,195,255,.09)';
  cx.fillRect(g.L, Y(oort), g.pw, g.T + g.ph - Y(oort));
  cx.strokeStyle = 'rgba(111,195,255,.28)';
  cx.beginPath(); cx.moveTo(g.L, Y(oort)); cx.lineTo(g.L + g.pw, Y(oort)); cx.stroke();
  cx.fillStyle = 'rgba(111,195,255,.5)'; cx.font = '9.5px system-ui';
  cx.textAlign = 'left';
  cx.fillText('outer Oort cloud', g.L + 5, Y(oort) - 4);

  // --- curves -------------------------------------------------------------
  const NS = 300;
  const labels = [];

  /* Sample times for one star's curve: the uniform screen-space grid, plus a
   * geometric ladder of offsets around that star's own perihelion.
   *
   * The uniform grid alone cannot draw a close encounter that happens far from
   * the present. The time axis is symlog, so a megayear out one grid step spans
   * thousands of years, while the approach itself is over in a few hundred:
   * Gliese 710's dip is 0.23 grid steps wide and HD 7977's is 0.03, so the
   * polyline stepped straight over both. Gliese 710 bottomed out at 1.24 ly
   * against a true 0.17, and HD 7977 at 2.69 against 0.07 - a curve that never
   * went near the perihelion dot drawn on top of it.
   *
   * Refining by looking at neighboring samples would not have found them
   * either: both endpoints sit at the same height on either side of the spike.
   * The perihelion time is in the data, so sample there directly. A ladder from
   * a year to 100 kyr resolves the V whatever its width, for 42 extra points. */
  const curveTimes = (i) => {
    const ts = [];
    for (let k = 0; k <= NS; k++) ts.push(invSymlog(k / NS * 2 - 1));
    const tp = S.d.tmin_myr[i];
    if (Number.isFinite(tp) && Math.abs(tp) < TSPAN) {
      ts.push(tp);
      for (let e = -6; e <= -1 + 1e-9; e += 0.25) {
        const dt = Math.pow(10, e);
        if (tp - dt >= -TSPAN) ts.push(tp - dt);
        if (tp + dt <= TSPAN) ts.push(tp + dt);
      }
      ts.sort((a, b) => a - b);
    }
    return ts;
  };
  for (const i of chartStars) {
    const isSel = i === S.sel, isHov = i === chartHover;
    const feat = S.d.flags[i] & 16;
    // Named curves lead; the rest are context and stay out of their way.
    const named = labelSet.has(i);
    const alpha = named ? 0.85 : (feat ? 0.42 : 0.26);
    cx.strokeStyle = (isSel || isHov) ? '#6fc3ff'
      : `rgba(${S.d.color[i*3]},${S.d.color[i*3+1]},${S.d.color[i*3+2]},${alpha})`;
    cx.lineWidth = (isSel || isHov) ? 2.2 : 1;
    cx.beginPath();
    let pen = false;
    for (const t of curveTimes(i)) {
      const ly = distAt(i, t) * LY;
      if (ly > D_HI) { pen = false; continue; }          // leave the plot, break
      const x = X(t), y = Y(ly);
      pen ? cx.lineTo(x, y) : cx.moveTo(x, y);
      pen = true;
    }
    cx.stroke();

    // Mark the perihelion the interface actually quotes - the Monte Carlo
    // median - not the nominal orbit's, which can sit somewhere else entirely.
    const tm = rankT(i), dm = rankD(i) * LY;
    if (isSel || isHov) {
      // whisker: the 16th-84th percentile of the closest approach under
      // resampled astrometry. This is the one number on the chart with an
      // honest error bar, so it is drawn rather than left in the text panel.
      const y1 = Y(Math.max(S.d.mc_dmin_lo[i] * LY, D_LO));
      const y2 = Y(Math.max(S.d.mc_dmin_hi[i] * LY, D_LO));
      const x0 = X(tm), lo = Math.min(y1, y2), hi = Math.max(y1, y2);
      cx.strokeStyle = 'rgba(111,195,255,.75)'; cx.lineWidth = 1.4;
      cx.beginPath();
      cx.moveTo(x0, lo); cx.lineTo(x0, hi);
      cx.moveTo(x0 - 4, lo); cx.lineTo(x0 + 4, lo);
      cx.moveTo(x0 - 4, hi); cx.lineTo(x0 + 4, hi);
      cx.stroke();
    }
    if (Math.abs(tm) <= TSPAN && dm < D_HI && (isSel || isHov || named || feat)) {
      if (isSel || isHov || named) {
        labels.push({ i, x: X(tm), y: Y(dm), sel: isSel || isHov });
      }
      cx.fillStyle = (isSel || isHov) ? '#6fc3ff'
        : (named ? 'rgba(210,225,245,.9)' : 'rgba(190,205,228,.45)');
      cx.beginPath(); cx.arc(X(tm), Y(dm), isSel ? 3.4 : 2, 0, 7); cx.fill();
    }
  }

  // --- labels, nudged apart so they stay readable -------------------------
  labels.sort((a, b) => a.y - b.y);
  const placed = [];
  cx.textAlign = 'center';
  for (const L of labels) {
    cx.font = L.sel ? '600 11px system-ui' : '10px system-ui';
    const text = chartLabel(L.i, L.sel);
    const half = cx.measureText(text).width / 2 + 5;
    // Keep the text inside the plot: a 27-character Gaia designation sitting on
    // a +3 Myr encounter used to run off the right edge and be unreadable.
    const x = Math.max(g.L + half, Math.min(g.L + g.pw - half, L.x));
    let y = L.y - 9, dir = -1;
    for (let guard = 0; guard < 60; guard++) {
      const clash = placed.some((p) =>
        Math.abs(p.y - y) < 11 && Math.abs(p.x - x) < p.half + half);
      if (!clash) break;
      y += dir * 11;
      if (y < g.T + 9) { dir = 1; y = L.y + 15; }
      if (y > g.T + g.ph - 3) { y = g.T + g.ph - 3; break; }
    }
    y = Math.max(g.T + 9, Math.min(g.T + g.ph - 3, y));
    placed.push({ x, y, half });
    // leader line back to the point, when the label had to move
    if (Math.abs(y - (L.y - 9)) > 6 || Math.abs(x - L.x) > 2) {
      cx.strokeStyle = 'rgba(160,185,215,.28)'; cx.lineWidth = 1;
      cx.beginPath(); cx.moveTo(L.x, L.y + (y > L.y ? 4 : -4));
      cx.lineTo(x, y + (y > L.y ? -8 : 3)); cx.stroke();
    }
    cx.fillStyle = L.sel ? '#6fc3ff' : 'rgba(205,220,242,.85)';
    cx.fillText(text, x, y);
  }

  // --- time cursor, with a grab handle -------------------------------------
  const cxp = X(S.t);
  const hot = timeDrag || hoverCursor;
  cx.fillStyle = hot ? 'rgba(111,195,255,.10)' : 'rgba(140,170,220,.05)';
  cx.fillRect(g.L, g.T - 4, g.pw, SCRUB_H + 2);

  cx.strokeStyle = hot ? 'rgba(150,215,255,.95)' : 'rgba(255,255,255,.55)';
  cx.lineWidth = hot ? 1.6 : 1;
  cx.beginPath(); cx.moveTo(cxp, g.T); cx.lineTo(cxp, g.T + g.ph); cx.stroke();

  cx.fillStyle = hot ? '#8fd4ff' : 'rgba(226,236,250,.85)';
  cx.beginPath();                                   // downward-pointing grip
  cx.moveTo(cxp - 6, g.T - 3); cx.lineTo(cxp + 6, g.T - 3);
  cx.lineTo(cxp, g.T + 8); cx.closePath(); cx.fill();

  if (hot) {
    const lab = fmtT(S.t);
    cx.font = '600 10px ui-monospace, monospace';
    const w = cx.measureText(lab).width + 10;
    const bx = Math.max(g.L, Math.min(g.L + g.pw - w, cxp - w / 2));
    cx.fillStyle = 'rgba(9,14,24,.92)';
    cx.fillRect(bx, g.T + 21, w, 15);
    cx.strokeStyle = 'rgba(111,195,255,.5)'; cx.lineWidth = 1;
    cx.strokeRect(bx + .5, g.T + 21.5, w - 1, 14);
    cx.fillStyle = '#8fd4ff'; cx.textAlign = 'center';
    cx.fillText(lab, bx + w / 2, g.T + 31.5);
  }
}

/* --------------------------------------------------------------------- list */
function listRows() {
  const q = S.query.toLowerCase();
  let idx = [];
  for (let i = 0; i < S.n; i++) {
    if (!S.visible[i]) continue;
    if (q && !S.d.names[i].toLowerCase().includes(q)) continue;
    idx.push(i);
  }
  if (S.tab === 'near') idx.sort((a, b) => distAt(a, S.t) - distAt(b, S.t));
  else {
    // Split on the same epoch the row displays. Using the nominal epoch here
    // while showing the Monte Carlo one put eight stars in the tab opposite to
    // the date printed beside them.
    const sign = S.tab === 'future' ? 1 : -1;
    idx = idx.filter((i) => sign * rankT(i) > 0);
    idx.sort((a, b) => rankD(a) - rankD(b));
  }
  return idx.slice(0, 160);
}

/* Rank on the Monte Carlo median rather than the nominal perihelion. For a star
 * 200 pc away the nominal value can land outside its own 68% interval, so
 * sorting by it promotes statistical flukes over real encounters. */
const rankD = (i) => (S.d.mc_dmin_med[i] > 0 ? S.d.mc_dmin_med[i] : S.d.dmin_pc[i]);
/* Epoch from the same integrated draw population as the distance. Pairing a
 * resampled distance with the nominal orbit's epoch mixed two models. */
const rankT = (i) => (S.d.mc_tmin_med ? S.d.mc_tmin_med[i] : S.d.tmin_myr[i]);
/* More than half this star's draws are still closing at the window edge, so it
 * has no perihelion inside the window - only a bound. */
const isCensored = (i) => !!(S.d.flags[i] & 256);
/* Which edge matters. At +5 Myr the star is still closing and its encounter is
 * later; at -5 Myr it was already receding when the window opened, so the
 * encounter happened earlier. Calling both "still approaching" - as this did -
 * is simply wrong at the past boundary. */
/* Returns 0 when the draws straddle both boundaries: some put the encounter
 * before the window, some after, and the measurements cannot say which. Those
 * get neutral wording rather than a side picked for them. */
const censEdge = (i) => {
  const t = rankT(i);
  return Math.abs(t) < TSPAN * 0.99 ? 0 : (t < 0 ? -1 : 1);
};
const censPhrase = (i) => {
  const e = censEdge(i);
  return e === 0 ? `closest approach lies outside the ±${TSPAN} Myr window · which side is unresolved`
       : e < 0 ? `already receding at −${TSPAN} Myr · minimum lies earlier`
       : `still approaching at +${TSPAN} Myr · minimum lies later`;
};
/* Terse form, for the detail panel's narrow value column and for tooltips. */
const censShort = (i) => {
  const e = censEdge(i);
  return e === 0 ? `outside ±${TSPAN} Myr · side unresolved`
       : e < 0 ? `earlier than −${TSPAN} Myr` : `later than +${TSPAN} Myr`;
};
/* How much of the doubt is the boundary rather than the measurement. */
const censPct = (i) => Math.round((S.d.mc_censored ? S.d.mc_censored[i] : 1) * 100);

const gradeColor = { A: '#7ee787', B: '#ffb454', C: '#ff6b6b' };
/* Built as DOM rather than an HTML string: the rows carry catalog names from
 * SIMBAD, which is external data, and the rows also need to be real focusable
 * controls rather than decorated divs. */
function renderList() {
  const rows = listRows();
  const el = $('list');
  if (!rows.length) {
    // A blank panel reads as a broken app. Say which control emptied it.
    const p = document.createElement('p');
    p.className = 'empty';
    p.textContent = S.query
      ? `No star matching “${S.query}”. Names come from SIMBAD, so try a `
        + `catalog designation — Gliese 710, Ross 248, HD 7977.`
      : (S.tab === 'near'
         ? 'No star passes the current filters. Widen the view radius, raise '
           + 'Max spread, or untick the quality filters.'
         : `No star has its closest approach ${S.tab === 'future' ? 'ahead' : 'behind'} `
           + 'us under the current filters.');
    el.replaceChildren(p);
    el.setAttribute('aria-busy', 'false');
    return;
  }
  const frag = document.createDocumentFragment();
  for (const i of rows) frag.appendChild(buildRow(i));
  el.replaceChildren(frag);
  // Roving tabindex: one stop for the whole list, arrows move within it. With
  // tabindex=0 on every row, tabbing past the list took 160 presses.
  const first = el.querySelector('.row.sel') || el.firstElementChild;
  for (const r of el.children) r.tabIndex = r === first ? 0 : -1;
  el.setAttribute('aria-busy', 'false');
}

/* One row of the star list, as real DOM: it carries a catalog name from
 * SIMBAD and needs to be a focusable control, not a decorated div. */
function buildRow(i) {
  const near = S.tab === 'near';
  const v = near ? fmtD(distAt(i, S.t))
                 : (isCensored(i) ? '\u2264 ' : '') + fmtD(rankD(i));
  const cens = isCensored(i);
  const sub = near
    ? (cens ? censPhrase(i)
            : `closest ${fmtD(rankD(i))} at ${fmtT(rankT(i))}`)
    : (cens ? `${censPhrase(i)} \u00b7 now ${fmtD(S.d.dist_now_pc[i])}`
            : `${fmtT(rankT(i))} \u00b7 now ${fmtD(S.d.dist_now_pc[i])}`
              + ` \u00b7 \u00b1${((S.d.mc_dmin_hi[i] - S.d.mc_dmin_lo[i]) / 2 * LY).toFixed(2)} ly`);
  const row = document.createElement('div');
  row.className = 'row' + (i === S.sel ? ' sel' : '');
  row.dataset.i = String(i);
  row.setAttribute('role', 'option');
  row.tabIndex = -1;                 // renderList promotes exactly one to 0
  row.setAttribute('aria-selected', String(i === S.sel));
  row.setAttribute('aria-label',
    `${S.d.names[i]}, ${v}, ${sub}, data quality grade ${S.d.grade[i]}`);

  const nm = document.createElement('div');
  nm.className = 'nm';
  const dot = document.createElement('span');
  dot.className = 'dot';
  dot.style.background = gradeColor[S.d.grade[i]] || '#888';
  nm.append(dot, document.createTextNode(S.d.names[i]));

  const vl = document.createElement('div');
  vl.className = 'vl';
  vl.textContent = v;

  const sb = document.createElement('div');
  sb.className = 'sub';
  sb.textContent = sub;

  row.append(nm, vl, sb);
  return row;
}

/* --------------------------------------------------------------------- info */
const FLAGS = [
  [1, 'white dwarf', ''],
  [2, 'unreliable radial velocity', 'bad'],
  [4, 'probable unresolved binary', 'warn'],
  [8, 'RV from literature', ''],
  [16, 'notable', 'good'],
  [32, 'binary RV contamination', 'warn'],
  [64, 'barycentric system velocity used', 'good'],
  [128, 'spectroscopic binary \u2014 RV is one orbital phase', 'warn'],
  [256, 'closest approach lies outside the \u00b15 Myr window', 'warn'],
  [512, 'Galaxy model matters more than the measurements', 'warn'],
  [1024, 'parallax zero-point matters more than the measurements', 'warn'],
];

/* Catalog names are external data; build the node rather than splicing
 * strings into innerHTML. */
function tipText(tip, title, sub) {
  tip.replaceChildren(document.createTextNode(title));
  const d = document.createElement('div');
  d.className = 't2';
  d.textContent = sub;
  tip.appendChild(d);
}
function showInfo(i) {
  if (i == null || i < 0 || !S.d.names[i]) { $('info').style.display = 'none'; return; }
  const d = S.d;
  $('info').style.display = 'block';
  $('i-name').textContent = d.names[i];
  $('i-sp').textContent = [d.sp_type[i], `grade ${d.grade[i]}`].filter(Boolean).join(' · ');
  const lo = d.mc_dmin_lo[i] * LY, hi = d.mc_dmin_hi[i] * LY;
  const rows = [
    ['distance now', fmtD(d.dist_now_pc[i])],
    ['distance at cursor', fmtD(distAt(i, S.t))],
  ];
  if (isCensored(i)) {
    rows.push(['closest within window', `\u2264 ${fmtD(rankD(i))}`],
              ['…minimum lies', censShort(i)],
              ['draws hitting the edge', `${censPct(i)}%`]);
  } else {
    rows.push(['closest approach', fmtD(rankD(i))],
              ['measurement range', `${lo.toFixed(2)} \u2013 ${hi.toFixed(2)} ly`],
              ['\u2026occurs', fmtT(rankT(i))]);
    // A minority of draws can still run into the boundary; the flag only fires
    // past half, so without this the reader never sees the near misses.
    if (censPct(i) > 0) {
      rows.push(['draws hitting the edge', `${censPct(i)}%`]);
    }
  }
  // How much of this prediction is the Galaxy or the astrometric zero-point
  // rather than the data. Stage 9 measures both per star; neither is part of
  // the interval above and they must not be added to it in quadrature.
  const sysRow = (key, label) => {
    if (!d[key]) return;
    const v = d[key][i] * LY;
    // "\u00b10.00 ly" told the reader nothing; for the nearby stars the answer
    // is that these do not move them at all, and that is worth saying.
    rows.push([label, v < 0.005 ? 'negligible' : `\u00b1${v.toFixed(2)} ly`]);
  };
  sysRow('sys_dmin_pc', 'Galaxy-model shift');
  sysRow('zp_dmin_pc', 'parallax zero-point');
  rows.push(['nominal astrometry', fmtD(d.dmin_pc[i])],
            ['straight-line model', `${fmtD(d.lin_dmin_pc[i])}`],
            ['radial velocity', `${d.rv_kms[i].toFixed(1)} km/s`]);
  if (d.abs_g[i] < 90) rows.push(['absolute mag', d.abs_g[i].toFixed(2)]);
  // Built as nodes rather than markup. Every value here is derived from numbers
  // today, but the catalog is external data and this panel is the one place it
  // is formatted for display, so it follows the same rule as the list rows.
  const dl = document.createDocumentFragment();
  for (const [k, v] of rows) {
    const dt = document.createElement('dt'); dt.textContent = k;
    const dd = document.createElement('dd'); dd.textContent = v;
    dl.append(dt, dd);
  }
  $('i-dl').replaceChildren(dl);
  const f = d.flags[i];
  const chips = FLAGS.filter(([b]) => f & b).map(([, l, c]) => [l, c]);
  const width = d.mc_dmin_hi[i] - d.mc_dmin_lo[i];
  const med = d.mc_dmin_med[i];
  if (width > 0.5 * Math.max(med, 1e-6)) {
    chips.push(['uncertainty rivals the prediction', 'warn']);
  }
  // Only meaningful when both numbers describe an encounter. For a censored
  // star the interval is a set of boundary values, so a nominal orbit falling
  // outside it says nothing about the estimator.
  if (!isCensored(i) && med > 0
      && (d.dmin_pc[i] < d.mc_dmin_lo[i] || d.dmin_pc[i] > d.mc_dmin_hi[i])) {
    chips.push(['nominal value outside its own interval', 'bad']);
  }
  $('i-flags').replaceChildren(...chips.map(([label, cls]) => {
    const el = document.createElement('span');
    el.className = `chip ${cls}`.trim();
    el.textContent = label;
    return el;
  }));
}

/* ------------------------------------------------------------------- events */
function bind() {
  let drag = null;
  canvas.addEventListener('pointerdown', (e) => {
    drag = { x: e.clientX, y: e.clientY, moved: 0 };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener('pointermove', (e) => {
    if (drag) {
      const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      drag.moved += Math.abs(dx) + Math.abs(dy);
      camAnim = null; markCamDirty();
      // Looking around a narrowed field should not whip past everything, so
      // the sky view turns in proportion to how much sky is on screen.
      const k = S.sky.on ? 0.006 * (S.sky.fov / 1.2) : 0.006;
      S.cam.yaw -= dx * k;
      S.cam.pitch = Math.max(-1.5, Math.min(1.5, S.cam.pitch + dy * k));
      drag.x = e.clientX; drag.y = e.clientY;
    } else {
      const r = canvas.getBoundingClientRect();
      const i = pick(e.clientX - r.left, e.clientY - r.top);
      if (S.hover !== i) invalidate();
      S.hover = i;
      const tip = $('tip');
      if (i >= 0) {
        tip.style.display = 'block';
        tip.style.left = `${e.clientX + 14}px`;
        tip.style.top = `${e.clientY + 12}px`;
        tipText(tip, S.d.names[i], isCensored(i)
          ? `${fmtD(distAt(i, S.t))} \u00b7 \u2264 ${fmtD(rankD(i))}, ${censShort(i)}`
          : `${fmtD(distAt(i, S.t))} \u00b7 closest ${fmtD(rankD(i))} ${fmtT(rankT(i))}`);
      } else tip.style.display = 'none';
    }
  });
  canvas.addEventListener('pointerup', (e) => {
    if (drag && drag.moved < 5) {
      const r = canvas.getBoundingClientRect();
      toggleSelect(pick(e.clientX - r.left, e.clientY - r.top));
    }
    drag = null;
  });
  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    if (S.sky.on) { setFov(S.sky.fov * Math.exp(e.deltaY * 0.0012)); return; }
    S.cam.dist = Math.max(DIST_MIN, Math.min(DIST_MAX, S.cam.dist * Math.exp(e.deltaY * 0.0012)));
    camAnim = null; markCamDirty();
  }, { passive: false });

  /* Chart pointer handling. The top strip and the cursor line itself scrub
   * time; anywhere else keeps the hover-a-curve / click-to-select behavior. */
  const chartPt = (e) => {
    const r = chart.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  };

  chart.addEventListener('pointerdown', (e) => {
    const [mx, my] = chartPt(e);
    if (!onScrubBand(mx, my) && !onCursor(mx, my)) return;
    timeDrag = true;
    setPlay(0);                       // grabbing the cursor stops playback
    chart.setPointerCapture(e.pointerId);
    $('tip').style.display = 'none';
    setT(chartT(mx));
    e.preventDefault();
  });

  chart.addEventListener('pointermove', (e) => {
    const [mx, my] = chartPt(e);

    if (timeDrag) { setT(chartT(mx)); return; }

    hoverCursor = onScrubBand(mx, my) || onCursor(mx, my);
    if (hoverCursor) {
      chart.style.cursor = 'ew-resize';
      chartHover = -1;
      $('tip').style.display = 'none';
      return;
    }

    let best = -1, bestD = 9;
    for (const i of chartStars) {
      // walk the drawn curve and keep the closest vertex in screen space
      for (let k = 0; k <= 60; k++) {
        const t = invSymlog(k / 60 * 2 - 1);
        const ly = distAt(i, t) * LY;
        if (ly > D_HI) continue;
        const d = Math.hypot(chartX(t) - mx, chartY(ly) - my);
        if (d < bestD) { bestD = d; best = i; }
      }
    }
    if (chartHover !== best) invalidate();
    chartHover = best;
    chart.style.cursor = best >= 0 ? 'pointer' : 'default';
    const tip = $('tip');
    if (best >= 0) {
      tip.style.display = 'block';
      tip.style.left = `${e.clientX + 14}px`;
      tip.style.top = `${e.clientY - 34}px`;
      tipText(tip, S.d.names[best], isCensored(best)
        ? `\u2264 ${fmtD(rankD(best))} \u00b7 ${censShort(best)}`
        : `closest ${fmtD(rankD(best))} \u00b7 ${fmtT(rankT(best))}`);
    } else tip.style.display = 'none';
  });

  const endDrag = (e) => {
    if (!timeDrag) return;
    timeDrag = false;
    if (chart.hasPointerCapture?.(e.pointerId)) chart.releasePointerCapture(e.pointerId);
  };
  chart.addEventListener('pointerup', endDrag);
  chart.addEventListener('pointercancel', endDrag);
  chart.addEventListener('pointerleave', () => {
    if (timeDrag) return;
    chartHover = -1; hoverCursor = false; $('tip').style.display = 'none';
  });
  // a click that was really a scrub must not also select a star
  chart.addEventListener('click', () => { if (!timeDrag && chartHover >= 0) toggleSelect(chartHover); });

  $('list').addEventListener('click', (e) => {
    const r = e.target.closest('.row');
    if (r) toggleSelect(+r.dataset.i, true);
  });
  $('list').addEventListener('keydown', (e) => {
    const r = e.target.closest?.('.row');
    if (!r) return;
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      toggleSelect(+r.dataset.i, true);
    } else if (e.key === 'ArrowDown' || e.key === 'ArrowUp'
               || e.key === 'Home' || e.key === 'End') {
      e.preventDefault();
      const sib = e.key === 'Home' ? r.parentElement.firstElementChild
                : e.key === 'End' ? r.parentElement.lastElementChild
                : e.key === 'ArrowDown' ? r.nextElementSibling : r.previousElementSibling;
      if (sib) {                       // move the single tab stop with the focus
        r.tabIndex = -1; sib.tabIndex = 0; sib.focus();
      }
    }
  });
  $('i-focus').addEventListener('click', () => { S.camDirty = false; focusStar(S.sel); });
  $('search').addEventListener('input', (e) => { S.query = e.target.value; renderList(); });
  document.querySelectorAll('.tab').forEach((t) => t.addEventListener('click', () => {
    document.querySelectorAll('.tab').forEach((x) => {
      x.classList.remove('on');
      x.setAttribute('aria-selected', 'false');
    });
    t.classList.add('on');
    t.setAttribute('aria-selected', 'true');
    S.tab = t.dataset.tab;
    renderList();
  }));
  $('fq').addEventListener('change', (e) => { S.filt.hideC = e.target.checked; refilter(); });
  $('func').addEventListener('input', (e) => {
    const v = +e.target.value / 100;
    S.filt.maxSpread = v;
    $('funcv').textContent = v >= 2.0 ? 'any' : `${Math.round(v * 100)}%`;
    refilter();
  });
  $('ff').addEventListener('change', (e) => { S.filt.onlyFeat = e.target.checked; refilter(); });
  $('frad').addEventListener('input', (e) => {
    S.filt.radius = uToRad(+e.target.value);
    $('fradv').textContent = fmtRadius(S.filt.radius);
  });
  // --- zoom + orientation --------------------------------------------------
  const setDist = (d) => {
    S.cam.dist = Math.max(DIST_MIN, Math.min(DIST_MAX, d));
    invalidate();
    camAnim = null;                     // a manual zoom cancels any fly-to
    markCamDirty();
  };
  // In the sky view there is no camera distance to change - you are already at
  // the Sun - so the same control narrows and widens the field of view.
  $('zoom').addEventListener('input', (e) => (S.sky.on
    ? setFov(uToFov(+e.target.value)) : setDist(uToDist(+e.target.value))));
  $('zin').addEventListener('click', () => (S.sky.on
    ? setFov(S.sky.fov / 1.3) : setDist(S.cam.dist / 1.4)));
  $('zout').addEventListener('click', () => (S.sky.on
    ? setFov(S.sky.fov * 1.3) : setDist(S.cam.dist * 1.4)));
  // yaw = PI looks from the anticenter toward the Galactic center
  const preset = (yaw, pitch, dist) => () => {
    S.preFocus = null;                  // an explicit view choice replaces it
    flyTo(yaw, pitch, dist);
  };
  $('vLevel').addEventListener('click', preset(Math.PI, 0, S.cam.dist));
  $('vTop').addEventListener('click', preset(Math.PI, 1.40, S.cam.dist));
  $('vReset').addEventListener('click', preset(0.6, 0.32, 46));
  $('vSky').addEventListener('click', () => setSky(!S.sky.on));
  // A platform that cannot do it should not offer it: iOS Safari fullscreens
  // video and nothing else. Removing the control beats a button that does
  // nothing when pressed.
  if (fsSupported()) {
    $('fsBtn').addEventListener('click', toggleFullscreen);
    for (const ev of ['fullscreenchange', 'webkitfullscreenchange']) {
      document.addEventListener(ev, syncFullscreen);
    }
  } else {
    $('fsBtn').remove();
  }
  // The HUD is rewritten every frame, so the handler lives on the container.
  for (const id of ['hud', 'skybar']) {
    $(id).addEventListener('click', (e) => {
      const b = e.target.closest('button.lnk');
      if (b) toggleSelect(+b.dataset.i, true);
    });
  }

  addEventListener('resize', () => {
    placeAboutBtn(); fitMobileChrome(); fitInfoCard(); pickChartStars(); invalidate();
  });
  $('menubtn').addEventListener('click', () => togglePanel('menu'));
  $('chartbtn').addEventListener('click', () => togglePanel('chart'));
  $('filterbtn').addEventListener('click', () => togglePanel('filters'));
  // Picking a star is a request to look at it, so get the sheet out of the way.
  $('list').addEventListener('click', () => { if (isNarrow()) setPanel('menu', false); });
  // So is touching the sky.
  canvas.addEventListener('pointerdown', () => { if (isNarrow()) setPanel('menu', false); });
  $('tourBtn').addEventListener('click', tourStart);
  $('tourexit').addEventListener('click', tourStop);
  $('tourplay').addEventListener('click', () => tourPause(!tour.paused));
  $('tournext').addEventListener('click', () => tourGo(tour.i + 1));
  $('tourprev').addEventListener('click', () => tourGo(Math.max(0, tour.i - 1)));
  $('aboutBtn').addEventListener('click', () => toggleAbout(true));
  $('aboutClose').addEventListener('click', () => toggleAbout(false));

  $('playBack').addEventListener('click', () => togglePlay(-1));
  $('playFwd').addEventListener('click', () => togglePlay(+1));
  $('reset').addEventListener('click', () => { setT(0); setPlay(0); });
  $('scrub').addEventListener('input', (e) => setT(invSymlog(+e.target.value)));
  $('infoclose').addEventListener('click', deselect);
  /* A control that has focus owns its own keys. Guarding only on INPUT meant
   * Space on a focused row selected the star *and* started playback, Space on a
   * button pressed it *and* started playback, and an arrow key on the speed
   * menu changed the speed *and* scrubbed time. */
  const editing = (el) => /^(INPUT|SELECT|TEXTAREA)$/.test(el.tagName)
                          || el.isContentEditable;
  const activatable = (el) => el.tagName === 'BUTTON' || el.tagName === 'A'
                          || el.getAttribute('role') === 'option';
  addEventListener('keydown', (e) => {
    if (editing(e.target)) return;
    // Escape leaves full screen before it does anything else, so one press
    // never both shrinks the window and stops the tour.
    if (e.code === 'Escape' && fsOn()) return;
    // Full screen is a property of the window, not of the scene, so it stays
    // available while the tour is driving - the button itself is inert then.
    if (e.code === 'KeyF' && !e.metaKey && !e.ctrlKey && !e.altKey) {
      if (fsSupported()) toggleFullscreen();
      return;
    }
    // While the tour is driving, the transport keys belong to it.
    if (tour.on) {
      if (e.code === 'Escape') { tourStop(); return; }
      if (e.code === 'Space') { e.preventDefault(); tourPause(!tour.paused); return; }
      if (e.code === 'ArrowRight') { tourGo(tour.i + 1); return; }
      if (e.code === 'ArrowLeft') { tourGo(Math.max(0, tour.i - 1)); return; }
      // Every other shortcut belongs to the tour too. J, L, K, E and the zoom
      // keys used to fall through and start playback, flip to the sky view or
      // change the zoom underneath a shot that was already animating - the
      // keyboard equivalent of the controls that are now inert.
      return;
    }
    if (e.code === 'Escape') {
      if (!$('about').hidden) toggleAbout(false);
      else if (document.body.classList.contains('filters-open')) setPanel('filters', false);
      else if (document.body.classList.contains('menu-open')) setPanel('menu', false);
      else if (document.body.classList.contains('chart-open')) setPanel('chart', false);
      else deselect();
    }
    if (activatable(e.target) && (e.code === 'Space' || e.code === 'Enter')) return;
    // through togglePlay so Space gets the same replay-at-the-end behavior
    if (e.code === 'Space') { e.preventDefault(); togglePlay(S.play || S.lastDir); }
    if (e.key === 'j' || e.key === 'J') togglePlay(-1);
    if (e.key === 'l' || e.key === 'L') togglePlay(+1);
    if (e.key === 'k' || e.key === 'K') setPlay(0);
    if (e.key === 'e' || e.key === 'E') setSky(!S.sky.on);
    if (e.key === '+' || e.key === '=') {
      if (S.sky.on) setFov(S.sky.fov / 1.3);
      else { S.cam.dist = Math.max(DIST_MIN, S.cam.dist / 1.4); camAnim = null; }
    }
    if (e.key === '-' || e.key === '_') {
      if (S.sky.on) setFov(S.sky.fov * 1.3);
      else { S.cam.dist = Math.min(DIST_MAX, S.cam.dist * 1.4); camAnim = null; }
    }
    if (e.code === 'ArrowLeft') setT(invSymlog(symlog(S.t) - 0.01));
    if (e.code === 'ArrowRight') setT(invSymlog(symlog(S.t) + 0.01));
  });
}
/* Playback direction is explicit. The old single button ping-ponged at the
 * ends, which made it impossible to say "run this forwards" and have it mean
 * that; now each button owns a direction, pressing the active one pauses, and
 * reaching either end of the window stops rather than silently reversing. */
// Within a step of the end, playing further is a no-op, so the button has to
// say something other than "play".
const atEnd = (d) => (d > 0 ? S.t >= TSPAN - 1e-9 : S.t <= -TSPAN + 1e-9);

let playBtnState = '';
function refreshPlayButtons() {
  // Cheap guard: setT runs every frame during playback, and rewriting four
  // DOM strings sixty times a second to say the same thing is pure churn.
  const key = `${S.play}|${atEnd(-1)}|${atEnd(1)}`;
  if (key === playBtnState) return;
  playBtnState = key;
  const d = S.play;
  for (const [id, dir, word] of [['playBack', -1, 'backwards into the past'],
                                 ['playFwd', 1, 'forwards into the future']]) {
    const b = $(id);
    const active = d === dir;
    b.classList.toggle('on', active);
    b.textContent = active ? '❚❚' : (dir < 0 ? '◀' : '▶');
    // The label has to track what the button will actually do. It showed
    // "Play backwards" while displaying a pause glyph, which is the tooltip
    // contradicting the icon, and screen readers got the bare glyph.
    const label = active ? 'Pause'
      : atEnd(dir) ? `Replay from ${dir > 0 ? '−' : '+'}${TSPAN} Myr`
      : `Play ${word}`;
    b.title = `${label} (${dir < 0 ? 'J' : 'L'})`;
    b.setAttribute('aria-label', label);
  }
}
function setPlay(d) {
  S.play = d;
  invalidate();
  if (d) S.lastDir = d;
  refreshPlayButtons();
}
function togglePlay(d) {
  if (S.play === d) { setPlay(0); return; }
  // Standing at the end of the window, "play" means replay: every media
  // control in the world restarts from the top rather than doing nothing, and
  // doing nothing is what this did - one frame of playback, immediately
  // clamped back off, with no hint as to why.
  if (atEnd(d)) setT(-d * TSPAN);
  setPlay(d);
}
function refilter() { recomputeVisible(); pickChartStars(); renderList(); invalidate(); }
function select(i, focus = false) {
  S.sel = (i == null || i < 0) ? -1 : i;
  i = S.sel;
  if (i >= 0) S.pulse = performance.now();
  invalidate();
  showInfo(i); pickChartStars(); renderList(); buildTrail(i);
  if (i >= 0 && focus) focusStar(i);
}

/* Clear the selection. Selecting a star can move the camera, so clearing it
 * puts the view back where it was - unless the user has since driven the camera
 * themselves, in which case their view is the one worth keeping. */
function deselect() {
  const back = S.preFocus;
  S.preFocus = null;
  select(-1);
  if (back && !S.camDirty && HAS_GL) {
    if (back.radius !== S.filt.radius) { setRadius(back.radius); refilter(); }
    flyTo(back.yaw, back.pitch, back.dist, back.target, 620);
    S.camDirty = false;
  }
}

/* Clicking the same star again clears it; clicking a different one moves on. */
function toggleSelect(i, focus = false) {
  if (i >= 0 && i === S.sel) deselect();
  else select(i, focus);
}
function setT(t) {
  // Math.min/max propagate NaN rather than clamping it, and a NaN epoch is
  // unrecoverable: every star's position becomes NaN, so the view empties, the
  // chart blanks, and no later input clears it because the arithmetic stays
  // NaN. No UI control can produce one today - the sliders all coerce - but
  // the cost of never finding out the hard way is one comparison.
  if (!Number.isFinite(t)) return;
  S.t = Math.max(-TSPAN, Math.min(TSPAN, t));
  invalidate();
  $('scrub').value = symlog(S.t);
  $('tnow').textContent = fmtT(S.t);
  refreshPlayButtons();      // scrubbing to an end changes what play means
}

/* --------------------------------------------------------- mobile chrome */
/* On a narrow screen the star list is an off-canvas sheet and the chart is
 * optional, so that the 3-D view - the thing worth looking at on a phone - gets
 * the whole viewport. Both are plain classes on <body>; the layout lives in CSS
 * and this only flips state and tells the canvases their size changed. */
const isNarrow = () => matchMedia('(max-width: 900px)').matches;

const PANEL_BTN = { menu: 'menubtn', chart: 'chartbtn', filters: 'filterbtn' };

function setPanel(name, open) {
  const btn = $(PANEL_BTN[name]);
  document.body.classList.toggle(`${name}-open`, open);
  btn.setAttribute('aria-expanded', String(open));
  // The chart's label budget is chosen from its pixel width, which is 0 while
  // it is display:none, so it has to be recomputed once it is on screen.
  if (name === 'chart' && open) requestAnimationFrame(() => { pickChartStars(); invalidate(); });
  invalidate();
}
const togglePanel = (name) =>
  setPanel(name, !document.body.classList.contains(`${name}-open`));

/* The dock's height changes when the chart is toggled, and the sheet has to
 * stop above it - a fixed inset left the last filter row hidden behind the
 * transport. Measured rather than assumed, same approach as the detail card. */
function fitMobileChrome() {
  const dock = $('dock');
  if (!dock) return;
  document.documentElement.style.setProperty(
    '--dockh', `${Math.round(dock.getBoundingClientRect().height)}px`);
}

/* The title bar has no room for the Methods button beside the menu toggle, and
 * hiding it would remove the one place the interface explains its own limits.
 * So it lives in the sheet on a phone and in the title bar on a desktop. */
/* Full screen.
 *
 * The whole stage goes fullscreen, not just the 3-D canvas. The list, the
 * chart, the caption and the orientation panel are part of the instrument, and
 * fullscreening the canvas alone would drop every one of them. #stage is
 * already `position: fixed; inset: 0`, so filling a screen needs no layout
 * change - only a background, since <body> is not painted behind a fullscreen
 * element.
 *
 * The button cannot own the state. Esc, the browser's own chrome and the OS
 * can all leave fullscreen without going through it, so the label follows the
 * `fullscreenchange` event and the click handler only ever asks. */
const fsEl = () => $('stage');
const fsOn = () => !!(document.fullscreenElement || document.webkitFullscreenElement);
const fsSupported = () => {
  const el = fsEl();
  return !!(el && (el.requestFullscreen || el.webkitRequestFullscreen));
};

function toggleFullscreen() {
  try {
    if (fsOn()) {
      (document.exitFullscreen || document.webkitExitFullscreen).call(document);
    } else {
      const el = fsEl();
      const req = el.requestFullscreen || el.webkitRequestFullscreen;
      // navigationUI:'hide' is a hint; browsers that do not know it ignore it.
      const p = req.call(el, { navigationUI: 'hide' });
      if (p && p.catch) p.catch(() => syncFullscreen());
    }
  } catch (e) {
    // Refused - an embedding page without allow="fullscreen", or a platform
    // that only does it for video. Leave the label telling the truth.
    syncFullscreen();
  }
}

function syncFullscreen() {
  const b = $('fsBtn');
  if (!b) return;
  const on = fsOn();
  b.setAttribute('aria-pressed', String(on));
  b.textContent = on ? '⛶ Exit full screen' : '⛶ Full screen';
  b.title = on ? 'Leave full screen (F or Esc)' : 'Fill the screen (F)';
  // Entering or leaving changes every measured dimension: the dock height, the
  // detail card's ceiling, and the chart's label budget, which is chosen from
  // its pixel width. The canvas resizes itself from clientWidth on the next
  // frame, so it only needs to be asked for one.
  fitMobileChrome(); fitInfoCard(); pickChartStars(); invalidate();
}

function placeAboutBtn() {
  const btn = $('aboutBtn');
  const want = isNarrow() ? $('filters') : $('title');
  if (btn && btn.parentElement !== want) want.appendChild(btn);
  // Full screen travels with Methods & limits: same reason, same destination.
  const fs = $('fsBtn');
  if (fs && fs.parentElement !== want) want.appendChild(fs);
  // Narrow layouts drop the whole orientation panel, which would strand the
  // From Earth toggle there with no way to reach it. Move it into the sheet
  // alongside Methods & limits rather than leaving the view inaccessible.
  const sky = $('vSky');
  const target = isNarrow() ? $('filters') : document.querySelector('.navbtns');
  if (sky && target && sky.parentElement !== target) target.appendChild(sky);
}

/* The detail card and the orientation panel share the right-hand column but are
 * positioned from opposite edges, so on a short window a long card ran into the
 * top of the panel. The CSS reserves a fixed guess; this replaces it with the
 * panel's measured extent, which stays correct if its contents ever change. */
function fitInfoCard() {
  const nav = $('nav');
  const shown = nav && nav.offsetParent !== null;
  // max-height is measured from the card's own top, so the reservation has to
  // cover the panel's extent, the card's 14px top offset and a 14px gap.
  const reserve = shown ? innerHeight - nav.getBoundingClientRect().top + 28 : 28;
  document.documentElement.style.setProperty('--navblock', `${Math.round(reserve)}px`);
}

/* The panel grows after boot - the star count and model lines are filled in by
 * script - so measuring it once is measuring the wrong thing. Watch it instead. */
function watchNav() {
  const nav = $('nav');
  if (typeof ResizeObserver === 'function') {
    if (nav) new ResizeObserver(fitInfoCard).observe(nav);
    new ResizeObserver(fitMobileChrome).observe($('dock'));
  } else {
    addEventListener('load', () => { fitInfoCard(); fitMobileChrome(); });
  }
}

/* Methods panel. The interface previously shipped none of this, so a reader had
 * no way to know what the numbers mean without finding the repository. Built
 * from meta.json so it cannot drift from the data. */
function toggleAbout(show) {
  const el = $('about');
  if (!show) { el.hidden = true; return; }
  const m = S.meta;
  const li = (xs) => xs.map((t) => `<li>${t}</li>`).join('');
  el.querySelector('#about-body').innerHTML = `
    <p><b>${m.n_stars_shipped ?? m.n_stars} stars</b> shown of
       ${m.n_stars_catalog ?? m.n_stars} screened
       (${m.display_cut ?? 'see repository'}).
       Source catalog <b>${m.catalog ?? 'Gaia DR3'}</b>, retrieved
       ${m.retrieved ?? 'unknown'}.</p>

    <h3>Dynamics</h3>
    <p>The Sun and every star are integrated as test particles through a Milky
       Way potential (${m.sources?.potential ?? ''}), not drifted in straight
       lines. Trajectories are compressed to degree-${m.cheb_degree} Chebyshev
       polynomials, accurate on the closest approach to
       ${(m.perihelion_err_max_pc * 1e3).toFixed(1)} milliparsec at worst.</p>

    ${m.encounter_rate ? `
    <h3>How often does this happen?</h3>
    <p>Counting encounters near the present, where the catalog is nearly
       complete, and correcting for what it is missing at each distance, gives
       <b>${m.encounter_rate.rate_within_1pc_per_myr} &plusmn;
       ${m.encounter_rate.stat_error} (stat) &plusmn;
       ${m.encounter_rate.sys_error} (sys) encounters per million years within
       1 parsec</b> — about one every 48,000 years. The published measurement is
       ${m.encounter_rate.published.value} &plusmn;
       ${m.encounter_rate.published.error}
       (${esc(m.encounter_rate.published.source)}), which this agrees with to
       ${m.encounter_rate.sigma_from_published}&sigma;. It is the one number
       here derived from the catalog as a whole rather than read off a single
       star, and it was arrived at independently.</p>
    <p class="warn">${esc(m.encounter_rate.caveat)}</p>` : ''}

    <h3>About those travel times</h3>
    <p>The familiar "Voyager would take seventy thousand years to reach Proxima"
       is distance divided by speed, with the star held still. Curvature is not
       what breaks it — integrating a coasting probe through the Galactic
       potential departs from a straight line by <b>0.0001% of the trip</b> over
       seventy-five thousand years, because the Sun and the probe fall through
       the Galaxy together. What breaks it is that the target moves. Proxima
       recedes at <b>32 km/s</b> after its closest approach, nearly twice
       Voyager&nbsp;1's 17 km/s, so a probe launched today never catches it —
       it falls about 1.1 light years short and loses ground from then on.
       Gliese 710 closes at only 14 km/s, slower than the probe, so an intercept
       does exist: integrated, about <b>594,000 years</b>, some 46% sooner than
       the naive figure because the star does most of the traveling.</p>

    <h3>Reading the chart</h3>
    <p>Distance from the Sun runs up the vertical axis, logarithmically. Time
       runs across it <i>symmetric-logarithmically</i>: spacing is linear within
       ${(1e4).toLocaleString()} years of now and logarithmic outside, so the
       near future is legible without hiding the far one. Only the closest few
       curves are labeled; hover any curve to name it. A selected or hovered
       star gets a vertical bar at its closest approach — that is the
       16th–84th percentile of the resampled draws, not a drawing flourish.</p>

    <h3>From Earth</h3>
    <p>The <b>From Earth</b> view stands at the Sun and draws each star at the
       brightness it would have in Earth's sky at that moment,
       m = M + 5&nbsp;log&#8321;&#8320;(d) &minus; 5, using the distance the
       trajectory is already reporting. Scrub time and the sky changes: Sirius
       leads it today at &minus;1.1, Gliese 710 reaches &minus;3.8 in 1.3 Myr,
       and HD 7977 touches roughly &minus;8 on its way past 2.8 Myr ago —
       brighter than Venus ever gets, though its perihelion is uncertain enough
       that the figure could be two magnitudes either way.</p>
    <p><b>These are Gaia G magnitudes</b>, not the visual scale the famous
       numbers are quoted on, because G is what Gaia measures for every star
       here. For a Sun-like star the two agree to a tenth; for a red one they
       do not. Gliese 710 is a K7 dwarf, so its &minus;3.8 in G is about
       &minus;3.2 visually, and the published figure — computed at a slightly
       more distant perihelion than this catalog finds — is about
       &minus;2.7. Same star, three defensible numbers, and the difference is
       filters and perihelion rather than disagreement about brightness.</p>
    <p class="warn">This is not the night sky. It contains only the stars this
       catalog tracks — the solar neighborhood and everything screened as an
       encounter candidate — so about 660 of them are naked-eye visible today
       against roughly 9,000 in the real sky. The distant giants that make up
       most of the constellations are not here, and no star is drawn without
       Gaia photometry. Read it as which of the Sun's <i>neighbors</i> are
       visible, and how that changes.</p>

    <h3 class="warn">What the intervals mean</h3>
    <p>${m.monte_carlo?.model ?? ''}.
       <b class="warn">${m.monte_carlo?.interval_meaning ?? ''}</b></p>

    <h3 class="warn">…and what they leave out</h3>
    <p>The interval says how far a prediction moves when the <i>measurements</i>
       are resampled. It says nothing about the Milky Way's mass being slightly
       different from the model assumed here. That was measured separately, by
       re-integrating the whole catalog under 21 perturbations of the disc,
       halo, bulge and solar parameters: for stars within 25 pc it moves nothing
       (a ten-thousandth of a parsec), but beyond 100 pc it is comparable to the
       measurement interval, and <b class="warn">for
       ${m.model_limited ? `${m.model_limited.limited} of the
       ${m.model_limited.candidates.toLocaleString()}` : 'a number of the'}
       candidates approaching within 5 pc it is larger than the interval
       entirely</b>.
       Those carry a flag, and every star's figure is the <i>Galaxy-model
       shift</i> line in its detail card. The Gaia parallax zero-point — a known
       offset of order 0.02–0.05 mas that this build does not correct for —
       was measured the same way and appears as the <i>parallax zero-point</i>
       line.</p>

    <h3>Known limits</h3>
    <ul>${li([
      'Stars with no measured radial velocity cannot be propagated and are absent.',
      'Star counts away from the present are a lower bound: the sample was selected on being near now or approaching.',
      'Stars are test particles — no encounters between stars, no bar, no spiral arms.',
      'Systems are candidate co-moving groups, not a binding-energy test.',
      'Some trajectories never reach a minimum inside the window; those are shown as bounds, not encounters.',
      'The Galaxy-model shift is an envelope over chosen variants, not a posterior; do not add it to the interval in quadrature.',
    ])}</ul>

    <h3>Sources</h3>
    <ul>${li([
      m.sources?.astrometry ?? '', m.sources?.radial_velocity ?? '',
      m.sources?.solar_motion ?? '',
    ].filter(Boolean))}</ul>

    <h3>Data</h3>
    <p>The full screened catalog, including every column behind these numbers,
       is <code>data/star_encounters.csv</code> in the repository, with column
       documentation in <code>data/COLUMNS.md</code>. The derived catalog is
       CC BY 4.0; the code is MIT.</p>

    <h3>Acknowledgments</h3>
    <p>This work has made use of data from the European Space Agency mission
       <b>Gaia</b>, processed by the Gaia Data Processing and Analysis
       Consortium (DPAC). Funding for the DPAC has been provided by national
       institutions, in particular the institutions participating in the Gaia
       Multilateral Agreement. It also uses the SIMBAD database and the VizieR
       catalog access tool, CDS, Strasbourg, France.</p>`;
  el.hidden = false;
  $('aboutClose').focus();
}

/* --------------------------------------------------------------------- loop */
/* Anything that changes what is on screen marks the frame dirty. Idle frames
 * then cost nothing instead of re-evaluating thousands of polynomials and
 * redrawing the whole chart sixty times a second for an unchanged picture. */
function invalidate() { S.dirty = true; }

let last = performance.now();
function frame(now) {
  const dt = Math.min(0.05, (now - last) / 1000); last = now;
  stepTimeTween(now);
  if (S.play && !timeDrag) {
    // Advance in YEARS, not in axis units. Sweeping the symlog axis at a
    // constant rate looked right on the chart and wrong everywhere else: the
    // real-time rate it implies runs from 0.06 to 31 Myr per axis unit, a
    // factor of 500, so the stars crawled near the present and then visibly
    // rocketed. Stars move at a fixed speed through space, so the only way
    // they look like they are moving at a fixed speed is for time itself to
    // advance at a fixed rate. The cursor now covers the compressed middle of
    // the axis quickly instead, which is the honest way round: the distortion
    // belongs to the axis, not to the physics.
    let t = S.t + dt * (+$('speed').value) * S.play;
    if (t >= TSPAN) { t = TSPAN; setPlay(0); }   // stop at the ends, do not bounce
    if (t <= -TSPAN) { t = -TSPAN; setPlay(0); }
    setT(t);
  }
  stepCamAnim(now);
  stepTourDrift(dt);

  const pulsing = (now - S.pulse) / 1000 < PULSE_SECS;
  const active = S.play || camAnim || timeDrag || pulsing || tTween || S.dirty;
  if (active) {
    S.dirty = false;
    render();
    drawOverlay();
    drawGizmo();
    drawChart();
    if (S.play && S.tab === 'near') renderList();
    if (S.sel >= 0) showInfo(S.sel);
    updateHud();
  }
  requestAnimationFrame(frame);
}
const esc = (s) => String(s).replace(/[&<>"']/g,
  (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

function updateHud() {
  if (S.sky.on) {
    // Naked eye is conventionally magnitude 6.5. The count is worth showing
    // next to the brightest star because it moves so much: the same sky holds
    // 660 of these stars today and under 200 at either end of the window.
    let n = 0;
    for (let i = 0; i < S.n; i++) {
      if (S.visible[i] && appMag(i, S.t) <= S.sky.magLim) n++;
    }
    const top = skyBrightest(1)[0];
    // Naming the brightest star is useless if you then have to hunt the sky
    // for it, so the name turns the view to face it.
    $('hud').innerHTML =
      (top ? `brightest <button type="button" class="lnk" id="hudtop" `
             + `data-i="${top.i}" title="Turn to face it">`
             + `${esc(S.d.names[top.i])}</button> mag ${top.m.toFixed(2)}`
           : 'nothing above the limit') +
      `<br>${n} naked-eye &nbsp;·&nbsp; ${fmtT(S.t)}`;
    $('zoom').value = fovToU(S.sky.fov);
    // Name the band. These are Gaia G, which for a red star runs most of a
    // magnitude brighter than the visual scale the familiar numbers are
    // quoted on, and an unlabeled "mag" invites a false comparison.
    $('zlabel').textContent =
      `${Math.round(S.sky.fov * 180 / Math.PI)}° field · to Gaia G ${S.sky.magLim}`;
    const bar = $('skybar');
    // Narrow layouts give the detail card most of the screen, and it already
    // names the star; two readouts would just collide at the card's bottom.
    bar.hidden = !isNarrow() || $('info').style.display === 'block';
    if (!bar.hidden) {
      bar.innerHTML = top
        ? `brightest <button type="button" class="lnk" id="hudtop2" `
          + `data-i="${top.i}">${esc(S.d.names[top.i])}</button> `
          + `mag ${top.m.toFixed(2)} · ${n} naked-eye`
        : `nothing above magnitude ${S.sky.magLim}`;
    }
    return;
  }
  $('skybar').hidden = true;
  let n = 0;
  for (let i = 0; i < S.n; i++) if (S.visible[i] && distAt(i, S.t) <= S.filt.radius) n++;
  $('hud').innerHTML =
    `${n} stars &lt; ${S.filt.radius} pc &nbsp;·&nbsp; ${fmtT(S.t)}<br>` +
    `Galactic potential · Chebyshev ${DEG}`;

  $('zoom').value = distToU(S.cam.dist);
  const across = viewAcross();
  const pc = across < 10 ? across.toFixed(1) : Math.round(across);
  $('zlabel').textContent = `${pc} pc · ${Math.round(across * LY)} ly across`;
}

/* --------------------------------------------------------------------- tour */
/* A narrated pass through the catalog. Two decisions shape the whole thing.
 *
 * The narration is spoken with the browser's own speech synthesiser rather than
 * shipped as audio: five minutes of voice track is several megabytes, would
 * have to be regenerated whenever a number changes, and would be one more thing
 * that can silently drift from the data. Captions are drawn either way, so the
 * tour is complete with the sound off and works where there is no voice at all.
 *
 * Beats advance when the SPEECH ENDS, not on a timer. Speech rate varies by
 * platform, voice and user setting, so any fixed schedule would drift out of
 * sync with the words within a minute. Where there is no speech to wait for,
 * the beat falls back to an estimate from its own word count. */

let tTween = null;
function tweenT(to, ms) {
  setPlay(0);
  tTween = { from: S.t, to, t0: performance.now(), ms };
}
function stepTimeTween(now) {
  if (!tTween) return;
  const u = Math.min(1, (now - tTween.t0) / tTween.ms);
  // Linear in years, eased at both ends, for the same reason playback is: the
  // tour is asking you to watch a star travel, and a tween through axis units
  // would have it drift and then bolt.
  const e = u < 0.5 ? 2 * u * u : 1 - Math.pow(-2 * u + 2, 2) / 2;
  setT(tTween.from + (tTween.to - tTween.from) * e);
  if (u >= 1) tTween = null;
}

/* A held shot of two points of light is a photograph: you cannot tell how far
 * apart they are, or which is in front. Turning slowly around the pair gives
 * the parallax that makes the geometry read. Slow enough that the next beat's
 * fly-to reads as a deliberate move rather than a correction, and off entirely
 * when the visitor has asked for reduced motion. */
const TOUR_DRIFT = 0.035;                 // rad/s, about 2 degrees a second
function stepTourDrift(dt) {
  if (!tour.on || tour.paused || camAnim || S.sky.on) return;
  if (REDUCED_MOTION.matches) return;
  S.cam.yaw += dt * TOUR_DRIFT;
  invalidate();
}

const idx = (name) => (S.d ? S.d.names.indexOf(name) : -1);
function look(name, ms = 1600) {
  const j = idx(name);
  if (j < 0) return;                  // a renamed star must not stall the tour
  select(j, false);                   // focusStar below does the flying
  focusStar(j, ms);
}
/* Point at a star without moving the camera. When the narration walks through
 * the members of one system, flying to each in turn would be seasickness; the
 * reticle and its name label moving is the whole signal. It has to be, for
 * Alpha Centauri A and B: 25 AU apart at 4.3 light years, there is no camera
 * position that separates them, and the changing label is the only honest way
 * to say which one is being talked about. */
function mark(name) {
  const j = idx(name);
  if (j >= 0) select(j, false);
}
const wide = (yaw, pitch, dist, ms = 2200) => flyTo(yaw, pitch, dist, [0, 0, 0], ms);

/* Watch a star come in. `focusStar` frames a star where it is *now*, so a shot
 * that flies to it and then scrubs a million years ends up pointed at empty
 * space while the star arrives somewhere else entirely. Re-frame on arrival. */
function arrive(name, toT, ms) {
  look(name, Math.min(2200, ms * 0.35));
  tweenT(toT, ms);
  const id = setTimeout(() => {
    const j = idx(name);
    if (tour.on && j >= 0) focusStar(j, 1800);
  }, ms + 150);
  tour.pending.push(id);
}

const TOUR = [
{ say: `Every star you can see is moving. Not slowly — tens of kilometers every second — but they are so far away, and we live so briefly, that the sky looks painted on. Only a few dozen sit close enough to matter to us: if we are ever to become more than a one-star species, these are the only doors within reach. This is what ten million years does to them.`,
  act() { deselect(); setSky(false); setT(0); setRadius(20); refilter(); wide(0.6, 0.32, 46); },
  // Narrow to the neighborhood on the sentence that says why it matters.
  // Dropping the view radius rather than flying makes the point by what it
  // removes: the field empties to the forty-eight stars inside five parsecs.
  cues: [null, null,
         () => { setRadius(5); refilter(); wide(0.75, 0.30, 16, 3200); },
         null] },

{ say: `This project collected twelve thousand two hundred and eighteen stars' positions and motions from the Gaia spacecraft, and carried every one of them through the gravity of the Milky Way. The Sun sits at the middle only because the map is drawn from where we stand — not at the middle of the Galaxy, whose center lies twenty-six thousand light years away, far outside the little neighborhood you see here. The rings around the Sun are measured in parsecs, and one parsec is about three and a quarter light years. They mark five, ten, twenty and forty parsecs out.`,
  act() { setRadius(40); refilter(); wide(1.5, 0.55, 95, 4000); } },

{ say: `Start close. This is Alpha Centauri — the nearest star system to the Sun right now, and the only one within five light years. It is a three-star system, held together by gravity.`,
  act() { setRadius(20); refilter(); look('Alpha Centauri A', 2400); } },

{ say: `Alpha Centauri A is a near twin of the Sun. Alpha Centauri B is cooler and smaller, and it circles A about as far out as Neptune is from the Sun — so close together that at this scale the two of them sit on a single point. Proxima is the third: a red dwarf, which means a small cool star far too faint to see without a telescope, and it sits about a fifth of a light year away from the other two — which still leaves it the closest star to Earth, at four point two five light years.`,
  // Fly into the system rather than viewing it from the Sun. At the previous
  // framing Proxima sits about eighteen pixels from the AB pair, which is not
  // a separation a viewer can see the marker move across; from inside, the
  // fifth of a light year the narration describes is most of the frame.
  act() {
    const j = idx('Alpha Centauri A');
    if (j < 0) return;
    const p = posAt(j, S.t);
    setRadius(20); refilter();
    flyTo(Math.atan2(p[1], p[0]) + 0.7, 0.25, 0.55, p, 2600);
    mark('Alpha Centauri A');
  },
  cues: [() => mark('Alpha Centauri A'),
         () => mark('Alpha Centauri B'),
         () => mark('Proxima Centauri')] },

{ say: `Those four point two five light years are not fixed. The whole Alpha Centauri system is drifting toward us. In about twenty-eight thousand years the gap closes to three point one light years — a quarter nearer than today — and then it starts to open again.`,
  // Pull back out of the system before scrubbing. This beat is about the
  // distance between two things, so both have to be in shot - and the tight
  // framing inherited from the previous beat is fatal here: over 28 kyr Alpha
  // Centauri travels 0.92 pc through a frame 0.46 pc wide and leaves the top
  // of the screen while the narration is still describing it. `arrive` frames
  // the Sun-star pair and re-frames once time stops moving.
  act() { arrive('Alpha Centauri A', 0.0278, 6000); } },

{ say: `But Alpha Centauri will not stay the nearest star for ever. In about thirty-three thousand years a dim red star called Ross 248, ten light years away today, will have moved close enough to take its place. Ross 248 stays nearest for only about nine thousand years, and then it carries on past us and away. A star called Gliese 445 is nearest after that, and it does no better — a few thousand years, and then it leaves too. And then the Alpha Centauri stars are the nearest once more, just as they are today. Over the next million years, the nearest star to the Sun changes about thirty times.`,
  // Three stars hold the record in turn here, so the camera has to hand over
  // with them; parking on Ross 248 for the whole beat contradicts the words.
  // Each `arrive` scrubs to the epoch the sentence is describing and re-frames
  // the Sun-star pair on landing, so you watch each holder come in and leave.
  // The sentence count is load-bearing: cues[i] fires as sentence i begins, so
  // this beat must stay at exactly six, in this order of subjects.
  act() {},
  cues: [null,
         () => arrive('Ross 248', 0.0365, 6000),
         () => arrive('Ross 248', 0.0463, 6000),      // watch it recede again
         () => look('Gliese 445', 2600),              // same epoch, new holder
         () => arrive('Alpha Centauri B', 0.0560, 5500),
         // Last sentence is a short one; keep the move inside it rather than
         // letting the next beat cut a half-finished scrub.
         () => { deselect(); setRadius(10); refilter();
                 wide(2.1, 0.38, 24, 3000); tweenT(0.5, 4200); }] },

{ say: `But none of that is really a close encounter. Alpha Centauri, Ross 248, Gliese 445 — not one of them gets inside three light years, and the tightest of the three, Ross 248, still stops at three point zero five. For a genuine near miss we have to look further out, and much further ahead.`,
  // The one beat that deliberately does *not* hand the camera over as it names
  // its stars: the sentence exists to dismiss all three at once, and marking
  // them in turn would argue against the words while the shot pulls away.
  act() { deselect(); tweenT(0, 3000); setRadius(25); refilter(); wide(0.6, 0.42, 55, 3000); } },

{ say: `This is Gliese 710, an orange dwarf in the constellation Serpens, about six tenths the mass of the Sun. Right now it is sixty-two light years away, and entirely unremarkable — far too faint to see with the naked eye.`,
  act() { look('Gliese 710', 2600); } },

{ say: `Except for one thing. It is heading almost exactly at us. Watch any other star for a year and it shifts a little against the ones behind it; Gliese 710 barely shifts at all. That is what a star looks like when it is coming straight toward you.`,
  act() {} },

{ say: `In one point three million years, Gliese 710 will pass about zero point one seven light years from the Sun.`,
  act() { arrive('Gliese 710', 1.296, 9000); } },

// The old wording here said "the refereed figure sits further out, near two
// tenths" — which implied one published number, and that this catalog was
// alone below it. Neither is true: the refereed values split into a group near
// 0.17 (de la Fuente Marcos 2020 on EDR3, essentially our answer) and a group
// near 0.20 (Bailer-Jones 2018/2022, Fernandez-Puig 2026). Stating the spread
// is both shorter and more honest than defending a disagreement.
{ say: `That is about eleven thousand times the distance from the Earth to the Sun. Different research teams get slightly different answers, somewhere between about zero point one seven and zero point two light years, and this project's answer is at the closer end of that range. Either way, Gliese 710 passes deep inside the Oort cloud — a huge shell of icy comets surrounding the solar system — and it spends tens of thousands of years in there, stirring them up.`,
  act() { setRadius(3); refilter(); wide(1.2, 0.3, 4.5, 3000); } },

// The old wording explained the difference between Gaia's green band and the
// visual scale, then quoted a number from each - the hardest thirty seconds in
// the tour, and the two-scale detour bought nothing a listener could use. Both
// figures now come from the band the readout on screen is already showing, so
// the narration and the display agree and only one scale has to be explained.
{ say: `Seen from Earth it would be the brightest star in the night sky, easily beating Sirius, which is the brightest star we see today. Astronomers measure brightness on a backwards scale: smaller numbers mean brighter, and the very brightest objects go below zero. Sirius reads about minus one, and Gliese 710 would reach nearly minus four — more than ten times brighter than Sirius — and it would shine a distinct orange.`,
  act() { setSky(true); setFov(1.0); look('Gliese 710', 2000); } },

{ say: `Gliese 710 may also be our best chance of ever reaching another star. Voyager One, the fastest thing we have ever launched, travels about seventeen kilometers every second. People often say it would take seventy thousand years to reach Proxima Centauri at that speed — but that answer pretends Proxima stays where it is. Proxima is moving away from us faster than Voyager can fly, so a probe launched today would never catch up with it at all. Gliese 710 is different: it travels slower than Voyager, and it comes almost all of the way to us by itself. When it is closest, a probe leaving Earth at Voyager's speed would take only about three thousand years to get there.`,
  act() {} },

{ say: `Now run it backwards.`,
  act() { setSky(false); deselect(); setRadius(20); refilter(); wide(0.6, 0.32, 46, 2500); tweenT(0, 4000); } },

{ say: `Seventy-nine thousand years ago, while our ancestors were spreading out of Africa, two very small stars swept past together, just under one light year away. The pair goes by a single name: Scholz's Star. One of them is a red dwarf, and the other is a brown dwarf — an object too small ever to shine properly as a star. They passed through the outer edge of the Oort cloud, and nobody would have noticed a thing, because both are far too faint to see.`,
  act() { arrive("Scholz's Star", -0.0789, 7000); } },

{ say: `Further back the record gets less certain, and a great deal more dramatic. Two point eight million years ago, HD 7977 came within about an eighth of a light year — closer than Gliese 710 will ever get.`,
  act() { arrive('HD 7977', -2.765, 9000); } },

{ say: `That word "about" matters here. HD 7977 is much further away than the others and much harder to measure, so its closest approach could really have been anywhere from around a sixteenth of a light year to a fifth — the far end of that range is three times the near end. That is not a mistake in the calculation. It is simply how much the measurements can tell us so far, and this project shows you that range rather than hiding it.`,
  act() { setRadius(3); refilter(); wide(2.4, 0.35, 4.5, 3000); } },

{ say: `At its closest, HD 7977 would have reached magnitude minus eight. Brighter than Venus has ever been. If anything had been looking up, it would have been the brightest thing in the night sky after the Moon, and it would have cast shadows.`,
  act() { setSky(true); setFov(1.2); look('HD 7977', 2000); } },

{ say: `None of this is unusual. This map cannot possibly show every star — some are too faint to have been spotted, and some have not been measured yet — so the count is done in two steps: add up the close passes we can see, then work out how many we must be missing. Put those together, and a star comes within one parsec of the Sun — three and a quarter light years — about twenty-one times every million years.`,
  act() { setSky(false); deselect(); setRadius(30); refilter(); tweenT(0, 4000); wide(2.2, 0.5, 70, 4000); } },

{ say: `That is one close pass every forty-eight thousand years. About six of them since our species appeared. Around fourteen hundred since the dinosaurs died out. Other astronomers, working from different data, counted nineteen point seven passes per million years; this project counted twenty point seven, and two separate attempts landing that close together is what it looks like when a result is real.`,
  act() { setPlay(1); } },

{ say: `The Sun does not really have neighbors. It has traffic. Alpha Centauri is not a fixed address — it is a car we happen to be driving alongside for a while.`,
  act() { setPlay(0); wide(4.0, 0.25, 80, 8000); } },

{ say: `Everything here is measured, not imagined. The Gaia spacecraft told us where these stars are and how fast they are moving, and every path you have seen was worked out by following the pull of the Milky Way's gravity, one small step at a time. Where a number is uncertain, this project tells you how uncertain, instead of hiding it. Have a look around.`,
  act() { deselect(); setT(0); setRadius(20); refilter(); wide(0.6, 0.32, 46, 3000); } },
];

const tour = { on: false, i: -1, paused: false, speech: 0,
               timer: null, ping: null, pending: [] };
function tourClearPending() {
  for (const id of tour.pending) clearTimeout(id);
  tour.pending.length = 0;
}

/* Voice choice matters more than anything else here. The browser's default is
 * whatever the OS considers its system voice, which on macOS is Samantha - a
 * 2009-era formant voice that makes five minutes of narration hard to sit
 * through. The modern neural voices are far better but are not always present,
 * so the pick is ranked and the user can override it. */

// British English female, in the order they are worth having. The first that
// exists on the machine wins. Chrome's "Google UK English Female" is the target
// and is what the narration was written for; the rest are the same register on
// Edge and macOS, ending with the enhanced and plain system voices.
const VOICE_WANTED = [
  'Google UK English Female',
  'Microsoft Libby Online (Natural) - English (United Kingdom)',
  'Microsoft Sonia Online (Natural) - English (United Kingdom)',
  'Microsoft Hazel - English (United Kingdom)',
  'Serena (Premium)', 'Serena (Enhanced)', 'Serena',
  'Kate (Enhanced)', 'Kate', 'Stephanie (Premium)', 'Stephanie',
  'Martha', 'Fiona',
];

// Character and alert voices. macOS lists them alongside the real ones and
// several sort early alphabetically, so any "first available English voice"
// fallback can genuinely land on Bubbles or Zarvox.
const NOVELTY = /^(albert|bad news|bahh|bells|boing|bubbles|cellos|good news|jester|junior|kathy|organ|ralph|superstar|trinoids|whisper|wobble|zarvox|fred|grandma|grandpa|deranged|hysterical|princess)/i;

// getVoices() is empty until the browser has loaded them, and on some it stays
// empty until something asks. Warm it early so the first beat is not the one
// that discovers there is no voice yet.
let voicesReady = [];
function usableVoices() {
  return (window.speechSynthesis?.getVoices() || voicesReady)
    .filter((v) => /^en/i.test(v.lang) && !NOVELTY.test(v.name));
}
function pickVoice() {
  const list = usableVoices();
  if (!list.length) return null;
  for (const want of VOICE_WANTED) {
    const hit = list.find((v) => v.name === want);
    if (hit) return hit;
  }
  // Nothing on the list is installed. Take the best British voice available,
  // preferring a neural one, and only then fall back to any English at all.
  const rank = (v) => (/premium|enhanced|neural|natural|siri/i.test(v.name) ? 2 : 0)
                    + (/^en-gb/i.test(v.lang) ? 1 : 0);
  return list.slice().sort((a, b) => rank(b) - rank(a))[0];
}
if (window.speechSynthesis) {
  const grab = () => { voicesReady = speechSynthesis.getVoices() || []; };
  grab();
  speechSynthesis.addEventListener?.('voiceschanged', grab);
}

function tourEl(id) { return $(id); }

/* Speak a beat one sentence at a time. Two reasons, both audible: a synthesiser
 * given a whole paragraph flattens its prosody and races the full stops, and
 * Chrome silently truncates utterances past about fifteen seconds. Short
 * utterances with a beat of silence between them are what makes it sound like
 * narration rather than a screen reader. */
function tourSpeak(text, done, cues) {
  const synth = window.speechSynthesis;
  const voice = pickVoice();
  const parts = (text.match(/[^.!?]+[.!?]*\s*/g) || [text]).map((s) => s.trim())
    .filter(Boolean);
  // A beat can name three stars in three sentences, and pointing at all three
  // at once says nothing. `cues[i]` runs as sentence i begins, so the marker
  // keeps step with the words.
  const fire = (i) => { try { cues?.[i]?.(); } catch (e) { /* never stall */ } };

  // No speech on this device: hold each sentence for as long as it takes to
  // read, so the captions - and the cues - still run at a sensible pace.
  if (!synth || !voice) {
    let k = 0;
    const token = ++tour.speech;
    const step = () => {
      if (token !== tour.speech) return;
      if (k >= parts.length) { done(); return; }
      fire(k);
      const ms = 500 + parts[k++].split(/\s+/).length * 385;
      tour.timer = setTimeout(step, ms);
    };
    step();
    return;
  }
  synth.cancel();
  const token = ++tour.speech;
  let k = 0;
  const next = () => {
    if (token !== tour.speech) return;               // superseded by a new beat
    if (k >= parts.length) { done(); return; }
    fire(k);
    const u = new SpeechSynthesisUtterance(parts[k++]);
    if (voice) { u.voice = voice; u.lang = voice.lang; }
    u.rate = 0.92; u.pitch = 0.98; u.volume = 1.0;
    let fired = false;
    const step = () => {
      if (fired || token !== tour.speech) return;
      fired = true;
      tour.timer = setTimeout(next, 260);            // a breath between sentences
    };
    u.onend = step;
    u.onerror = (e) => {
      if (e.error === 'interrupted' || e.error === 'canceled') return;
      step();
    };
    // Belt and braces: if the utterance never reports back, move on anyway.
    clearTimeout(tour.timer);
    tour.timer = setTimeout(step, 2500 + u.text.split(/\s+/).length * 700);
    synth.speak(u);
  };
  next();
}

function tourShow() {
  const b = TOUR[tour.i];
  tourEl('tourtext').textContent = b ? b.say : '';
  tourEl('tourstep').textContent = `${tour.i + 1} / ${TOUR.length}`;
  tourEl('tourplay').textContent = tour.paused ? '▶' : '❚❚';
  tourEl('tourplay').title = tour.paused ? 'Resume' : 'Pause';
}

function tourGo(i) {
  clearTimeout(tour.timer);
  tourClearPending();
  if (window.speechSynthesis) window.speechSynthesis.cancel();
  if (i >= TOUR.length) { tourStop(); return; }
  tour.i = i;
  const b = TOUR[i];
  try { b.act(); } catch (e) { /* a beat must never strand the tour */ }
  tourShow();
  if (tour.paused) return;
  tourSpeak(b.say, () => { if (tour.on && !tour.paused) tourGo(tour.i + 1); }, b.cues);
}

/* While the tour is driving, everything except the caption is inert.
 *
 * Dimming the rail and the orientation panel was not enough: they stayed
 * clickable, and so did the transport. Pressing Today, dragging the scrubber,
 * choosing a camera preset or picking a list row all fight the shot the
 * narration is talking over - the tour keeps animating, the user's action
 * lands on top of it, and the result looks broken rather than interrupted.
 *
 * `inert` removes a subtree from the pointer, focus and accessibility trees in
 * one attribute, and one line puts it back, so the two states cannot drift.
 * The panels stay visible: the epoch readout and the moving scrubber are part
 * of what the tour is showing. Only #tourbar stays live, which is also what
 * makes Exit Tour reachable. */
const TOUR_INERT = ['title', 'rail', 'nav', 'dock'];
function setTourInert(on) {
  for (const id of TOUR_INERT) { const el = $(id); if (el) el.inert = on; }
}

function tourStart() {
  if (tour.on) return;
  tour.on = true; tour.paused = false;
  document.body.classList.add('touring');
  setTourInert(true);
  tourEl('tourbar').hidden = false;
  if (isNarrow()) { setPanel('menu', false); setPanel('chart', false); }
  toggleAbout(false);
  fitMobileChrome();               // the dock just lost the chart; re-measure
  // Chrome pauses synthesis after ~15s of a backgrounded or long utterance;
  // a periodic resume is the documented way to keep it talking.
  tour.ping = setInterval(() => {
    const s = window.speechSynthesis;
    if (s && s.speaking && !s.paused && !tour.paused) { s.pause(); s.resume(); }
  }, 8000);
  tourGo(0);
}

function tourStop() {
  if (!tour.on) return;
  tour.on = false; tour.speech++;          // orphan any in-flight utterance
  clearTimeout(tour.timer); clearInterval(tour.ping); tourClearPending();
  if (window.speechSynthesis) window.speechSynthesis.cancel();
  document.body.classList.remove('touring');
  setTourInert(false);                     // the controls come back
  tourEl('tourbar').hidden = true;
  setPlay(0); tTween = null;
  fitMobileChrome(); pickChartStars(); invalidate();
}

function tourPause(p) {
  tour.paused = p;
  if (p) {
    clearTimeout(tour.timer); tourClearPending();
    if (window.speechSynthesis) window.speechSynthesis.cancel();
    setPlay(0); tTween = null;
  } else {
    tourGo(tour.i);          // replay the current beat rather than half of it
  }
  tourShow();
}

/* --------------------------------------------------------------------- boot */
async function grab(path, how) {
  let r;
  try {
    // Always revalidate. These three files carry the science, and a browser
    // holding a cached copy from an earlier deploy shows stale numbers with no
    // sign that anything is wrong - the page looks perfect and is out of date.
    // 'no-cache' revalidates rather than refetching, so an unchanged 1.9 MB
    // trajectory blob still costs only a 304.
    r = await fetch(path, { cache: 'no-cache' });
  } catch (e) {
    throw new Error(`could not load ${path} — serve this directory over HTTP `
      + `(python3 -m http.server 8777 --directory web) rather than opening the file directly`);
  }
  if (!r.ok) throw new Error(`${path} — HTTP ${r.status}. `
    + `Is the whole data/ directory present? cheb.bin is ~1.9 MB and is required.`);
  return how(r);
}

async function boot() {
  const [meta, stars, cheb] = await Promise.all([
    grab('data/meta.json', (r) => r.json()),
    grab('data/stars.json', (r) => r.json()),
    grab('data/cheb.bin', (r) => r.arrayBuffer()),
  ]);
  // Adopt the build's own parameters *before* validating against them. Checking
  // the blob against the fallback degree would reject a perfectly good bundle
  // fitted at any degree other than 12.
  if (meta.t_span_myr) { TSPAN = meta.t_span_myr; symK = Math.log10(1 + TSPAN / T0); }
  if (meta.cheb_degree != null) { DEG = meta.cheb_degree; NC = DEG + 1; }
  const wantBytes = stars.names.length * 3 * NC * 4;
  if (cheb.byteLength !== wantBytes) {
    throw new Error(`data/cheb.bin is ${cheb.byteLength} bytes but `
      + `${stars.names.length} stars x 3 x ${NC} float32 (degree ${DEG}) needs `
      + `${wantBytes}. The file is truncated or does not match stars.json.`);
  }
  S.meta = meta; S.d = stars; S.n = stars.names.length;
  S.coeffs = new Float32Array(cheb);
  $('nstars').textContent = S.n.toLocaleString();
  $('stage').appendChild(sunEl);

  initGL();
  if (!HAS_GL) {
    document.getElementById('nav').style.display = 'none';
    canvas.style.display = 'none';
    overlay.style.display = 'none';
    const t = document.querySelector('#title p');
    t.innerHTML = `<b style="color:var(--warn)">WebGL2 unavailable</b> — the 3-D view is off. `
      + `The distance-time chart, search and encounter tables below all still work.`;
  }
  recomputeVisible();
  placeAboutBtn();
  syncFullscreen();          // a reload inside full screen must not say "Full screen"
  fitMobileChrome();
  fitInfoCard();
  watchNav();
  pickChartStars();
  renderList();
  bind();
  setT(0);
  setPlay(0);
  $('loading').style.display = 'none';
  requestAnimationFrame(frame);
}
boot().catch((e) => { $('loading').textContent = 'error: ' + e.message; console.error(e); });
