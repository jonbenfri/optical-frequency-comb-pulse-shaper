'use strict';

const C0 = 299_792_458.0;
const MAX_LINES = 16;
const DEFAULT_LINE_COUNT = 8;
const DEFAULT_CENTER_WAVELENGTH_NM = 1550.0;
const PLOT_COLORS = {
  intensity: '#0f9f8e',
  target: '#f59e0b',
  snapshot: '#ef4444',
  spectrum: '#8a5cf6',
};

const state = {
  lineCount: DEFAULT_LINE_COUNT,
  amplitudes: Array(MAX_LINES).fill(1.0),
  phasesDeg: Array(MAX_LINES).fill(0.0),
  targetKind: null,
  targetType: null,
  targetLabel: '',
  targetRmse: null,
  presetFamily: 'manual',
  batch: false,
};

const coefficientCache = new Map();

function $(id) {
  return document.getElementById(id);
}

function lineIndices(n = state.lineCount) {
  return Array.from({ length: n }, (_, k) => k);
}

function centeredIndices(n = state.lineCount) {
  return lineIndices(n).map(k => k - (n - 1) / 2);
}

function activeAmplitudes(n = state.lineCount) {
  return state.amplitudes.slice(0, n);
}

function activePhases(n = state.lineCount) {
  return state.phasesDeg.slice(0, n);
}

function clamp(value, lo, hi) {
  return Math.min(hi, Math.max(lo, value));
}

function cleanNumber(value, decimals = 6) {
  let rounded = Number(Number(value).toFixed(decimals));
  if (Math.abs(rounded) < Math.pow(10, -decimals)) rounded = 0;
  for (const special of [-1, -0.5, 0, 0.5, 1]) {
    if (Math.abs(rounded - special) < 5 * Math.pow(10, -decimals)) rounded = special;
  }
  return rounded;
}

function formatAmp(value) {
  return cleanNumber(value, 5).toFixed(3);
}

function formatPhase(value) {
  return cleanNumber(value, 4).toFixed(1);
}

function parseNumber(value, fallback = 0) {
  const parsed = Number(String(value).replace('°', '').replaceAll(',', '').trim());
  return Number.isFinite(parsed) ? parsed : fallback;
}

function wrapPhaseToSliderRange(phaseDeg) {
  return ((phaseDeg + 180) % 360 + 360) % 360 - 180;
}

function mean(values) {
  return values.reduce((a, b) => a + b, 0) / values.length;
}

function normalizeToMax(values) {
  const maxAbs = Math.max(...values.map(v => Math.abs(v)));
  if (maxAbs === 0 || !Number.isFinite(maxAbs)) return values.slice();
  return values.map(v => v / maxAbs);
}

function rmseNormalized(values, target) {
  const y = normalizeToMax(values);
  const t = normalizeToMax(target);
  const mse = mean(y.map((v, i) => (v - t[i]) ** 2));
  return Math.sqrt(mse);
}

function combFrequencies(spacingGhz, centerWavelengthNm = DEFAULT_CENTER_WAVELENGTH_NM, nLines = state.lineCount) {
  const spacingHz = spacingGhz * 1e9;
  const centerFreqHz = C0 / (centerWavelengthNm * 1e-9);
  const offsetsHz = centeredIndices(nLines).map(idx => idx * spacingHz);
  const freqsHz = offsetsHz.map(offset => centerFreqHz + offset);
  const wavelengthsNm = freqsHz.map(freq => C0 / freq * 1e9);
  const periodPs = 1e12 / spacingHz;
  return { spacingHz, offsetsHz, freqsHz, wavelengthsNm, periodPs };
}

function targetProfile(kind, u) {
  const floor = 0.05;
  switch (kind) {
    case 'intensity_saw': return floor + (1 - floor) * u;
    case 'intensity_reverse_saw': return 1 - (1 - floor) * u;
    case 'intensity_triangle': return floor + (1 - floor) * (1 - 2 * Math.abs(u - 0.5));
    case 'intensity_square': return u < 0.5 ? 1 : floor;
    default: throw new Error(`Unknown target kind: ${kind}`);
  }
}

function fitRealFieldCoefficients(samples, nLines = state.lineCount) {
  const m = samples.length;
  const coeff = [];
  for (const k of lineIndices(nLines)) {
    let cosNum = 0;
    let cosDen = 0;
    let sinNum = 0;
    let sinDen = 0;
    for (let j = 0; j < m; j++) {
      const u = j / m;
      const theta = 2 * Math.PI * k * u;
      const cos = Math.cos(theta);
      const minusSin = -Math.sin(theta);
      cosNum += samples[j] * cos;
      cosDen += cos * cos;
      sinNum += samples[j] * minusSin;
      sinDen += minusSin * minusSin;
    }
    const a = cosDen > 0 ? cosNum / cosDen : 0;
    const b = sinDen > 1e-14 ? sinNum / sinDen : 0;
    coeff.push({ re: a, im: b });
  }
  const maxAbs = Math.max(...coeff.map(c => Math.hypot(c.re, c.im)));
  if (maxAbs > 0) {
    coeff.forEach(c => { c.re /= maxAbs; c.im /= maxAbs; });
  }
  return coeff;
}

function coefficientsForTarget(kind, type, nLines = state.lineCount) {
  const cacheKey = `${nLines}:${type}:${kind}`;
  if (coefficientCache.has(cacheKey)) return coefficientCache.get(cacheKey);

  const m = 4096;
  let samples = Array.from({ length: m }, (_, j) => targetProfile(kind, j / m));

  if (type === 'intensity') {
    samples = samples.map(v => Math.sqrt(Math.max(0, v)));
  }

  const coeff = fitRealFieldCoefficients(samples, nLines);
  coefficientCache.set(cacheKey, coeff);
  return coeff;
}

function coeffToMagnitudePhase(coeff, nLines = state.lineCount) {
  const magnitudes = coeff.map(c => cleanNumber(Math.hypot(c.re, c.im), 6));
  const phasesDeg = coeff.map((c, i) => {
    if (magnitudes[i] < 1e-8) return 0;
    return cleanNumber(Math.atan2(c.im, c.re) * 180 / Math.PI, 4);
  });
  while (magnitudes.length < nLines) magnitudes.push(0);
  while (phasesDeg.length < nLines) phasesDeg.push(0);
  return { magnitudes, phasesDeg };
}

function makeGaussianLinePreset(sigmaLines = 1.35, nLines = state.lineCount) {
  let magnitudes = centeredIndices(nLines).map(x => Math.exp(-0.5 * (x / sigmaLines) ** 2));
  const maxVal = Math.max(...magnitudes);
  magnitudes = magnitudes.map(v => cleanNumber(v / maxVal, 6));
  return { magnitudes, phasesDeg: Array(nLines).fill(0) };
}

function makeDoublePulsePreset(delayFraction = 0.35, nLines = state.lineCount) {
  const coeff = lineIndices(nLines).map(k => {
    const theta = -2 * Math.PI * k * delayFraction;
    return { re: 1 + Math.cos(theta), im: Math.sin(theta) };
  });
  const maxAbs = Math.max(...coeff.map(c => Math.hypot(c.re, c.im)));
  coeff.forEach(c => { c.re /= maxAbs; c.im /= maxAbs; });
  return coeffToMagnitudePhase(coeff, nLines);
}

function synthesizeEnvelope(magnitudes, phasesDeg, spacingGhz, nPeriods, centerWavelengthNm = DEFAULT_CENTER_WAVELENGTH_NM, samplesPerPeriod = 520, nLines = state.lineCount) {
  const { spacingHz, offsetsHz, freqsHz, wavelengthsNm, periodPs } = combFrequencies(spacingGhz, centerWavelengthNm, nLines);
  const nSamples = Math.max(900, Math.round(samplesPerPeriod * nPeriods));
  const tMin = -0.5 * nPeriods * periodPs;
  const tMax = 0.5 * nPeriods * periodPs;
  const tPs = Array.from({ length: nSamples }, (_, j) => tMin + (tMax - tMin) * j / (nSamples - 1));
  const centeredRe = Array(nSamples).fill(0);
  const centeredIm = Array(nSamples).fill(0);
  const lineRe = Array(nSamples).fill(0);
  const lineIm = Array(nSamples).fill(0);
  const indices = lineIndices(nLines);

  for (let k = 0; k < nLines; k++) {
    const amp = magnitudes[k] ?? 0;
    const phi = (phasesDeg[k] ?? 0) * Math.PI / 180;
    for (let j = 0; j < nSamples; j++) {
      const tS = tPs[j] * 1e-12;
      const thetaCentered = 2 * Math.PI * offsetsHz[k] * tS + phi;
      centeredRe[j] += amp * Math.cos(thetaCentered);
      centeredIm[j] += amp * Math.sin(thetaCentered);
      const thetaLine = 2 * Math.PI * indices[k] * spacingHz * tS + phi;
      lineRe[j] += amp * Math.cos(thetaLine);
      lineIm[j] += amp * Math.sin(thetaLine);
    }
  }

  const intensity = centeredRe.map((re, j) => re * re + centeredIm[j] * centeredIm[j]);
  return { tPs, centeredRe, centeredIm, lineRe, lineIm, intensity, offsetsHz, freqsHz, wavelengthsNm, periodPs };
}

function targetValuesForTime(tPs, periodPs, kind) {
  return tPs.map(t => {
    let u = (t / periodPs) % 1;
    if (u < 0) u += 1;
    return targetProfile(kind, u);
  });
}

function getGlobalControls() {
  const lineCount = clamp(parseNumber($('lineCountSelect').value, state.lineCount), 4, MAX_LINES);
  state.lineCount = lineCount;
  return {
    centerWavelengthNm: clamp(parseNumber($('centerWavelengthInput').value, DEFAULT_CENTER_WAVELENGTH_NM), 1200, 1700),
    spacingGhz: parseNumber($('spacingSelect').value, 200),
    lineCount,
    nPeriods: parseNumber($('periodsRange').value, 3),
    phasorTOverT: parseNumber($('phasorRange').value, 0),
    normalizePlots: $('normalizeCheckbox').checked,
    showTarget: $('targetCheckbox').checked,
  };
}


function presetFamilyLabel(family) {
  switch (family) {
    case 'intensity': return 'intensity target';
    case 'amplitude': return 'amplitude shaping';
    case 'phase': return 'phase/time';
    case 'random': return 'random';
    case 'manual':
    default: return 'manual';
  }
}

function updatePresetBadge(family = state.presetFamily) {
  const badge = $('presetBadge');
  if (!badge) return;
  const cleanFamily = ['intensity', 'amplitude', 'phase', 'random', 'manual'].includes(family) ? family : 'manual';
  badge.className = `focus-badge preset-badge family-${cleanFamily}`;
  badge.textContent = presetFamilyLabel(cleanFamily);
}

function markManualEdit() {
  state.targetKind = null;
  state.targetType = null;
  state.targetLabel = '';
  state.targetRmse = null;
  state.presetFamily = 'manual';
  updatePresetBadge('manual');
}

function applyStateToControls() {
  state.batch = true;
  $('lineCountSelect').value = String(state.lineCount);
  for (let k = 0; k < MAX_LINES; k++) {
    const row = $(`lineRow${k}`);
    if (row) row.classList.toggle('line-hidden', k >= state.lineCount);
    const amp = cleanNumber(clamp(state.amplitudes[k], 0, 5), 6);
    const phase = cleanNumber(state.phasesDeg[k], 4);
    const ampRange = $(`ampRange${k}`);
    const ampNumber = $(`ampNumber${k}`);
    const phaseRange = $(`phaseRange${k}`);
    const phaseNumber = $(`phaseNumber${k}`);
    if (!ampRange) continue;
    ampRange.value = String(clamp(amp, 0, 1.5));
    ampNumber.value = formatAmp(amp);
    phaseRange.value = String(wrapPhaseToSliderRange(phase));
    phaseNumber.value = formatPhase(phase);
  }
  state.batch = false;
}

function readLineControls() {
  for (let k = 0; k < state.lineCount; k++) {
    state.amplitudes[k] = clamp(parseNumber($(`ampNumber${k}`).value, state.amplitudes[k]), 0, 5);
    state.phasesDeg[k] = parseNumber($(`phaseNumber${k}`).value, state.phasesDeg[k]);
  }
}

function updateSummary(lineCount, centerWavelengthNm, spacingGhz, periodPs, wavelengthsNm) {
  const dlambdas = [];
  for (let i = 1; i < wavelengthsNm.length; i++) dlambdas.push(Math.abs(wavelengthsNm[i] - wavelengthsNm[i - 1]));
  $('lineCountText').textContent = String(lineCount);
  $('centerWavelengthText').textContent = `${centerWavelengthNm.toFixed(1)} nm`;
  $('spacingText').textContent = `${spacingGhz} GHz`;
  $('periodText').textContent = `${periodPs.toFixed(3)} ps`;
  $('lambdaSpacingText').textContent = dlambdas.length ? `${mean(dlambdas).toFixed(4)} nm` : '—';
  $('rmseText').textContent = state.targetRmse === null ? '—' : `normalized RMSE ≈ ${state.targetRmse.toFixed(3)}`;
}

function currentTargetForPlot(synth, controls) {
  state.targetRmse = null;
  if (!controls.showTarget || state.targetKind === null) return null;

  const targetRaw = targetValuesForTime(synth.tPs, synth.periodPs, state.targetKind);
  let compareRaw;
  if (state.targetType === 'intensity') compareRaw = synth.intensity;
  else return null;

  state.targetRmse = rmseNormalized(compareRaw, targetRaw);
  return controls.normalizePlots ? normalizeToMax(targetRaw) : targetRaw;
}

function makeSpectrumPlot(synth, controls) {
  const amps = activeAmplitudes(controls.lineCount);
  const phases = activePhases(controls.lineCount);
  const offsetsGhz = synth.offsetsHz.map(v => v / 1e9);
  const text = amps.map((amp, k) => `${synth.wavelengthsNm[k].toFixed(2)} nm<br>φ=${phases[k].toFixed(1)}°`);
  const yMax = Math.max(1.55, 1.1 * Math.max(...amps, 1));
  const trace = {
    type: 'bar',
    x: offsetsGhz,
    y: amps,
    text,
    textposition: 'outside',
    hovertemplate: 'Offset %{x:.1f} GHz<br>Amplitude %{y:.3f}<br>%{text}<extra></extra>',
    marker: { color: PLOT_COLORS.spectrum, line: { width: 1 } },
  };
  const layout = {
    margin: { l: 58, r: 16, t: 20, b: 62 },
    xaxis: { title: `Frequency offset from ${controls.centerWavelengthNm.toFixed(1)}-nm carrier (GHz)`, zeroline: true },
    yaxis: { title: 'Line amplitude', range: [0, yMax] },
    bargap: 0.35,
    showlegend: false,
  };
  Plotly.react('spectrumPlot', [trace], layout, plotConfig());
}

function makePhasorPlot(synth, controls) {
  const amps = activeAmplitudes(controls.lineCount);
  const phases = activePhases(controls.lineCount);
  const offsetsGhz = synth.offsetsHz.map(v => v / 1e9);
  const phasorTimePs = controls.phasorTOverT * synth.periodPs;
  const phasorTimeS = phasorTimePs * 1e-12;
  const maxMag = Math.max(1, ...amps);
  const traces = [];

  traces.push({
    type: 'scatter3d',
    mode: 'lines',
    x: [Math.min(...offsetsGhz) - controls.spacingGhz, Math.max(...offsetsGhz) + controls.spacingGhz],
    y: [0, 0],
    z: [0, 0],
    line: { width: 4 },
    hoverinfo: 'skip',
    showlegend: false,
  });

  for (let k = 0; k < controls.lineCount; k++) {
    const theta = phases[k] * Math.PI / 180 + 2 * Math.PI * synth.offsetsHz[k] * phasorTimeS;
    const y = amps[k] * Math.cos(theta);
    const z = amps[k] * Math.sin(theta);
    traces.push({
      type: 'scatter3d',
      mode: 'lines+markers+text',
      x: [offsetsGhz[k], offsetsGhz[k]],
      y: [0, y],
      z: [0, z],
      text: ['', String(k)],
      textposition: 'top center',
      marker: { size: [2, 5] },
      line: { width: 6 },
      hovertemplate: `Line ${k}<br>Offset ${offsetsGhz[k].toFixed(1)} GHz<br>Amplitude ${amps[k].toFixed(3)}<br>Phase ${phases[k].toFixed(1)}°<extra></extra>`,
      showlegend: false,
    });
  }

  const circleX = [];
  const circleY = [];
  const circleZ = [];
  for (let j = 0; j <= 160; j++) {
    const theta = 2 * Math.PI * j / 160;
    circleX.push(0);
    circleY.push(maxMag * Math.cos(theta));
    circleZ.push(maxMag * Math.sin(theta));
  }
  traces.push({
    type: 'scatter3d',
    mode: 'lines',
    x: circleX,
    y: circleY,
    z: circleZ,
    line: { width: 2, dash: 'dot' },
    hoverinfo: 'skip',
    showlegend: false,
  });

  const layout = {
    margin: { l: 18, r: 34, t: 38, b: 28 },
    scene: {
      xaxis: {
        title: { text: 'Freq. offset<br>(GHz)', font: { size: 11 } },
        tickfont: { size: 10 },
      },
      yaxis: {
        title: { text: 'In-phase', font: { size: 11 } },
        tickfont: { size: 10 },
        range: [-1.25 * maxMag, 1.25 * maxMag],
      },
      zaxis: {
        title: { text: 'Quadrature', font: { size: 11 } },
        tickfont: { size: 10 },
        range: [-1.25 * maxMag, 1.25 * maxMag],
      },
      camera: { eye: { x: 1.55, y: -1.75, z: 1.15 } },
    },
    title: {
      text: `Phasors at t = ${phasorTimePs.toFixed(3)} ps = ${controls.phasorTOverT.toFixed(2)}T`,
      font: { size: 13 },
      y: 0.97,
    },
  };
  Plotly.react('phasorPlot', traces, layout, plotConfig());
}

function makeTimePlot(synth, controls) {
  const intensityPlot = controls.normalizePlots ? normalizeToMax(synth.intensity) : synth.intensity;
  const traces = [{
    type: 'scatter',
    mode: 'lines',
    x: synth.tPs,
    y: intensityPlot,
    name: controls.normalizePlots ? 'Normalized intensity |E(t)|²' : 'Intensity |E(t)|²',
    line: { width: 3, color: PLOT_COLORS.intensity },
    hovertemplate: 't = %{x:.3f} ps<br>%{y:.4g}<extra></extra>',
  }];

  const targetPlot = currentTargetForPlot(synth, controls);
  if (targetPlot) {
    traces.push({
      type: 'scatter',
      mode: 'lines',
      x: synth.tPs,
      y: targetPlot,
      name: `Target: ${state.targetLabel}`,
      line: { width: 2.5, dash: 'dot', color: PLOT_COLORS.target },
      hovertemplate: 't = %{x:.3f} ps<br>target %{y:.4g}<extra></extra>',
    });
  }

  const phasorTimePs = controls.phasorTOverT * synth.periodPs;
  const tMin = Math.min(...synth.tPs);
  let displayTime = tMin + (((phasorTimePs - tMin) % synth.periodPs + synth.periodPs) % synth.periodPs);
  traces.push({
    type: 'scatter',
    mode: 'lines',
    x: [displayTime, displayTime],
    y: controls.normalizePlots ? [-0.08, 1.08] : [0, Math.max(...synth.intensity) * 1.05],
    name: 'phasor snapshot time',
    line: { width: 1.5, dash: 'dashdot', color: PLOT_COLORS.snapshot },
    hoverinfo: 'skip',
  });

  const yRange = controls.normalizePlots ? [-0.10, 1.12] : undefined;
  const rmseText = state.targetRmse === null ? '' : `target RMSE ≈ ${state.targetRmse.toFixed(3)}`;
  const layout = {
    margin: { l: 66, r: 16, t: 20, b: 78 },
    xaxis: { title: `Time relative to envelope, T = ${synth.periodPs.toFixed(3)} ps` },
    yaxis: { title: controls.normalizePlots ? 'Normalized intensity' : 'Intensity |E(t)|²', range: yRange },
    legend: { orientation: 'h', y: -0.25, x: 0, yanchor: 'top' },
    annotations: rmseText ? [{
      text: rmseText,
      xref: 'paper', yref: 'paper', x: 1, y: 1.06,
      xanchor: 'right', showarrow: false,
      font: { size: 12 },
    }] : [],
  };
  Plotly.react('timePlot', traces, layout, plotConfig());
}

function plotConfig() {
  return {
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
  };
}

let pendingPlot = null;
function scheduleUpdate() {
  if (state.batch) return;
  if (pendingPlot !== null) cancelAnimationFrame(pendingPlot);
  pendingPlot = requestAnimationFrame(() => {
    pendingPlot = null;
    updatePlots();
  });
}

function updatePlots() {
  if (typeof Plotly === 'undefined') {
    const status = $('plotStatus');
    if (status) status.textContent = 'Plotly did not load. Check the network connection or CDN access.';
    return;
  }
  const status = $('plotStatus');
  if (status) status.textContent = '';
  readLineControls();
  const controls = getGlobalControls();
  $('periodsOutput').value = controls.nPeriods;
  $('phasorOutput').value = controls.phasorTOverT.toFixed(2);
  $('centerWavelengthInput').value = controls.centerWavelengthNm.toFixed(1);
  $('lineCountSelect').value = String(controls.lineCount);

  const synth = synthesizeEnvelope(activeAmplitudes(controls.lineCount), activePhases(controls.lineCount), controls.spacingGhz, controls.nPeriods, controls.centerWavelengthNm, 520, controls.lineCount);
  makeSpectrumPlot(synth, controls);
  makePhasorPlot(synth, controls);
  makeTimePlot(synth, controls);
  updateSummary(controls.lineCount, controls.centerWavelengthNm, controls.spacingGhz, synth.periodPs, synth.wavelengthsNm);
}

function setValues({ magnitudes, phasesDeg, note = '', targetKind = null, targetType = null, targetLabel = '', family = 'manual' }) {
  for (let k = 0; k < MAX_LINES; k++) {
    state.amplitudes[k] = k < magnitudes.length ? clamp(cleanNumber(magnitudes[k], 6), 0, 5) : 0;
    state.phasesDeg[k] = k < phasesDeg.length ? cleanNumber(phasesDeg[k], 4) : 0;
  }
  state.targetKind = targetKind;
  state.targetType = targetType;
  state.targetLabel = targetLabel;
  state.targetRmse = null;
  state.presetFamily = family;
  updatePresetBadge(family);
  if (note) $('presetNote').innerHTML = note;
  applyStateToControls();
  updatePlots();
}

function amplitudePreset(magnitudes, label, extraNote = '') {
  return {
    magnitudes: magnitudes.map(v => cleanNumber(v, 6)),
    phasesDeg: Array(magnitudes.length).fill(0),
    note: ('<strong>' + label + ':</strong> flat spectral phase with shaped line amplitudes. ' + extraNote).trim(),
    family: 'amplitude',
  };
}

function makeLinearAmplitudeRamp(direction = 'up', nLines = state.lineCount) {
  const lo = 0.18;
  const hi = 1.0;
  let magnitudes = lineIndices(nLines).map(k => lo + (hi - lo) * k / Math.max(nLines - 1, 1));
  if (direction === 'down') magnitudes = magnitudes.slice().reverse();
  return magnitudes;
}

function makeCenterWeightedAmplitudes(sigmaLines = null, nLines = state.lineCount) {
  const sigma = sigmaLines ?? Math.max(1.05, 0.22 * nLines);
  let magnitudes = centeredIndices(nLines).map(x => Math.exp(-0.5 * (x / sigma) ** 2));
  const maxVal = Math.max(...magnitudes);
  return magnitudes.map(v => v / maxVal);
}

function makeEdgeWeightedAmplitudes(power = 1.2, nLines = state.lineCount) {
  const centers = centeredIndices(nLines);
  const maxAbs = Math.max(...centers.map(Math.abs), 1);
  return centers.map(x => 0.18 + 0.82 * (Math.abs(x) / maxAbs) ** power);
}

function makeCenterNotchAmplitudes(nLines = state.lineCount) {
  const centers = centeredIndices(nLines);
  const maxAbs = Math.max(...centers.map(Math.abs), 1);
  return centers.map(x => {
    const d = Math.abs(x) / maxAbs;
    return 0.12 + 0.88 * (1 - Math.exp(-4.2 * d ** 2));
  });
}

function makeAlternatingAmplitudePattern(nLines = state.lineCount) {
  return lineIndices(nLines).map(k => (k % 2 === 0 ? 1.0 : 0.22));
}

function intensityTargetPreset(kind, label) {
  const coeff = coefficientsForTarget(kind, 'intensity', state.lineCount);
  const { magnitudes, phasesDeg } = coeffToMagnitudePhase(coeff, state.lineCount);
  return {
    magnitudes,
    phasesDeg,
    note: `<strong>Intensity target: ${label}.</strong> Compare the solid intensity curve with the dotted target. Finite ${state.lineCount}-line bandwidth causes rounded edges and ripple.`,
    targetKind: kind,
    targetType: 'intensity',
    targetLabel: `${label} intensity target`,
    family: 'intensity',
  };
}

const presetFunctions = {
  'Transform-limited pulse train': () => ({
    magnitudes: Array(state.lineCount).fill(1),
    phasesDeg: Array(state.lineCount).fill(0),
    note: `<strong>Transform-limited:</strong> equal amplitudes and flat spectral phase. Produces the shortest pulse train for this flat ${state.lineCount}-line spectrum.`,
    family: 'amplitude',
  }),
  'Time shift: linear phase ramp': () => ({
    magnitudes: Array(state.lineCount).fill(1),
    phasesDeg: lineIndices(state.lineCount).map(k => -360 * k * 0.25),
    note: '<strong>Linear phase ramp:</strong> shifts the waveform in time while preserving the intensity shape.',
    family: 'phase',
  }),
  'Intensity target: sawtooth': () => intensityTargetPreset('intensity_saw', 'sawtooth'),
  'Intensity target: reverse sawtooth': () => intensityTargetPreset('intensity_reverse_saw', 'reverse sawtooth'),
  'Intensity target: triangle': () => intensityTargetPreset('intensity_triangle', 'triangle'),
  'Intensity target: square': () => intensityTargetPreset('intensity_square', 'square'),
  'Amplitude ramp: low → high frequency': () => amplitudePreset(
    makeLinearAmplitudeRamp('up', state.lineCount),
    'Amplitude ramp, low to high frequency',
    'This emphasizes the high-frequency side of the comb and gives an asymmetric spectral envelope.'
  ),
  'Amplitude ramp: high → low frequency': () => amplitudePreset(
    makeLinearAmplitudeRamp('down', state.lineCount),
    'Amplitude ramp, high to low frequency',
    'This emphasizes the low-frequency side of the comb and gives the opposite spectral tilt.'
  ),
  'Amplitude chirp: center weighted': () => amplitudePreset(
    makeCenterWeightedAmplitudes(null, state.lineCount),
    'Center-weighted amplitude chirp',
    'Most power sits near the carrier; the time waveform broadens and sidelobes are reduced.'
  ),
  'Amplitude chirp: edge weighted': () => amplitudePreset(
    makeEdgeWeightedAmplitudes(1.2, state.lineCount),
    'Edge-weighted amplitude chirp',
    'Most power sits at the comb edges, which tends to sharpen features and increase ripple.'
  ),
  'Amplitude notch: weak center lines': () => amplitudePreset(
    makeCenterNotchAmplitudes(state.lineCount),
    'Center-notched amplitude chirp',
    'The center lines are suppressed, producing a split-spectrum example with a more structured intensity envelope.'
  ),
  'Amplitude comb: alternating strong / weak': () => amplitudePreset(
    makeAlternatingAmplitudePattern(state.lineCount),
    'Alternating strong/weak amplitudes',
    'Every other line is suppressed, so the intensity envelope visibly changes periodic structure.'
  ),
  'Gaussian-like line amplitudes': () => {
    const sigma = Math.max(1.35, 0.24 * state.lineCount);
    const { magnitudes, phasesDeg } = makeGaussianLinePreset(sigma, state.lineCount);
    return {
      magnitudes,
      phasesDeg,
      note: `<strong>Gaussian-like spectrum:</strong> tapered line amplitudes reduce sidelobes compared with a flat ${state.lineCount}-line spectrum.`,
      family: 'amplitude',
    };
  },
  'Double-pulse-like intensity': () => {
    const { magnitudes, phasesDeg } = makeDoublePulsePreset(0.35, state.lineCount);
    return {
      magnitudes,
      phasesDeg,
      note: '<strong>Double-pulse-like intensity:</strong> spectral interference factor approximates two pulses separated by about 0.35T.',
      family: 'phase',
    };
  },
  'Alternating 0 / π phase': () => ({
    magnitudes: Array(state.lineCount).fill(1),
    phasesDeg: lineIndices(state.lineCount).map(k => k % 2 === 0 ? 0 : 180),
    note: '<strong>Alternating 0/π phase:</strong> flips the sign of alternating spectral lines, producing a strongly reshaped periodic envelope.',
    family: 'phase',
  }),
  'Random amplitudes and phases': () => ({
    magnitudes: lineIndices(state.lineCount).map(() => cleanNumber(0.45 + Math.random() * 0.75, 4)),
    phasesDeg: lineIndices(state.lineCount).map(() => cleanNumber(-180 + Math.random() * 360, 3)),
    note: '<strong>Random amplitudes and phases:</strong> a quick way to explore arbitrary line-by-line settings.',
    family: 'random',
  }),
};

function buildLineControls() {
  const table = $('lineTable');
  table.innerHTML = `
    <div class="line-header">
      <span class="column-header">Line</span>
      <span class="column-header">Amplitude slider</span>
      <span class="column-header">Amplitude</span>
      <span class="column-header phase-header">Phase slider</span>
      <span class="column-header phase-number-header">Phase (deg)</span>
    </div>
  `;

  for (let k = 0; k < MAX_LINES; k++) {
    const row = document.createElement('div');
    row.className = 'line-row';
    row.id = `lineRow${k}`;
    row.innerHTML = `
      <span class="line-index">${k}</span>
      <input id="ampRange${k}" type="range" min="0" max="1.5" step="0.01" value="1" aria-label="Amplitude slider for line ${k}" />
      <input id="ampNumber${k}" type="number" min="0" max="5" step="0.001" value="1.000" aria-label="Amplitude value for line ${k}" />
      <input id="phaseRange${k}" class="phase-slider" type="range" min="-180" max="180" step="1" value="0" aria-label="Phase slider for line ${k}" />
      <input id="phaseNumber${k}" class="phase-number" type="number" step="0.1" value="0.0" aria-label="Phase value for line ${k} in degrees" />
    `;
    table.appendChild(row);

    const ampRange = $(`ampRange${k}`);
    const ampNumber = $(`ampNumber${k}`);
    const phaseRange = $(`phaseRange${k}`);
    const phaseNumber = $(`phaseNumber${k}`);

    ampRange.addEventListener('input', () => {
      ampNumber.value = formatAmp(ampRange.value);
      markManualEdit();
      scheduleUpdate();
    });
    ampNumber.addEventListener('change', () => {
      const amp = clamp(parseNumber(ampNumber.value, state.amplitudes[k]), 0, 5);
      ampNumber.value = formatAmp(amp);
      ampRange.value = String(clamp(amp, 0, 1.5));
      markManualEdit();
      scheduleUpdate();
    });

    phaseRange.addEventListener('input', () => {
      phaseNumber.value = formatPhase(phaseRange.value);
      markManualEdit();
      scheduleUpdate();
    });
    phaseNumber.addEventListener('change', () => {
      const phase = parseNumber(phaseNumber.value, state.phasesDeg[k]);
      phaseNumber.value = formatPhase(phase);
      phaseRange.value = String(wrapPhaseToSliderRange(phase));
      markManualEdit();
      scheduleUpdate();
    });
  }
}


function presetFamilyForName(name) {
  if (name.startsWith('Intensity target:')) return 'intensity';
  if (name.startsWith('Amplitude ') || name.startsWith('Gaussian-like') || name.startsWith('Transform-limited')) return 'amplitude';
  if (name.includes('phase') || name.includes('Time shift') || name.startsWith('Double-pulse')) return 'phase';
  if (name.startsWith('Random')) return 'random';
  return 'manual';
}

function buildPresetMenu() {
  const select = $('presetSelect');
  for (const name of Object.keys(presetFunctions)) {
    const option = document.createElement('option');
    option.value = name;
    option.textContent = name;
    select.appendChild(option);
  }
}

function attachGlobalEvents() {
  $('centerWavelengthInput').addEventListener('change', scheduleUpdate);
  $('spacingSelect').addEventListener('change', scheduleUpdate);
  $('lineCountSelect').addEventListener('change', () => {
    state.lineCount = clamp(parseNumber($('lineCountSelect').value, state.lineCount), 4, MAX_LINES);
    applyStateToControls();
    scheduleUpdate();
  });
  $('periodsRange').addEventListener('input', () => {
    $('periodsOutput').value = $('periodsRange').value;
    scheduleUpdate();
  });
  $('phasorRange').addEventListener('input', () => {
    $('phasorOutput').value = Number($('phasorRange').value).toFixed(2);
    scheduleUpdate();
  });
  $('normalizeCheckbox').addEventListener('change', scheduleUpdate);
  $('targetCheckbox').addEventListener('change', scheduleUpdate);

  $('presetSelect').addEventListener('change', () => {
    updatePresetBadge(presetFamilyForName($('presetSelect').value));
  });

  $('applyPresetButton').addEventListener('click', () => {
    const preset = presetFunctions[$('presetSelect').value]();
    setValues(preset);
  });

  $('resetButton').addEventListener('click', () => {
    $('presetSelect').value = 'Transform-limited pulse train';
    setValues(presetFunctions['Transform-limited pulse train']());
  });

  $('exportButton').addEventListener('click', () => {
    readLineControls();
    const payload = {
      description: `${state.lineCount}-line optical pulse shaping settings`,
      lineCount: state.lineCount,
      centerWavelengthNm: clamp(parseNumber($('centerWavelengthInput').value, DEFAULT_CENTER_WAVELENGTH_NM), 1200, 1700),
      spacingGhz: parseNumber($('spacingSelect').value, 200),
      amplitudes: activeAmplitudes(state.lineCount).map(v => cleanNumber(v, 6)),
      phasesDeg: activePhases(state.lineCount).map(v => cleanNumber(v, 4)),
      targetKind: state.targetKind,
      targetType: state.targetType,
      targetLabel: state.targetLabel,
      presetFamily: state.presetFamily,
    };
    const blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = 'pulse-shaping-settings.json';
    link.click();
    URL.revokeObjectURL(url);
  });
}

function init() {
  buildLineControls();
  buildPresetMenu();
  attachGlobalEvents();
  applyStateToControls();
  setValues(presetFunctions['Transform-limited pulse train']());
}

window.addEventListener('DOMContentLoaded', init);
