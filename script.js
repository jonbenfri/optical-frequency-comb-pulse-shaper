'use strict';

const C0 = 299_792_458.0;
const N_LINES = 8;
const CENTER_WAVELENGTH_NM = 1550.0;
const CENTER_WAVELENGTH_M = CENTER_WAVELENGTH_NM * 1e-9;
const CENTER_FREQ_HZ = C0 / CENTER_WAVELENGTH_M;
const lineIndices = Array.from({ length: N_LINES }, (_, k) => k);
const centeredIndices = lineIndices.map(k => k - (N_LINES - 1) / 2);

const state = {
  amplitudes: Array(N_LINES).fill(1.0),
  phasesDeg: Array(N_LINES).fill(0.0),
  targetKind: null,
  targetType: null,
  targetLabel: '',
  targetRmse: null,
  batch: false,
};

const coefficientCache = new Map();

function $(id) {
  return document.getElementById(id);
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

function combFrequencies(spacingGhz) {
  const spacingHz = spacingGhz * 1e9;
  const offsetsHz = centeredIndices.map(idx => idx * spacingHz);
  const freqsHz = offsetsHz.map(offset => CENTER_FREQ_HZ + offset);
  const wavelengthsNm = freqsHz.map(freq => C0 / freq * 1e9);
  const periodPs = 1e12 / spacingHz;
  return { spacingHz, offsetsHz, freqsHz, wavelengthsNm, periodPs };
}

function targetProfile(kind, u) {
  const floor = 0.05;
  switch (kind) {
    case 'field_saw': return 2 * u - 1;
    case 'field_reverse_saw': return 1 - 2 * u;
    case 'field_triangle': return 1 - 4 * Math.abs(u - 0.5);
    case 'field_square': return u < 0.5 ? 1 : -1;
    case 'intensity_saw': return floor + (1 - floor) * u;
    case 'intensity_reverse_saw': return 1 - (1 - floor) * u;
    case 'intensity_triangle': return floor + (1 - floor) * (1 - 2 * Math.abs(u - 0.5));
    case 'intensity_square': return u < 0.5 ? 1 : floor;
    default: throw new Error(`Unknown target kind: ${kind}`);
  }
}

function fitRealFieldCoefficients(samples) {
  // Fit y(u) with Re{sum_k C_k exp(i 2πku)}, k = 0..7.
  // Continuous orthogonality is approximated by dense uniform quadrature.
  const m = samples.length;
  const coeff = [];
  for (const k of lineIndices) {
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

function coefficientsForTarget(kind, type) {
  const cacheKey = `${type}:${kind}`;
  if (coefficientCache.has(cacheKey)) return coefficientCache.get(cacheKey);

  const m = 4096;
  let samples = Array.from({ length: m }, (_, j) => targetProfile(kind, j / m));

  if (type === 'field') {
    const avg = mean(samples);
    samples = normalizeToMax(samples.map(v => v - avg));
  } else if (type === 'intensity') {
    samples = samples.map(v => Math.sqrt(Math.max(0, v)));
  }

  const coeff = fitRealFieldCoefficients(samples);
  coefficientCache.set(cacheKey, coeff);
  return coeff;
}

function coeffToMagnitudePhase(coeff) {
  const magnitudes = coeff.map(c => cleanNumber(Math.hypot(c.re, c.im), 6));
  const phasesDeg = coeff.map((c, i) => {
    if (magnitudes[i] < 1e-8) return 0;
    return cleanNumber(Math.atan2(c.im, c.re) * 180 / Math.PI, 4);
  });
  return { magnitudes, phasesDeg };
}

function makeGaussianLinePreset(sigmaLines = 1.35) {
  let magnitudes = centeredIndices.map(x => Math.exp(-0.5 * (x / sigmaLines) ** 2));
  const maxVal = Math.max(...magnitudes);
  magnitudes = magnitudes.map(v => cleanNumber(v / maxVal, 6));
  return { magnitudes, phasesDeg: Array(N_LINES).fill(0) };
}

function makeDoublePulsePreset(delayFraction = 0.35) {
  const coeff = lineIndices.map(k => {
    const theta = -2 * Math.PI * k * delayFraction;
    return { re: 1 + Math.cos(theta), im: Math.sin(theta) };
  });
  const maxAbs = Math.max(...coeff.map(c => Math.hypot(c.re, c.im)));
  coeff.forEach(c => { c.re /= maxAbs; c.im /= maxAbs; });
  return coeffToMagnitudePhase(coeff);
}

function synthesizeEnvelope(magnitudes, phasesDeg, spacingGhz, nPeriods, samplesPerPeriod = 620) {
  const { spacingHz, offsetsHz, freqsHz, wavelengthsNm, periodPs } = combFrequencies(spacingGhz);
  const nSamples = Math.max(900, Math.round(samplesPerPeriod * nPeriods));
  const tMin = -0.5 * nPeriods * periodPs;
  const tMax = 0.5 * nPeriods * periodPs;
  const tPs = Array.from({ length: nSamples }, (_, j) => tMin + (tMax - tMin) * j / (nSamples - 1));
  const centeredRe = Array(nSamples).fill(0);
  const centeredIm = Array(nSamples).fill(0);
  const lineRe = Array(nSamples).fill(0);
  const lineIm = Array(nSamples).fill(0);

  for (let k = 0; k < N_LINES; k++) {
    const amp = magnitudes[k];
    const phi = phasesDeg[k] * Math.PI / 180;
    for (let j = 0; j < nSamples; j++) {
      const tS = tPs[j] * 1e-12;
      const thetaCentered = 2 * Math.PI * offsetsHz[k] * tS + phi;
      centeredRe[j] += amp * Math.cos(thetaCentered);
      centeredIm[j] += amp * Math.sin(thetaCentered);
      const thetaLine = 2 * Math.PI * lineIndices[k] * spacingHz * tS + phi;
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
  return {
    spacingGhz: parseNumber($('spacingSelect').value, 200),
    nPeriods: parseNumber($('periodsRange').value, 3),
    phasorTOverT: parseNumber($('phasorRange').value, 0),
    normalizePlots: $('normalizeCheckbox').checked,
    showField: $('fieldCheckbox').checked,
    showTarget: $('targetCheckbox').checked,
  };
}

function applyStateToControls() {
  state.batch = true;
  for (let k = 0; k < N_LINES; k++) {
    const amp = cleanNumber(clamp(state.amplitudes[k], 0, 5), 6);
    const phase = cleanNumber(state.phasesDeg[k], 4);
    $(`ampRange${k}`).value = String(clamp(amp, 0, 1.5));
    $(`ampNumber${k}`).value = formatAmp(amp);
    $(`phaseRange${k}`).value = String(wrapPhaseToSliderRange(phase));
    $(`phaseNumber${k}`).value = formatPhase(phase);
  }
  state.batch = false;
}

function readLineControls() {
  for (let k = 0; k < N_LINES; k++) {
    state.amplitudes[k] = clamp(parseNumber($(`ampNumber${k}`).value, state.amplitudes[k]), 0, 5);
    state.phasesDeg[k] = parseNumber($(`phaseNumber${k}`).value, state.phasesDeg[k]);
  }
}

function updateSummary(spacingGhz, periodPs, wavelengthsNm) {
  const dlambdas = [];
  for (let i = 1; i < wavelengthsNm.length; i++) dlambdas.push(Math.abs(wavelengthsNm[i] - wavelengthsNm[i - 1]));
  $('centerWavelengthText').textContent = `${CENTER_WAVELENGTH_NM.toFixed(1)} nm`;
  $('spacingText').textContent = `${spacingGhz} GHz`;
  $('periodText').textContent = `${periodPs.toFixed(3)} ps`;
  $('lambdaSpacingText').textContent = `${mean(dlambdas).toFixed(4)} nm`;
  $('rmseText').textContent = state.targetRmse === null ? '—' : `normalized RMSE ≈ ${state.targetRmse.toFixed(3)}`;
}

function currentTargetForPlot(synth, controls) {
  state.targetRmse = null;
  if (!controls.showTarget || state.targetKind === null) return null;

  const targetRaw = targetValuesForTime(synth.tPs, synth.periodPs, state.targetKind);
  let compareRaw;
  if (state.targetType === 'field') compareRaw = synth.lineRe;
  else if (state.targetType === 'intensity') compareRaw = synth.intensity;
  else return null;

  state.targetRmse = rmseNormalized(compareRaw, targetRaw);
  return controls.normalizePlots ? normalizeToMax(targetRaw) : targetRaw;
}

function makeSpectrumPlot(synth, controls) {
  const offsetsGhz = synth.offsetsHz.map(v => v / 1e9);
  const phases = state.phasesDeg;
  const text = state.amplitudes.map((amp, k) => `${synth.wavelengthsNm[k].toFixed(2)} nm<br>φ=${phases[k].toFixed(1)}°`);
  const yMax = Math.max(1.55, 1.1 * Math.max(...state.amplitudes));
  const trace = {
    type: 'bar',
    x: offsetsGhz,
    y: state.amplitudes,
    text,
    textposition: 'outside',
    hovertemplate: 'Offset %{x:.1f} GHz<br>Amplitude %{y:.3f}<br>%{text}<extra></extra>',
    marker: { line: { width: 1 } },
  };
  const layout = {
    margin: { l: 58, r: 16, t: 20, b: 62 },
    xaxis: { title: 'Frequency offset from 1550-nm carrier (GHz)', zeroline: true },
    yaxis: { title: 'Line amplitude', range: [0, yMax] },
    bargap: 0.35,
    showlegend: false,
  };
  Plotly.react('spectrumPlot', [trace], layout, plotConfig());
}

function makePhasorPlot(synth, controls) {
  const offsetsGhz = synth.offsetsHz.map(v => v / 1e9);
  const phasorTimePs = controls.phasorTOverT * synth.periodPs;
  const phasorTimeS = phasorTimePs * 1e-12;
  const maxMag = Math.max(1, ...state.amplitudes);
  const traces = [];

  traces.push({
    type: 'scatter3d',
    mode: 'lines',
    x: [Math.min(...offsetsGhz) - controls.spacingGhz, Math.max(...offsetsGhz) + controls.spacingGhz],
    y: [0, 0],
    z: [0, 0],
    line: { width: 4 },
    hoverinfo: 'skip',
    name: 'frequency axis',
    showlegend: false,
  });

  for (let k = 0; k < N_LINES; k++) {
    const theta = state.phasesDeg[k] * Math.PI / 180 + 2 * Math.PI * synth.offsetsHz[k] * phasorTimeS;
    const y = state.amplitudes[k] * Math.cos(theta);
    const z = state.amplitudes[k] * Math.sin(theta);
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
      name: `Line ${k}`,
      hovertemplate: `Line ${k}<br>Offset ${offsetsGhz[k].toFixed(1)} GHz<br>Amplitude ${state.amplitudes[k].toFixed(3)}<br>Phase ${state.phasesDeg[k].toFixed(1)}°<extra></extra>`,
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
    name: 'unit phase circle',
    showlegend: false,
  });

  const layout = {
    margin: { l: 0, r: 0, t: 24, b: 0 },
    scene: {
      xaxis: { title: 'Frequency offset (GHz)' },
      yaxis: { title: 'In-phase', range: [-1.25 * maxMag, 1.25 * maxMag] },
      zaxis: { title: 'Quadrature', range: [-1.25 * maxMag, 1.25 * maxMag] },
      camera: { eye: { x: 1.5, y: -1.7, z: 1.1 } },
    },
    title: { text: `Phasors at t = ${phasorTimePs.toFixed(3)} ps = ${controls.phasorTOverT.toFixed(2)}T`, font: { size: 14 } },
  };
  Plotly.react('phasorPlot', traces, layout, plotConfig());
}

function makeTimePlot(synth, controls) {
  const intensityPlot = controls.normalizePlots ? normalizeToMax(synth.intensity) : synth.intensity;
  const fieldPlot = controls.normalizePlots ? normalizeToMax(synth.lineRe) : synth.lineRe;
  const traces = [{
    type: 'scatter',
    mode: 'lines',
    x: synth.tPs,
    y: intensityPlot,
    name: controls.normalizePlots ? 'Normalized intensity |E(t)|²' : 'Intensity |E(t)|²',
    line: { width: 3 },
    hovertemplate: 't = %{x:.3f} ps<br>%{y:.4g}<extra></extra>',
  }];

  if (controls.showField) {
    traces.push({
      type: 'scatter',
      mode: 'lines',
      x: synth.tPs,
      y: fieldPlot,
      name: controls.normalizePlots ? 'Normalized Re{field}' : 'Re{field}',
      line: { width: 2, dash: 'dash' },
      hovertemplate: 't = %{x:.3f} ps<br>%{y:.4g}<extra></extra>',
    });
  }

  const targetPlot = currentTargetForPlot(synth, controls);
  if (targetPlot) {
    traces.push({
      type: 'scatter',
      mode: 'lines',
      x: synth.tPs,
      y: targetPlot,
      name: `Target: ${state.targetLabel}`,
      line: { width: 2.5, dash: 'dot' },
      hovertemplate: 't = %{x:.3f} ps<br>target %{y:.4g}<extra></extra>',
    });
  }

  const phasorTimePs = controls.phasorTOverT * synth.periodPs;
  const tMin = Math.min(...synth.tPs);
  const tMax = Math.max(...synth.tPs);
  let displayTime = tMin + (((phasorTimePs - tMin) % synth.periodPs + synth.periodPs) % synth.periodPs);
  traces.push({
    type: 'scatter',
    mode: 'lines',
    x: [displayTime, displayTime],
    y: controls.normalizePlots ? [-1.08, 1.08] : [0, Math.max(...synth.intensity) * 1.05],
    name: 'phasor snapshot time',
    line: { width: 1.5, dash: 'dashdot' },
    hoverinfo: 'skip',
  });

  const title = state.targetRmse === null
    ? `Resultant periodic time-domain signal, T = ${synth.periodPs.toFixed(3)} ps`
    : `Resultant periodic time-domain signal, T = ${synth.periodPs.toFixed(3)} ps | target RMSE ≈ ${state.targetRmse.toFixed(3)}`;

  const yRange = controls.normalizePlots ? [-1.15, 1.15] : undefined;
  const layout = {
    margin: { l: 70, r: 20, t: 48, b: 64 },
    title: { text: title, font: { size: 16 } },
    xaxis: { title: 'Time relative to optical carrier envelope (ps)' },
    yaxis: { title: 'Amplitude / intensity', range: yRange },
    legend: { orientation: 'h', y: 1.12, x: 0 },
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
  readLineControls();
  const controls = getGlobalControls();
  $('periodsOutput').value = controls.nPeriods;
  $('phasorOutput').value = controls.phasorTOverT.toFixed(2);

  const synth = synthesizeEnvelope(state.amplitudes, state.phasesDeg, controls.spacingGhz, controls.nPeriods);
  makeSpectrumPlot(synth, controls);
  makePhasorPlot(synth, controls);
  makeTimePlot(synth, controls);
  updateSummary(controls.spacingGhz, synth.periodPs, synth.wavelengthsNm);
}

function setValues({ magnitudes, phasesDeg, note = '', targetKind = null, targetType = null, targetLabel = '' }) {
  state.amplitudes = magnitudes.map(v => clamp(cleanNumber(v, 6), 0, 5));
  state.phasesDeg = phasesDeg.map(v => cleanNumber(v, 4));
  state.targetKind = targetKind;
  state.targetType = targetType;
  state.targetLabel = targetLabel;
  state.targetRmse = null;
  if (note) $('presetNote').innerHTML = note;
  applyStateToControls();
  updatePlots();
}

function fieldTargetPreset(kind, label) {
  const coeff = coefficientsForTarget(kind, 'field');
  const { magnitudes, phasesDeg } = coeffToMagnitudePhase(coeff);
  return {
    magnitudes,
    phasesDeg,
    note: `<strong>Field target: ${label}.</strong> Compare the dashed Re{field} curve with the dotted target. The solid intensity is |E|², so it will not itself look like the signed field target.`,
    targetKind: kind,
    targetType: 'field',
    targetLabel: `${label} field target`,
  };
}

function intensityTargetPreset(kind, label) {
  const coeff = coefficientsForTarget(kind, 'intensity');
  const { magnitudes, phasesDeg } = coeffToMagnitudePhase(coeff);
  return {
    magnitudes,
    phasesDeg,
    note: `<strong>Intensity target: ${label}.</strong> Compare the solid intensity curve with the dotted target. Finite 8-line bandwidth causes rounded edges and ripple.`,
    targetKind: kind,
    targetType: 'intensity',
    targetLabel: `${label} intensity target`,
  };
}

const presetFunctions = {
  'Transform-limited pulse train': () => ({
    magnitudes: Array(N_LINES).fill(1),
    phasesDeg: Array(N_LINES).fill(0),
    note: '<strong>Transform-limited:</strong> equal amplitudes and flat spectral phase. Produces the shortest pulse train for this flat 8-line spectrum.',
  }),
  'Time shift: linear phase ramp': () => ({
    magnitudes: Array(N_LINES).fill(1),
    phasesDeg: lineIndices.map(k => -360 * k * 0.25),
    note: '<strong>Linear phase ramp:</strong> shifts the waveform in time while preserving the intensity shape.',
  }),
  'Chirped: quadratic spectral phase': () => {
    const norm = Math.max(...centeredIndices.map(Math.abs));
    let phasesDeg = centeredIndices.map(x => 170 * (x / norm) ** 2);
    const avg = mean(phasesDeg);
    phasesDeg = phasesDeg.map(v => v - avg);
    return {
      magnitudes: Array(N_LINES).fill(1),
      phasesDeg,
      note: '<strong>Chirped field:</strong> quadratic spectral phase, analogous to group-delay dispersion.',
    };
  },
  'Chirped: strong quadratic spectral phase': () => {
    const norm = Math.max(...centeredIndices.map(Math.abs));
    let phasesDeg = centeredIndices.map(x => 380 * (x / norm) ** 2);
    const avg = mean(phasesDeg);
    phasesDeg = phasesDeg.map(v => v - avg);
    return {
      magnitudes: Array(N_LINES).fill(1),
      phasesDeg,
      note: '<strong>Strong chirped field:</strong> larger quadratic phase excursion for stronger temporal spreading.',
    };
  },
  'Chirped: cubic / asymmetric phase': () => {
    const norm = Math.max(...centeredIndices.map(Math.abs));
    return {
      magnitudes: Array(N_LINES).fill(1),
      phasesDeg: centeredIndices.map(x => 260 * (x / norm) ** 3),
      note: '<strong>Cubic / asymmetric chirp:</strong> third-order spectral phase that makes the field more asymmetric.',
    };
  },
  'Field target: sawtooth': () => fieldTargetPreset('field_saw', 'sawtooth'),
  'Field target: reverse sawtooth': () => fieldTargetPreset('field_reverse_saw', 'reverse sawtooth'),
  'Field target: triangle': () => fieldTargetPreset('field_triangle', 'triangle'),
  'Field target: square': () => fieldTargetPreset('field_square', 'square'),
  'Intensity target: sawtooth': () => intensityTargetPreset('intensity_saw', 'sawtooth'),
  'Intensity target: reverse sawtooth': () => intensityTargetPreset('intensity_reverse_saw', 'reverse sawtooth'),
  'Intensity target: triangle': () => intensityTargetPreset('intensity_triangle', 'triangle'),
  'Intensity target: square': () => intensityTargetPreset('intensity_square', 'square'),
  'Gaussian-like line amplitudes': () => {
    const { magnitudes, phasesDeg } = makeGaussianLinePreset(1.35);
    return {
      magnitudes,
      phasesDeg,
      note: '<strong>Gaussian-like spectrum:</strong> tapered line amplitudes reduce sidelobes compared with a flat spectrum.',
    };
  },
  'Double-pulse-like field': () => {
    const { magnitudes, phasesDeg } = makeDoublePulsePreset(0.35);
    return {
      magnitudes,
      phasesDeg,
      note: '<strong>Double-pulse-like field:</strong> spectral interference factor approximates two field pulses separated by about 0.35T.',
    };
  },
  'Alternating 0 / π phase': () => ({
    magnitudes: Array(N_LINES).fill(1),
    phasesDeg: lineIndices.map(k => k % 2 === 0 ? 0 : 180),
    note: '<strong>Alternating 0/π phase:</strong> flips the sign of alternating spectral lines, producing a strongly reshaped periodic field.',
  }),
  'Random amplitudes and phases': () => ({
    magnitudes: lineIndices.map(() => cleanNumber(0.45 + Math.random() * 0.75, 4)),
    phasesDeg: lineIndices.map(() => cleanNumber(-180 + Math.random() * 360, 3)),
    note: '<strong>Random amplitudes and phases:</strong> a quick way to explore arbitrary line-by-line settings.',
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

  for (let k = 0; k < N_LINES; k++) {
    const row = document.createElement('div');
    row.className = 'line-row';
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

    ampRange.addEventListener('input', () => { ampNumber.value = formatAmp(ampRange.value); });
    ampRange.addEventListener('change', () => { state.targetKind = null; scheduleUpdate(); });
    ampNumber.addEventListener('change', () => {
      const amp = clamp(parseNumber(ampNumber.value, state.amplitudes[k]), 0, 5);
      ampNumber.value = formatAmp(amp);
      ampRange.value = String(clamp(amp, 0, 1.5));
      state.targetKind = null;
      scheduleUpdate();
    });

    phaseRange.addEventListener('input', () => { phaseNumber.value = formatPhase(phaseRange.value); });
    phaseRange.addEventListener('change', () => { state.targetKind = null; scheduleUpdate(); });
    phaseNumber.addEventListener('change', () => {
      const phase = parseNumber(phaseNumber.value, state.phasesDeg[k]);
      phaseNumber.value = formatPhase(phase);
      phaseRange.value = String(wrapPhaseToSliderRange(phase));
      state.targetKind = null;
      scheduleUpdate();
    });
  }
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
  $('spacingSelect').addEventListener('change', scheduleUpdate);
  $('periodsRange').addEventListener('input', () => { $('periodsOutput').value = $('periodsRange').value; });
  $('periodsRange').addEventListener('change', scheduleUpdate);
  $('phasorRange').addEventListener('input', () => { $('phasorOutput').value = Number($('phasorRange').value).toFixed(2); });
  $('phasorRange').addEventListener('change', scheduleUpdate);
  $('normalizeCheckbox').addEventListener('change', scheduleUpdate);
  $('fieldCheckbox').addEventListener('change', scheduleUpdate);
  $('targetCheckbox').addEventListener('change', scheduleUpdate);

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
      description: '8-line optical pulse shaping settings',
      centerWavelengthNm: CENTER_WAVELENGTH_NM,
      spacingGhz: parseNumber($('spacingSelect').value, 200),
      amplitudes: state.amplitudes.map(v => cleanNumber(v, 6)),
      phasesDeg: state.phasesDeg.map(v => cleanNumber(v, 4)),
      targetKind: state.targetKind,
      targetType: state.targetType,
      targetLabel: state.targetLabel,
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
  setValues(presetFunctions['Transform-limited pulse train']());
}

window.addEventListener('DOMContentLoaded', init);
