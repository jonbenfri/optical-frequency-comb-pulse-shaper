# Interactive Optical Pulse Shaping

**Live demo:** https://jonbenfri.github.io/optical-frequency-comb-pulse-shaper/

Static GitHub Pages app for demonstrating line-by-line optical pulse shaping with up to 16 equal-frequency-spaced CW laser inputs.

## Features

- Selectable active line count: `4`, `8`, `12`, or `16`
- Laser lines centered around an editable center wavelength, default `1550.0 nm`
- Comb spacing selector: `20`, `50`, `100`, and `200 GHz`
- Editable amplitude and phase controls for every active line
- Plotly-based input spectrum, 3D phasor view, and time-domain intensity envelope
- Intensity-focused presets for sawtooth-like, reverse-sawtooth-like, triangle-like, and square-like targets
- Amplitude-shaping presets: ramped, center-weighted, edge-weighted, notched, and alternating strong/weak line amplitudes
- Compact horizontal desktop layout for seeing controls and plots at once
- JSON export for the current line settings

## Deploy with GitHub Pages

The simplest deployment is to put these files at the root of a repository and enable GitHub Pages from the repository settings.

1. Copy `index.html`, `styles.css`, `script.js`, `.nojekyll`, and optionally `.github/workflows/pages.yml` into your repository.
2. Commit and push.
3. In GitHub, open **Settings → Pages**.
4. Choose either branch deployment from the root directory, or use the included GitHub Actions workflow.

No Python backend is required. Plotly is loaded from a CDN in `index.html`.

## Notes

The app synthesizes the complex optical envelope using equal frequency spacing. The plotted waveform is the intensity envelope `|E(t)|²`, not the raw optical carrier field. With more active lines, discontinuous targets such as square and sawtooth waveforms are approximated more faithfully, though finite-bandwidth smoothing and ringing remain visible.

## Interaction notes

The amplitude, phase, time-window, and phasor-time sliders update the plots live while dragging. The JavaScript uses `requestAnimationFrame` to coalesce rapid slider events, so the app redraws at most once per browser animation frame instead of attempting a full redraw for every raw input event.

## v7 layout note

The interpretation notes now sit directly below the waveform plot inside the plot column. The preset panel and waveform panel have colored visual emphasis, and the waveform plot uses consistent colors for synthesized intensity, preset target, and phasor snapshot time.
