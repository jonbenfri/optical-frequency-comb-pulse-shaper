# Interactive Optical Pulse Shaping

Static GitHub Pages app for demonstrating line-by-line optical pulse shaping with 8 equal-frequency-spaced CW laser inputs.

## Features

- 8 laser lines centered around an editable center wavelength, default `1550.0 nm`
- Comb spacing selector: `20`, `50`, `100`, and `200 GHz`
- Editable amplitude and phase controls for every line
- Plotly-based input spectrum, 3D phasor view, and time-domain intensity envelope
- Intensity-focused presets for sawtooth-like, reverse-sawtooth-like, triangle-like, and square-like targets
- Chirped spectral-phase presets
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

The app synthesizes the complex optical envelope using equal frequency spacing. The plotted waveform is the intensity envelope `|E(t)|²`, not the raw optical carrier field. With only eight lines, discontinuous targets such as square and sawtooth waveforms will show finite-bandwidth smoothing and ringing.
