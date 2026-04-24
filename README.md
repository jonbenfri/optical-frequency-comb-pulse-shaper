# Interactive Optical Pulse Shaping

A static GitHub Pages web app for demonstrating line-by-line optical pulse shaping with eight equal-frequency-spaced CW laser lines near 1550 nm.

The app lets you adjust each line's amplitude and phase, apply presets, and view:

- the input comb spectrum,
- a 3D phasor representation of the line phases,
- the synthesized periodic time-domain field/intensity envelope.

## Files

```text
index.html
styles.css
script.js
.nojekyll
.github/workflows/pages.yml   # optional GitHub Actions deployment workflow
```

The site uses Plotly from a CDN, so the browser needs internet access when loading the page.

## Run locally

Because this is a static site, you can open `index.html` directly in a browser. A local server is also fine:

```bash
python -m http.server 8000
```

Then open `http://localhost:8000`.

## Deploy with GitHub Pages, simple branch method

1. Create a new repository, for example `pulse-shaping-demo`.
2. Copy these files into the repository root.
3. Commit and push to `main`.
4. In GitHub, go to **Settings → Pages**.
5. Under **Build and deployment**, choose **Deploy from a branch**.
6. Select `main` and `/(root)`, then save.

The site will appear at:

```text
https://<your-username>.github.io/<repo-name>/
```

## Deploy with GitHub Actions instead

This package also includes `.github/workflows/pages.yml`. To use it:

1. In **Settings → Pages**, set **Source** to **GitHub Actions**.
2. Push to `main`.
3. The workflow uploads the repository root as a static Pages artifact and deploys it.

The branch method is simpler for this no-build site. The Actions method is useful if you later add a build step.

## Notes on the presets

- Field presets target `Re{E(t)}` in a line-index/Fourier-series reference frame.
- Intensity presets target the measurable `|E(t)|²`.
- With only eight lines, discontinuous targets such as square and sawtooth waves show finite-bandwidth ringing and rounded edges.
