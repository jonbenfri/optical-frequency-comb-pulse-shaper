# Changelog

## v9

- Added optimized finite-line presets for triangle and square intensity targets.
- Removed the redundant `Alternating 0 / π phase` preset because it is equivalent to a half-period linear phase ramp.
- Kept version-specific notes out of the README.
- Added this standalone changelog.
- Bumped static asset cache strings to `v=10`.

## v8

- Improved sawtooth and reverse-sawtooth intensity presets using finite-line spectral-factor coefficients.
- Kept the targets discontinuous; no smoothing was added.

## v7

- Moved interpretation notes below the waveform plot.
- Added visual emphasis to the presets and waveform panels.
- Added consistent Plotly trace colors for synthesized intensity, target, and phasor snapshot time.

## v6 and earlier

- Added selectable line counts up to 16.
- Added editable center wavelength.
- Added compact horizontal layout for GitHub Pages.
- Added live slider updates using `requestAnimationFrame`.
- Added amplitude-shaping presets and removed phase-chirp presets.
