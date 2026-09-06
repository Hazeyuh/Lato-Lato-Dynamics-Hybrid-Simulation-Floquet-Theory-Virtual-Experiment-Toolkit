# Lato-lato Simulation Code

This directory contains the MATLAB code used for the hybrid flexible-cord model, the small-angle Floquet calculation, finite-time startup scans, and the interactive virtual-experiment interface.

## Layout

- `+lato/`: simulation, drive, theory, scan, and plotting functions.
- `configs/`: parameter sets for the full model and the collision-suppressed startup reference.
- `gui/`: interactive interface for fixed-frequency, continuous-sweep, and step-sweep simulations.
- `scripts/`: three reproducible entry points.
- `setup_lato.m`: adds the package, configuration, GUI, and script directories to the MATLAB path.

## Requirements

Use MATLAB R2021a or later. The core calculations use standard MATLAB functions, including `ode45`.

## Run

Start MATLAB in this directory and run:

```matlab
setup_lato
make_F1_mechanism_timeline
make_F2_startup_phase_diagram
```

`make_F7_startup_cross_maps` performs a substantially larger scan. Set `formalMode = false` in that script for a quick check before a full run.

To open the graphical interface, run:

```matlab
lato_gui
```

The GUI uses the same `lato.simulate_hybrid` solver as the scripts. It computes a run before playback; pause and stop affect playback only, and changed controls apply to the next run.
The supplementary package also includes a reference screenshot at [`../F9_gui_virtual_experiment.png`](../F9_gui_virtual_experiment.png).

Scripts create `results/` and `figures/` as needed. All paths are resolved relative to this directory.
