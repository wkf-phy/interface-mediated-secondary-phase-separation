# Figure S2

This folder contains the MATLAB scripts and data used to reproduce Figure S2.

## Folder Structure

```text
FigS2/
├── plot_FigS2.m
├── critical_profile.mat
│
└── critical_profile_scan/
    ├── Scan_critical_profile.m
    ├── critical_profile.mat
    └── ccode/
        ├── get_steady_droplet.m
        ├── get_steady_ternary.m
        ├── plot_cartesian_phase_diagram.m
        └── prepare_fft.m
```

## Fig. S2

`plot_FigS2.m` reproduces the size dependence of the critical interfacial density shown in Fig. S2.

### Files

- `plot_FigS2.m` — Main plotting script for Fig. S2.
- `critical_profile.mat` — Precomputed critical concentration profiles used to generate Fig. S2.

### Run

Open the `FigS2` folder in MATLAB and run:

```matlab
plot_FigS2
```

## Critical Profile Scan

The `critical_profile_scan` folder contains the code used to calculate the critical concentration profiles stored in `critical_profile.mat`.

`Scan_critical_profile.m` generates a pre-relaxed two-dimensional concentration field for each scaffold concentration and uses radial BVP continuation to determine the critical profile.

### Files

- `Scan_critical_profile.m` — Main script for scanning and extracting the critical concentration profiles.
- `critical_profile.mat` — Output of the critical-profile scan.
- `ccode/get_steady_ternary.m` — Generates the pre-relaxed ternary concentration field used to initialize the radial BVP.
- `ccode/get_steady_droplet.m` — Generates the pre-existing scaffold droplet used during initialization.
- `ccode/prepare_fft.m` — Constructs the Fourier-space wave numbers used by the spectral solver.
- `ccode/plot_cartesian_phase_diagram.m` — Helper function for plotting the ternary phase diagram.

### Run

Open the `critical_profile_scan` folder in MATLAB and run:

```matlab
Scan_critical_profile
```

The resulting `critical_profile.mat` contains the critical concentration profiles used by `plot_FigS2.m`.

## Requirements

The scripts were developed and tested in MATLAB. Some functions may require the Image Processing Toolbox.