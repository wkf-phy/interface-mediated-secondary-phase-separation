# Code for *Interface-mediated secondary phase separation*

This repository contains the MATLAB code and numerical data associated with the manuscript:

**Interface-mediated secondary phase separation**

The code is provided to reproduce the numerical results and figures presented in the main text and Supplementary Information.

## Overview

The numerical calculations are based on a ternary phase-field model describing a scaffold, a client, and a solvent. The repository includes codes for Cahn–Hilliard dynamics, phase-boundary analysis, and generation of the figures presented in the manuscript.

The code is organized primarily by figure. Each figure folder contains the scripts and data required to reproduce the corresponding figure. Additional subfolders are provided when intermediate numerical calculations are required to generate the plotting data.

## Repository Structure

The repository contains folders corresponding to figures in the main text and Supplementary Information.

In general, each figure folder contains:

- MATLAB scripts (`.m`) used to reproduce the figure.
- Precomputed numerical data (`.mat`) used by the plotting scripts.
- Supporting MATLAB functions when required.
- Additional numerical scripts for regenerating the precomputed data when applicable.

## Numerical Methods

The repository includes implementations of:

- Cahn–Hilliard simulations for the dynamics of the ternary mixture.
- Radially symmetric boundary-value problem (BVP) calculations for equilibrium concentration profiles.
- Continuation calculations for determining instability boundaries and critical profiles.
- Analysis of interfacial enrichment and droplet morphology.

## System Requirements

### Software requirements

The code was developed and tested using:

- MATLAB R2025b
- Windows 11

Additional toolbox requirements depend on the scripts being run:

- Image Processing Toolbox: required for image-based morphology analysis
  using functions such as `bwlabel`.
- Parallel Computing Toolbox: required for parallel execution of the
  parameter-scan scripts that use `parfor`.

Boundary-value problem calculations use MATLAB's built-in `bvp4c` solver.

Other MATLAB versions and operating systems have not been explicitly tested.

### Hardware requirements

The code runs on a standard desktop computer using a CPU.
No GPU or other non-standard hardware is required.

## Installation Guide

1. Install MATLAB and the toolboxes required for the selected scripts.
2. Download this repository using **Code > Download ZIP** on GitHub,
   and extract the archive.
3. Open MATLAB and set the current folder to the directory containing
   the script you wish to run.

The repository contains MATLAB source files and does not require
compilation or a separate software installation step.

Once MATLAB and the required toolboxes are installed and the repository
has been downloaded, setup consists only of extracting the files and
selecting the appropriate working folder. The time required depends
on download speed and file extraction performance.

Keep the supplied data files and supporting function folders in their
original locations.

## Demo

Two examples are provided for Figure 1: plotting precomputed results
and rerunning the numerical simulation.

### 1. Plot precomputed results

In MATLAB, change the current folder to `Fig1` and run:

    plot_Fig1

**Input:** `Fig1_dynamics_snapshots.mat`, included in the same folder.

**Expected output:** A figure showing concentration maps and corresponding
concentration profiles at four selected times.

**Expected runtime:** Approximately 10–20 seconds on the authors'
standard desktop computer.

This example visualizes precomputed simulation data; it does not rerun
the dynamics.

### 2. Run the dynamics simulation

In MATLAB, change the current folder to `Fig1/dynamic` and run:

    run_dynamics

The script uses the supplied initial-condition file
`phi1_phi3_equilibrium.mat` and supporting functions in the `ccode`
subfolder.

**Expected behavior:** The script evolves the ternary concentration
fields from the supplied initial condition using the parameters
specified in the script.

**Expected runtime:** Approximately 2–3 minutes on the authors'
standard desktop computer.

Actual runtimes may vary with hardware and system load.

## Instructions for Use

The code is organized by figure. To run a particular calculation or
plotting script:

1. Set the MATLAB current folder to the directory containing the script.
2. Ensure that the required input data and supporting functions are
   available in their original locations.
3. Run the script from the MATLAB Editor or Command Window.

Plotting scripts use the supplied precomputed `.mat` files.
Simulation and numerical-analysis scripts are provided in the
corresponding subfolders where applicable.

To investigate other parameter settings, modify the relevant parameters
in the simulation script, such as interaction strengths, average
compositions, domain size, grid resolution, or simulation duration.
Initial conditions should remain consistent with the selected model
parameters and discretization.

Changing parameters may alter the resulting dynamics and computational
cost. The supplied parameter settings should be retained when reproducing
the corresponding results of the manuscript.

## Data and Code Availability

The MATLAB source code and numerical data accompanying
*Interface-mediated secondary phase separation* are available at:

https://github.com/wkf-phy/interface-mediated-secondary-phase-separation

The repository includes precomputed data for plotting and numerical
scripts for the calculations provided in the figure-specific folders.
