# Code for *Interface-mediated secondary phase separation*

This repository contains the MATLAB code and numerical data associated with the manuscript:

**Interface-mediated secondary phase separation**

The code is provided to reproduce the numerical results and figures presented in the main text and Supplementary Information.

## Overview

The numerical calculations are based on a ternary phase-field model describing a scaffold, a client, and a solvent. The repository includes codes for equilibrium calculations, Cahn–Hilliard dynamics, phase-boundary analysis, and generation of the figures presented in the manuscript.

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

## Requirements

The code was developed and tested in MATLAB.

Some scripts may require MATLAB toolboxes, including:

- Image Processing Toolbox

Numerical BVP calculations use MATLAB's `bvp4c` solver.

## Data and Code Availability

All MATLAB scripts and numerical data required to reproduce the computational results presented in the manuscript are included in this repository.