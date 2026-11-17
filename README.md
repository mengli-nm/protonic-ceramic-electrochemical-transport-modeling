# Protonic Ceramic Electrochemical Transport Modeling

This repository provides computational modeling tools for protonic ceramic materials. The current MATLAB implementation accompanies *Predicting Faradaic efficiency under high steam by reconstructing defect chemistry in protonic ceramic electrolysis cells* (Zhang et al., *Joule*).

## Overview

The current MATLAB routines implement a one-dimensional, steady-state electrochemical transport model for proton-conducting ceramic electrolytes. The model couples the transport of protonic defects, oxygen vacancies, electrons, and holes with defect equilibria, electrostatic potential, and electrode boundary conditions to calculate carrier distributions and Faradaic efficiency (FE).

The supplied calculation uses BaZr₀.₄Ce₀.₄Y₀.₁Yb₀.₁O₃₋δ (BZCYYb4411) as a representative electrolyte. It uses one reconstructed set of equilibrium constants and diffusivities as model input. The defect-chemistry reconstruction procedure and associated equations are described in the manuscript and Supplementary Information. This release demonstrates the electrochemical transport model using one precomputed reconstruction parameter set; it does not execute the inverse reconstruction workflow.

## Requirements

- MATLAB with Optimization Toolbox (`fsolve`)
- No external MATLAB packages are required for the included plotting scripts

## Run the transport model

Open MATLAB, set the working directory to the `code` folder, and run:

```matlab
run_example_4411
```

The script calculates steady-state transport at 600 °C and −1 A cm⁻² for oxygen-electrode steam fractions of 50%, 70%, and 90%. Material-specific model parameters are loaded from `example_parameters_4411.mat`.

The calculation saves `example_results_4411.mat`, three condition-specific `transport_*.mat` files, and six PNG plots in the `code` folder:

| Output | Quantity |
| --- | --- |
| `example_electronic_carriers.png` | Electron and hole concentration profiles |
| `example_ionic_carriers.png` | Proton and oxygen-vacancy concentration profiles |
| `example_electronic_currents.png` | Electron and hole partial-current profiles |
| `example_ionic_currents.png` | Proton and oxygen-vacancy partial-current profiles |
| `example_potential.png` | Electrostatic-potential profiles |
| `example_gas_pressures.png` | Internal gas partial-pressure profiles |

These plots illustrate the spatial distributions for the supplied parameter set. The manuscript's ensemble FE distributions incorporate multiple retained parameter sets rather than this single example.

## Plot ensemble FE results

To visualize the supplied ensemble FE distributions, run the following command from the `code` folder:

```matlab
plot_summary_4411
```

The script reads `figure_data/figure_summary_4411.mat` and writes two plots to the `generated` folder:

| Output | Quantity |
| --- | --- |
| `FE_All_4411.png` | FE distributions at 50%, 70%, and 90% H₂O across current densities |
| `RP_70_4411.png` | FE sensitivity to electrode polarization resistance at 70% H₂O |

The `figure_summary_4411.mat` file contains precomputed FE distributions for the representative electrolyte. Plotting these distributions is separate from the single-parameter forward simulation above.

## File organization

- `code/` contains the transport solver, boundary-condition and charge-balance functions, the example input parameters, and the execution and plotting scripts
- `figure_data/` contains precomputed FE distribution arrays used to illustrate the ensemble plots with `plot_summary_4411.m`

## Model parameters

`example_parameters_4411.mat` contains the MATLAB variable `P`, an 8 × 1 vector with the following ordering:

| Index | Parameter |
| --- | --- |
| 1 | Proton diffusivity, `Dp` |
| 2 | Oxygen-vacancy diffusivity, `Dv` |
| 3 | Hole diffusivity, `Dh` |
| 4 | Electron diffusivity, `De` |
| 5 | Hydration equilibrium parameter, `Kp,H2O` |
| 6 | Oxygen-incorporation equilibrium parameter, `Kp,O2` |
| 7 | Hydrogen-incorporation equilibrium parameter, `Kp,H2` |
| 8 | Electron–hole equilibrium parameter, `Keh` |

The transport solver also specifies material and cell properties and electrode boundary conditions in its source code. Application to other electrolyte compositions requires their corresponding reconstructed parameters and material-specific settings. See the manuscript and Supplementary Information for definitions and modeling assumptions.