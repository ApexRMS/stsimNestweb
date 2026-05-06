# stsimNestweb

### Nesting bird habitat suitability modeling for ST-Sim

**stsimNestweb** is an add-on package for [SyncroSim](https://syncrosim.com/) that extends [ST-Sim](https://syncrosim.com/packages/) to model habitat suitability for aspen-associated nesting bird species. It applies pre-trained generalized linear mixed models (GLMMs) to landscape simulation outputs, producing spatially explicit habitat suitability maps and tabular summaries across stochastic iterations and timesteps.

[![GitHub release](https://img.shields.io/github/v/release/ApexRMS/stsimNestweb)](https://github.com/ApexRMS/stsimNestweb/releases)

---

## Requirements

- [SyncroSim](https://syncrosim.com/download/) >= 3.1.0
- [ST-Sim](https://syncrosim.com/packages/) package

---

## Installation

1. Clone or download this repository to your local machine.
2. Open SyncroSim Studio and go to **File > Local Packages > Install from Folder...** and select the `src` folder of this repository.
3. SyncroSim will prompt you to set up the required conda environment (`stsimnestweb_env`) during installation.

---

## Overview

stsimNestweb integrates with an existing ST-Sim scenario to predict per-pixel habitat suitability for one or more nesting bird species at each simulation timestep. The package uses species-specific GLMMs (stored as external `.RData` files) fit outside of SyncroSim, and propagates model parameter uncertainty by sampling fixed-effect coefficients each iteration.

### How it works

For each iteration and timestep, the package:

1. Samples GLMM coefficients from their uncertainty distributions to propagate model parameter uncertainty across stochastic runs.
2. Loads landscape rasters produced by ST-Sim:
   - Aspen cover (%) and diameter (cm) from stock & flow spatial outputs
   - State class and Time Since Transition (TST) from ST-Sim spatial outputs
3. Derives spatial covariates from those rasters:
   - Distance to most recent harvest (from TST)
   - Distance to forest edge (from state class reclassification)
4. Predicts habitat suitability (0–1) for each species using the GLMM.
5. Masks out habitat in user-defined invalid stratum–state class combinations.
6. Writes spatial rasters and/or tabular habitat amount summaries.

Optionally, a second transformer averages habitat suitability rasters across all iterations to produce mean maps.

---

## Inputs

All inputs are configured at the scenario level in SyncroSim Studio under the **ST-Sim Nestweb** section.

| Datasheet | Description |
|---|---|
| **Habitat Models** | Path to a pre-trained GLMM (`.RData` file) for each species |
| **Sites** | A raster classifying each pixel by site type |
| **Invalid Habitat** | Stratum–state class–species combinations to mask out (set to 0) |
| **Output Options** | Flags and timestep intervals controlling which outputs are generated |

**Project-level inputs** (shared across scenarios):

| Datasheet | Description |
|---|---|
| **Species** | Lookup table of species names |
| **Sites** | Lookup table of site type names and IDs |

---

## Outputs

| Output | Type | Description |
|---|---|---|
| **Habitat Amount** | Tabular | Habitat suitability summed to hectares, by stratum, secondary stratum, site, and species |
| **Habitat** | Spatial raster | Per-species habitat suitability (0–1) at each output timestep |
| **Habitat Change** | Spatial raster | Change in suitability relative to the first output timestep |
| **Habitat Average** | Spatial raster | Mean habitat suitability across all iterations |
| **Habitat Change Average** | Spatial raster | Mean habitat change across all iterations |

Tabular and spatial outputs can be configured to be written at different timestep intervals via **Output Options**.

---

## Pipeline

stsimNestweb adds two transformers to the ST-Sim pipeline:

| Transformer | Description |
|---|---|
| **Model Habitat** | Core habitat modeling — runs for each iteration and timestep |
| **Summarize Habitat** | Post-processing — averages spatial outputs across iterations (local only) |

---

## Habitat models

Each species requires a pre-trained GLMM saved as an `.RData` file containing a `glmmTMB` model object. The model must include the following fixed-effect covariates:

| Covariate | Description |
|---|---|
| `Perc_At` | Aspen cover (proportion) |
| `Median_DBH` | Aspen diameter (cm) |
| `edge_near` | Distance to forest edge (m) |
| `dist_to_cut` | Distance to most recent harvest (m) |
| `Num_Trees` | Mean number of trees |
| `Num_2BI` | Mean number of two-by-inch stems |
| `Mean_decay` | Mean snag decay class |
| `Site` | Site type (factor) |

---

## Contact

Developed by [ApexRMS](https://apexrms.com/). For questions or issues, please open a [GitHub issue](https://github.com/ApexRMS/stsimNestweb/issues).
