# SWAT-SVR-Roubidoux

## Overview

This repository contains R scripts implementing hybrid SWAT+–SVR streamflow prediction models developed for the Roubidoux Creek Watershed (RCW), Missouri, USA. The modeling framework combines physically based SWAT+ simulations with Support Vector Regression (SVR) to improve streamflow prediction performance under both monthly and daily time scales.

Two hybrid modeling approaches are provided:

* **SVR_mon_final1.R** – Monthly SWAT+–SVR hybrid model
* **SVR_daily_final1.R** – Daily two-stage SWAT+–SVR hybrid model

The workflow was developed as part of a hydrologic modeling study evaluating machine-learning-assisted improvements to SWAT+ simulations in a karst-dominated watershed.

---

## Model Framework

### Monthly SWAT+–SVR Model

The monthly model uses SWAT+-simulated streamflow together with hydrologically relevant predictor variables to estimate monthly observed streamflow.

Key predictor variables include:

* SWAT+-simulated monthly streamflow
* Estimated monthly baseflow proxy
* Monthly precipitation
* Antecedent Precipitation Index (API)
* Rolling three-month precipitation accumulation (RR3)
* Lagged hydrologic variables
* Drainage area
* Seasonal sinusoidal predictors

---

### Daily Two-Stage SWAT+–SVR Model

The daily model employs a two-stage framework:

#### Stage 1: Baseflow Proxy Estimation

Observed daily streamflow is separated into baseflow and quickflow components using the Eckhardt digital filter. An intermediate SVR model is then trained to estimate a daily baseflow proxy from:

* SWAT+-simulated streamflow
* Daily precipitation
* Antecedent precipitation
* Antecedent Precipitation Index (API)
* Drainage area
* Seasonal predictors

#### Stage 2: Streamflow Prediction

The predicted baseflow proxy is incorporated into a second SVR model used to estimate daily streamflow.

Predictors include:

* Current SWAT+ streamflow
* Lagged SWAT+ streamflow
* Predicted baseflow proxy
* Lagged predicted baseflow proxy
* Daily precipitation
* Antecedent precipitation
* Antecedent Precipitation Index (API)
* Drainage area
* Seasonal sinusoidal predictors

This design prevents information leakage by ensuring that observed baseflow from the validation period is never used during prediction.

---

## Study Area

Roubidoux Creek Watershed (RCW), Missouri, USA

Monitoring stations included:

| USGS Site No. | Station Description                                 |
| ------------- | --------------------------------------------------- |
| 06928420      | Roubidoux Creek at Polla Rd below Fort Leonard Wood |
| 06928380      | Upper Smith Branch below Engineering Ponds          |
| 06928330      | Hurd Hollow Tributary (Outfall 14)                  |
| 06928320      | Musgrave Hollow Tributary                           |
| 06928300      | Roubidoux Creek above Fort Leonard Wood             |

---

## Required Inputs

### SWAT+ Outputs

Daily or monthly SWAT+ channel output files:

* `channel_sd_day.txt`
* `channel_sd_mon.txt`

### Observed Streamflow

Observed streamflow records with:

* Site ID
* Date
* Streamflow (m³ s⁻¹)
* Calibration/validation period designation

### Precipitation Data

Daily precipitation time series mapped to each monitoring station.

---

## Required R Packages

```r
dplyr
tidyr
lubridate
ggplot2
e1071
hydroGOF
FlowScreen
patchwork
scales
tibble
```

Install missing packages using:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "lubridate",
  "ggplot2",
  "e1071",
  "hydroGOF",
  "FlowScreen",
  "patchwork",
  "scales",
  "tibble"
))
```

---

## Model Outputs

The scripts generate:

* Streamflow predictions
* Baseflow proxy estimates
* Calibration and validation performance statistics
* Site-specific performance metrics
* Plotting datasets
* Publication-quality hydrographs

Performance metrics include:

* Nash–Sutcliffe Efficiency (NSE)
* Kling–Gupta Efficiency (KGE)
* Log-transformed KGE (KGElf)
* Percent Bias (PBIAS)
* Root Mean Square Error (RMSE)
* Coefficient of Determination (R²)

---

## Reproducibility

Random seeds are fixed to ensure reproducibility of model results.

The scripts were developed and tested using R (version 4.x) on Windows operating systems.

---

## Citation

If you use this repository, please cite the associated manuscript:

Lifeng Yuan, ...

---
