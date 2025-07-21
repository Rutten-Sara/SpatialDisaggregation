# A Bayesian Geoadditive Model for Spatial Disaggregation
Sara Rutten, Thomas Neyens, Elisa Duarte, Christel Faes

## About this repository
This repository contains Rcodes used to generate the results from the simulation study and PBC data application from the paper "A Bayesian Geoadditive Model for Spatial Disaggregation".

## Data
The dataset on mortality in Belgium can be requested through Statbel. Following datasets are provided in the **data** folder:

| Dataset | Description | Downloaded from |
| --- | --- | --- |
| PBCshp.RData | Dataset containing incidence on PBC in Newcastle upon Tyne, UK | Available in SDALGCP R package or downloaded from https://github.com/olatunjijohnson/SDALGCP|
| pop_den.RData | Dataset containing population density in Newcastle upon Tyne, UK at grid level | Available in SDALGCP R package or downloaded from https://github.com/olatunjijohnson/SDALGCP |

## Rcode

The main Rcode used in this paper is structured in two files:

| Function | Description |
| --- | ---|
| simulation.R | Simulation study code |
| data_application_PBC.R | Main code for data application on PBC |

The code on the Belgium data application is not provided, as the dataset is not publicly available. The code of this data application is available on request.

The folder **functions** contains all functions necessary to fit our proposed method. These functions are called by the simulation.R and data_application_PBC.R main files. This folder contains the following functions:

| Function | Description |
| --- | ---|
| Krig_fit.R | Main function to fit the proposed model |
| Krig_approx.R | Function to fit the proposed model when using a Taylor approximation (SSDAM). This function is called by Krig_fit.R |
| Krig_exact.R | Function to fit the proposed model when using the exact method (SSDEM). This function is called by Krig_fit.R |
| cov.spatial.R | Function defining the spatial spline basis functions. This function is called by Krig_approx.R and Krig_exact.R |
| help_functions.R | Additional functions for data preparation and prediction |

The **functions** folder contains an additional subfolder, called **Other functions**, which are the functions to fit the SDALGCP method (downloaded directly from https://github.com/olatunjijohnson/SDALGCP), ATP Poisson Kriging (downloaded directly from https://github.com/DavidPayares/ATA-Poisson-Cokriging/tree/main) and additional functions to process the results from the disaggregation package.
