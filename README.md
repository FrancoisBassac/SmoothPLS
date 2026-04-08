# SmoothPLS 

**SmoothPLS** is an R package designed for **Hybrid Functional Data Analysis**. It implements an advanced PLS regression framework capable of handling both **Categorical Functional Data (CFD)** and **Scalar Functional Data (SFD)**.

This package is developed as part of a PhD research project in collaboration with **Decathlon** and **INRIA**.

## Key Features
* **Unified Framework**: Integrate categorical states and continuous signals in a single PLS model.
* **S3 Implementation**: (In progress) Clean R interface for model fitting and prediction.
* **High Precision**: Segment-based integration for precise handling of state transitions.

## Installation
Currently in development. Install the stable v0.1.2 using:
```R
# install.packages("devtools")
devtools::install_github("F-BSC/SmoothPLS")
```


<!-- badges: start -->
[![R-CMD-check](https://github.com/F-BSC/SmoothPLS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/F-BSC/SmoothPLS/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
