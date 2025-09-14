# vdisp_fit: A Stellar Velocity Dispersion Fitting Pipeline
---
## Abstract

Stellar velocity dispersion (σ) is a fundamental physical parameter that quantifies the random motions of stars within a galaxy and provides critical insights into its dynamical mass, gravitational potential, and evolutionary history. Accurate measurement of σ is essential for understanding galaxy formation, the scaling relations between galaxies and their central black holes, and the mass-to-light ratio of stellar populations.  

This project, **vdisp_fit**, presents a robust and automated Python pipeline for measuring stellar velocity dispersions from galaxy spectra, particularly from large spectroscopic surveys such as the Sloan Digital Sky Survey (SDSS). The pipeline employs direct spectral fitting by convolving high-resolution stellar templates from the MILES library with a Gaussian broadening function, combined with an additive polynomial to account for continuum variations. A stable χ² minimization routine is used to separate non-linear broadening from linear components, ensuring accurate and efficient fitting.  

**vdisp_fit** offers both a command-line interface and a library function, enabling flexible batch processing, robust error estimation, instrumental correction, and visualization of results. By providing an end-to-end automated solution, this pipeline facilitates the study of galaxy kinematics and contributes to the broader understanding of galaxy structure and evolution.

---

**vdisp_fit** is a robust and automated Python pipeline designed to measure the stellar velocity dispersion (σ) of galaxies from spectroscopic surveys, such as the Sloan Digital Sky Survey (SDSS).  
It employs a direct spectral fitting technique by convolving a high-resolution stellar template from the MILES library with a Gaussian broadening function to model the observed galaxy spectrum.

The core methodology is built on a stable χ² minimization routine that separates the non-linear broadening search from the linear components of the fit, ensuring a robust and efficient solution.

---

## Key Features

- **Robust Fitting**: Utilizes a stable additive model with an analytical linear least-squares solution, avoiding issues with local minima.  
- **Flexible Input**: Can process a single galaxy FITS file or a batch of files from an entire directory.  
- **Instrumental Correction**: Corrects the fitted velocity dispersion for instrumental resolution differences between the galaxy and template spectra.  
- **Error Estimation**: Provides a robust 1-σ uncertainty on the final velocity dispersion measurement.  
- **Command-Line Tool**: Can be run directly from the terminal with a simple command structure.  
- **Plotting Capabilities**: Generates detailed single-fit plots and summary plots for batch processing to visualize results.  
- **Single-File Module**: Entire pipeline contained in a single Python script (`vdisp_fit.py`) for easy distribution and use.

---

## Methodology

The pipeline's methodology follows a four-step process:

### 1. Data Pre-processing
- **Redshift Correction**: Galaxy spectra are shifted to their rest-frame using the SDSS-provided redshift.  
- **Logarithmic Rebinning**: Spectra are rebinned to a logarithmic wavelength scale, ensuring that kinematic broadening corresponds to a constant pixel shift.  
- **Masking**: Strong emission lines and bad pixels are masked out to prevent contamination of the stellar continuum fit.

### 2. The Fitting Model

The observed galaxy spectrum *G_obs* is modeled as a broadened stellar template *T* combined with an additive polynomial *Pₖ* to account for continuum variations:

$$
G_{model}(\ln \lambda) = c_0 \cdot [T(\ln \lambda) * B(v, \sigma)] + \sum_{k=1}^{N} c_k P_k(\ln \lambda)
$$

- *B* is a Gaussian broadening function.  
- *c₀,...,c_N* are linear coefficients.

### 3. χ² Minimization

- The pipeline iterates through a grid of trial σ values.  
- For each σ, the optimal linear coefficients are solved analytically, and the corresponding χ² value is calculated.  
- The σ that yields the minimum χ² is identified as the best-fit velocity dispersion.

### 4. Instrumental Correction and Error Estimation

- The measured dispersion is corrected for instrumental broadening:

  $$
  \sigma_{gal} = \sqrt{\sigma_{fit}^2 - \sigma_{inst}^2}
  $$

- The 1-σ error is estimated from the curvature of the χ² curve near its minimum.

---

## Requirements

Install the required Python libraries using pip:

```bash
pip install numpy matplotlib scipy astropy
````

---

## Usage

### As a Command-Line Tool

Run the pipeline from the terminal:

```bash
python vdisp_fit.py <input_path> <template_folder> [options]
```

#### Arguments

* `<input_path>`: Path to a single galaxy FITS file or directory containing multiple FITS files.
* `<template_folder>`: Path to the directory containing the MILES stellar template FITS files.

#### Options

* `-o`, `--output`: Name of the CSV file to save results (e.g., `results.csv`).
* `-p`, `--plot`: Display plots of the fits.
* `--poly_degree`: Degree of the additive polynomial (default: 8).
* `--velscale`: Velocity scale in km/s per pixel (default: 60.0).
* `--fit_range`: Wavelength range for fitting in Å (e.g., `4200 6500`).
* `--instrumental_sigma`: Instrumental dispersion difference in km/s (default: 25.0).

#### Example

```bash
python vdisp_fit.py ./galaxy_spectra/ ./miles_templates/ -o results.csv --plot
```

---

### As a Library Function

Import and use `vdisp_fit` in another Python script:

```python
from vdisp_fit import vdisp_fit

results = vdisp_fit(
    input_path='./galaxy_spectra/galaxy_1.fits',
    template_folder='./miles_templates/',
    output_csv='single_fit_result.csv',
    plot=True
)

print(results)
```

---

## Acknowledgements

This pipeline was developed as a project report at IUCAA.
It builds upon foundational work in galaxy kinematics, including the pioneering work of Faber & Jackson (1976) and modern techniques like the **pPXF method** developed by Cappellari & Emsellem (2004).
It utilizes publicly available data from the Sloan Digital Sky Survey (SDSS) and the MILES stellar library.

---

**Author**: Jashanpreet Singh Dingra
**Supervisor**: Dr. Gulab Chand Dewangan
**Date**: September 1, 2025
