
# vdisp_fit: A Stellar Velocity Dispersion Fitting Pipeline
---

## Abstract

```diff
+ Stellar velocity dispersion (σ) is a fundamental parameter quantifying random stellar motions in galaxies.
+ Accurate σ measurement is crucial for understanding galaxy dynamics, mass-to-light ratios, and black hole scaling relations.
````

This project, **vdisp\_fit**, presents a robust Python pipeline to measure stellar velocity dispersions from galaxy spectra, particularly from large spectroscopic surveys like **SDSS**.

The pipeline employs **direct spectral fitting** by convolving high-resolution **MILES stellar templates** with a **Gaussian broadening function**, combined with an additive polynomial to account for continuum variations. A stable **χ² minimization** routine separates non-linear broadening from linear components, ensuring accurate and efficient fitting.

**vdisp\_fit** supports both **command-line interface** and **library usage**, batch processing, error estimation, instrumental correction, and visualization.

---

## Key Features

```diff
+ Robust Fitting             : Stable additive model avoids local minima
+ Flexible Input             : Single FITS file or batch directory
+ Instrumental Correction    : Corrects for resolution differences
+ Error Estimation           : Provides 1-σ uncertainty
+ Command-Line Tool          : Simple terminal interface
+ Plotting                   : Generates single-fit and summary plots
+ Single-File Module         : Entire pipeline in `vdisp_fit.py`
```

---

## Methodology

### 1. Data Pre-processing

```python
# Shift galaxy spectra to rest-frame using SDSS redshift
# Rebin spectra to logarithmic wavelength scale for constant velocity pixels
# Mask strong emission lines and bad pixels
```

### 2. Fitting Model

```python
# Observed galaxy spectrum G_obs is modeled as:
G_model = c0 * convolve(T, B(v, sigma)) + sum(c_k * P_k)
# B is Gaussian broadening
# c0,...,cN are linear coefficients
# P_k is additive polynomial for continuum variations
```

### 3. χ² Minimization

```diff
+ Iterate over trial σ values
+ Solve linear coefficients analytically for each σ
+ Compute χ² for each trial
+ Minimum χ² → best-fit velocity dispersion
```

### 4. Instrumental Correction & Error Estimation

```python
sigma_gal = sqrt(sigma_fit**2 - sigma_inst**2)
# 1-σ error estimated from curvature of χ² near minimum
```

---

## Requirements

```bash
pip install numpy matplotlib scipy astropy
```

---

## Usage

### As a Command-Line Tool

```bash
python vdisp_fit.py <input_path> <template_folder> [options]
```

#### Arguments

```diff
+ <input_path>       : Path to a single galaxy FITS file or directory of FITS files
+ <template_folder>  : Directory containing MILES stellar templates
```

#### Options

```diff
+ -o, --output           : CSV file name for results (e.g., results.csv)
+ -p, --plot             : Display plots of fits
+ --poly_degree          : Degree of additive polynomial (default: 8)
+ --velscale             : Velocity scale in km/s per pixel (default: 60.0)
+ --fit_range            : Wavelength range for fitting in Å (e.g., 4200 6500)
+ --instrumental_sigma   : Instrumental dispersion difference in km/s (default: 25.0)
```

#### Example

```bash
python vdisp_fit.py ./galaxy_spectra/ ./miles_templates/ -o results.csv --plot
```

---

### As a Library Function

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

This pipeline was developed as a project report at **IUCAA**.
It builds upon foundational work in galaxy kinematics, including:

```diff
+ Faber & Jackson (1976)
+ pPXF method by Cappellari & Emsellem (2004)
```

It uses publicly available data from **SDSS** and the **MILES stellar library**.

---

**Author**: Jashanpreet Singh Dingra
**Supervisor**: Dr. Gulab Chand Dewangan
**Date**: September 1, 2025
