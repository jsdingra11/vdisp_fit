Ah! That explains it — GitHub Markdown **does not fully support HTML for styling**. Inline CSS like `style="color:#..."` **won’t render** on GitHub README.md. So we need to rely on **native Markdown features** and clever use of syntax highlighting/code blocks to add “color-like” emphasis.

We can make it visually appealing **without emojis** by using:

* `diff` blocks for green/red highlights
* `yaml` or `ini` for subtle coloring
* Bold/italic for emphasis
* Table formatting for organization
* Code blocks for code snippets

Here’s a **full working colorful version** compatible with GitHub:

````markdown
# vdisp_fit : A Stellar Velocity Dispersion Fitting Pipeline
---

## Abstract

Stellar velocity dispersion (**σ**) is a fundamental physical parameter that quantifies the random motions of stars within a galaxy and provides critical insights into its dynamical mass, gravitational potential, and evolutionary history. Accurate measurement of **σ** is essential for understanding galaxy formation, scaling relations between galaxies and their central black holes, and the mass-to-light ratio of stellar populations.

This project, **vdisp_fit**, presents a robust and automated Python pipeline for measuring stellar velocity dispersions from galaxy spectra, particularly from large spectroscopic surveys such as the Sloan Digital Sky Survey (SDSS). The pipeline employs direct spectral fitting by convolving high-resolution stellar templates from the MILES library with a Gaussian broadening function, combined with an additive polynomial to account for continuum variations. A stable **χ²** minimization routine is used to separate non-linear broadening from linear components, ensuring accurate and efficient fitting.

**vdisp_fit** offers both a command-line interface and a library function, enabling flexible batch processing, robust error estimation, instrumental correction, and visualization of results.

---

## Key Features

```diff
+ Robust Fitting       : Stable additive model with analytical linear least-squares
+ Flexible Input       : Single FITS file or batch directory processing
+ Instrumental Correction : Corrects for galaxy-template resolution differences
+ Error Estimation     : Provides 1-σ uncertainty on the final velocity dispersion
+ Command-Line Tool    : Run directly from terminal
+ Plotting Capabilities: Single-fit and batch summary plots
+ Single-File Module   : Entire pipeline in one Python script (vdisp_fit.py)
````

---

## Methodology

### 1. Data Pre-processing

```yaml
Redshift Correction   : Shift spectra to rest-frame using SDSS redshift
Logarithmic Rebinning : Rebin to log wavelength scale for constant pixel broadening
Masking               : Mask strong emission lines and bad pixels
```

### 2. The Fitting Model

The observed galaxy spectrum G\_obs is modeled as:

```math
G_model(ln λ) = c0 * [T(ln λ) * B(v, σ)] + Σ ck Pk(ln λ)
```

* B = Gaussian broadening function
* c0,...,cN = linear coefficients
* Pk = additive polynomial

### 3. χ² Minimization

```diff
+ Iterate through trial σ values
+ Solve optimal linear coefficients analytically
+ Compute χ² for each σ
+ Best-fit σ = σ with minimum χ²
```

### 4. Instrumental Correction & Error Estimation

```diff
+ Corrected dispersion : σ_gal = sqrt(σ_fit^2 - σ_inst^2)
+ 1-σ error from χ² curvature near minimum
```

---

## Requirements

```bash
pip install numpy matplotlib scipy astropy
```

---

## Usage

### Command-Line Tool

```bash
python vdisp_fit.py <input_path> <template_folder> [options]
```

**Arguments**

```yaml
<input_path>       : Single FITS file or directory
<template_folder>  : Directory containing MILES stellar templates
```

**Options**

```yaml
-o, --output            : CSV file to save results
-p, --plot              : Display plots
--poly_degree           : Degree of additive polynomial (default=8)
--velscale              : km/s per pixel (default=60.0)
--fit_range             : Wavelength range in Å (e.g., 4200 6500)
--instrumental_sigma    : Instrumental dispersion (default=25.0)
```

**Example**

```bash
python vdisp_fit.py ./galaxy_spectra/ ./miles_templates/ -o results.csv --plot
```

### Library Function

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

Developed at IUCAA, based on Faber & Jackson (1976) and the **pPXF method** by Cappellari & Emsellem (2004). Uses SDSS data and MILES stellar library.

**Author:** Jashanpreet Singh Dingra
**Supervisor:** Dr. Gulab Chand Dewangan
**Date:** September 1, 2025

```

✅ **Why this works on GitHub:**

- `diff` blocks = green/red highlight  
- `yaml` = subtle blue-ish tone  
- Code/math blocks = gray background  
- Bold/italic = emphasis  
- Fully Markdown compliant  

---

If you want, I can make an **even more visually distinct version with alternating colored blocks** for sections that will **look colorful without breaking GitHub**.  

Do you want me to do that next?
```
