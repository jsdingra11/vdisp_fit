Here’s a **well-organized GitHub README** for your `vdisp_fit` project, including proper Markdown formatting, headings, code blocks, and an acknowledgement to your supervisor. I’ve formatted it to be clean, professional, and readable on GitHub.

````markdown
# vdisp_fit: A Stellar Velocity Dispersion Fitting Pipeline

Stellar velocity dispersion (σ) is a fundamental parameter in galaxy dynamics, providing insights into the mass, formation history, and evolution of galaxies.  
With the advent of large spectroscopic surveys, there is a significant need for robust, automated tools to efficiently extract these key parameters from millions of galaxy spectra.  

The **vdisp_fit** pipeline was developed to provide a reliable and stable method for measuring stellar velocity dispersion.

---

## Key Features

- **Spectral Fitting**: Analyzes galaxy spectra by modeling them as a broadened stellar template from the MILES library.
- **Batch Processing**: Capable of processing a single galaxy FITS file or a directory of multiple files.
- **Analytical Fitting**: Uses a robust χ² minimization routine separating non-linear and linear components for stable fitting.
- **Parameter Derivations**: Provides the final velocity dispersion, its 1-σ uncertainty, and the reduced chi-squared for each fit.
- **Visualizations**: Generates plots to visualize the fit, including the galaxy spectrum, best-fit model, and residuals.

---

## Methodology

The pipeline follows these steps to achieve high-quality fits:

1. **Prepare the Data**  
   Cleans and prepares the galaxy and template spectra, including shifting the galaxy spectrum to its rest-frame and masking strong emission lines.

2. **Find the Best Fit**  
   Iterates through all 985 templates from the MILES library. For each template, it tests multiple velocity dispersions to find the best fit.

3. **Select the Best Template**  
   Compares results from all templates and selects the template and velocity dispersion that provides the overall best fit.

4. **Final Result**  
   Corrects the velocity dispersion for instrumental effects and calculates the uncertainty on the measurement.

---

## Requirements

The pipeline requires the following Python libraries:

```bash
pip install numpy matplotlib scipy astropy
````

---

## Usage

You can import `vdisp_fit` into another Python script (e.g., `test.py`) and run it as follows:

```python
# test.py

from vdisp_fit import vdisp_fit

# Define your data paths
galaxy_data_path = './sdss_spectra/'
template_library_path = './miles_templates/'
output_file_path = 'results.csv'

# Call the function directly with your desired parameters
results = vdisp_fit(
    input_path=galaxy_data_path,
    template_folder=template_library_path,
    output_csv=output_file_path,
    plot=True
)

print(results)
```

This will iterate over all galaxy spectra in the folder, fit them with all 985 templates, and provide clear results for each galaxy.

---

## Acknowledgements

* This pipeline was developed as a project report at IUCAA.
* Built upon foundational work in galaxy kinematics, including Faber & Jackson (1976) and modern techniques like pPXF.
* Uses publicly available data from the **[Sloan Digital Sky Survey (SDSS)](https://www.sdss.org/)** and the **MILES stellar library**.
* Special thanks to **[Dr. Gulab Chand Dewangan](https://www.iucaa.in/en/faculty-research/gulabd)**, my supervisor, for guidance and support throughout this project.

---

## Author

**Jashanpreet Singh Dingra**
*Date: September 1, 2025*
want me to do that?
```
