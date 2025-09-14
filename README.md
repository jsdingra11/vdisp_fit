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
*(Make sure to keep Mile Template, vdisp_fit.py in the same working folder.)* <br>
This will iterate over all galaxy spectra in the folder, fit them with all 985 templates, and provide clear results for each galaxy.

---

## Acknowledgements

This project would not have been possible without the invaluable guidance, encouragement, and mentorship of **[Dr. Gulab Chand Dewangan](https://www.iucaa.in/en/faculty-research/gulabd)**.  
I am deeply grateful for his patience in answering my questions, for inspiring me to explore complex problems in galaxy kinematics, and for providing unwavering support throughout the development of this pipeline.  
His expertise and dedication have been a constant source of learning and motivation, and I am truly humbled to have had the opportunity to work under his guidance.  

I also acknowledge the use of publicly available data from the **[Sloan Digital Sky Survey (SDSS)](https://www.sdss.org/)** and the **MILES stellar library**, as well as the foundational work in galaxy kinematics by Faber & Jackson (1976) and modern spectral fitting techniques like pPXF.

---

## Author

**Jashanpreet Singh Dingra** <br>
**Date: September 1, 2025**
