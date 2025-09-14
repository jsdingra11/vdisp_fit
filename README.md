<h1 style="color:#1f618d;">vdisp_fit : A Stellar Velocity Dispersion Fitting Pipeline</h1>
<hr style="border:1px solid #1f618d;">

<h2 style="color:#2874a6;">Abstract</h2>

<p style="color:#17202a;">
Stellar velocity dispersion (<b>&sigma;</b>) is a fundamental physical parameter that quantifies the random motions of stars within a galaxy and provides critical insights into its dynamical mass, gravitational potential, and evolutionary history. Accurate measurement of <b>&sigma;</b> is essential for understanding galaxy formation, scaling relations between galaxies and their central black holes, and the mass-to-light ratio of stellar populations.
</p>

<p style="color:#17202a;">
This project, <b>vdisp_fit</b>, presents a robust and automated Python pipeline for measuring stellar velocity dispersions from galaxy spectra, particularly from large spectroscopic surveys such as the Sloan Digital Sky Survey (SDSS). The pipeline employs direct spectral fitting by convolving high-resolution stellar templates from the MILES library with a Gaussian broadening function, combined with an additive polynomial to account for continuum variations. A stable <b>&chi;<sup>2</sup></b> minimization routine is used to separate non-linear broadening from linear components, ensuring accurate and efficient fitting.
</p>

<p style="color:#17202a;">
<b>vdisp_fit</b> offers both a command-line interface and a library function, enabling flexible batch processing, robust error estimation, instrumental correction, and visualization of results. By providing an end-to-end automated solution, this pipeline facilitates the study of galaxy kinematics and contributes to the broader understanding of galaxy structure and evolution.
</p>

<hr style="border:1px solid #1f618d;">

<h2 style="color:#2874a6;">Key Features</h2>

<ul style="color:#1b2631;">
<li><b style="color:#117864;">Robust Fitting:</b> Utilizes a stable additive model with an analytical linear least-squares solution, avoiding issues with local minima.</li>
<li><b style="color:#117864;">Flexible Input:</b> Can process a single galaxy FITS file or a batch of files from an entire directory.</li>
<li><b style="color:#117864;">Instrumental Correction:</b> Corrects the fitted velocity dispersion for instrumental resolution differences between the galaxy and template spectra.</li>
<li><b style="color:#117864;">Error Estimation:</b> Provides a robust 1-&sigma; uncertainty on the final velocity dispersion measurement.</li>
<li><b style="color:#117864;">Command-Line Tool:</b> Can be run directly from the terminal with a simple command structure.</li>
<li><b style="color:#117864;">Plotting Capabilities:</b> Generates detailed single-fit plots and summary plots for batch processing to visualize results.</li>
<li><b style="color:#117864;">Single-File Module:</b> Entire pipeline contained in a single Python script (<code>vdisp_fit.py</code>) for easy distribution and use.</li>
</ul>

<hr style="border:1px solid #1f618d;">

<h2 style="color:#2874a6;">Methodology</h2>

<h3 style="color:#239b56;">1. Data Pre-processing</h3>
<ul style="color:#1b2631;">
<li><b>Redshift Correction:</b> Galaxy spectra are shifted to their rest-frame using the SDSS-provided redshift.</li>
<li><b>Logarithmic Rebinning:</b> Spectra are rebinned to a logarithmic wavelength scale, ensuring that kinematic broadening corresponds to a constant pixel shift.</li>
<li><b>Masking:</b> Strong emission lines and bad pixels are masked out to prevent contamination of the stellar continuum fit.</li>
</ul>

<h3 style="color:#239b56;">2. The Fitting Model</h3>
<p style="color:#1b2631;">
The observed galaxy spectrum <i>G_obs</i> is modeled as a broadened stellar template <i>T</i> combined with an additive polynomial <i>P<sub>k</sub></i> to account for continuum variations:
</p>

<div style="color:#7d3c98; background-color:#f4f6f7; padding:10px; border-radius:5px;">
$$
G_{model}(\ln \lambda) = c_0 \cdot [T(\ln \lambda) * B(v, \sigma)] + \sum_{k=1}^{N} c_k P_k(\ln \lambda)
$$
</div>

<ul style="color:#1b2631;">
<li><b>B</b> is a Gaussian broadening function.</li>
<li><b>c<sub>0</sub>,...,c<sub>N</sub></b> are linear coefficients.</li>
</ul>

<h3 style="color:#239b56;">3. &chi;<sup>2</sup> Minimization</h3>
<ul style="color:#1b2631;">
<li>The pipeline iterates through a grid of trial &sigma; values.</li>
<li>For each &sigma;, the optimal linear coefficients are solved analytically, and the corresponding &chi;<sup>2</sup> value is calculated.</li>
<li>The &sigma; that yields the minimum &chi;<sup>2</sup> is identified as the best-fit velocity dispersion.</li>
</ul>

<h3 style="color:#239b56;">4. Instrumental Correction and Error Estimation</h3>
<ul style="color:#1b2631;">
<li>The measured dispersion is corrected for instrumental broadening:
<div style="color:#7d3c98; background-color:#f4f6f7; padding:10px; border-radius:5px;">
$$
\sigma_{gal} = \sqrt{\sigma_{fit}^2 - \sigma_{inst}^2}
$$
</div>
</li>
<li>The 1-&sigma; error is estimated from the curvature of the &chi;<sup>2</sup> curve near its minimum.</li>
</ul>

<hr style="border:1px solid #1f618d;">

<h2 style="color:#2874a6;">Requirements</h2>

<pre style="background-color:#f4f6f7; padding:10px; border-radius:5px; color:#1b2631;">
pip install numpy matplotlib scipy astropy
</pre>

<hr style="border:1px solid #1f618d;">

<h2 style="color:#2874a6;">Usage</h2>

<h3 style="color:#239b56;">As a Command-Line Tool</h3>

<pre style="background-color:#f4f6f7; padding:10px; border-radius:5px; color:#1b2631;">
python vdisp_fit.py &lt;input_path&gt; &lt;template_folder&gt; [options]
</pre>

<p style="color:#1b2631;"><b>Arguments:</b></p>
<ul style="color:#1b2631;">
<li><code>&lt;input_path&gt;</code>: Path to a single galaxy FITS file or directory containing multiple FITS files.</li>
<li><code>&lt;template_folder&gt;</code>: Path to the directory containing the MILES stellar template FITS files.</li>
</ul>

<p style="color:#1b2631;"><b>Options:</b></p>
<ul style="color:#1b2631;">
<li><code>-o, --output</code>: Name of the CSV file to save results (e.g., <code>results.csv</code>).</li>
<li><code>-p, --plot</code>: Display plots of the fits.</li>
<li><code>--poly_degree</code>: Degree of the additive polynomial (default: 8).</li>
<li><code>--velscale</code>: Velocity scale in km/s per pixel (default: 60.0).</li>
<li><code>--fit_range</code>: Wavelength range for fitting in Å (e.g., <code>4200 6500</code>).</li>
<li><code>--instrumental_sigma</code>: Instrumental dispersion difference in km/s (default: 25.0).</li>
</ul>

<p style="color:#1b2631;"><b>Example:</b></p>
<pre style="background-color:#f4f6f7; padding:10px; border-radius:5px; color:#1b2631;">
python vdisp_fit.py ./galaxy_spectra/ ./miles_templates/ -o results.csv --plot
</pre>

<h3 style="color:#239b56;">As a Library Function</h3>

<pre style="background-color:#f4f6f7; padding:10px; border-radius:5px; color:#1b2631;">
from vdisp_fit import vdisp_fit

results = vdisp_fit(
    input_path='./galaxy_spectra/galaxy_1.fits',
    template_folder='./miles_templates/',
    output_csv='single_fit_result.csv',
    plot=True
)

print(results)
</pre>

<hr style="border:1px solid #1f618d;">

<h2 style="color:#2874a6;">Acknowledgements</h2>

<p style="color:#1b2631;">
This pipeline was developed as a project report at IUCAA.  
It builds upon foundational work in galaxy kinematics, including the pioneering work of Faber & Jackson (1976) and modern techniques like the <b>pPXF method</b> developed by Cappellari & Emsellem (2004).  
It utilizes publicly available data from the Sloan Digital Sky Survey (SDSS) and the MILES stellar library.
</p>

<p style="color:#1b2631;">
<b>Author:</b> Jashanpreet Singh Dingra <br>
<b>Supervisor:</b> Dr. Gulab Chand Dewangan <br>
<b>Date:</b> September 1, 2025
</p>
