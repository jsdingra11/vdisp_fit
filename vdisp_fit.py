# -----------------------------------------------------------------------------
#         vdisp_fit: Velocity Dispersion Fitting Module
# -----------------------------------------------------------------------------
#
# Author: Jashanpreet Singh Dingra 
# Date: September 1, 2025
#
# Description:
# This single-file module provides a robust pipeline to measure the velocity
# dispersion of galaxies. It is designed to be both a callable library function
# (`vdisp_fit`) and a standalone command-line tool.
#
# The core of the module uses a professional fitting technique that separates the
# linear (continuum + template scaling) and non-linear (broadening) parts of the
# fit. For each trial sigma, the optimal template scaling and additive polynomial
# are solved for analytically using a linear least-squares fit, which is extremely
# robust and avoids local minima issues.
#
# Key Features:
# - Flexible input: Handles single files or entire directories.
# - Robust fitting: Uses a stable additive model with linear least-squares.
# - Error Estimation: Calculates 1-sigma uncertainties on the dispersion.
# - Plotting: Can generate detailed single-fit plots or batch summary plots.
# -----------------------------------------------------------------------------

import argparse
import csv
from pathlib import Path
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from astropy.io import fits
from astropy.constants import c
from numpy.polynomial import legendre

# --- Constants ---
C_KMS = c.to('km/s').value

# -----------------------------------------------------------------------------
# Core Algorithmic Functions
# -----------------------------------------------------------------------------

def load_spectrum(filepath, is_sdss=False):
    """ Loads a spectrum from a FITS file, supporting SDSS and MILES formats. """
    try:
        with fits.open(filepath) as hdul:
            if is_sdss:
                data = hdul[1].data
                flux, ivar = data['flux'], data['ivar']
                wave = 10**data['loglam']
                redshift = hdul[2].data['Z'][0]
                return wave, flux, redshift, ivar
            else:
                flux = hdul[0].data.flatten()
                header = hdul[0].header
                crval, cdelt = header['CRVAL1'], header['CDELT1']
                wave = crval + (np.arange(len(flux)) - (header['CRPIX1'] - 1)) * cdelt
                return wave, flux
    except Exception as e:
        print(f"ERROR: Could not load FITS file {filepath}: {e}")
        return (None, None, None, None) if is_sdss else (None, None)

def log_rebin(wave, flux, velscale):
    """ Rebins a spectrum to a logarithmic wavelength scale. """
    if wave is None or flux is None: return None, None
    ln_wave_range = np.log(wave[[0, -1]])
    n_pixels = int(np.ceil((ln_wave_range[1] - ln_wave_range[0]) / (velscale / C_KMS)))
    new_ln_wave = np.linspace(ln_wave_range[0], ln_wave_range[1], n_pixels)
    new_wave = np.exp(new_ln_wave)
    new_flux = np.interp(new_wave, wave, flux)
    return new_wave, new_flux

def broaden_spectrum(flux, sigma_kms, velscale):
    """ Broadens a spectrum by a Gaussian kernel. """
    if sigma_kms <= 0: return flux
    sigma_pixels = sigma_kms / velscale
    kernel_size = int(np.ceil(sigma_pixels * 10))
    if kernel_size % 2 == 0: kernel_size += 1
    x = np.arange(kernel_size) - kernel_size // 2
    kernel = np.exp(-0.5 * (x / sigma_pixels)**2)
    kernel /= kernel.sum()
    return signal.convolve(flux, kernel, mode='same')

def _estimate_sigma_error(sigma_grid, chi2_grid):
    """ Estimates the 1-sigma error on sigma by fitting a parabola to the chi^2 curve. """
    try:
        min_idx = np.argmin(chi2_grid)
        fit_slice = slice(max(0, min_idx - 5), min_idx + 6)
        x, y = sigma_grid[fit_slice], chi2_grid[fit_slice]
        p, V = np.polyfit(x, y, 2, cov=True)
        if p[0] <= 0: return None
        return 1 / np.sqrt(p[0])
    except (np.linalg.LinAlgError, ValueError):
        return None

def _fit_single_template(gal_flux, template_flux, velscale, weights, mask, poly_degree):
    """
    Core fitting routine for a single template using a robust additive model.
    Model: Galaxy ~ c0 * Broadened_Template + Additive_Polynomial
    """
    sigma_grid_coarse = np.arange(50, 451, 10)
    chi2_values = []
    
    x_coords = np.linspace(-1, 1, len(gal_flux))
    leg_basis = np.vstack([legendre.legval(x_coords, [0]*k + [1]) for k in range(poly_degree + 1)]).T

    w = np.sqrt(weights[mask])
    w_gal_flux = gal_flux[mask] * w

    for sigma_val in sigma_grid_coarse:
        broadened_template = broaden_spectrum(template_flux, sigma_val, velscale)
        A = np.hstack([broadened_template[mask, np.newaxis], leg_basis[mask]])
        w_A = A * w[:, np.newaxis]
        try:
            _, res, _, _ = np.linalg.lstsq(w_A, w_gal_flux, rcond=None)
            chi2 = res[0] if len(res) > 0 else np.inf
        except np.linalg.LinAlgError: chi2 = np.inf
        chi2_values.append(chi2)
            
    min_idx_coarse = np.argmin(chi2_values)
    best_sigma_coarse = sigma_grid_coarse[min_idx_coarse]
    
    sigma_grid_fine = np.linspace(best_sigma_coarse - 15, best_sigma_coarse + 15, 31)
    min_chi2, best_sigma, best_coeffs = np.inf, best_sigma_coarse, None
    fine_chi2_values = []

    for sigma_val in sigma_grid_fine:
        broadened_template = broaden_spectrum(template_flux, sigma_val, velscale)
        A = np.hstack([broadened_template[mask, np.newaxis], leg_basis[mask]])
        w_A = A * w[:, np.newaxis]
        try:
            coeffs, res, _, _ = np.linalg.lstsq(w_A, w_gal_flux, rcond=None)
            chi2 = res[0] if len(res) > 0 else np.inf
        except np.linalg.LinAlgError: chi2 = np.inf
        fine_chi2_values.append(chi2)
        if chi2 < min_chi2:
            min_chi2, best_sigma, best_coeffs = chi2, sigma_val, coeffs

    sigma_err = _estimate_sigma_error(sigma_grid_fine, np.array(fine_chi2_values))
    return best_sigma, sigma_err, min_chi2, best_coeffs

# -----------------------------------------------------------------------------
# Plotting Functions
# -----------------------------------------------------------------------------

def plot_single_fit(fit_data):
    """ Displays a detailed plot for a single galaxy fit. """
    res = fit_data['result']
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    err_str = f" \u00B1 {res['sigma_err_kms']:.2f}" if res['sigma_err_kms'] is not None else ""
    title = (f"Galaxy: {res['galaxy_file']} | "
             f"$\sigma_{{true}} = {res['best_sigma_kms']:.2f}{err_str}$ km/s | "
             f"$\chi_v^2 = {res['min_chi2_reduced']:.3f}$\nTemplate: {res['best_template_file']}")
    fig.suptitle(title, fontsize=16)

    ax1.plot(fit_data['wave'], fit_data['flux'], 'k-', label='Galaxy Spectrum', lw=1.5, ds='steps-mid')
    ax1.plot(fit_data['wave'], fit_data['model'], 'r-', label='Best-fit Model', lw=1.5, alpha=0.8)
    ax1.plot(fit_data['wave'][~fit_data['mask']], fit_data['flux'][~fit_data['mask']], 'x', color='limegreen', ms=6, label='Masked Pixels')
    ax1.set_ylabel("Flux"), ax1.legend(), ax1.set_title("Fit Comparison"), ax1.set_ylim(bottom=0)

    residuals = (fit_data['flux'] - fit_data['model']) * np.sqrt(fit_data['ivar'])
    ax2.plot(fit_data['wave'][fit_data['mask']], residuals[fit_data['mask']], 'g-', lw=1, ds='steps-mid')
    ax2.axhline(0, color='grey', linestyle='--'), ax2.set_xlabel("Rest-frame Wavelength ($\AA$)")
    ax2.set_ylabel("Residuals (Flux/Error)"), ax2.set_ylim(-5, 5)
    plt.tight_layout(rect=[0, 0.03, 1, 0.94]), plt.show()

def plot_batch_summary(all_fit_data):
    """ Displays a summary plot with multiple subplots for a batch job. """
    n_plots = len(all_fit_data)
    if n_plots == 0: return
    
    n_cols = int(np.ceil(np.sqrt(n_plots)))
    n_rows = int(np.ceil(n_plots / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows), squeeze=False)
    fig.suptitle("Batch Fitting Summary", fontsize=20)
    
    for i, fit_data in enumerate(all_fit_data):
        ax = axes.flat[i]
        res = fit_data['result']
        err_str = f"$\pm${res['sigma_err_kms']:.1f}" if res['sigma_err_kms'] else ""
        ax.plot(fit_data['wave'], fit_data['flux'], 'k-', lw=1, ds='steps-mid')
        ax.plot(fit_data['wave'], fit_data['model'], 'r-', lw=1, alpha=0.8)
        ax.set_title(f"{res['galaxy_file']}\n$\sigma={res['best_sigma_kms']:.1f}{err_str}$ km/s, $\chi_v^2={res['min_chi2_reduced']:.2f}$")
        ax.tick_params(axis='x', labelsize=8), ax.tick_params(axis='y', labelsize=8)

    for j in range(i + 1, len(axes.flat)):
        axes.flat[j].set_visible(False)
    plt.tight_layout(rect=[0, 0, 1, 0.96]), plt.show()

# -----------------------------------------------------------------------------
# Main Library Function
# -----------------------------------------------------------------------------

def vdisp_fit(input_path, template_folder, output_csv=None, plot=False,
              poly_degree=8, velscale=60.0, fit_range=(4200, 6500),
              instrumental_sigma=25.0):
    """
    Measures velocity dispersion for one or more galaxies.

    Args:
        input_path (str or Path): Path to a galaxy FITS file or a directory of FITS files.
        template_folder (str or Path): Path to the folder with template FITS files.
        output_csv (str or Path, optional): If provided, saves results to this CSV file. Defaults to None.
        plot (bool, optional): If True, shows plots of the fits. Defaults to False.
        poly_degree (int, optional): Degree of the additive polynomial. Defaults to 8.
        velscale (float, optional): Velocity scale in km/s per pixel. Defaults to 60.0.
        fit_range (tuple, optional): Wavelength range (min, max) for fitting. Defaults to (4200, 6500).
        instrumental_sigma (float, optional): Instrumental dispersion difference in km/s. Defaults to 25.0.

    Returns:
        list: A list of dictionaries, where each dictionary contains the fitting result for one galaxy.
    """
    input_path, template_folder = Path(input_path), Path(template_folder)
    
    # --- Find Input Files ---
    if input_path.is_dir():
        galaxy_files = sorted(list(input_path.glob('*.fits')))
    elif input_path.is_file() and input_path.suffix.lower() == '.fits':
        galaxy_files = [input_path]
    else:
        print(f"FATAL ERROR: No FITS files found at '{input_path}'.")
        return []

    template_files = sorted(list(template_folder.glob('*.fits')))
    if not template_files:
        print(f"FATAL ERROR: No template FITS files found in '{template_folder}'.")
        return []

    # --- Setup ---
    all_results = []
    all_plot_data = []
    emission_lines = {
        'H-beta': (4850, 4875), '[OIII]': (4950, 5020), '[NI]': (5190, 5210),
        'HeI': (5865, 5885), '[OI]': (6290, 6310), '[NII]': (6535, 6595),
    }

    print(f"Found {len(galaxy_files)} galaxy spectra and {len(template_files)} templates.")
    
    # --- Main Processing Loop ---
    for i, galaxy_file in enumerate(galaxy_files):
        print(f"\n--- Processing ({i+1}/{len(galaxy_files)}): {galaxy_file.name} ---")
        
        # 1. Load and preprocess galaxy data
        gal_wave_obs, gal_flux, gal_z, gal_ivar = load_spectrum(galaxy_file, is_sdss=True)
        if gal_wave_obs is None: continue

        gal_wave_rest = gal_wave_obs / (1 + gal_z)
        gal_wave_log, gal_flux_log = log_rebin(gal_wave_rest, gal_flux, velscale)
        _, gal_ivar_log = log_rebin(gal_wave_rest, gal_ivar, velscale)

        fit_idx = np.where((gal_wave_log >= fit_range[0]) & (gal_wave_log <= fit_range[1]))
        gal_wave_fit, gal_flux_fit, gal_ivar_fit = gal_wave_log[fit_idx], gal_flux_log[fit_idx], gal_ivar_log[fit_idx]
        
        mask = gal_ivar_fit > 0
        for line_range in emission_lines.values():
            em_mask = (gal_wave_fit >= line_range[0]) & (gal_wave_fit <= line_range[1])
            mask[em_mask] = False

        # 2. Iterate through templates to find the best fit
        min_chi2_reduced, best_template_file, best_fit_params = np.inf, None, None
        overall_best_sigma_fit, overall_sigma_err = None, None

        for j, template_file in enumerate(template_files):
            print(f"\r  Fitting templates: {j+1}/{len(template_files)}", end="")
            template_wave, template_flux = load_spectrum(Path(template_file))
            if template_wave is None: continue
            
            _, template_flux_log = log_rebin(template_wave, template_flux, velscale)
            template_interp = np.interp(gal_wave_fit, np.exp(np.linspace(np.log(template_wave[0]), np.log(template_wave[-1]), len(template_flux_log))), template_flux_log)
            
            sigma, sigma_err, chi2, coeffs = _fit_single_template(gal_flux_fit, template_interp, velscale, gal_ivar_fit, mask, poly_degree)
            if sigma is None: continue

            dof = np.sum(mask) - (poly_degree + 1 + 1)
            current_chi2_reduced = chi2 / dof if dof > 0 else np.inf

            if current_chi2_reduced < min_chi2_reduced:
                min_chi2_reduced, overall_best_sigma_fit, overall_sigma_err = current_chi2_reduced, sigma, sigma_err
                best_template_file, best_fit_params = template_file, coeffs

        print("\n  Fit complete.")
        if best_template_file is None:
            print("  -> WARNING: Fit failed for this galaxy.")
            continue
        
        # 3. Final calculations
        sigma_corr_sq = instrumental_sigma**2
        final_sigma = np.sqrt(overall_best_sigma_fit**2 - sigma_corr_sq) if overall_best_sigma_fit**2 > sigma_corr_sq else 0.0

        result = {
            'galaxy_file': galaxy_file.name, 'best_sigma_kms': final_sigma,
            'sigma_err_kms': overall_sigma_err, 'min_chi2_reduced': min_chi2_reduced,
            'best_template_file': Path(best_template_file).name
        }
        all_results.append(result)
        
        err_str = f" \u00B1 {result['sigma_err_kms']:.2f}" if result['sigma_err_kms'] is not None else ""
        print(f"  -> RESULT: Sigma = {result['best_sigma_kms']:.2f}{err_str} km/s, Chi^2_red = {result['min_chi2_reduced']:.3f}")

        # 4. Store data for plotting if requested
        if plot:
            best_template_wave, best_template_flux = load_spectrum(Path(best_template_file))
            _, best_template_flux_log = log_rebin(best_template_wave, best_template_flux, velscale)
            best_template_interp = np.interp(gal_wave_fit, np.exp(np.linspace(np.log(best_template_wave[0]), np.log(best_template_wave[-1]), len(best_template_flux_log))), best_template_flux_log)
            best_broadened_template = broaden_spectrum(best_template_interp, overall_best_sigma_fit, velscale)
            x_coords_plot = np.linspace(-1, 1, len(gal_wave_fit))
            additive_poly = legendre.legval(x_coords_plot, best_fit_params[1:])
            best_fit_model = best_fit_params[0] * best_broadened_template + additive_poly
            all_plot_data.append({
                'result': result, 'wave': gal_wave_fit, 'flux': gal_flux_fit,
                'model': best_fit_model, 'mask': mask, 'ivar': gal_ivar_fit
            })

    # --- Final Output ---
    if output_csv:
        try:
            with open(output_csv, 'w', newline='') as csvfile:
                fieldnames = ['galaxy_file', 'best_sigma_kms', 'sigma_err_kms', 'min_chi2_reduced', 'best_template_file']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_results)
            print(f"\nSuccessfully saved results to {output_csv}")
        except IOError as e:
            print(f"\nERROR: Could not write to CSV file {output_csv}: {e}")
            
    if plot:
        if len(all_plot_data) == 1:
            plot_single_fit(all_plot_data[0])
        elif len(all_plot_data) > 1:
            plot_batch_summary(all_plot_data)

    return all_results

# -----------------------------------------------------------------------------
# Command-Line Execution
# -----------------------------------------------------------------------------

def _parse_cli_args():
    """ Parses command-line arguments for standalone execution. """
    parser = argparse.ArgumentParser(description="Fit galaxy velocity dispersion.")
    parser.add_argument("input_path", type=Path, help="Path to a galaxy FITS file or a directory of FITS files.")
    parser.add_argument("template_folder", type=Path, help="Path to the folder with template FITS files.")
    parser.add_argument("-o", "--output", type=Path, default="velocity_dispersion_results.csv", help="Output CSV file name.")
    parser.add_argument("-p", "--plot", action="store_true", help="Show plots of the fits.")
    parser.add_argument("--poly_degree", type=int, default=8, help="Degree of the additive polynomial.")
    parser.add_argument("--velscale", type=float, default=60.0, help="Velocity scale in km/s per pixel.")
    parser.add_argument("--fit_range", type=float, nargs=2, default=[4200, 6500], help="Wavelength range for fitting (min max).")
    parser.add_argument("--instrumental_sigma", type=float, default=25.0, help="Instrumental dispersion difference in km/s.")
    return parser.parse_args()

if __name__ == '__main__':
    """ This block runs when the script is executed directly from the command line. """
    args = _parse_cli_args()
    vdisp_fit(
        input_path=args.input_path,
        template_folder=args.template_folder,
        output_csv=args.output,
        plot=args.plot,
        poly_degree=args.poly_degree,
        velscale=args.velscale,
        fit_range=tuple(args.fit_range),
        instrumental_sigma=args.instrumental_sigma
    )