#!/usr/bin/env python3
# ==============================================================================
# 🧠 FORENSIC SCANNER V1.5_beta
# Pulsed D-T + N-16 + 5D Layered + Delayed Low-Energy X-rays (Fe-59 etc.)
# ==============================================================================

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.ndimage import gaussian_filter1d
from matplotlib.animation import ArtistAnimation
import scipy.optimize as opt
import datetime

try:
    import openmc
except Exception as e:
    print("ERROR: openmc Python package not found.")
    raise

try:
    from tqdm import tqdm
    TQDM = True
except:
    TQDM = False

detector_volume = 1.0 

# ==========================================
# 0. SMART THERMAL SCATTERING + MATERIALS (unchanged)
# ==========================================
def apply_smart_thermal_scattering(mat, name):
    nuclides = [n.split('_')[0] for n in mat.get_nuclide_atom_densities()]
    name_lower = name.lower()
    if 'H' in nuclides:
        if any(x in name_lower for x in ['poly', 'plastic', 'rubber', 'pe', 'pvc', 'tnt', 'rdx', 'tatp', 'nq', 'cardboard']):
            mat.add_s_alpha_beta('c_H_in_CH2')
        else:
            mat.add_s_alpha_beta('c_H_in_H2O')
    if 'C' in nuclides and 'H' not in nuclides:
        mat.add_s_alpha_beta('c_C_in_graphite')

# (All micro_dict, air, csg_materials, material_dict definitions keep exactly the same as before in 1.4)

# ==========================================
# 1. ROIs & CONFIGURATION
# ==========================================
fast_rois = { ... }      # unchanged
thermal_rois = { ... }   # unchanged (includes new 4.1 and 9.1 MeV)

# NEW: Delayed low-energy X-ray ROIs
delayed_low_rois = {
    'Fe_Xray_cluster': (5.0e3, 8.0e3),      # Fe Kα/Kβ pile-up (main EC signature)
    'Fe59_1099keV': 1099.2e3,
    'Fe59_1292keV': 1291.6e3,
    'Fe59_192keV': 192.3e3,
    'Fe59_142keV': 142.4e3
}

CONFIG = {
    'sigma_certainty': 4,
    'particles_per_batch': 500000, 
    'batches': 100,                
    'energy_min_eV': 1.0,
    'energy_max_eV': 1.4e7,
    'energy_bins': 2000,
    'time_bins_ns': [0, 10, 20, 30, 40, 50, 100, 200, 300, 500, 1000],
    'use_delayed_low': False   # will be set by CLI
}

# ==========================================
# 2. GEOMETRY & BEAM (unchanged)
# ==========================================
def configure_beam(beam_mode="fan"): ...   # your clean pulsed version
def build_geometry_and_settings(...): ...  # unchanged

# ==========================================
# 3. TALLIES — now with optional low-energy
# ==========================================
def build_tallies(detector_cell, use_layers=False, use_delayed_low=False):
    # ... all your original tallies (scanner_signal, pgaa_signal, n16_production, 3d_mesh_flux, layered_flux_5d) ...

    tallies_list = [tally_n, tally_g, rxn_tally, tally_mesh]
    if use_layers:
        # ( 5D layered tally code unchanged)
        ...

    # NEW: Delayed low-energy X-ray tally
    if use_delayed_low:
        print("[DELAYED LOW] Adding 1–30 keV photon tally for activation X-rays")
        low_e_filter = openmc.EnergyFilter([1e3, 30e3])
        tally_low_delayed = openmc.Tally(name='delayed_low_energy_photons')
        tally_low_delayed.filters = [
            low_e_filter,
            openmc.ParticleFilter('photon'),
            openmc.TimeFilter([0.0, 1e9])   # wide gate → captures delayed emissions
        ]
        tally_low_delayed.scores = ['flux']
        tallies_list.append(tally_low_delayed)

    return openmc.Tallies(tallies_list), energies, gamma_energies, neutron_edges, photon_edges

# ==========================================
# 4. Delayed low-energy analysis function
# ==========================================
def delayed_low_energy_analysis(mean_low_spectrum, energies_low, threshold_sigma=3.0):
    """Net count extraction + Fe-59 / EC inference."""
    roi_net = {}
    background_mean = np.mean(mean_low_spectrum)
    background_std = np.std(mean_low_spectrum)

    for name, e_val in delayed_low_rois.items():
        if isinstance(e_val, tuple):  # range (X-ray cluster)
            mask = (energies_low >= e_val[0]) & (energies_low <= e_val[1])
        else:  # single peak
            mask = (energies_low >= e_val * 0.98) & (energies_low <= e_val * 1.02)
        if mask.sum() > 0:
            net = np.sum(mean_low_spectrum[mask]) - background_mean * mask.sum()
            roi_net[name] = max(net, 0)

    warnings = []
    if roi_net.get('Fe_Xray_cluster', 0) > threshold_sigma * background_std * 3000:
        warnings.append("Strong Fe K X-ray cluster detected → possible Fe-59 activation (EC branch)")
    if (roi_net.get('Fe59_1099keV', 0) > threshold_sigma * background_std * 1000 and
        roi_net.get('Fe59_1292keV', 0) > threshold_sigma * background_std * 700):
        warnings.append("Fe-59 1099 + 1292 keV pair detected → high-confidence iron activation")

    return roi_net, warnings

# ==========================================
# 5. RUN SCAN (updated to extract low-energy if enabled)
# ==========================================
def run_scan(material: openmc.Material, label: str, use_layers=False, use_delayed_low=False):
    # ... keep existing run_scan code (pulsed source, tallies, N-16 injection, 5D extraction) ...

    mean_low = None
    if use_delayed_low:
        try:
            t_low = sp.get_tally(name='delayed_low_energy_photons')
            mean_low = t_low.get_values(scores=['flux'], value='mean').flatten()
            print(f"    [LOW-E] Low-energy (1-30 keV) spectrum extracted")
        except Exception as e:
            print(f"[WARNING] Low-energy tally extraction failed: {e}")

    # ... rest of your existing return ...
    return mean_n, mean_g, layered_5d, mean_low   # added mean_low

# ==========================================
# 6. SOLVER — hook in low-energy analysis
# ==========================================
def solve_mystery_target(mystery_target_name, Matrix_A=None, labels=None, save_prefix="Matrix_A_Macro", mode="macro", use_layers=False, use_delayed_low=False):
    # ... existing solve code up to run_scan ...
    mean_n, mean_g, layered_5d, mean_low = run_scan(target_mat, f"Target: {mystery_target_name}", use_layers=use_layers, use_delayed_low=use_delayed_low)

    # === NEW: Low-energy analysis ===
    if use_delayed_low and mean_low is not None:
        energies_low = np.linspace(1e3, 30e3, len(mean_low))   # simple linear bins matching tally
        roi_counts, warnings = delayed_low_energy_analysis(mean_low, energies_low)
        print("\n[DELAYED LOW-ENERGY ANALYSIS]")
        for name, counts in roi_counts.items():
            print(f" ➤ {name:18}: {counts:.2e} net counts")
        for w in warnings:
            print(f"⚠️ {w}")

    # ... rest of solve_mystery_target here (percentages, threat detection, plot_pftna_spectra) ...

    plot_pftna_spectra(mean_g, gamma_energies, mystery_target_name, mean_low=mean_low if use_delayed_low else None)

# ==========================================
# 7. VISUALIZATION — with optional low-energy inset
# ==========================================
def plot_pftna_spectra(mean_g, gamma_energies, mystery_target_name, mean_low=None):
    # ... existing dual-window plot code from 1.4 go here ...

    # NEW: Low-energy inset (if data available)
    if mean_low is not None:
        energies_low = np.linspace(1e3, 30e3, len(mean_low))
        ax_inset = fig.add_axes([0.65, 0.65, 0.25, 0.25])
        ax_inset.semilogy(energies_low / 1e3, mean_low + 1e-10, 'b-', lw=1)
        ax_inset.set_xlim(1, 20)
        ax_inset.set_title('Delayed 1–20 keV', fontsize=8)
        ax_inset.tick_params(labelsize=6)
        ax_inset.axvspan(5, 8, alpha=0.15, color='orange')   # Fe X-ray cluster
        ax_inset.text(6.5, max(mean_low)*0.3, 'Fe K X-rays', fontsize=7, color='darkorange')

    plt.tight_layout()
    plt.savefig(f"Scan_Report_{mystery_target_name}.png", dpi=300)
    plt.close()
    print(f"[PLOT] Spectroscopic diagnostic saved to: Scan_Report_{mystery_target_name}.png")

# ==========================================
# 8. CLI (with new flag)
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description="🧠 Forensic Scanner V2.0_beta\nPulsed D-T + 5D Layered + Delayed Low-Energy X-rays",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--build-matrix', action='store_true')
    parser.add_argument('--solve', type=str, default=None)
    parser.add_argument('--check-geometry', type=str, choices=['cone', 'fan', 'isotropic'])
    parser.add_argument('--mode', type=str, choices=['macro', 'micro'], default='macro')
    parser.add_argument('--use-csg', action='store_true')
    parser.add_argument('--use-layers', action='store_true')
    parser.add_argument('--use-delayed-low', action='store_true', help='Enable 1–30 keV delayed X-ray tally + Fe-59 monitoring')
    parser.add_argument('--time-bins-ns', type=float, nargs='*', default=None)

    args = parser.parse_args()

    if args.time_bins_ns is not None:
        CONFIG['time_bins_ns'] = args.time_bins_ns
    CONFIG['use_delayed_low'] = args.use_delayed_low

    # ... geometry and tallies building (now passing use_delayed_low) ...
    tallies, energies, gamma_energies, neutron_edges, photon_edges = build_tallies(
        detector_cell, use_layers=args.use_layers, use_delayed_low=args.use_delayed_low)

    # ... rest of main() unchanged except passing use_delayed_low to solve_mystery_target ...

    if args.solve:
        solve_mystery_target(
            args.solve, save_prefix=prefix, mode=args.mode,
            use_layers=args.use_layers, use_delayed_low=args.use_delayed_low
        )

if __name__ == "__main__":
    main()
