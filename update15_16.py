#!/usr/bin/env python3
# ==============================================================================
# 🧠 FORENSIC SCANNER V1.6
# Pulsed D-T + N-16 + 5D Layered + Expanded Delayed Low-Energy X-rays (Fe/Cu/Zn/Ga)
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
# 0. SMART THERMAL SCATTERING ENGINE
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

# ==========================================
# 1. MATERIALS (unchanged from 1.5)
# ==========================================
micro_elements = { ... }  # ← your full micro_elements dict (kept exactly as in 1.5)
# (All micro_dict, air, csg_materials, material_dict code is identical as in 1.5 — omitted here for brevity but fully present in the real file)

# ==========================================
# 2. ROIs & CONFIGURATION — EXPANDED LOW-ENERGY
# ==========================================
fast_rois = { ... }      # unchanged from 1.5
thermal_rois = { ... }   # unchanged from 1.5 (the  4.1 / 9.1 MeV peaks stay)

# new delayed rois follow below to be added:
delayed_low_rois = {
    'Fe_Xray_cluster': (5.0e3, 8.0e3),      # Fe Kα/Kβ pile-up
    'Cu_Xray_cluster': (7.8e3, 8.3e3),      # Cu Kα ~8.05 keV (from Cu-64 EC)
    'Zn_Xray_cluster': (8.4e3, 8.9e3),      # Zn Kα ~8.64 keV
    'Ga_Xray_cluster': (9.1e3, 9.6e3),      # Ga Kα ~9.25 keV
    'Fe59_1099keV': 1099.2e3,
    'Fe59_1292keV': 1291.6e3,
    'Fe59_192keV': 192.3e3,
    'Fe59_142keV': 142.4e3
}

#use the default config paramaters agaain
CONFIG = {
    'sigma_certainty': 4,
    'particles_per_batch': 500000, 
    'batches': 100,                
    'energy_min_eV': 1.0,
    'energy_max_eV': 1.4e7,
    'energy_bins': 2000,
    'time_bins_ns': [0, 10, 20, 30, 40, 50, 100, 200, 300, 500, 1000],
    'use_delayed_low': False
}

# ==========================================
# 3-4. GEOMETRY, TALLIES, ANALYSIS (unchanged except low-energy bins)
# ==========================================
# (configure_beam, build_geometry_and_settings, build_tallies, delayed_low_energy_analysis are identical to 1.5 except the expanded rois above)

# ==========================================
# 5. RUN SCAN & SOLVER — updated for expanded low-energy
# ==========================================
def run_scan(...): ...   # unchanged from 1.5 (already returns mean_low)

#only merge this next part to the solve target function:
def solve_mystery_target(...):
    # ... existing code ...
    if use_delayed_low and mean_low is not None:
        energies_low = low_energies[:-1]
        roi_counts, warnings = delayed_low_energy_analysis(mean_low, energies_low)
        print("\n[DELAYED LOW-ENERGY ANALYSIS — EXPANDED]")
        for name, counts in roi_counts.items():
            print(f" ➤ {name:18}: {counts:.2e} net counts")
        for w in warnings:
            print(f"⚠️ {w}")

    # ... rest of solve unchanged from 1.5 ...

# ==========================================
# 6. VISUALIZATION — bigger & prettier low-energy inset
# ==========================================
def plot_pftna_spectra(mean_g, gamma_energies, mystery_target_name, mean_low=None):
    # ... The full dual-window plot code unchanged from 1.5 ...
    # add the below part for the new low xray function
    if mean_low is not None:
        energies_low = low_energies[:-1]
        ax_inset = fig.add_axes([0.62, 0.62, 0.32, 0.32])   # bigger inset
        ax_inset.semilogy(energies_low / 1e3, mean_low + 1e-10, 'b-', lw=1.2)
        ax_inset.set_xlim(1, 30)
        ax_inset.set_title('Delayed 1–30 keV X-rays', fontsize=9)
        ax_inset.tick_params(labelsize=7)
        # Shade all clusters
        ax_inset.axvspan(5.0, 8.0, alpha=0.18, color='orange', label='Fe K X-rays')
        ax_inset.axvspan(7.8, 8.3, alpha=0.18, color='red', label='Cu K X-rays')
        ax_inset.axvspan(8.4, 8.9, alpha=0.18, color='purple', label='Zn K X-rays')
        ax_inset.axvspan(9.1, 9.6, alpha=0.18, color='cyan', label='Ga K X-rays')
        ax_inset.legend(fontsize=6, loc='upper right')

    plt.tight_layout()
    plt.savefig(f"Scan_Report_{mystery_target_name}.png", dpi=300)
    plt.close()

# ==========================================
# 7. CLI (unchanged)
# ==========================================
def main():
    # ... exactly as in 1.5 but with the expanded delayed_low_rois already in the global scope ...
    # (parser, geometry, tallies, solve calls — all identical from 1.5 )

if __name__ == "__main__":
    main()
