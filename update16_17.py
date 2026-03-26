#!/usr/bin/env python3
# ==============================================================================
# 🧠 FORENSIC SCANNER V1.7
# Pulsed D-T + 5D Layered + Delayed Low-Energy X-rays + Bayesian NNFL Layer
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
from scipy.stats import dirichlet
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

# (All the materials, micro_dict, csg_materials, material_dict are identical to 1.6)

# ==========================================
# 1. ROIs & CONFIG (unchanged + new Bayesian flag)
# ==========================================
fast_rois = { ... }      # unchanged
thermal_rois = { ... }   # unchanged
delayed_low_rois = { ... }  # unchanged from 1.6

CONFIG = {
    'sigma_certainty': 4,
    'particles_per_batch': 500000,
    'batches': 100,
    'energy_min_eV': 1.0,
    'energy_max_eV': 1.4e7,
    'energy_bins': 2000,
    'time_bins_ns': [0, 10, 20, 30, 40, 50, 100, 200, 300, 500, 1000],
    'use_delayed_low': False,
    'use_bayesian_nnfl': False   # ← NEW SWITCH
}

# ==========================================
# 2. GEOMETRY + TALLIES (unchanged)
# ==========================================
# (configure_beam, build_geometry_and_settings, build_tallies, delayed_low_energy_analysis unchanged)

# ==========================================
# 3. NEW: Bayesian NNFL Layer
# ==========================================
def bayesian_nnfl_analysis(x_nnls, Matrix_A, labels, prior_alpha=1.0):
    """Dirichlet posterior update for compositional attribution (NNFL-style)."""
    # Prior: uniform Dirichlet (or weak informative)
    alpha_prior = np.full_like(x_nnls, prior_alpha)
    # Likelihood scaling from NNLS solution
    alpha_post = alpha_prior + x_nnls * 1000.0   # scaling factor for sharpness
    # Draw posterior
    post = dirichlet(alpha_post)
    mean_post = post.mean()
    # 95% credible intervals
    lower, upper = post.interval(0.95)
    
    print("\n[NNFL BAYESIAN POSTERIOR]")
    print("Material                  Posterior %    95% Credible Interval")
    print("-" * 70)
    for i, label in enumerate(labels):
        if mean_post[i] > 0.001:  # only show meaningful contributions
            print(f"{label:24} {mean_post[i]*100:6.1f}%       [{lower[i]*100:5.1f} – {upper[i]*100:5.1f}]%")
    
    # Threat confidence
    threat_idx = [i for i, l in enumerate(labels) if l in ['TNT', 'RDX', 'TATP_Ghost', 'Threat_NQ']]
    threat_prob = sum(mean_post[i] for i in threat_idx) if threat_idx else 0.0
    if threat_prob > 0.80:
        print(f"\n🔴 HIGH-CONFIDENCE THREAT: {threat_prob*100:.1f}% posterior probability")
    return mean_post, lower, upper

# ==========================================
# 4. RUN SCAN & SOLVER (now with Bayesian)
# ==========================================
def run_scan(...): ...   # unchanged

def solve_mystery_target(mystery_target_name, Matrix_A=None, labels=None, save_prefix="Matrix_A_Macro", mode="macro", use_layers=False, use_delayed_low=False, use_bayesian_nnfl=False):
    # ... existing NNLS code ...
    x, residue = opt.nnls(Matrix_A, b)

    # === NEW BAYESIAN NNFL LAYER ===
    if use_bayesian_nnfl:
        mean_post, lower, upper = bayesian_nnfl_analysis(x, Matrix_A, labels)

    # ... rest of percentage/threat logic unchanged ...

    print(f"------------------------------------------------------------\n 📉 Mathematical Residual (Unexplained Noise): {residue:.6f}")
    print("\n🔴 STATUS: RED ALARM - EXPLOSIVE SIGNATURE UNMIXED" if threat_detected else "\n🟢 STATUS: GREEN - No threat detected")
    print("============================================================")

    plot_pftna_spectra(mean_g, gamma_energies, mystery_target_name, mean_low=mean_low if use_delayed_low else None)
    if use_layers and layered_5d is not None: plot_layered_flux(...)

# ==========================================
# 5. CLI — new flag added
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description="🧠 Forensic Scanner V1.7\nPulsed D-T + 5D + Low-Energy X-rays + Bayesian NNFL",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # ... all previous arguments ...
    parser.add_argument('--use-bayesian-nnfl', action='store_true', help='Enable Bayesian NNFL posterior attribution (Dirichlet)')

    args = parser.parse_args()
    CONFIG['use_delayed_low'] = args.use_delayed_low
    CONFIG['use_bayesian_nnfl'] = args.use_bayesian_nnfl   # ← activated here

    # ... geometry and tallies building unchanged ...

    if args.solve:
        solve_mystery_target(
            args.solve, save_prefix=prefix, mode=args.mode,
            use_layers=args.use_layers,
            use_delayed_low=args.use_delayed_low,
            use_bayesian_nnfl=args.use_bayesian_nnfl   # ← passed through
        )

if __name__ == "__main__":
    main()
