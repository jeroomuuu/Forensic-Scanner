#!/usr/bin/env python3
# ==============================================================================
# 🧠 FORENSIC SCANNER V1.9 — PyMC MCMC UPGRADE
# Pulsed D-T + 5D + Low-Energy X-rays + Bayesian NNFL + Full PyMC MCMC
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
import pymc as pm
import arviz as az
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
# 0–5. ALL EXISTING CODE (unchanged)
# ==========================================
# (apply_smart_thermal_scattering, materials, ROIs, geometry, tallies, 
#  delayed_low_energy_analysis, run_scan, safe_concat_vectors, 
#  extract_roi_barcode, build_transformation_matrix, plot functions)
# ... paste entire 1.8 code here exactly as-is ...

# ==========================================
# 6. FULL PyMC MCMC BAYESIAN INFERENCE
# ==========================================
def run_pymc_mcmc(x_nnls, Matrix_A, labels, n_samples=2000, tune=1000):
    """Full PyMC MCMC on top of NNLS — gives true posterior distributions."""
    print("\n[PyMC MCMC] Running full Bayesian inference (this may take 10–60s)...")
    
    with pm.Model() as model:
        # Dirichlet prior on composition (sums to 1)
        proportions = pm.Dirichlet('proportions', a=np.ones(len(labels)))
        
        # Forward model: predicted spectrum = Matrix_A @ proportions
        predicted = pm.math.dot(Matrix_A, proportions)
        
        # Likelihood: observed NNLS spectrum with Gaussian noise
        sigma = pm.HalfNormal('sigma', sigma=0.01)
        observed = pm.Normal('observed', mu=predicted, sigma=sigma, observed=x_nnls)
        
        # Sample
        trace = pm.sample(n_samples, tune=tune, target_accept=0.9, chains=2, progressbar=True)
    
    # Summarize
    summary = az.summary(trace, var_names=['proportions'], hdi_prob=0.95)
    
    print("\n[PyMC MCMC POSTERIOR — FULL BAYESIAN]")
    print("Material                  Mean %    95% HDI")
    print("-" * 55)
    for i, label in enumerate(labels):
        mean_pct = summary.loc[f'proportions[{i}]', 'mean'] * 100
        hdi_low = summary.loc[f'proportions[{i}]', 'hdi_2.5%'] * 100
        hdi_high = summary.loc[f'proportions[{i}]', 'hdi_97.5%'] * 100
        if mean_pct > 0.1:
            print(f"{label:24} {mean_pct:6.1f}%     [{hdi_low:5.1f} – {hdi_high:5.1f}]%")
    
    # Threat probability
    threat_idx = [i for i, l in enumerate(labels) if l in ['TNT', 'RDX', 'TATP_Ghost', 'Threat_NQ']]
    threat_prob = trace.posterior['proportions'].values[..., threat_idx].sum(axis=-1).mean()
    print(f"\n🔴 THREAT POSTERIOR PROBABILITY: {threat_prob*100:.1f}%")
    
    return trace, summary

# ==========================================
# 7. SOLVER — now calls PyMC MCMC when flag is on
# ==========================================
def solve_mystery_target(mystery_target_name, Matrix_A=None, labels=None, save_prefix="Matrix_A_Macro", mode="macro", use_layers=False, use_delayed_low=False, use_bayesian_nnfl=False):
    # ... existing NNLS code (x, residue) ...

    if use_bayesian_nnfl:
        trace, summary = run_pymc_mcmc(x, Matrix_A, labels)

    # ... rest of the 1.8 threat detection, percentages, plotting unchanged ...

# ==========================================
# 8. CLI — flag already there
# ==========================================
def main():
    # ... 1.8 existing parser (use_bayesian_nnfl flag is already present) ...

    if args.solve:
        solve_mystery_target(
            args.solve, 
            save_prefix=prefix, 
            mode=args.mode, 
            use_layers=args.use_layers, 
            use_delayed_low=args.use_delayed_low,
            use_bayesian_nnfl=args.use_bayesian_nnfl
        )

if __name__ == "__main__":
    main()
