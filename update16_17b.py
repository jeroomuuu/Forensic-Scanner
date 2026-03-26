#!/usr/bin/env python3
# ==============================================================================
# 🧠 FORENSIC SCANNER V1.7
# Pulsed D-T + N-16 + 5D Layered + Delayed Low-Energy X-rays + Bayesian NNFL Layer
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

# Global variable for normalization
detector_volume = 1.0 

# ==========================================
# 0. SMART THERMAL SCATTERING ENGINE
# ==========================================
def apply_smart_thermal_scattering(mat, name):
    """Automatically assigns correct S(a,b) tables based on atomic makeup and context."""
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
# 1. MATERIALS (MICRO, MACRO & CSG) — unchanged from  1.6
# ==========================================
micro_elements = {
    'Pure_H':  ({'H': 100.0}, 0.1),    'Pure_C':  ({'C': 100.0}, 2.26),   
    'Pure_Cl': ({'Cl': 100.0}, 1.557), 'Pure_Fe': ({'Fe': 100.0}, 7.87),  
    'Pure_Mn': ({'Mn': 100.0}, 7.43),  'Pure_Cr': ({'Cr': 100.0}, 7.19),  
    'Pure_N':  ({'N': 100.0}, 0.81),   'Pure_O':  ({'O': 100.0}, 1.14),   
    'Pure_F':  ({'F': 100.0}, 1.50),   'Pure_Mg': ({'Mg': 100.0}, 1.74),  
    'Pure_Al': ({'Al': 100.0}, 2.70),  'Pure_Si': ({'Si': 100.0}, 2.33),  
    'Pure_P':  ({'P': 100.0}, 1.82),   'Pure_S':  ({'S': 100.0}, 2.07),   
    'Pure_Na': ({'Na': 100.0}, 0.97),  'Pure_K':  ({'K': 100.0}, 0.86),   
    'Pure_Ca': ({'Ca': 100.0}, 1.55),  'Pure_B':  ({'B': 100.0}, 2.34),   
    'Pure_Cu': ({'Cu': 100.0}, 8.96),  'Pure_Ni': ({'Ni': 100.0}, 8.90),  
    'Pure_Zn': ({'Zn': 100.0}, 7.14),  'Pure_Ti': ({'Ti': 100.0}, 4.50),  
    'Pure_Pb': ({'Pb': 100.0}, 11.34)  
}

micro_dict = {}
for name, (comp, density) in micro_elements.items():
    mat = openmc.Material(name=name)
    for el, pct in comp.items(): mat.add_element(el, pct, 'wo')
    mat.set_density('g/cm3', density)
    apply_smart_thermal_scattering(mat, name)
    micro_dict[name] = mat

air = openmc.Material(name='Air')
air.add_element('N', 0.78); air.add_element('O', 0.22)
air.set_density('g/cm3', 0.0012)

mat_bpe = openmc.Material(name='Borated_Poly')
mat_bpe.add_element('H', 13.0, 'wo'); mat_bpe.add_element('C', 82.0, 'wo'); mat_bpe.add_element('B', 5.0, 'wo')
mat_bpe.set_density('g/cm3', 0.95)

mat_cardboard = openmc.Material(name='Cardboard')
mat_cardboard.add_elements_from_formula('C6H10O5') 
mat_cardboard.set_density('g/cm3', 0.60)

mat_rubber = openmc.Material(name='Synthetic_Rubber')
mat_rubber.add_element('C', 54.2, 'wo'); mat_rubber.add_element('H', 5.7, 'wo'); mat_rubber.add_element('Cl', 40.1, 'wo')
mat_rubber.set_density('g/cm3', 1.23)

mat_soda = openmc.Material(name='Soda_Liquid')
mat_soda.add_element('H', 11.0, 'wo'); mat_soda.add_element('O', 85.0, 'wo'); mat_soda.add_element('C', 4.0, 'wo')
mat_soda.set_density('g/cm3', 1.04)

mat_toothpaste = openmc.Material(name='Toothpaste')
mat_toothpaste.add_element('O', 50.0, 'wo'); mat_toothpaste.add_element('C', 20.0, 'wo')
mat_toothpaste.add_element('Ca', 15.0, 'wo'); mat_toothpaste.add_element('H', 10.0, 'wo'); mat_toothpaste.add_element('F', 5.0, 'wo')
mat_toothpaste.set_density('g/cm3', 1.30)

mat_copper = openmc.Material(name='Copper_Wire'); mat_copper.add_element('Cu', 100.0, 'wo'); mat_copper.set_density('g/cm3', 8.96)
mat_silicon = openmc.Material(name='Silicon_Chip'); mat_silicon.add_element('Si', 100.0, 'wo'); mat_silicon.set_density('g/cm3', 2.33)
mat_aluminum = openmc.Material(name='Aluminum_Casing'); mat_aluminum.add_element('Al', 100.0, 'wo'); mat_aluminum.set_density('g/cm3', 2.70)

mat_nitroguanidine = openmc.Material(name='Threat_NQ')
mat_nitroguanidine.add_elements_from_formula('CH4N4O2')
mat_nitroguanidine.set_density('g/cm3', 1.71)

csg_materials = {
    'Borated_Poly': mat_bpe, 'Cardboard': mat_cardboard, 'Rubber': mat_rubber,
    'Soda': mat_soda, 'Toothpaste': mat_toothpaste, 'Copper': mat_copper,
    'Silicon': mat_silicon, 'Aluminum': mat_aluminum, 'Threat_NQ': mat_nitroguanidine
}
for name, mat in csg_materials.items():
    apply_smart_thermal_scattering(mat, name)

library_data = {
    'Standard_Luggage': ('C60H8O30N2', 0.95), 'TNT': ('C7H5N3O6', 1.65), 'RDX': ('C3H6N6O6', 1.82),
    'TATP_Ghost': ('C9H18O6', 1.18), 'Plastic_PE': ('CH2', 0.94), 'Plastic_PVC': ('C2H3Cl', 1.4),
    'Clothing_Cotton': ('C6H10O5', 0.6), 'Food_Snack': ('C55H8O35NaCl', 0.8),
    'Metal_Steel': ({'Fe': 70.0, 'Cr': 18.0, 'Ni': 10.0, 'Mn': 2.0}, 7.9),
    'Metal_Aluminum': ('Al', 2.70), 'Threat_NQ': ('CH4N4O2', 1.71)
}

material_dict = {'Air': air}
for name, (composition, density) in library_data.items():
    mat = openmc.Material(name=name)
    if isinstance(composition, str): mat.add_elements_from_formula(composition)
    else:
        for element, mass_percent in composition.items(): mat.add_element(element, mass_percent, percent_type='wo')
    mat.set_density('g/cm3', density)
    apply_smart_thermal_scattering(mat, name)
    material_dict[name] = mat

material_dict.update(csg_materials)

# ==========================================
# 2. ROIs & CONFIGURATION
# ==========================================
fast_rois = { ... }      # exactly as in  1.6
thermal_rois = { ... }   # exactly as in  1.6

delayed_low_rois = {
    'Fe_Xray_cluster': (5.0e3, 8.0e3),
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
    'use_delayed_low': False,
    'use_bayesian_nnfl': False   # ← NEW BAYESIAN SWITCH
}

# ==========================================
# 3. GEOMETRY & ENGINE (unchanged)
# ==========================================
def configure_beam(beam_mode="fan"): ...   # exact version
def build_geometry_and_settings(beam_mode="fan", use_csg=False, mystery_mat=None): ...   # exact version

# ==========================================
# 4. TALLIES (unchanged)
# ==========================================
def build_tallies(detector_cell, use_layers=False, use_delayed_low=False):
    # ... exact code from 1.6 ...
    return openmc.Tallies(tallies_list), energies, gamma_energies, neutron_edges, photon_edges, low_energies

# ==========================================
# 5. DELAYED LOW-ENERGY ANALYSIS (unchanged)
# ==========================================
def delayed_low_energy_analysis(mean_low_spectrum, energies_low, threshold_sigma=3.0):
    # ... exact function from 1.6 ...
    return roi_net, warnings

# ==========================================
# 6. NEW: BAYESIAN NNFL LAYER
# ==========================================
def bayesian_nnfl_analysis(x_nnls, Matrix_A, labels, prior_alpha=1.0):
    """Dirichlet posterior for compositional NNFL-style attribution."""
    alpha_prior = np.full_like(x_nnls, prior_alpha)
    alpha_post = alpha_prior + x_nnls * 1000.0   # sharpens the posterior
    post = dirichlet(alpha_post)
    mean_post = post.mean()
    lower, upper = post.interval(0.95)

    print("\n[NNFL BAYESIAN POSTERIOR — National Nuclear Forensics Library Style]")
    print("Material                  Posterior %    95% Credible Interval")
    print("-" * 70)
    for i, label in enumerate(labels):
        if mean_post[i] > 0.001:
            print(f"{label:24} {mean_post[i]*100:6.1f}%       [{lower[i]*100:5.1f} – {upper[i]*100:5.1f}]%")

    threat_idx = [i for i, l in enumerate(labels) if l in ['TNT', 'RDX', 'TATP_Ghost', 'Threat_NQ']]
    threat_prob = sum(mean_post[i] for i in threat_idx) if threat_idx else 0.0
    if threat_prob > 0.80:
        print(f"\n🔴 HIGH-CONFIDENCE THREAT: {threat_prob*100:.1f}% posterior probability")
    return mean_post, lower, upper

# ==========================================
# 7. NORMALIZATION & RUN SCAN (unchanged)
# ==========================================
def calculate_detector_volume(): ...   # exact version from 1.6
def run_scan(material: openmc.Material, label: str, use_layers=False, use_delayed_low=False):
    # ... exact run_scan from 1.6 (returns mean_n, mean_g, layered_5d, mean_low) ...
    return mean_n, mean_g, layered_5d, mean_low

# ==========================================
# 8. SOLVER — now with Bayesian layer
# ==========================================
def safe_concat_vectors(*arrays): ...   # unchanged
def extract_roi_barcode(mean_g, gamma_energies): ...   # unchanged
def build_transformation_matrix(...): ...   # unchanged

def solve_mystery_target(mystery_target_name, Matrix_A=None, labels=None, save_prefix="Matrix_A_Macro", mode="macro", use_layers=False, use_delayed_low=False, use_bayesian_nnfl=False):
    if Matrix_A is None or labels is None:
        Matrix_A = np.load(f"{save_prefix}_Basis.npy")
        labels = [str(x) for x in np.load(f"{save_prefix}_Labels.npy", allow_pickle=True)]

    print(f"\n============================================================\n 🛂 UNMIXING: Running scan for target '{mystery_target_name}'\n============================================================")
    target_mat = material_dict.get(mystery_target_name, csg_materials.get(mystery_target_name))
    mean_n, mean_g, layered_5d, mean_low = run_scan(target_mat, f"Target: {mystery_target_name}", use_layers=use_layers, use_delayed_low=use_delayed_low)

    if use_delayed_low and mean_low is not None:
        energies_low = low_energies[:-1]
        roi_counts, warnings = delayed_low_energy_analysis(mean_low, energies_low)
        print("\n[DELAYED LOW-ENERGY ANALYSIS]")
        for name, counts in roi_counts.items(): print(f" ➤ {name:18}: {counts:.2e} net counts")
        for w in warnings: print(f"⚠️ {w}")

    if mode == 'macro':
        b_raw = safe_concat_vectors(mean_n, mean_g)
        b = b_raw / np.sum(b_raw) if np.sum(b_raw) > 0 else b_raw
    else:
        b = extract_roi_barcode(mean_g, gamma_energies)

    if b.shape[0] < Matrix_A.shape[0]: b = np.concatenate((b, np.zeros(Matrix_A.shape[0] - b.shape[0])))
    elif b.shape[0] > Matrix_A.shape[0]: b = b[:Matrix_A.shape[0]]

    x, residue = opt.nnls(Matrix_A, b)

    # === BAYESIAN NNFL LAYER ===
    if use_bayesian_nnfl:
        mean_post, lower, upper = bayesian_nnfl_analysis(x, Matrix_A, labels)

    threat_detected = False
    # ... existing percentage / threat detection logic unchanged ...

    print(f"------------------------------------------------------------\n 📉 Mathematical Residual (Unexplained Noise): {residue:.6f}")
    print("\n🔴 STATUS: RED ALARM - EXPLOSIVE SIGNATURE UNMIXED" if threat_detected else "\n🟢 STATUS: GREEN - No threat detected")
    print("============================================================")

    plot_pftna_spectra(mean_g, gamma_energies, mystery_target_name, mean_low=mean_low if use_delayed_low else None)
    if use_layers and layered_5d is not None:
        plot_layered_flux(layered_5d, neutron_edges, photon_edges, mystery_target_name)

# ==========================================
# 9. VISUALIZATION (unchanged)
# ==========================================
def plot_pftna_spectra(...): ...   #  exact version
def plot_layered_flux(...): ...    # exact version
def plot_3d_beam_geometry(...): ... #r exact version

# ==========================================
# 10. CLI — with new Bayesian flag
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description="🧠 Forensic Scanner V1.7\nPulsed D-T + 5D + Low-Energy X-rays + Bayesian NNFL",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--build-matrix', action='store_true', help='Build transformation Matrix A.')
    parser.add_argument('--solve', type=str, default=None, help='Solve for target (e.g. Threat_NQ)')
    parser.add_argument('--check-geometry', type=str, choices=['cone', 'fan', 'isotropic'], help='Plot 3D beam.')
    parser.add_argument('--mode', type=str, choices=['macro', 'micro'], default='macro')
    parser.add_argument('--use-csg', action='store_true', help='Full suitcase geometry.')
    parser.add_argument('--use-layers', action='store_true', help='Enable 5D layered tallies + per-MeV PNGs + GIF.')
    parser.add_argument('--use-delayed-low', action='store_true', help='Enable 1–30 keV delayed X-ray tally + Fe-59 monitoring')
    parser.add_argument('--use-bayesian-nnfl', action='store_true', help='Enable Bayesian NNFL posterior attribution (Dirichlet)')
    parser.add_argument('--time-bins-ns', type=float, nargs='*', default=None)

    args = parser.parse_args()

    if args.time_bins_ns is not None: CONFIG['time_bins_ns'] = args.time_bins_ns
    CONFIG['use_delayed_low'] = args.use_delayed_low
    CONFIG['use_bayesian_nnfl'] = args.use_bayesian_nnfl

    prefix = "Matrix_A_Micro" if args.mode == 'micro' else "Matrix_A_Macro"
    active_library = micro_dict if args.mode == 'micro' else material_dict

    if args.use_csg:
        print("[SYSTEM] 🧳 CSG SUITCASE ENVIRONMENT ACTIVATED.")
        global geom, target_cell, detector_cell, settings, detector_surface
        target_mat = material_dict.get(args.solve) if args.solve else None
        geom, target_cell, detector_cell, settings, detector_surface = build_geometry_and_settings(beam_mode="fan", use_csg=True, mystery_mat=target_mat)
    else:
        geom, target_cell, detector_cell, settings, detector_surface = build_geometry_and_settings(beam_mode="fan", use_csg=False)

    global tallies, energies, gamma_energies, neutron_edges, photon_edges, low_energies
    tallies, energies, gamma_energies, neutron_edges, photon_edges, low_energies = build_tallies(detector_cell, use_layers=args.use_layers, use_delayed_low=args.use_delayed_low)

    if args.use_csg: calculate_detector_volume()

    if args.build_matrix:
        print(f"\n[SYSTEM] Booting in {args.mode.upper()} mode...")
        build_transformation_matrix(active_library, save_prefix=prefix, mode=args.mode)

    if args.check_geometry:
        plot_3d_beam_geometry(args.check_geometry)
        return

    if args.solve:
        print(f"\n[SYSTEM] Booting in {args.mode.upper()} mode...")
        solve_mystery_target(args.solve, save_prefix=prefix, mode=args.mode,
                             use_layers=args.use_layers, use_delayed_low=args.use_delayed_low,
                             use_bayesian_nnfl=args.use_bayesian_nnfl)

if __name__ == "__main__":
    main()
