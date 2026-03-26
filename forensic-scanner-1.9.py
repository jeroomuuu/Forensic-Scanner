#!/usr/bin/env python3
# ==============================================================================
# 🧠 FORENSIC SCANNER V1.9 ALPHA
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
# 1. MATERIALS (MICRO, MACRO & CSG)
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
fast_rois = {
    'Carbon_Inelastic': 4.44e6, 'Oxygen_Prompt': 6.13e6, 'Silicon_Prompt': 1.78e6, 
    'Aluminum_Prompt': 2.21e6, 'Sulfur_Fast': 2.23e6, 'Calcium_Inelastic': 4.42e6,
    'Nitrogen_Fast_8MeV': 8.61e6, 'Nitrogen_Fast_10MeV': 10.53e6, 'Iron_Inelastic': 0.84e6,
    'Copper_Inelastic': 0.96e6, 'Chlorine_Inelastic': 1.22e6, 'Lead_Fast': 2.61e6,
    'N16_Decay_Main': 6.128e6, 'N16_Decay_Sec': 7.115e6, 'N16_Decay_Tert': 2.741e6
}

thermal_rois = {
    'Hydrogen_Capture': 2.22e6, 'Nitrogen_Prompt1': 5.27e6, 'Nitrogen_Prompt2': 10.83e6, 
    'Fluorine_Prompt': 1.35e6, 'Potassium_Prompt': 1.53e6, 'Chlorine_Major': 6.11e6, 
    'Nitrogen_Cascade_9MeV': 9.1e6, 'Carbon_Oxygen_4MeV': 4.1e6, 'Chlorine_Secondary': 6.62e6, 
    'Iron_High': 7.63e6, 'Sodium_Prompt': 6.39e6, 'Calcium_Prompt': 6.42e6, 'Aluminum_Capture': 7.72e6
}

delayed_low_rois = {
    'Fe_Xray_cluster': (5.0e3, 8.0e3),      
    'Cu_Xray_cluster': (7.8e3, 8.3e3),      
    'Zn_Xray_cluster': (8.4e3, 8.9e3),      
    'Ga_Xray_cluster': (9.1e3, 9.6e3),      
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
    'use_bayesian_nnfl': False   
}

# ==========================================
# 3. GEOMETRY & ENGINE
# ==========================================
def configure_beam(beam_mode="fan"):
    source_space = openmc.stats.Point((0.0, 0.0, -30.0))
    angle = openmc.stats.PolarAzimuthal(
        mu=openmc.stats.Uniform(0.90, 1.0), phi=openmc.stats.Uniform(np.pi/2 - 0.05, np.pi/2 + 0.05)
    ) if beam_mode == "fan" else openmc.stats.Isotropic()
    
    sources = []
    for i in range(1, 21):
        t_center = i * 1e-9; energy_val = float(i * 1e6) if i < 20 else 19.9e6
        time_dist = openmc.stats.Uniform(t_center - 0.5e-9, t_center + 0.5e-9)
        s = openmc.IndependentSource(space=source_space, angle=angle, energy=openmc.stats.Discrete([energy_val], [1.0]), time=time_dist)
        s.particle = 'neutron'; sources.append(s)
    for i in range(19, 0, -1):
        t_center = (40 - i) * 1e-9; energy_val = float(i * 1e6)
        time_dist = openmc.stats.Uniform(t_center - 0.5e-9, t_center + 0.5e-9)
        s = openmc.IndependentSource(space=source_space, angle=angle, energy=openmc.stats.Discrete([energy_val], [1.0]), time=time_dist)
        s.particle = 'neutron'; sources.append(s)

    for s in sources: s.strength = 1.0 / len(sources)
    return sources

def build_geometry_and_settings(beam_mode="fan", use_csg=False, mystery_mat=None):
    sources = configure_beam(beam_mode)
    settings = openmc.Settings(); settings.run_mode = 'fixed source'; settings.source = sources
    settings.batches = CONFIG['batches']; settings.particles = CONFIG['particles_per_batch']
    settings.photon_transport = True; settings.output = {'summary': False, 'tallies': False}

    if use_csg and mystery_mat is not None:
        collimator_block = openmc.model.RectangularParallelepiped(-15,15,-15,15,-25,-15)
        collimator_slit   = openmc.model.RectangularParallelepiped(-15,15,-2,2,-25,-15)
        bpe_cell = openmc.Cell(name='Collimator_BPE', fill=csg_materials['Borated_Poly'], region=-collimator_block & ~(-collimator_slit))

        outer_box = openmc.model.RectangularParallelepiped(-25,25,-25,25,-5,5)
        mid_box   = openmc.model.RectangularParallelepiped(-24,24,-24,24,-4,4)
        inner_void= openmc.model.RectangularParallelepiped(-23.8,23.8,-23.8,23.8,-3.8,3.8)
        rubber_cell = openmc.Cell(name='Rubber_Shell', fill=csg_materials['Rubber'], region=-outer_box & +mid_box)
        cardboard_cell = openmc.Cell(name='Cardboard_Lining', fill=csg_materials['Cardboard'], region=-mid_box & +inner_void)

        can_outer = openmc.ZCylinder(r=3.2, x0=10, y0=10); can_inner = openmc.ZCylinder(r=3.0, x0=10, y0=10)
        can_top = openmc.ZPlane(z0=3); can_bottom = openmc.ZPlane(z0=-3)
        can_wall_cell = openmc.Cell(fill=csg_materials['Aluminum'], region=-can_outer & +can_inner & -can_top & +can_bottom)
        can_liq_cell = openmc.Cell(fill=csg_materials['Soda'], region=-can_inner & -can_top & +can_bottom)
        region_can_total = -can_outer & -can_top & +can_bottom

        detector_bpe = openmc.model.RectangularParallelepiped(-35,35,-35,35,-35,35)
        detector_bpe_inner = openmc.model.RectangularParallelepiped(-32,32,-32,32,-32,32)
        bpe_collar_cell = openmc.Cell(name='Detector_BPE_Collar', fill=csg_materials['Borated_Poly'], region=-detector_bpe & +detector_bpe_inner)

        nq_box = openmc.model.RectangularParallelepiped(-2,2,-2,2,-1,1)
        nq_cell = openmc.Cell(fill=mystery_mat, region=-nq_box)

        wire_cyl = openmc.XCylinder(r=0.5, y0=-15, z0=0)
        wire_left = openmc.XPlane(x0=-10); wire_right = openmc.XPlane(x0=10)
        wire_cell = openmc.Cell(fill=csg_materials['Copper'], region=-wire_cyl & +wire_left & -wire_right)
        laptop_box = openmc.model.RectangularParallelepiped(-20,0,-20,0,-1,1)
        laptop_cell = openmc.Cell(fill=csg_materials['Silicon'], region=-laptop_box)
        paste_cyl = openmc.YCylinder(r=2.0, x0=-15, z0=0)
        paste_top = openmc.YPlane(y0=15); paste_bottom = openmc.YPlane(y0=5)
        paste_cell = openmc.Cell(fill=csg_materials['Toothpaste'], region=-paste_cyl & -paste_top & +paste_bottom)

        air_inside_region = -inner_void & ~region_can_total & ~(-nq_box) & ~(-wire_cyl & +wire_left & -wire_right) & ~(-laptop_box) & ~(-paste_cyl & -paste_top & +paste_bottom)
        air_inside_cell = openmc.Cell(name='Suitcase_Air', fill=air, region=air_inside_region)

        detector_surface = openmc.Sphere(r=70.0); outer_boundary = openmc.Sphere(r=75.0, boundary_type='vacuum')
        detector_cell = openmc.Cell(name='Detector_Volume', fill=air, region=+outer_box & ~(-collimator_block) & ~(-detector_bpe & +detector_bpe_inner) & -detector_surface)
        collimator_air_cell = openmc.Cell(name='Collimator_Air', fill=air, region=-collimator_slit)
        world_cell = openmc.Cell(name='World', fill=air, region=+detector_surface & -outer_boundary)

        geom = openmc.Geometry([bpe_cell, collimator_air_cell, rubber_cell, cardboard_cell, can_wall_cell, can_liq_cell, nq_cell, wire_cell, laptop_cell, paste_cell, air_inside_cell, detector_cell, world_cell, bpe_collar_cell])
        return geom, None, detector_cell, settings, detector_surface
    else:
        target_region = openmc.model.RectangularParallelepiped(-10.0, 10.0, -10.0, 10.0, -1.0, 1.0)
        detector_surface = openmc.Sphere(r=20.0); outer_boundary = openmc.Sphere(r=40.0, boundary_type='vacuum')
        target_cell = openmc.Cell(name='Target_Sheet', region=-target_region)
        detector_cell = openmc.Cell(name='Detector_Volume', fill=air, region=+target_region & -detector_surface)
        world_cell = openmc.Cell(name='World', fill=air, region=+detector_surface & -outer_boundary)
        geom = openmc.Geometry([target_cell, detector_cell, world_cell])
        return geom, target_cell, detector_cell, settings, detector_surface

# ==========================================
# 4. TALLIES
# ==========================================
def build_tallies(detector_cell, use_layers=False, use_delayed_low=False):
    energies = np.logspace(np.log10(CONFIG['energy_min_eV']), np.log10(CONFIG['energy_max_eV']), CONFIG['energy_bins'])
    neutron_filter = openmc.ParticleFilter('neutron'); energy_filter = openmc.EnergyFilter(energies)
    cell_filter = openmc.CellFilter([detector_cell])

    tally_n = openmc.Tally(name='scanner_signal'); tally_n.filters = [cell_filter, neutron_filter, energy_filter]; tally_n.scores = ['flux']

    gamma_energies = np.linspace(1e6, 12e6, 2000)
    photon_filter = openmc.ParticleFilter('photon'); gamma_energy_filter = openmc.EnergyFilter(gamma_energies)
    time_filter = openmc.TimeFilter([t * 1e-9 for t in CONFIG['time_bins_ns']]) if CONFIG.get('time_bins_ns') else None

    tally_g = openmc.Tally(name='pgaa_signal')
    tally_g.filters = [cell_filter, photon_filter, gamma_energy_filter] + ([time_filter] if time_filter else []); tally_g.scores = ['flux']

    rxn_tally = openmc.Tally(name='n16_production'); rxn_tally.scores = ['(n,p)']; rxn_tally.nuclides = ['O16']

    mesh = openmc.RegularMesh(); mesh.dimension = [40, 40, 20]; mesh.lower_left = [-25, -25, -5]; mesh.upper_right = [25, 25, 5]
    tally_mesh = openmc.Tally(name='3d_mesh_flux'); tally_mesh.filters = [openmc.MeshFilter(mesh), openmc.ParticleFilter('photon')]; tally_mesh.scores = ['flux']

    tallies_list = [tally_n, tally_g, rxn_tally, tally_mesh]
    neutron_edges = photon_edges = low_energies = None

    if use_layers:
        neutron_edges = np.arange(1, 15, 1) * 1e6
        photon_edges = np.arange(0, 14.1e6, 100e3)
        layered_tally = openmc.Tally(name='layered_flux_5d')
        layered_tally.filters = [openmc.MeshFilter(mesh), openmc.ParticleFilter('photon'), openmc.EnergyFilter(photon_edges), openmc.TimeFilter([t * 1e-9 for t in CONFIG['time_bins_ns']])]
        layered_tally.scores = ['flux']; tallies_list.append(layered_tally)

    if use_delayed_low:
        print("[DELAYED LOW] Adding 1–30 keV photon tally for activation X-rays")
        low_energies = np.linspace(1e3, 30e3, 1000)
        low_e_filter = openmc.EnergyFilter(low_energies)
        tally_low_delayed = openmc.Tally(name='delayed_low_energy_photons')
        tally_low_delayed.filters = [low_e_filter, openmc.ParticleFilter('photon'), openmc.TimeFilter([0.0, 1e9])]
        tally_low_delayed.scores = ['flux']; tallies_list.append(tally_low_delayed)

    return openmc.Tallies(tallies_list), energies, gamma_energies, neutron_edges, photon_edges, low_energies

# ==========================================
# 5. DELAYED LOW-ENERGY ANALYSIS
# ==========================================
def delayed_low_energy_analysis(mean_low_spectrum, energies_low, threshold_sigma=3.0):
    roi_net = {}
    background_mean = np.mean(mean_low_spectrum)
    background_std = np.std(mean_low_spectrum)

    for name, e_val in delayed_low_rois.items():
        if isinstance(e_val, tuple):
            mask = (energies_low >= e_val[0]) & (energies_low <= e_val[1])
        else:
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
# 6. BAYESIAN NNFL LAYER
# ==========================================
def bayesian_nnfl_analysis(x_nnls, Matrix_A, labels, prior_alpha=1.0):
    """Dirichlet posterior update for compositional attribution (NNFL-style)."""
    # 1. Normalize x_nnls to 1.0 so the confidence scalar always works mathematically
    x_norm = x_nnls / (np.sum(x_nnls) + 1e-12) if np.sum(x_nnls) > 0 else np.zeros_like(x_nnls)
    
    alpha_prior = np.full_like(x_nnls, prior_alpha)
    # 2. Apply the scalar to the normalized array
    alpha_post = alpha_prior + (x_norm * 1000.0)   
    post = dirichlet(alpha_post)
    mean_post = post.mean()
    
    samples = post.rvs(size=10000)
    lower = np.percentile(samples, 2.5, axis=0)
    upper = np.percentile(samples, 97.5, axis=0)
    # ... (the rest of the function stays exactly the same)
    # ---------------------------------------------------------------------------
    
    print("\n[NNFL BAYESIAN POSTERIOR]")
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
# 7. NORMALIZATION & RUN SCAN
# ==========================================
def calculate_detector_volume():
    global detector_volume
    print("\n[SYSTEM] Pre-calculating stochastic volume for physics normalization...")

    vol_calc = openmc.VolumeCalculation([detector_cell], samples=100000, lower_left=[-80,-80,-80], upper_right=[80,80,80])
    settings_vol = openmc.Settings(); settings_vol.volume_calculations = [vol_calc]; settings_vol.run_mode = 'volume'
        
    geom.export_to_xml(); settings_vol.export_to_xml()
    openmc.Materials(geom.get_all_materials().values()).export_to_xml()   
    openmc.Tallies().export_to_xml(); openmc.calculate_volumes(output=False)
    
    try:
        vol_stat = openmc.VolumeCalculation.from_hdf5('volume_1.h5')
        detector_volume = vol_stat.volumes[detector_cell.id].n
        print(f"[INFO] 📐 Detector Volume successfully mapped: {detector_volume:.2f} cm³")
    except Exception as e:
        print(f"[WARNING] Volume calculation failed ({e}). Using default 1.0."); detector_volume = 1.0

def run_scan(material: openmc.Material, label: str, use_layers=False, use_delayed_low=False):
    print(f"\n[SCANNING] Interrogating: {label} ...")
    if target_cell is not None:
        target_cell.fill = material; openmc.Materials([material, air]).export_to_xml()
    else:
        all_mats = list(dict.fromkeys([air] + list(csg_materials.values()) + list(material_dict.values())))
        openmc.Materials(all_mats).export_to_xml()
        
    geom.export_to_xml(); settings.export_to_xml(); tallies.export_to_xml()

    t0 = time.time(); openmc.run(output=False); elapsed = time.time() - t0
    sp = openmc.StatePoint(f'statepoint.{settings.batches}.h5')

    mean_n = mean_g = np.array([])
    layered_5d = mean_low = None

    try: t_n = sp.get_tally(name='scanner_signal'); mean_n = t_n.get_values(scores=['flux'], value='mean').flatten() / detector_volume
    except: pass
    try: t_g = sp.get_tally(name='pgaa_signal'); mean_g = t_g.get_values(scores=['flux'], value='mean').flatten() / detector_volume
    except: pass
    
    if use_delayed_low:
        try:
            t_low = sp.get_tally(name='delayed_low_energy_photons')
            mean_low = t_low.get_values(scores=['flux'], value='mean').flatten()
            print(f"    [LOW-E] Low-energy (1-30 keV) spectrum extracted")
        except Exception as e: print(f"[WARNING] Low-energy tally extraction failed: {e}")

    try:
        t_rxn = sp.get_tally(name='n16_production')
        n16_rate = t_rxn.get_values(scores=['(n,p)'], value='mean').sum()
        if n16_rate > 0 and len(mean_g) > 0:
            def inject_gamma(energy_eV, branch_ratio):
                idx = np.abs(gamma_energies[:-1] - energy_eV).argmin()
                injected = (n16_rate * branch_ratio * 0.05) / detector_volume
                n_energy_bins = len(gamma_energies) - 1
                if len(mean_g) > n_energy_bins:
                    mean_g[((len(mean_g) // n_energy_bins - 1) * n_energy_bins) + idx] += injected
                else: mean_g[idx] += injected
            inject_gamma(6.128e6, 0.672); inject_gamma(7.115e6, 0.049); inject_gamma(2.741e6, 0.008)
            print(f"    [N-16] Injected delayed gammas (rate: {n16_rate:.2e})")
    except: pass

    if use_layers:
        try:
            t_layered = sp.get_tally(name='layered_flux_5d')
            layered_5d = t_layered.get_values(scores=['flux'], value='mean').reshape((40, 40, 20, len(CONFIG['time_bins_ns'])-1, len(photon_edges)-1))
            print(f"    [5D] Layered flux extracted → shape {layered_5d.shape}")
            np.savez_compressed(f"5D_flux_full_{label}.npz", flux=layered_5d)
            for i in range(layered_5d.shape[3]): np.savez_compressed(f"5D_flux_neutron_{int(neutron_edges[i]/1e6)}MeV_{label}.npz", flux=layered_5d[:,:,:,i,:])
        except Exception as e: print(f"[WARNING] 5D layered extraction failed: {e}")

    sp.close()
    print(f"[INFO] OpenMC run done in {elapsed:.1f}s")
    return mean_n, mean_g, layered_5d, mean_low

# ==========================================
# 8. SOLVER & ROI BARCODE
# ==========================================
def safe_concat_vectors(*arrays):
    cleaned = [np.zeros(0) if arr is None or len(arr) == 0 else arr.flatten() for arr in arrays]
    return np.concatenate(cleaned) if cleaned else np.zeros(0)

def extract_roi_barcode(mean_g, gamma_energies):
    n_energy_bins = len(gamma_energies) - 1
    if len(mean_g) % n_energy_bins == 0 and len(mean_g) != n_energy_bins:
        n_time_bins = len(mean_g) // n_energy_bins
        mean_g_2d = mean_g.reshape((n_time_bins, n_energy_bins))
    else:
        n_time_bins = 1
        mean_g_2d = mean_g.reshape((1, len(mean_g)))

    e_mids = gamma_energies[:-1]
    masked_spectrum = np.zeros_like(mean_g_2d)

    def apply_compton_subtraction(spectrum_slice, mask_boolean):
        idx = np.where(mask_boolean)[0]
        if len(idx) < 3: return spectrum_slice[mask_boolean]
        baseline = np.linspace(spectrum_slice[idx[0]], spectrum_slice[idx[-1]], len(idx))
        return np.maximum(spectrum_slice[idx] - baseline, 0.0)

    for roi_dict in [fast_rois, thermal_rois] if n_time_bins >= 2 else [{**fast_rois, **thermal_rois}]:
        window_idx = 0 if roi_dict is fast_rois else -1   
        for name, energy in roi_dict.items():
            mask = (e_mids >= energy * 0.98) & (e_mids <= energy * 1.02)
            masked_spectrum[window_idx, mask] = apply_compton_subtraction(mean_g_2d[window_idx, :], mask)
    return masked_spectrum.flatten()

def build_transformation_matrix(material_dictionary, save_prefix="Matrix_A", mode="macro"):
    print(f"\n============================================================\n 🧠 INITIALIZING MATRIX A: {mode.upper()} BASIS GENERATION\n============================================================")
    matrix_columns, material_names = [], []
    iterator = tqdm(list(material_dictionary.items()), desc="Building basis vectors") if TQDM else material_dictionary.items()

    for name, mat in iterator:
        if name.startswith("Mystery_"): continue 
        mean_n, mean_g, _, _ = run_scan(mat, f"Library Generation: {name}") 
        if mode == 'macro':
            combined = safe_concat_vectors(mean_n, mean_g)
            norm_vector = combined / np.sum(combined) if np.sum(combined) > 0 else combined.copy()
        else:
            norm_vector = extract_roi_barcode(mean_g, gamma_energies).copy()
        matrix_columns.append(norm_vector)
        material_names.append(name)

    max_len = max(v.shape[0] for v in matrix_columns)
    padded = [np.concatenate((v, np.zeros(max_len - v.shape[0]))) if v.shape[0] < max_len else v for v in matrix_columns]
    Matrix_A = np.column_stack(padded)
    
    np.save(f"{save_prefix}_Basis.npy", Matrix_A)
    np.save(f"{save_prefix}_Labels.npy", np.array(material_names))
    print(f"\n ✅ MATRIX A COMPILED AND SAVED TO DISK\n    ➤ Dimensions: {Matrix_A.shape[0]} Bins x {Matrix_A.shape[1]} Materials")
    return Matrix_A, material_names

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
        print("\n[DELAYED LOW-ENERGY ANALYSIS — EXPANDED]")
        for name, counts in roi_counts.items():
            print(f" ➤ {name:18}: {counts:.2e} net counts")
        for w in warnings:
            print(f"⚠️ {w}")

    if mode == 'macro':
        b_raw = safe_concat_vectors(mean_n, mean_g)
        b = b_raw / np.sum(b_raw) if np.sum(b_raw) > 0 else b_raw
    else: b = extract_roi_barcode(mean_g, gamma_energies)

    if b.shape[0] < Matrix_A.shape[0]: b = np.concatenate((b, np.zeros(Matrix_A.shape[0] - b.shape[0])))
    elif b.shape[0] > Matrix_A.shape[0]: b = b[:Matrix_A.shape[0]]

    x, residue = opt.nnls(Matrix_A, b)
    threat_detected = False

    # === BAYESIAN NNFL LAYER ===
    if use_bayesian_nnfl:
        mean_post, lower, upper = bayesian_nnfl_analysis(x, Matrix_A, labels)

    if mode == 'micro':
        actual_masses = np.array([x[i] * (800.0 * micro_dict[l].density) for i, l in enumerate(labels)])
        percentages = (actual_masses / np.sum(actual_masses)) * 100.0 if np.sum(actual_masses) > 0 else np.zeros_like(x)
    else: percentages = (x / np.sum(x)) * 100.0 if np.sum(x) > 0 else np.zeros_like(x)

    print("\n 📊 DETECTED MATERIAL COMPOSITION:\n------------------------------------------------------------")
    percent_dict = {}
    for i, label in enumerate(labels):
        pct = percentages[i]
        percent_dict[label] = pct
        if pct > 0.1: print(f"    ➤ {label}: {pct:.1f}%")
        if mode == 'macro' and label in ['TNT', 'RDX', 'PETN', 'TATP', 'Threat_NQ'] and pct > 3.0: 
            threat_detected = True

    if mode == 'micro':
        n_pct = percent_dict.get('Pure_N', 0.0); o_pct = percent_dict.get('Pure_O', 0.0)
        # 1. Standard Nitrogen Explosive Trigger
    if n_pct > 0.5 or (n_pct > 0.2 and o_pct > 1.0):
        threat_detected = True
        print("\n⚠️ STOICHIOMETRIC ANOMALY: Nitrogen density exceeds background clutter limits!")

        # 2. NEW: The Ghost Explosive Trigger (TATP/Peroxides)
        # b[0] = Carbon 4.44 MeV ROI | b[1] = Oxygen 6.13 MeV ROI
        if b[0] > 0:
            o_c_ratio = b[1] / b[0]
            # If the Oxygen/Carbon ratio is unusually high, and there is NO Nitrogen...
            if o_c_ratio > 0.8 and n_pct < 0.2:
                threat_detected = True
                print(f"\n👻 GHOST ANOMALY: Abnormal Fast O/C ratio ({o_c_ratio:.2f}) with ZERO Nitrogen! Peroxide threat (TATP) det>

    print(f"------------------------------------------------------------\n 📉 Mathematical Residual (Unexplained Noise): {residue:.6f}")
    print("\n🔴 STATUS: RED ALARM - EXPLOSIVE SIGNATURE UNMIXED" if threat_detected else "\n🟢 STATUS: GREEN - No threat detected")
    print("============================================================")

    plot_pftna_spectra(mean_g, gamma_energies, mystery_target_name, mean_low=mean_low if use_delayed_low else None)
    if use_layers and layered_5d is not None: plot_layered_flux(layered_5d, neutron_edges, photon_edges, mystery_target_name)

# ==========================================
# 9. VISUALIZATION
# ==========================================
def plot_pftna_spectra(mean_g, gamma_energies, mystery_target_name, mean_low=None):
    print("\n[PLOT] Generating Dual-Window Spectroscopic Dashboard...")
    e_mids_MeV = gamma_energies[:-1] / 1e6
    n_energy_bins = len(e_mids_MeV)
    
    if detector_volume > 1.0: mean_g = mean_g * detector_volume
    
    if len(mean_g) % n_energy_bins == 0 and len(mean_g) != n_energy_bins:
        n_time_bins = len(mean_g) // n_energy_bins
        mean_g_2d = mean_g.reshape((n_time_bins, n_energy_bins))
        fast_spectrum   = gaussian_filter1d(mean_g_2d[0, :], sigma=1.5)
        thermal_spectrum = gaussian_filter1d(mean_g_2d[-1, :], sigma=4.0) if n_time_bins > 1 else np.zeros_like(fast_spectrum)
    else:
        fast_spectrum = gaussian_filter1d(mean_g, sigma=1.5)
        thermal_spectrum = np.zeros_like(fast_spectrum)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(f"PFTNA Scan Diagnostics: {mystery_target_name}", fontsize=18, fontweight='bold')

    ax1.plot(e_mids_MeV, fast_spectrum, color='crimson', lw=1.5, label='Fast Window (0 - 100ns)')
    ax1.fill_between(e_mids_MeV, fast_spectrum, color='crimson', alpha=0.15)
    ax1.set_title("Fast Inelastic Scattering Window", fontsize=14); ax1.grid(True, linestyle='--', alpha=0.5)
    for name, eV in fast_rois.items():
        MeV = eV / 1e6
        ax1.axvspan(MeV * 0.98, MeV * 1.02, color='gold', alpha=0.3); ax1.axvline(MeV, color='goldenrod', linestyle='--', lw=1)
        ax1.text(MeV, max(fast_spectrum)*0.8, name.replace('_', '\n'), rotation=90, va='top', ha='center', fontsize=9, color='darkred')

    ax2.plot(e_mids_MeV, thermal_spectrum, color='dodgerblue', lw=1.5, label='Thermal Window (100ns - 1ms)')
    ax2.fill_between(e_mids_MeV, thermal_spectrum, color='dodgerblue', alpha=0.15)
    ax2.set_title("Thermal Capture Window", fontsize=14); ax2.grid(True, linestyle='--', alpha=0.5)
    for name, eV in thermal_rois.items():
        MeV = eV / 1e6
        ax2.axvspan(MeV * 0.98, MeV * 1.02, color='limegreen', alpha=0.3); ax2.axvline(MeV, color='forestgreen', linestyle='--', lw=1)
        y_pos = max(thermal_spectrum)*0.8 if max(thermal_spectrum) > 0 else 1
        ax2.text(MeV, y_pos, name.replace('_', '\n'), rotation=90, va='top', ha='center', fontsize=9, color='midnightblue')

    if mean_low is not None:
        energies_low = low_energies[:-1]
        ax_inset = fig.add_axes([0.62, 0.62, 0.32, 0.32])   
        ax_inset.semilogy(energies_low / 1e3, mean_low + 1e-10, 'b-', lw=1.2)
        ax_inset.set_xlim(1, 30)
        ax_inset.set_title('Delayed 1–30 keV X-rays', fontsize=9)
        ax_inset.tick_params(labelsize=7)
        ax_inset.axvspan(5.0, 8.0, alpha=0.18, color='orange', label='Fe K X-rays')
        ax_inset.axvspan(7.8, 8.3, alpha=0.18, color='red', label='Cu K X-rays')
        ax_inset.axvspan(8.4, 8.9, alpha=0.18, color='purple', label='Zn K X-rays')
        ax_inset.axvspan(9.1, 9.6, alpha=0.18, color='cyan', label='Ga K X-rays')
        ax_inset.legend(fontsize=6, loc='upper right')

    plt.tight_layout()
    plt.savefig(f"Scan_Report_{mystery_target_name}.png", dpi=300)
    plt.close()
    print(f"[PLOT] Spectroscopic diagnostic saved to: Scan_Report_{mystery_target_name}.png")

def plot_layered_flux(layered_5d, neutron_edges, photon_edges, mystery_target_name):
    if layered_5d is None: return
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"\n[PLOT] Generating 14× per-MeV PNGs + timestamped GIF...")

    mid_z = layered_5d.shape[2] // 2
    n_groups = layered_5d.shape[3]; neutron_edges_MeV = neutron_edges / 1e6

    for i in range(n_groups):
        flux_slice = layered_5d[:, :, mid_z, i, :].sum(axis=2)
        flux_slice = gaussian_filter1d(gaussian_filter1d(flux_slice, sigma=1, axis=0), sigma=1, axis=1)
        plt.figure(figsize=(9, 7)); plt.imshow(flux_slice.T, origin='lower', cmap='viridis', aspect='auto')
        plt.title(f"Neutron {neutron_edges_MeV[i]:.0f} MeV → Total Photon Flux (mid-z slice)")
        plt.xlabel("X mesh"); plt.ylabel("Y mesh"); plt.colorbar(label="Photon flux (arb. units)")
        plt.savefig(f"Layered_Neutron_{neutron_edges_MeV[i]:.0f}MeV_{mystery_target_name}.png", dpi=300); plt.close()

    frames = []; fig, ax = plt.subplots(figsize=(10, 8))
    for i in range(n_groups):
        flux_slice = layered_5d[:, :, mid_z, i, :].sum(axis=2)
        im = ax.imshow(flux_slice.T, origin='lower', cmap='viridis', animated=True)
        ax.set_title(f"Neutron energy bin: {neutron_edges_MeV[i]:.0f} MeV"); frames.append([im])

    total_flux = layered_5d[:, :, mid_z, :, :].sum(axis=(2, 3))
    im = ax.imshow(total_flux.T, origin='lower', cmap='viridis', animated=True)
    ax.set_title("TOTAL photon flux (all neutron energies)"); frames.append([im])

    ani = ArtistAnimation(fig, frames, interval=700, blit=True)
    gif_path = f"neutron_energy_sweep_{timestamp}_{mystery_target_name}.gif"
    ani.save(gif_path, writer='pillow', fps=2); plt.close()
    print(f"[PLOT] GIF saved → {gif_path}  (14 frames + total-flux frame)")

def plot_3d_beam_geometry(beam_mode="fan", n_particles=3000):
    print(f"\n[PLOT] Generating 3D Geometry Diagnostic for '{beam_mode.upper()}' beam...")
    if beam_mode == "cone": mu, phi, color = np.random.uniform(0.95, 1.0, n_particles), np.random.uniform(0, 2*np.pi, n_particles), "crimson"
    elif beam_mode == "fan": mu, phi, color = np.random.uniform(0.90, 1.0, n_particles), np.random.uniform(np.pi/2 - 0.05, np.pi/2 + 0.05, n_particles), "darkorange"
    else: mu, phi, color = np.random.uniform(-1.0, 1.0, n_particles), np.random.uniform(0, 2*np.pi, n_particles), "purple"

    r = np.random.uniform(0, 40, n_particles); theta = np.arccos(mu)
    x = r * np.sin(theta) * np.cos(phi); y = r * np.sin(theta) * np.sin(phi); z = -30.0 + r * np.cos(theta)

    fig = plt.figure(figsize=(10, 8)); ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, color=color, s=2, alpha=0.15); ax.scatter([0], [0], [-30], color='black', s=150, marker='*')

    tx = [-25, 25, 25, -25, -25, 25, 25, -25]; ty = [-25, -25, 25, 25, -25, -25, 25, 25]; tz = [-5, -5, -5, -5, 5, 5, 5, 5]
    faces = [[0,1,2,3], [4,5,6,7], [0,1,5,4], [2,3,7,6], [1,2,6,5], [4,7,3,0]]
    poly3d = [[list(zip(tx, ty, tz))[ix] for ix in face] for face in faces]
    ax.add_collection3d(Poly3DCollection(poly3d, facecolors='cyan', linewidths=1, edgecolors='blue', alpha=0.2))

    ax.set_xlim([-30, 30]); ax.set_ylim([-30, 30]); ax.set_zlim([-35, 15]); ax.view_init(elev=15, azim=-60)
    plt.savefig(f"Geometry_Check_{beam_mode}.png", dpi=300); plt.close()
    print(f"[PLOT] Saved 3D geometry check to: Geometry_Check_{beam_mode}.png")

# ==========================================
# 10. CLI
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description="🧠 Forensic Scanner V1.8 Alpha\nPulsed D-T + N-16 + 5D Layered + Delayed Low-Energy X-rays + Bayesian NNFL",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--build-matrix', action='store_true', help='Build transformation Matrix A.')
    parser.add_argument('--solve', type=str, default=None, help='Solve for target (e.g. Threat_NQ)')
    parser.add_argument('--check-geometry', type=str, choices=['cone', 'fan', 'isotropic'], help='Plot 3D beam.')
    parser.add_argument('--mode', type=str, choices=['macro', 'micro'], default='macro')
    parser.add_argument('--use-csg', action='store_true', help='Full suitcase geometry.')
    parser.add_argument('--use-layers', action='store_true', help='Enable 5D layered tallies + per-MeV PNGs + GIF.')
    parser.add_argument('--use-delayed-low', action='store_true', help='Enable 1–30 keV delayed X-ray tally + Activation monitoring')
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
