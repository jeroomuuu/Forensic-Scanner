# Forensic Scanner V1.4 (tempu): Monte Carlo PFTNA Engine

<img width="1000" height="auto" alt="BigIntro" src="https://github.com/jeroomuuu/Forensic-Scanner/blob/main/Scan_Report_Visual.png" />

A high-performance, dual-mode (Macro/Micro) Pulsed Fast Thermal Neutron Analysis (PFTNA) simulation engine. 

This system uses OpenMC to fire 14 MeV neutrons into complex 3D CSG geometries (like shielded luggage) 
tracks subatomic particle collisions, and uses a Non-Negative Least Squares (NNLS) solver to mathematically 
unmix and identify hidden explosive threats from the resulting gamma-ray spectra.


## Installation & Environment Setup

This project requires a specific scientific Python stack and the official OpenMC physics engine. 
It is highly recommended to deploy this on a Linux machine (Ubuntu) using Conda.


### Step 1.1: Install Miniconda

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b

### Step 1.2: Build the enviroment

    source ~/miniconda3/bin/activate
    conda create -n physics_wsl python=3.11 -y
    conda activate physics_wsl
    conda install -c conda-forge openmc pandas scipy matplotlib tqdm -y

### Step 1.3: Get the Nuclear Data

Get the OpenMC endf-b-viii-1 pack & optional depletion data
look for the cross_sections.xml location after extracting
    
    wget https://anl.box.com/shared/static/6qr7jezzihkj9p9esl5jn19qgpujyjyz.xz
    tar -xJf 6qr7jezzihkj9p9esl5jn19qgpujyjyz.xz

might as well get the depletion data - you might need it later

    wget -O thermal.xml https://anl.box.com/shared/static/q6ev8pl7xct179ke7kq148smde8gzni6.xml
    wget -O fast.xml https://anl.box.com/shared/static/n0pkqe66uotskoljr93szvjyvtvycgze.xml

###### if you cannot download from the exact link, go find it at https://openmc.org/data/


### Step 1.4: Tell the enviroment where to find it & tell your startup bash as well

    export OPENMC_CROSS_SECTIONS=/path-to-file/cross_sections.xml
    echo 'export OPENMC_CROSS_SECTIONS="/path-to-file/cross_sections.xml"' >> ~/.bashrc


## Usage & Startup Syntax

The script operates in two phases: Matrix Generation (building the physics library) and Solving (scanning a mystery target)
You can run these in either macro (full spectrum) or micro (ROI-Barcoded) modes


### Step 2.1: Build the Transformation Matrix

    python forensic-scanner-1.0.py --build-matrix --mode micro
    python forensic-scanner-1.0.py --build-matrix --mode macro

### Step 2.2: Run a detection scan

    python forensic-scanner-1.0.py --solve Threat_NQ --use-csg --mode micro
    python forensic-scanner-1.0.py --solve Threat_NQ --use-csg --mode macro

### Step 2.3: Options to check beam output visually


    python forensic-scanner-1.0.py --check-geometry fan
    python forensic-scanner-1.0.py --check-geometry cone


## Tweaking Settings

If you need to change the resolution, particle count, ROIs; 
open forensic-scanner-1.0.py and modify the CONFIG dictionary

    'sigma_certainty': 4.0,           # Threshold for triggering Stoichiometric Anomaly warnings
    'particles_per_batch': 500000,    # THE CPU SLIDER: Change to 50,000,000 for Enterprise HPC runs
    'batches': 100,                   # Total batches to run
    'energy_min_eV': 1.0,
    'energy_max_eV': 1.4e7,
    'energy_bins': 2000,              # Neutron tally resolution
    'time_bins_ns': [0, 10, 20, 30, 40, 50, 100, 200, 300, 500, 1000]

## To increase Gamma-Ray Resolution

If you want true per-keV resolution on your final plots and matrices, 
find the gamma_energies variable inside the build_tallies() function,
and increase the resolution bins.

    gamma_energies = np.linspace(1e6, 12e6, 2000)



###### Note: If you change the settings, you must run --build-matrix again so your library matches your new resolution!
