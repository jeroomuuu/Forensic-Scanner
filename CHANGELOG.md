#### Current Versions

V0.0 -> V0.9:

- Moved to private repo

V1.0:

- First stable build and working wiki instructions

V1.1:

- Tried some stuff -> Didn't work that well
- ( > 20MeV beams and such )

V1.2:

- Added missing delays and better image generation + missing yellow bands

V1.3:

- Delayed gamma now works - doesnt need much
- 14 crisp per-MeV PNGs (no more 6–7 MeV clumping)
- 14 neutron energy groups (1 MeV bins 1–14 MeV)
- 141 photon bins (100 keV resolution 0–14 MeV)
- 40×40×20 spatial mesh
- One beautiful timestamped GIF sweep + total-flux frame
- Full 5D .npz files saved automatically
- Time bins act as neutron-energy proxy (thanks to pulsed source)
- Timestamped output

V1.4:

- real cascade gamma from thermal neutron capture on nitrogen-14¹⁴N(n,γ)¹⁵N
- Several lower-energy cascade gammas from the de-excitation of higher excited states in ¹⁵N
- Threat_NQ is nitrogen-rich (CH₄N₄O₂) and the borated-poly C/M now stands out clearly in the delayed/thermal bin
- Stuff I forgot about

V1.5:

- Beautiful orange-shaded Fe K X-ray cluster in the new low-energy inset
- Pulsed D-T, 5D layered tallies, GIFs, per-MeV PNGs, N-16 injection
- Delayed low-energy X-ray monitoring
- All global variables and unpacking fixed (no more NameError hell)


#### Roadmap

V1.6: [Proposed Updates]

- Tune & expand low-energy ROIs

V1.7: [Proposed Updates]

- Add Cu, Zn, Ga X-ray clusters + more Fe-59 lines

V1.8: [Proposed Updates]

- Bayesian / NNFL layer, probabilistic attribution instead of just NNLS percentages

V1.9: [Proposed Updates]

- DAGMC CAD import for real suitcase STL files in instead of CSG boxes

V2.0: [Proposed Updates]

- Cluster watchdog script for Vast.ai API (one VM per element, auto-assemble etc etc).

