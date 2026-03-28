#### CHANGELOG

V0.01 -> V0.99:

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

V1.6:

- Tuned & expanded low-energy ROIs (x-ray and photon)
- Something something - i forgot

V1.7:

- Added Cu, Zn, Ga X-ray clusters + more Fe-59 lines
- Finetuned a tiny bit more

V1.8:

- Bayesian_nnfl_analysis() function present and called correctly
- Dirichlet posterior + 95% credible intervals printing
- Threat probability warning for >80%
- No breakage to any previous features

V1.9: Problems fixed

- Better math
- Solver now also works with other target compositions, like "TATP_Ghost"
- Bugfixes all the way to 2.0

V2.0.1: - Stable - OUT NOW

- Small bugfix from V2.0
- 5d Mappings of Neutron Flux
- 2d layered neutron png output per MeV channel
- ready for future ml tensors

###### Run the micro and macro builds with high particle / batch counts on a cpu cluster ( 1T+ per vector) , then run the solver on the Tensor driven GPU.
###### Connect the GPU machine via datalink to the larger matrix sets on the cpu cluster. This way the ground truth has almost a net zero error rate to start.


---

### Roadmap

#### [Proposed Updates]

- Stochastic Noise Injection
- Full PyMC MCMC upgrade ->
- Probabilistic inference layer (PyMC, Pyro, TensorFlow Probability, etc.)
- DAGMC CAD import for real suitcase STL files in instead of CSG boxes
- Cluster watchdog script for Vast.ai API (one VM per element, auto-assemble etc etc).

- Models with:
- Collision Tracking list-mode digitizer
- Posterior distributions instead of point estimates
- Credible intervals on every material percentage
- Proper uncertainty propagation
- Threat probability with uncertainty (e.g. “92% ± 4% chance of explosive”)
