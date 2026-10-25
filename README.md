# Differential Burst dynamics of Slow and Fast gamma rhythms in Macaque V1
- Computed burst durations of stimulus-induced slow and fast gamma rhythms in LFP signals recorded from area V1 of two awake adult female monkeys
- Employed MP, OMP-GEAR, and HT to estimate burst durations
- Measured functional connectivity (PLV and WPLI) across all pair of electrodes
- Used WC network model (firing-rate) to simulate bursty slow and fast gamma rhythms
- For further details refer to: https://doi.org/10.1162/IMAG.a.1404
## Results
- Slow gamma rhythm exhibited significantly longer burst durations and longer latencies as compared to fast gamma.
- Slow gamma exhibited higher long-range synchrony compared to fast gamma.
- Results were consistent across different burst estimators.
- Simply changing the firing-rate time-constant of the corresponding inhibitory interneuronal population in the noisy WC model, leads to both slower and longer bursts.
## Requirements
- MATLAB R2024b or higher (MathWorks, RRID: SCR_001622)
- Chronux toolbox (Mitra and Bokil, 2008) (http://chronux.org/, RRID: SCR_005547)
- Fieldtrip Toolbox ((Oostenveld et al., 2011), RRID: SCR_004849)
## Code Structure
- Burst_detection_methods: Codes to estimate burst dynamics of slow and fast gamma, using MP, OMP-GEAR and HT.
- LFP_analysis: Codes for plotting the figures and complete analysis of estimated burst duration and latencies.
- Model_simulation: Used to simulate non-linear self oscillating WC model to generate slow and fast gamma rhythms.
- Model_analysis: Analysis of simulated LFP dataset (from model).
## Data availability
- The entire LFP dataset (decimated version) of the two monkeys- Monkey-1 (alpaH) and Monekey-2 (kesariH) is available in the data folder.
- Download the data and add to your path, and you can either use the pre-computed burst durations (default) or run the algorithms again.
- Modify the parent_file_path variable defined in the codes accordingly.
## References
  - https://github.com/supratimray/GammaLengthProjectCodes.git
  - https://www.fuw.edu.pl/~durka/software/mp/
  - https://github.com/avianand8/GammaLengthProjectCodes.git
  - https://github.com/RKrishnakumaran/NeOLab-WC-models.git

