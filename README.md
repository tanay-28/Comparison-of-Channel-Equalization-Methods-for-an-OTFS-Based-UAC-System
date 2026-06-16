# Comparison of Channel Equalization Methods for an OTFS-Based Underwater Acoustic Communication System

MATLAB implementation of an OTFS-based underwater acoustic (UWA) communication system with multiple equalization and detection methods, as presented in:

> **T. Agrawal, R. Kadlimatti, and N. Gupta**, "Comparison of Channel Equalization Methods for an OTFS-Based Underwater Acoustic Communication System," *2025 IEEE Asia Pacific Conference on Wireless and Mobile (APWiMob)*, Bali, Indonesia, Nov. 2025. DOI: 10.1109/APWiMob67231.2025.11269161

---

## Overview

The underwater acoustic (UWA) channel is characterized by large delay spreads and high Doppler shifts, making it doubly dispersive and challenging for conventional OFDM-based systems. OTFS modulation represents the time-varying channel in a sparse two-dimensional delay-Doppler (DD) domain, making it well-suited for UWA environments.

This work implements the statistically equivalent UWA channel model from Qarabaqi and Stojanovic [5], which accounts for both large-scale and small-scale perturbations as well as motion-induced Doppler effects. Four equalization and detection methods of increasing complexity are implemented and compared across three underwater channel scenarios:

- **Single-tap MMSE equalizer** — low complexity, suitable for frequency-flat channels
- **Time-domain LMMSE equalizer** — block-wise equalization in the delay-time domain
- **Maximal Ratio Combining (MRC) detector** — iterative delay-time domain detector
- **Message Passing Algorithm (MPA) detector** — factor graph-based iterative detector

Key finding: the MRC delay-time detector most effectively equalizes channels with large numbers of delay and Doppler coefficients, while the time-domain LMMSE equalizer is preferable in sparser channel conditions due to its lower complexity.

---

## Simulation Scenarios

Three underwater channel scenarios of increasing depth are evaluated:

| Parameter | Scenario 1 | Scenario 2 | Scenario 3 |
|-----------|-----------|-----------|-----------|
| Depth (m) | 50 | 100 | 1000 |
| TX height (m) | 4 | 40 | 40 |
| RX height (m) | 2 | 20 | 20 |
| Channel distance (m) | 100 | 1000 | 1000 |
| Spreading factor | 1.7 | 1.7 | 1.7 |
| Speed of sound in water (m/s) | 1500 | 1500 | 1500 |
| Speed of sound at bottom (m/s) | 1200 | 1200 | 1200 |

---

## System Parameters

| Parameter | Value |
|-----------|-------|
| Symbols per frame (N) | 16 |
| Subcarriers (M) | 64 or 128 |
| Modulation | 4-QAM, 16-QAM |
| Bandwidth (B) | 6400 Hz |
| Minimum frequency | 10 kHz |
| Detector iterations | 15 |
| Total simulation time | 204.8 s (4 × coherence time) |

---

## Channel Model

The statistically equivalent UWA channel model consists of four components:

**1. Large-scale geometry** — nominal path lengths determined by surface height, TX/RX depth, and channel distance via `mpgeometry.m`.

**2. Large-scale uncertainty** — slow AR(1) variations in surface height, TX height, RX height, and channel distance modeled over the total simulation time.

**3. Small-scale scattering** — intrapath fading modeled by a statistically equivalent AR(1) process using the Takagi factorization for efficient computation.

**4. Motion-induced Doppler** — three motion components:
- Drifting (transmitter and receiver)
- Vehicular motion
- Surface wave motion

The discrete delay-time channel model is:

```
r(q) = Σ H_s(l,q) · s(q-l) + z(q)
```

where `1 ≤ l ≤ taps` and `z ~ CN(0, σ²)`.

---

## Repository Structure

```
├── main_script.m                          # Main simulation script
├── Generate_2D_data_grid.m                # DD grid data placement
├── Gen_time_domain_channel_2.m            # Time-domain channel matrix generation
├── Gen_DD_and_DT_channel_matrices.m       # DD and delay-time channel matrices
├── Gen_DT_and_DD_channel_vectors_2.m      # Channel vector generation
├── Generate_time_frequency_channel_ZP_2.m # Time-frequency channel with ZP
├── mpgeometry.m                           # Underwater channel geometry (multipath)
├── absorption.m                           # Acoustic absorption coefficient
├── find_takagi_factor.m                   # Takagi factorization for small-scale model
├── MRC_delay_time_detector_21.m           # MRC iterative delay-time detector
├── MPA_detector_21.m                      # Message passing algorithm detector
├── TF_single_tap_equalizer_21.m           # Single-tap time-frequency equalizer
└── Block_LMMSE_detector_21.m             # Block time-domain LMMSE equalizer
```

---

## Running the Simulation

1. Clone the repository and open MATLAB
2. Navigate to this folder
3. Set the desired scenario parameters at the top of the main script:

```matlab
h0  = 100;   % surface depth [m]
ht0 = 40;    % TX height [m]
hr0 = 20;    % RX height [m]
d0  = 1000;  % channel distance [m]
M   = 64;    % number of subcarriers (64 or 128)
M_mod = 4;   % QAM order (4 or 16)
```

4. Run the main script. Set `plotting_key = 1` to visualize the delay-time channel:

```matlab
plotting_key = 1;
```

The script generates:
- **BER vs SNR** curves comparing all four equalizers
- **Delay-time channel plot** (when `plotting_key = 1`)

---

## Results Summary

**Scenario 1 (50 m depth)** — sparse DD channel. MPA and LMMSE perform best. BER ≈ 10⁻⁴ achievable with 4-QAM.

**Scenario 2 (100 m depth)** — dense DD channel with many delay-Doppler components. MRC delay-time detector outperforms all others. MPA fails to converge in 15 iterations.

**Scenario 3 (1 km depth)** — sparse DD channel similar to Scenario 1. LMMSE and MPA recover performance. MRC effective but over-complex for sparse channels.

General observation: increasing M from 64 to 128 improves BER across all scenarios due to finer Doppler resolution.

---

## Citation

If you use this code, please cite:

```bibtex
@inproceedings{agrawal2025comparison,
  title={Comparison of Channel Equalization Methods for an OTFS-Based Underwater Acoustic Communication System},
  author={Agrawal, Tanay and Kadlimatti, Ravi and Gupta, Naveen},
  booktitle={2025 IEEE Asia Pacific Conference on Wireless and Mobile (APWiMob)},
  pages={193--197},
  year={2025},
  address={Bali, Indonesia},
  doi={10.1109/APWiMob67231.2025.11269161}
}
```

---

## References

[1] P. Qarabaqi and M. Stojanovic, "Statistical Characterization and Computationally Efficient Modeling of a Class of Underwater Acoustic Communication Channels," *IEEE Journal of Oceanic Engineering*, vol. 38, no. 4, pp. 701–717, 2013.

[2] T. Thaj and E. Viterbo, "Low Complexity Iterative Rake Decision Feedback Equalizer for Zero-Padded OTFS Systems," *IEEE Transactions on Vehicular Technology*, vol. 69, no. 12, pp. 15606–15622, 2020.

[3] P. Raviteja, K. T. Phan, Y. Hong and E. Viterbo, "Interference Cancellation and Iterative Detection for Orthogonal Time Frequency Space Modulation," *IEEE Transactions on Wireless Communications*, vol. 17, no. 10, pp. 6501–6515, 2018.

[4] R. Hadani et al., "Orthogonal Time Frequency Space Modulation," *2017 IEEE WCNC*, San Francisco, CA, 2017.
