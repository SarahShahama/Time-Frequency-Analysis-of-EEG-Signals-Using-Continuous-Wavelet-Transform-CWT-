# Time-Frequency-Analysis-of-EEG-Signals-Using-Continuous-Wavelet-Transform-CWT-
# Project Overview
This project explores time–frequency dynamics of EEG signals from the motor cortex using wavelet transform. The main goal is to analyze how neural oscillations in the C3 (left motor cortex) and C4 (right motor cortex) channels change during rest (baseline) and movement conditions.

The study computes:

Band Power (delta, theta, alpha, beta)

Wavelet Coherence between C3 and C4

Phase Locking Value (PLV) across time windows

Statistical tests to compare baseline vs movement

These analyses help in understanding event-related desynchronization (ERD) and functional connectivity of the motor cortex.

Dataset

The data is from the PhysioNet EEG Motor Movement/Imagery Dataset:

Source:https://physionet.org/content/eegmmidb/1.0.0/

Subject: S001

Files used:

S001R01.edf – baseline (resting, eyes open)

S001R03.edf – movement task (imagined or actual hand movement)

Channels analyzed:

C3 (CH1_IDX = 9) – Left motor cortex

C4 (CH2_IDX = 13) – Right motor cortex

Note: EDF files are downloaded automatically using the provided R code.

Methodology

Data Loading:

EDF files are read using edfReader in R.

Signals are extracted for C3 and C4 and preprocessed (mean removal).

Time–Frequency Analysis:

Continuous Wavelet Transform (CWT) is applied using WaveletComp to visualize power fluctuations over time.

Band Power Calculation:

Mean-square power is computed for delta, theta, alpha, beta using windowed RMS.

Sliding windows: 2-second length, 1-second step.

Functional Connectivity:

Wavelet Coherence: Measures synchronization between C3 and C4 over time and frequency.

Phase Locking Value (PLV): Captures phase synchrony between channels using Hilbert transform.

Statistical Analysis:

Paired t-test and permutation tests compare PLV values between baseline and movement conditions.

Results

Alpha Band Power: Decreases during movement → event-related desynchronization (ERD).

Wavelet Coherence: Strong at baseline, weaker during movement → functional decoupling of hemispheres.

Phase Locking Value (PLV): Stable in both conditions; movement slightly stabilizes synchrony.

Statistical Tests: No significant difference in PLV (t-test p = 0.356, permutation test p = 0.373).

Inference: Movement triggers local desynchronization while maintaining stable inter-hemispheric coordination.

Installation & Usage

Clone the repository:

git clone https://github.com/<your-username>/EEG-Wavelet-Analysis.git
cd EEG-Wavelet-Analysis


Install required R packages:

install.packages(c("edfReader", "WaveletComp", "signal", "ggplot2", "tidyr"))


Run the analysis:

Open EEG_analysis.R in RStudio.

The script will automatically download the dataset, preprocess signals, and compute:

CWT scalograms

Band power

Wavelet coherence

Phase Locking Value

Results are saved as eeg_analysis_results_S001.rds.

Plotting:

Band power and PLV plots are generated using ggplot2.

Coherence and scalograms are plotted with WaveletComp.

Learning Outcomes

Through this project, I learned:

EEG signal processing and time–frequency analysis

Implementation of band power, coherence, and PLV

Application of wavelet transform for neural oscillations

Handling real EEG data in R and troubleshooting complex analyses

Interpretation of neural dynamics in motor cortex during rest and movement

Challenges Faced

Steep learning curve for EEG, frequency bands, and connectivity measures

Debugging Hilbert transform for PLV calculation

Interpreting time–frequency plots and coherence maps

References

Goldberger AL et al., “PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals,” Circulation, 2000.

Torrence C, Compo GP, “A Practical Guide to Wavelet Analysis,” Bulletin of the American Meteorological Society, 1998.

Cohen MX, Analyzing Neural Time Series Data: Theory and Practice, MIT Press, 2014.
