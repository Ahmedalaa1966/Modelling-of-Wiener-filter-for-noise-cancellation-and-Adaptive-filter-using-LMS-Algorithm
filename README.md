# Modelling of Wiener Filter for Noise Cancellation and Adaptive Filter Using LMS Algorithm

[cite_start]This project implements optimal and adaptive signal processing techniques in MATLAB to isolate signals and predict stochastic processes[cite: 13, 14, 42]. [cite_start]The first part builds an FIR Wiener filter (orders 6 and 12) to successfully perform noise cancellation on a sinusoidal signal corrupted by autoregressive noise[cite: 14, 35]. [cite_start]The second part utilizes the Least Mean Squares (LMS) algorithm to construct a 2-tap adaptive linear predictor that tracks a second-order autoregressive process, evaluating coefficient convergence across different step sizes ($\mu = 0.02$ and $0.04$) and verifying them against analytical stability bounds[cite: 43, 49, 57, 59, 60].

---

## 📁 Repository Structure

### 📁 `Matlab_codes/`
[cite_start]This directory contains the core implementation scripts for both tasks[cite: 70]. The code is structured for modularity and easy execution, handling everything from synthetic signal generation to statistical calculation:
* [cite_start]**Wiener Filter Design:** Scripts compute the empirical autocorrelation matrix and cross-correlation vector to solve the Wiener-Hopf equations directly, applying the calculated optimal weights to filter out the autoregressive noise[cite: 34].
* [cite_start]**LMS Adaptive Learning:** Implements the real-time stochastic gradient descent loop, updating the two predictor taps at every iteration based on the instantaneous error signal[cite: 49, 57].

### 📁 `sim/`
This directory contains the raw simulation environment data and output artifacts generated during execution:
* **Saved Workspaces:** Mat-files containing variables, generated stochastic processes, and calculated filter weights, allowing you to reload the exact simulation states without re-running the entire script.
* [cite_start]**Figure Exports:** High-resolution exports of the time-domain waveforms, error performance trends, and coefficient tracking plots that are featured in the final report[cite: 71].

### 📁 `Report/`
[cite_start]Contains the comprehensive final project report document[cite: 71].

---

## 📄 Project Report Contents

[cite_start]For a deep technical dive, please refer to the PDF file inside the **`Report/`** directory[cite: 71]. The document includes:
* [cite_start]**Mathematical Analysis:** Step-by-step theoretical derivations for the Wiener-Hopf equations, autocorrelation matrices, cross-correlation vectors, and the analytical eigenvalue bounds ($\lambda_{max}$) dictating LMS stability[cite: 60].
* [cite_start]**Simulation Results:** Comprehensive performance graphs comparing the original, corrupted, and filtered waveforms[cite: 36, 71].
* [cite_start]**Convergence Curves:** Detailed plots tracking the real-time adaptation and learning trajectory of the filter coefficients ($w_n(1)$ and $w_n(2)$) under different step sizes ($\mu = 0.02$ and $\mu = 0.04$)[cite: 59].

---

## 🛠️ Requirements & Tools
* [cite_start]MATLAB / GNU Octave [cite: 13]
* Signal Processing Toolbox
