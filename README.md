# Repo Name

## Project Overview

This repository contains the implementation and analysis of controlled rectifier circuits for battery charging applications and various load conditions. The project is divided into two root parts:

### Part I: Design of Battery Charger using Controlled Rectifiers

Design and simulate battery charger circuits using thyristor-based controlled rectifiers with analytical MATLAB functions.

### Part II: Analysis and Simulation of Controlled Rectifiers for Various Load Conditions

Comprehensive Simulink simulation and analysis of controlled rectifiers under different load types (R, RL, and highly inductive loads).

---

## Repository Structure

```
matlab/
├── battery_charger/                     # Part I: Battery Charger Design
│   └── matlab/                          # MATLAB analytical functions
│       ├── half_wave_charger.m          # Half-wave rectifier implementation
│       ├── full_wave_ct_charger.m       # Center-tapped full-wave rectifier
│       └── full_wave_bridge_charger.m   # Bridge full-wave rectifier
│
├── load_analysis/                       # Part II: Load Analysis
│   ├── simulink/                        # Simulink models
│   │   ├── half_wave_rectifier/         # Half-wave configurations
│   │   ├── full_wave_ct/                # Center-tapped configurations
│   │   └── full_wave_bridge/            # Bridge configurations
│   └── matlab/                          # Analytical calculation scripts
│
├── report/                              # LaTeX report files
│   ├── root.tex                         # main report document
│   └── figures/                         # Generated figures and plots
│
└── docs/                                # Additional documentation
```

---

## Project Objectives

- Simulate half-wave and full-wave controlled rectifiers using thyristors (SCRs)
- Analyze the influence of firing angle on average and RMS output voltage/current
- Compare performance of different rectifier configurations
- Investigate EMF load performance
- Develop reusable MATLAB functions for battery charger design
- Simulate controlled rectifiers with various load types
- Study the effect of different loads (R, RL, highly inductive) on waveforms
- Analyze firing angle influence on output characteristics
- Determine the role of free-wheeling diodes
- Compare half-wave and full-wave rectifier configurations

---

## Requirements

### Software
- **MATLAB** R2020a or later
- **Simulink** with Power Systems Toolbox
- **LaTeX Distribution** (TeX Live, MiKTeX, or MacTeX)
- **TeXstudio** (recommended LaTeX editor)

### MATLAB Toolboxes
- Simulink
- Simscape Electrical (Power Systems)

---

## Getting Started

### 1. Clone the Repository
```bash
git clone <repository-url>
cd matlab
```

### 2. Part I: Battery Charger Design
Navigate to `battery_charger/matlab/` and run the MATLAB functions:

```matlab
% Example: Half-wave rectifier
Vrms = 230;           % Supply voltage (RMS)
f = 50;               % Frequency (Hz)
Vbat = 12;            % Battery voltage (V)
Rbat = 0.1;           % Internal resistance (Ohm)
capacity = 50;        % Battery capacity (Ah)

% Call the function
half_wave_charger(Vrms, f, Vbat, Rbat, capacity);
```

### 3. Part II: Load Analysis
Open Simulink models in `load_analysis/simulink/`:

```matlab
% Open Simulink
simulink
% Navigate to models and run simulations
```

### 4. Generate Report
Navigate to `report/` and compile the LaTeX document:

```bash
cd report
pdflatex root.tex
bibtex root
pdflatex root.tex
pdflatex root.tex
```

Or use TeXstudio to open `root.tex` and compile.

---

## Key Specifications

### Part I: Battery Charger Function Inputs
- Supply voltage (Vrms) and frequency f [Hz]
- Battery voltage [V]
- Battery equivalent internal resistance [Ω]
- Battery capacity (Ah or Wh)
- **Optional:** Charging time [sec] and initial State of Charge (SoC) [%]
- **Optional:** Thyristor parameters (forward drop, leakage, switching times)

### Part I: Function Outputs
- Figure: Firing angle [°] vs. charging time [sec]
- **Optional:** Final State of Charge (SoC) [%]
- **Optional:** Average power losses [W]

### Part II: Simulation Parameters
- Supply: 230 V, 50 Hz single-phase
- Load types: R, RL, highly inductive
- Firing angle: 0° - 180° (adjustable)
- Measurements: Load voltage/current, thyristor voltage/current

---

## References

- Power Electronics Course Materials - Dr.-Ing. Moustafa Adly
- MATLAB/Simulink Documentation
- Thyristor-Based Power Converters Literature

---

## Author
- **Author Name:** Mohammed Azab
- **Email:** mohammed@azab.io

---