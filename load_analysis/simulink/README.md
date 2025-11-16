# Part II: Load Analysis - Simulink Models

This directory contains Simulink models for analyzing controlled rectifiers under various load conditions.

## Directory Structure

```
simulink/
├── half_wave_rectifier/     # Half-wave controlled rectifier models
├── full_wave_ct/            # Full-wave center-tapped models
├── full_wave_bridge/        # Full-wave bridge models
└── README.md               # This file
```

## Overview

Part II focuses on simulating and analyzing controlled rectifiers with three different load types:
1. **Resistive Load (R)** - Simple resistive loads like heaters
2. **RL Load** - Resistive-inductive loads like DC motors
3. **Highly Inductive Load** - Loads with very high inductance (L >> R)

## Supply Specifications

All simulations should use:
- **Voltage:** 230 V RMS, single-phase
- **Frequency:** 50 Hz
- **Waveform:** Sinusoidal

## Load Specifications

### Suggested Load Parameters

| Load Type | Resistance (Ω) | Inductance (H) | Time Constant L/R (ms) |
|-----------|----------------|----------------|------------------------|
| Resistive | 10            | 0              | 0                      |
| RL Load   | 10            | 0.05           | 5                      |
| Highly Inductive | 5     | 0.5            | 100                    |

*Note: Adjust these values based on desired behavior and convergence*

## Required Simulink Models

### 1. Half-Wave Controlled Rectifier

**Models to Build:**
- `hw_resistive_load.slx` - Half-wave with resistive load
- `hw_RL_load.slx` - Half-wave with RL load
- `hw_highly_inductive.slx` - Half-wave with highly inductive load
- `hw_with_FWD.slx` - Half-wave with free-wheeling diode (RL load)

**Key Components:**
- AC Voltage Source (Simscape > Electrical > Sources)
- Thyristor/SCR (Simscape > Electrical > Semiconductors & Converters)
- Resistor, Inductor (Simscape > Electrical > Passive)
- Pulse Generator for firing circuit
- Voltage and Current Sensors
- Scopes for visualization
- To Workspace blocks for data export

**Measurements Needed:**
- Load voltage waveform
- Load current waveform
- Thyristor voltage
- Thyristor current
- Average output voltage
- RMS output voltage
- Average current
- RMS current

### 2. Full-Wave Center-Tapped Rectifier

**Models to Build:**
- `fw_ct_resistive_load.slx`
- `fw_ct_RL_load.slx`
- `fw_ct_highly_inductive.slx`
- `fw_ct_with_FWD.slx`

**Key Components:**
- Same as half-wave, plus:
- Center-tapped transformer (or two separate voltage sources)
- Two thyristors
- Two pulse generators (180° phase shift)

**Implementation Notes:**
- T1 fires at α (positive half-cycle)
- T2 fires at π + α (negative half-cycle)
- Each SCR conducts for approximately half the cycle

### 3. Full-Wave Bridge Rectifier

**Models to Build:**
- `fw_bridge_resistive_load.slx`
- `fw_bridge_RL_load.slx`
- `fw_bridge_highly_inductive.slx`
- `fw_bridge_with_FWD.slx`

**Key Components:**
- Four thyristors in bridge configuration
- T1-T3 fire together
- T2-T4 fire together
- 180° phase shift between pairs

## Firing Circuit Design

### Requirements
1. **Synchronization:** Must be synchronized with AC source zero-crossing
2. **Adjustable Firing Angle:** 0° to 180°
3. **Pulse Width:** Wide enough for thyristor to latch (1-10 ms typical)
4. **Pulse Amplitude:** Sufficient gate current (check thyristor datasheet)

### Implementation Approaches

#### Method 1: Manual Pulse Generator
```
- Use standard Pulse Generator block
- Set Period = 1/f (0.02 sec for 50 Hz)
- Set Phase Delay = α / (360 × f) seconds
- Set Pulse Width appropriate for load
```

#### Method 2: Zero-Crossing Detection + Delay
```
- Detect zero-crossings of AC source
- Add delay equal to α / (360 × f)
- Generate pulse of appropriate width
```

#### Method 3: Comparator-Based
```
- Generate reference ramp synchronized with AC
- Compare with DC level representing desired α
- Output pulse when ramp exceeds threshold
```

## Simulation Settings

### Solver Configuration
- **Solver Type:** Variable-step (ode23tb recommended for stiff systems)
- **Max Step Size:** Auto or 1e-5
- **Relative Tolerance:** 1e-3
- **Absolute Tolerance:** 1e-6
- **Simulation Time:** At least 5 cycles (0.1 sec for 50 Hz)

### Initial Conditions
- Set initial capacitor voltages and inductor currents to zero (unless studying transient startup)
- Allow at least 2-3 cycles for steady-state

## Measurements and Analysis

### For Each Configuration and Load Type

**Waveforms to Capture:**
1. Input voltage (AC source)
2. Load voltage
3. Load current
4. Thyristor voltage
5. Thyristor current

**Metrics to Calculate:**
- Average output voltage (Vdc)
- RMS output voltage (Vrms)
- Average output current (Idc)
- RMS output current (Irms)
- Form Factor = Vrms / Vdc
- Ripple Factor = sqrt(FF² - 1)
- Extinction angle β (for discontinuous conduction)
- Conduction angle = β - α

**Firing Angles to Test:**
- 0°, 30°, 60°, 90°, 120°, 150°

## Key Observations to Document

### Resistive Load
- Current and voltage in phase
- Always discontinuous conduction
- Simple behavior, matches analytical predictions

### RL Load
- Current lags voltage
- Extinction angle β > π possible
- May have continuous or discontinuous conduction depending on α and L/R
- Free-wheeling diode eliminates negative voltage

### Highly Inductive Load
- Nearly constant (continuous) current
- Small current ripple
- Voltage can be negative while current is positive (α > 90°)
- Demonstrates four-quadrant operation

## Comparison Tasks

### Half-Wave vs Full-Wave
- Compare average voltage
- Compare ripple frequency and amplitude
- Compare efficiency
- Document advantages/disadvantages

### Center-Tapped vs Bridge
- Same electrical performance (approximately)
- Different component requirements
- Different transformer needs

### With vs Without Free-Wheeling Diode
- Effect on output voltage waveform
- Impact on average voltage
- Elimination of negative voltage excursions
- Loss of regenerative capability

## Tips for Successful Simulation

1. **Start Simple**
   - Begin with resistive load and α = 0°
   - Gradually add complexity (inductance, firing angle)
   
2. **Verify Synchronization**
   - Ensure firing pulses align with expected phase
   - Check zero-crossing detection is accurate
   
3. **Check Steady-State**
   - Allow sufficient simulation time
   - Look for periodic steady-state
   
4. **Validate Results**
   - Compare average voltage with theoretical equations
   - Check physical plausibility (e.g., positive power for α < 90°)
   
5. **Troubleshooting**
   - If simulation is slow: increase max step size
   - If results are inaccurate: decrease max step size
   - If convergence issues: try different solver (ode15s, ode23t)
   - If thyristor won't turn on: check gate pulse amplitude and width
   
6. **Data Export**
   - Use To Workspace blocks to export data
   - Process in MATLAB for detailed analysis
   - Generate publication-quality plots

## MATLAB Post-Processing

Create companion MATLAB scripts in `../matlab/` to:
- Import data from Simulink
- Calculate average and RMS values
- Generate comparison plots
- Validate against analytical calculations
- Create tables for report

## Expected Deliverables

For each configuration:
1. Simulink model (.slx file)
2. Screenshots of circuit diagrams
3. Waveform plots for different α and load types
4. Performance metric tables
5. Comparison plots

## Due Date
Part II is due on **December 6, 2025**

## Getting Started Checklist

- [ ] Open Simulink: `simulink` in MATLAB command window
- [ ] Create new model: File > New > Model
- [ ] Add Simscape Electrical blocks from library browser
- [ ] Build half-wave resistive load model first (simplest)
- [ ] Test with α = 0° (uncontrolled rectifier)
- [ ] Add firing circuit with variable α
- [ ] Verify against analytical calculations
- [ ] Expand to other load types
- [ ] Build full-wave configurations
- [ ] Perform systematic comparison

## Resources

- MATLAB/Simulink Documentation: [mathworks.com/help/simulink](https://www.mathworks.com/help/simulink/)
- Simscape Electrical Documentation: [mathworks.com/help/sps](https://www.mathworks.com/help/sps/)
- Thyristor Examples: Search MATLAB File Exchange for "thyristor rectifier"

## Questions?

Document any questions or issues in the project report. Include troubleshooting steps attempted and how problems were resolved.
