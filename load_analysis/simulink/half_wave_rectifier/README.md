# Half-Wave Controlled Rectifier Models

This directory contains Simulink models for half-wave controlled rectifiers with various load types.

## Models to Create

### 1. `hw_resistive_load.slx`
**Description:** Half-wave SCR rectifier with purely resistive load

**Circuit Components:**
- AC Voltage Source: 230V RMS, 50 Hz
- Thyristor (SCR)
- Resistor: 10 Ω (suggested)
- Firing pulse generator
- Voltage and current measurement blocks

**Expected Behavior:**
- Current follows voltage when SCR is on
- Discontinuous current (always)
- Current and voltage in phase
- Zero current when voltage is negative or SCR is off

**Test Cases:**
- α = 0°: Maximum output (uncontrolled rectifier)
- α = 30°: Moderate control
- α = 60°: Significant reduction
- α = 90°: Half the maximum output
- α = 120°: Heavily controlled

---

### 2. `hw_RL_load.slx`
**Description:** Half-wave SCR rectifier with resistive-inductive load

**Circuit Components:**
- AC Voltage Source: 230V RMS, 50 Hz
- Thyristor (SCR)
- Resistor: 10 Ω
- Inductor: 50 mH (suggested, adjust as needed)
- Firing pulse generator
- Measurement blocks

**Expected Behavior:**
- Current lags voltage
- Extinction angle β > 180° possible
- May have discontinuous or continuous conduction
- Smoother current waveform than resistive load

**Observations to Make:**
- Determine extinction angle for different α
- Identify boundary between DCM and CCM
- Compare with resistive load case

**Test Cases:**
- α = 0° to 120° in steps of 30°
- Note whether conduction is continuous or discontinuous

---

### 3. `hw_highly_inductive.slx`
**Description:** Half-wave SCR rectifier with highly inductive load

**Circuit Components:**
- AC Voltage Source: 230V RMS, 50 Hz
- Thyristor (SCR)
- Resistor: 5 Ω
- Inductor: 500 mH (L >> R condition)
- Firing pulse generator
- Measurement blocks

**Expected Behavior:**
- Nearly constant current (acts like current source)
- Always continuous conduction
- Small current ripple
- **Important:** Voltage can be negative while current is positive!

**Special Focus:**
- Observe behavior for α > 90°
- Document negative voltage with positive current
- Explain four-quadrant operation (Quadrant II)
- Calculate average voltage (should be negative for α > 90°)

**Test Cases:**
- α = 0°, 30°, 60° (normal rectification)
- α = 90° (boundary case, Vdc ≈ 0)
- α = 120°, 150° (regenerative operation, negative Vdc)

---

### 4. `hw_with_FWD.slx`
**Description:** Half-wave SCR rectifier with free-wheeling diode and RL load

**Circuit Components:**
- AC Voltage Source: 230V RMS, 50 Hz
- Thyristor (SCR)
- Free-wheeling Diode (parallel with load)
- Resistor: 10 Ω
- Inductor: 50 mH
- Firing pulse generator
- Measurement blocks

**Expected Behavior:**
- Output voltage never goes negative (clamped to zero by FWD)
- When SCR turns off, current free-wheels through diode
- Higher average voltage compared to without FWD
- Improved load voltage waveform

**Comparison Task:**
- Create identical model without FWD
- Compare waveforms side-by-side
- Calculate difference in average voltage
- Document advantages and disadvantages

**Test Cases:**
- α = 30°, 60°, 90° with and without FWD
- Measure improvement in Vdc

---

## Firing Circuit Implementation

### Basic Pulse Generator Setup

For half-wave rectifier:
```
Period: 0.02 sec (for 50 Hz)
Amplitude: 5 (gate signal)
Pulse Width: 10% (or 0.002 sec)
Phase Delay: α / (360 × 50) seconds

Examples:
- α = 0°   → Delay = 0 sec
- α = 30°  → Delay = 0.00167 sec
- α = 60°  → Delay = 0.00333 sec
- α = 90°  → Delay = 0.005 sec
```

### Alternative: Variable Firing Angle

Use MATLAB variable for firing angle:
```matlab
alpha_deg = 60;  % in workspace
alpha_delay = alpha_deg / (360 * 50);  % calculate delay
% Reference this in Pulse Generator Phase Delay parameter
```

---

## Measurements Required

### Voltage Measurements
1. **Supply Voltage** - Verify it's sinusoidal
2. **Load Voltage** - Main output
3. **Thyristor Voltage** - Shows blocking and conduction states

### Current Measurements
1. **Load Current** - Main output current
2. **Thyristor Current** - Same as load current when conducting

### Calculated Metrics
Using MATLAB scripts or Simulink measurement blocks:
- Vdc (average output voltage)
- Vrms (RMS output voltage)
- Idc (average load current)
- Irms (RMS load current)
- Form Factor = Vrms / Vdc
- Ripple Factor = sqrt(FF² - 1)

---

## Scope Configuration

### Scope 1: Voltage Waveforms
- Input AC voltage (reference)
- Load voltage
- Thyristor voltage

### Scope 2: Current Waveforms
- Load current
- Thyristor current

### Scope 3: Combined View
- Load voltage and current on same axes
- Helps visualize phase relationship

---

## Data Export

Add "To Workspace" blocks for:
- Load voltage array: `Vo_hw_R` (resistive), `Vo_hw_RL` (RL), etc.
- Load current array: `Io_hw_R`, `Io_hw_RL`, etc.
- Time array: `tout`

This allows post-processing in MATLAB scripts.

---

## Troubleshooting

### SCR Won't Turn On
- Check gate pulse amplitude (should be > 0)
- Check gate pulse timing (should be during positive voltage)
- Verify pulse width is sufficient

### Simulation is Slow
- Increase max step size
- Reduce simulation time (fewer cycles)
- Use ode23tb solver

### Current is Too Large/Small
- Check load resistance value
- Verify voltage source amplitude
- Check units (V, A, Ω)

### Unexpected Waveforms
- Allow more cycles for steady-state
- Check initial conditions
- Verify component connections

---

## Report Requirements

For each model, document:
1. Circuit diagram (screenshot from Simulink)
2. Component values used
3. Waveforms for multiple firing angles
4. Performance metrics table
5. Observations and analysis

---

## Tips

1. **Build incrementally** - Start with α = 0° and verify, then add control
2. **Compare with theory** - Check that Vdc matches analytical calculation
3. **Use consistent colors** - Same signals same colors across all scopes
4. **Label everything** - Use meaningful signal names and scope titles
5. **Save variations** - Save different α values as separate model versions or use parameter sweep

---

## Next Steps

After completing half-wave models:
1. Validate results against Part I MATLAB functions (if applicable)
2. Create comparison plots
3. Move to full-wave configurations
4. Prepare findings for report
