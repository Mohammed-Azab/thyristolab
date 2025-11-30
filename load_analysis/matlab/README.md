# Part I: Battery Charger Design - MATLAB Functions

This directory contains MATLAB functions for analyzing controlled rectifiers used in battery charging applications.

## Files

### Main Functions

1. **`half_wave_charger.m`**
   - Implements half-wave controlled rectifier analysis
   - Single thyristor configuration
   - Calculates firing angle vs. charging time relationship
   - Optional SoC tracking and power loss calculation

2. **`full_wave_ct_charger.m`**
   - Implements full-wave center-tapped rectifier analysis
   - Two thyristors with center-tapped transformer
   - Approximately 2x performance of half-wave
   - Same optional features as half-wave

3. **`full_wave_bridge_charger.m`**
   - Implements full-wave bridge rectifier analysis
   - Four thyristors in bridge configuration
   - No center-tap required
   - Similar performance to center-tapped

## Usage Examples

### Basic Usage - Half-Wave Rectifier
```matlab
% Define parameters
Vrms = 230;        % Supply voltage (RMS) [V]
f = 50;            % Frequency [Hz]
Vbat = 12;         % Battery voltage [V]
Rbat = 0.1;        % Internal resistance [Ohm]
capacity = 50;     % Battery capacity [Ah]

% Run analysis
[alpha, t_charge] = half_wave_charger(Vrms, f, Vbat, Rbat, capacity);
```

### With Optional Parameters
```matlab
[alpha, t_charge, SoC_final, P_loss] = half_wave_charger(Vrms, f, Vbat, Rbat, capacity, ...
    't_charge', 3600, ...      % Charging time: 1 hour
    'SoC_init', 20, ...        % Initial SoC: 20%
    'Vt', 1.5, ...             % Thyristor forward drop: 1.5V
    'Ileak', 1e-6);            % Leakage current: 1μA
```

### Comparing All Configurations
```matlab
compare_rectifier_configs(230, 50, 12, 0.1, 50);
```

## Function Inputs

### Required Parameters (all functions)
- `Vrms` - Supply voltage RMS value [V]
- `f` - Supply frequency [Hz]
- `Vbat` - Battery nominal voltage [V]
- `Rbat` - Battery internal resistance [Ω]
- `capacity` - Battery capacity [Ah]

### Optional 1 Parameters (Name-Value pairs)
- `'t_charge'` - Charging time [sec] (if known)
- `'SoC_init'` - Initial State of Charge [%] (default: 20)

### Optional 2 Parameters (Name-Value pairs)
- `'Vt'` - Thyristor forward voltage drop [V] (default: 0)
- `'Ileak'` - Thyristor reverse leakage current [A] (default: 0)
- `'t_rise'` - Voltage/current rise time [sec] (default: 0)
- `'t_fall'` - Voltage/current fall time [sec] (default: 0)

## Function Outputs

### Standard Outputs
- `alpha_deg` - Array of firing angles [degrees]
- `charging_time_hours` - Corresponding charging times [hours]

### Optional 1 Outputs
- `SoC_final` - Final State of Charge [%] (if t_charge provided)

### Optional 2 Outputs
- `P_loss_avg` - Average power losses [W] (if thyristor params provided)

## Implementation Notes

### Key Equations

**Half-Wave Rectifier:**
```
Vdc = (Vm / 2π) × (1 + cos(α))
Idc = (Vdc - Vbat) / Rbat
t_charge = (Capacity × ΔSoC / 100) / Idc
```

**Full-Wave Rectifiers (both CT and Bridge):**
```
Vdc = (Vm / π) × (1 + cos(α))
Idc = (Vdc - Vbat) / Rbat
t_charge = (Capacity × ΔSoC / 100) / Idc
```

Where:
- Vm = √2 × Vrms (peak voltage)
- α = firing angle in radians
- ΔSoC = desired change in State of Charge

### Full-Wave Center-Tapped Rectifier - Detailed Equations

#### Transformer and Voltage Equations
```
Vm = √2 × Vrms                    % Peak voltage per half-winding
V_abs(θ) = Vm × |sin(θ)|          % Rectified voltage magnitude
v_conv(θ) = V_abs(θ) - Vt         % Converter output (after thyristor drop)
```

#### Conduction Conditions
```
Gate enabled:  θ_mod ≥ α          % Within each half-cycle (0 to π)
Conduction:    v_conv > Vbat       % Voltage must exceed battery clamp
```

#### Current Equations
```
i_t(θ) = (v_conv(θ) - Vbat) / Rbat    % Instantaneous charging current
I_avg = (1/2π) ∫[0 to 2π] i_t(θ) dθ   % Average current
I_rms = √[(1/2π) ∫[0 to 2π] i_t²(θ) dθ]  % RMS current
```

#### Output Voltage Equations
```
v_out(θ) = v_conv(θ)  when conducting, else 0
V_avg = (1/2π) ∫[0 to 2π] v_out(θ) dθ
V_rms = √[(1/2π) ∫[0 to 2π] v_out²(θ) dθ]
```

For ideal resistive load:
```
Vdc_ideal = (Vm/π) × (1 + cos(α))
```

#### Power Loss Equations

**Battery Internal Losses:**
```
P_batt = I_rms² × Rbat            % I²R losses in battery resistance
```

**Thyristor Conduction Losses:**
```
P_thyristor = Vt × I_avg + Rth × I_rms²

Components:
  - Vt × I_avg:     Average loss due to forward voltage drop
  - Rth × I_rms²:   I²R losses in thyristor on-state resistance
```

**Thyristor Leakage Losses** (when Ileak > 0):
```
During blocking periods:
P_leak = (1/2π) ∫[blocking] V_blocking(θ) × Ileak dθ
```

**Thyristor Switching Losses** (when t_rise, t_fall > 0):
```
E_on = (1/6) × V_block × I_peak × t_rise
E_off = (1/6) × V_block × I_peak × t_fall
P_switching = f × (E_on + E_off)
```

**Total Power Loss:**
```
P_total = P_batt + P_thyristor + P_leak + P_switching
```

#### Battery Charging Equations

**Capacity Conversion:**
```
Q_C = Capacity [Ah] × 3600        % Total charge in Coulombs
If Capacity in Wh: Capacity_Ah = Capacity_Wh / Vbat
```

**Charging Time:**
```
ΔSoC = (SoC_target - SoC_init) / 100    % Fractional SoC change
t_charge = (Q_C × ΔSoC) / I_avg         % Time in seconds
```

**State of Charge Update:**
```
SoC_final = SoC_init + 100 × (I_avg × t_charge) / Q_C
```

#### Circuit Characteristics
- **Number of Thyristors:** 2 (one per half-winding)
- **Transformer:** Center-tapped secondary required
- **Conduction:** Each thyristor conducts for (π - α) per cycle
- **Frequency:** Output frequency = 2 × supply frequency
- **Voltage Drop:** Single Vt per half-cycle (one thyristor conducts)

## Testing

Suggested test cases:
1. α = 0° (maximum output)
2. α = 30° (typical operation)
3. α = 60° (moderate control)
4. α = 90° (zero average voltage for highly inductive loads)
5. α = 120° (reduced output)

## Expected Results

- Full-wave rectifiers should charge approximately 2× faster than half-wave
- Increasing firing angle increases charging time (non-linear relationship)
- With thyristor drops, actual performance will be slightly reduced
- Results should match theoretical calculations within acceptable tolerance

## Tips

1. Start with ideal thyristors (Vt = 0, Ileak = 0) to validate basic equations
2. Add non-idealities incrementally
3. Compare simulation results with hand calculations
4. Use small firing angle steps (e.g., 5°) for smooth plots
5. Validate that charging time → ∞ as α approaches values where Vdc ≈ Vbat
