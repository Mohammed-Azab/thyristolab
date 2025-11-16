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

4. **`compare_rectifier_configs.m`**
   - Utility function to compare all three configurations
   - Generates comparative plots and tables
   - Helps in configuration selection

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

### Optional Parameters (Name-Value pairs)
- `'t_charge'` - Charging time [sec] (if known)
- `'SoC_init'` - Initial State of Charge [%] (default: 20)
- `'Vt'` - Thyristor forward voltage drop [V] (default: 0)
- `'Ileak'` - Thyristor reverse leakage current [A] (default: 0)
- `'t_rise'` - Voltage/current rise time [sec] (default: 0)
- `'t_fall'` - Voltage/current fall time [sec] (default: 0)

## Function Outputs

### Standard Outputs
- `alpha_deg` - Array of firing angles [degrees]
- `charging_time_hours` - Corresponding charging times [hours]

### Optional Outputs
- `SoC_final` - Final State of Charge [%] (if t_charge provided)
- `P_loss_avg` - Average power losses [W] (if thyristor params provided)

## Implementation Notes

### TODO List for Implementation
Each function file contains detailed TODO comments indicating what needs to be implemented:

1. ✅ Function structure and documentation (provided)
2. ⬜ Calculate average output voltage for each firing angle
3. ⬜ Calculate average charging current
4. ⬜ Determine charging time based on current and capacity
5. ⬜ Implement SoC tracking (optional)
6. ⬜ Calculate power losses (optional)
7. ⬜ Generate plots and results

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

## Due Date
Part I is due on **November 22, 2025**
