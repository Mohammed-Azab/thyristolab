# Part II: Load Analysis - MATLAB Scripts

This directory contains MATLAB scripts for analytical calculations and post-processing of Simulink simulation results.

## Purpose

These scripts serve to:
1. Perform analytical calculations for validation
2. Process Simulink simulation data
3. Generate comparison plots and tables
4. Calculate performance metrics
5. Validate simulation results against theory

## Suggested Scripts to Create

### 1. `analytical_half_wave.m`
Calculate theoretical values for half-wave rectifier with different loads

```matlab
function [Vdc, Vrms, Idc, Irms] = analytical_half_wave(Vs, f, alpha, R, L)
% Analytical calculations for half-wave controlled rectifier
% TODO: Implement equations from theory
end
```

### 2. `analytical_full_wave.m`
Calculate theoretical values for full-wave rectifiers

```matlab
function [Vdc, Vrms, Idc, Irms] = analytical_full_wave(Vs, f, alpha, R, L)
% Analytical calculations for full-wave controlled rectifier
% TODO: Implement equations from theory
end
```

### 3. `process_simulation_data.m`
Process exported Simulink data to calculate metrics

```matlab
% Load simulation data from workspace
% Calculate Vdc, Vrms, Idc, Irms
% Compare with analytical values
% TODO: Implement data processing
```

### 4. `plot_waveforms.m`
Generate publication-quality plots from simulation data

```matlab
% Create formatted plots of voltage and current waveforms
% TODO: Implement plotting functions
```

### 5. `compare_configurations.m`
Compare different rectifier configurations

```matlab
% Compare half-wave vs full-wave
% Compare center-tapped vs bridge
% Compare with vs without FWD
% TODO: Implement comparison analysis
```

### 6. `firing_angle_sweep.m`
Automate simulation with varying firing angles

```matlab
% Run simulations for multiple alpha values
% Collect results
% Generate plots of Vdc vs alpha, etc.
% TODO: Implement parameter sweep
```

## Key Equations for Implementation

### Half-Wave Rectifier

**Resistive Load:**
```matlab
Vm = sqrt(2) * Vrms;
omega = 2 * pi * f;
alpha_rad = deg2rad(alpha);

% Average voltage
Vdc = (Vm / (2*pi)) * (1 + cos(alpha_rad));

% RMS voltage
Vrms_out = (Vm/2) * sqrt((1/pi) * ((pi - alpha_rad) + sin(2*alpha_rad)/2));

% Currents (divide by R)
Idc = Vdc / R;
Irms = Vrms_out / R;
```

**RL Load (more complex):**
- Requires solving transcendental equation for extinction angle β
- See theory section in report for detailed equations

### Full-Wave Rectifier

**Resistive Load:**
```matlab
% Average voltage (note: factor of pi instead of 2*pi)
Vdc = (Vm / pi) * (1 + cos(alpha_rad));

% RMS voltage
Vrms_out = (Vm/sqrt(2)) * sqrt((1/pi) * ((pi - alpha_rad) + sin(2*alpha_rad)/2));
```

## Example: Complete Analytical Script

```matlab
%% analytical_comparison.m
% Compare analytical predictions with simulation results

clear; clc; close all;

% Parameters
Vrms_supply = 230;  % V
f = 50;             % Hz
R = 10;             % Ohm
L_RL = 0.05;        % H (for RL load)
alpha_values = 0:30:150;  % degrees

% Initialize results
Vdc_hw_theory = zeros(size(alpha_values));
Vdc_fw_theory = zeros(size(alpha_values));

% Calculate for each firing angle
for i = 1:length(alpha_values)
    alpha = alpha_values(i);
    
    % Half-wave analytical
    [Vdc_hw_theory(i), ~, ~, ~] = analytical_half_wave(Vrms_supply, f, alpha, R, 0);
    
    % Full-wave analytical
    [Vdc_fw_theory(i), ~, ~, ~] = analytical_full_wave(Vrms_supply, f, alpha, R, 0);
end

% Plot results
figure;
plot(alpha_values, Vdc_hw_theory, 'b-o', 'LineWidth', 2, 'DisplayName', 'Half-Wave');
hold on;
plot(alpha_values, Vdc_fw_theory, 'r-s', 'LineWidth', 2, 'DisplayName', 'Full-Wave');
grid on;
xlabel('Firing Angle \alpha (degrees)');
ylabel('Average Output Voltage V_{dc} (V)');
title('Analytical: Average Voltage vs Firing Angle');
legend('Location', 'best');

% TODO: Load simulation data and compare
% TODO: Calculate percentage error
% TODO: Generate comparison table
```

## Data Import from Simulink

After running Simulink simulation with "To Workspace" blocks:

```matlab
% Data should be in workspace as arrays
% Example variable names: Vo_hw_R, Io_hw_R, tout

% Calculate average
Vdc_sim = mean(Vo_hw_R);

% Calculate RMS
Vrms_sim = rms(Vo_hw_R);

% For periodic signals, calculate over one period
T = 1/f;  % Period
N_samples_per_period = round(T / (tout(2) - tout(1)));
% Use last period for steady-state
Vdc_steady = mean(Vo_hw_R(end-N_samples_per_period:end));
```

## Validation Procedure

1. **Run analytical script** - Get theoretical values
2. **Run Simulink simulation** - Get simulated values
3. **Compare results** - Calculate error
4. **Investigate discrepancies** - If error > 5%, check:
   - Simulation reached steady-state?
   - Correct parameter values?
   - Numerical accuracy settings?

## Performance Metrics to Calculate

```matlab
% Form Factor
FF = Vrms / Vdc;

% Ripple Factor
RF = sqrt(FF^2 - 1);

% Rectification Efficiency
eta_rect = (Vdc / Vrms_supply)^2;

% Transformer Utilization Factor (TUF)
TUF = Pdc / (Vs_rated * Is_rated);

% Power Factor (more complex, requires harmonic analysis)
```

## Tips

1. **Modular Functions** - Create reusable functions for common calculations
2. **Clear Documentation** - Comment your code thoroughly
3. **Error Checking** - Validate inputs, handle edge cases
4. **Plotting Standards** - Use consistent colors, labels, and formats
5. **Save Results** - Export data and figures for report

## MATLAB Live Scripts

Consider creating MATLAB Live Scripts (.mlx) for:
- Interactive parameter exploration
- Combined code, results, and explanations
- Easy inclusion in report

## TODO Checklist

- [ ] Create analytical calculation functions
- [ ] Implement data import from Simulink
- [ ] Calculate performance metrics
- [ ] Generate comparison plots
- [ ] Create validation script
- [ ] Document all functions
- [ ] Test with sample data
- [ ] Prepare for report integration

## Due Date
Part II is due on **December 6, 2025**
