# Battery Charger Figures

This directory contains automatically generated plots from MATLAB simulations of various rectifier configurations.

## Directory Structure

- `half_wave/` - Plots from half-wave controlled rectifier simulations
- `full_wave_bridge/` - Plots from full-wave bridge rectifier simulations  
- `full_wave_ct/` - Plots from full-wave center-tap rectifier simulations

## Plot Types

Each rectifier configuration generates the following plots:

1. **Voltages** - Source, battery, output, and thyristor voltages
2. **Currents** - Battery and thyristor currents with instantaneous power losses
3. **Charging Time vs Alpha** - Relationship between firing angle and charging duration
4. **Power Losses vs Alpha** - Battery, thyristor, blocking, and switching losses
5. **SoC vs Time** - Battery state of charge profiles for different firing angles

## File Formats

Plots are saved in two formats:
- `.png` - High-resolution raster images for reports and documentation
- `.fig` - MATLAB figure files for further editing and analysis

## Usage

To enable plot saving, set `savePlots = true` in `params.m`. Plots will be automatically saved when running the charger functions with `enablePlots = true`.

## Note

Generated plot files (`.png` and `.fig`) are not tracked by Git to keep the repository lightweight.
