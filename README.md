# ThyristoLab - Thyristor Rectifier Analysis & Design Toolkit

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Simulink](https://img.shields.io/badge/Simulink-Power%20Systems-blue.svg)](https://www.mathworks.com/products/simulink.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **A comprehensive MATLAB/Simulink toolkit for analyzing, designing, and simulating thyristor-based controlled rectifier circuits**

ThyristoLab provides engineers and researchers with powerful analytical tools and pre-built simulation models for controlled rectifier analysis. Whether you're designing battery chargers, DC motor drives, or power supplies, this toolkit offers validated functions and customizable Simulink models to accelerate your development process.

---

## Features

###  Battery Charger Design Tools
- **Analytical MATLAB Functions** for rapid prototyping and optimization
- Three rectifier topologies: half-wave, full-wave center-tapped, and bridge
- Firing angle optimization for charging profiles
- State of Charge (SoC) tracking and power loss analysis
- Comparative performance analysis across configurations

###  Load Analysis & Simulation
- **Pre-built Simulink Models** for detailed waveform analysis
- Support for multiple load types: resistive (R), inductive (RL), and highly inductive
- Study of continuous and discontinuous conduction modes
- Free-wheeling diode analysis and four-quadrant operation
- Comprehensive firing circuit implementations

### Analysis Capabilities
- Average and RMS voltage/current calculations
- Form factor and ripple factor analysis
- Extinction angle determination
- Power factor and efficiency metrics
- Theoretical vs. simulation validation

---

## Repository Structure

```
thyristolab/
├── battery_charger/              # Battery charging analysis tools
│   └── matlab/                   # Analytical functions
│       ├── half_wave_charger.m
│       ├── full_wave_ct_charger.m
│       ├── full_wave_bridge_charger.m
│       └── compare_rectifier_configs.m
│
├── load_analysis/                # Load behavior analysis
│   ├── simulink/                 # Simulation models
│   │   ├── half_wave_rectifier/
│   │   ├── full_wave_ct/
│   │   └── full_wave_bridge/
│   └── matlab/                   # Post-processing scripts
│
├── report/                       # Documentation templates
│   ├── root.tex                  # LaTeX report template
│   ├── chapters/                 # Individual sections
│   └── figures/                  # Generated plots
│
└── docs/                         # Additional documentation
```

---

## Quick Start

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/thyristolab.git
   cd thyristolab
   ```

2. **Requirements**
   - MATLAB R2020a or later
   - Simulink with Simscape Electrical (Power Systems Toolbox)
   - (Optional) LaTeX distribution for documentation

### Basic Usage

#### Battery Charger Analysis
```matlab
% Navigate to battery charger tools
cd battery_charger/matlab

% Define system parameters
Vrms = 230;        % Supply voltage (RMS) [V]
f = 50;            % Frequency [Hz]
Vbat = 12;         % Battery voltage [V]
Rbat = 0.1;        % Internal resistance [Ω]
capacity = 50;     % Battery capacity [Ah]

% Analyze half-wave rectifier
[alpha, charging_time] = half_wave_charger(Vrms, f, Vbat, Rbat, capacity);

% Compare all configurations
compare_rectifier_configs(Vrms, f, Vbat, Rbat, capacity);
```

#### Load Analysis (Simulink)
```matlab
% Open Simulink
simulink

% Navigate to load_analysis/simulink/
% Open desired model (e.g., half_wave_rectifier/)
% Run simulation and analyze waveforms
```

---

## Documentation

### Battery Charger Functions

All charger functions support the following interface:

**Inputs:**
- `Vrms` - Supply voltage RMS value [V]
- `f` - Supply frequency [Hz]
- `Vbat` - Battery nominal voltage [V]
- `Rbat` - Battery internal resistance [Ω]
- `capacity` - Battery capacity [Ah]

**Optional Parameters (Name-Value pairs):**
- `'t_charge'` - Charging time [sec]
- `'SoC_init'` - Initial State of Charge [%] (default: 20)
- `'Vt'` - Thyristor forward voltage drop [V] (default: 0)
- `'Ileak'` - Thyristor reverse leakage current [A] (default: 0)

**Outputs:**
- `alpha_deg` - Array of firing angles [degrees]
- `charging_time_hours` - Corresponding charging times [hours]
- `SoC_final` - Final State of Charge [%] (optional)
- `P_loss_avg` - Average power losses [W] (optional)

### Simulink Models

Pre-configured models for:
- **Half-wave rectifiers** with R, RL, and highly inductive loads
- **Full-wave center-tapped** configurations
- **Full-wave bridge** configurations
- **Free-wheeling diode** variants

Each model includes:
- Configurable firing angle control
- Voltage and current measurement points
- Scope displays for real-time visualization
- Data export capabilities for post-processing

Detailed setup instructions available in `load_analysis/simulink/README.md`

---

## Use Cases

### Power Electronics Design
- Battery charging system optimization
- DC motor drive analysis
- Industrial power supply design
- Renewable energy converter development

### Education & Research
- Power electronics course projects
- Rectifier behavior demonstration
- Load type impact studies
- Firing angle control strategies

### Validation & Testing
- Compare theoretical calculations with simulation
- Verify hardware designs before implementation
- Analyze edge cases and failure modes
- Optimize control parameters

---

## Technical Background

### Supported Rectifier Topologies

#### 1. Half-Wave Controlled Rectifier
- Single thyristor configuration
- Average voltage: `Vdc = (Vm/2π)(1 + cos α)`
- Lower cost, higher ripple
- Suitable for low-power applications

#### 2. Full-Wave Center-Tapped Rectifier
- Two thyristors with center-tapped transformer
- Average voltage: `Vdc = (Vm/π)(1 + cos α)`
- Better performance, requires special transformer
- Common in medium-power applications

#### 3. Full-Wave Bridge Rectifier
- Four thyristors in bridge configuration
- Average voltage: `Vdc = (Vm/π)(1 + cos α)`
- Standard transformer, higher component count
- Industry standard for general applications

### Load Types

- **Resistive (R)**: Simple behavior, discontinuous current
- **RL Load**: Phase lag, possible continuous conduction
- **Highly Inductive**: Nearly constant current, four-quadrant operation

---

## Advanced Features

### State of Charge Tracking
Monitor battery SoC during charging process:
```matlab
[alpha, time, SoC_final] = half_wave_charger(Vrms, f, Vbat, Rbat, capacity, ...
    't_charge', 3600, 'SoC_init', 20);
fprintf('Final SoC: %.1f%%\n', SoC_final);
```

### Thyristor Loss Analysis
Include non-ideal thyristor characteristics:
```matlab
[alpha, time, ~, P_loss] = full_wave_bridge_charger(Vrms, f, Vbat, Rbat, capacity, ...
    'Vt', 1.5, 'Ileak', 1e-6);
fprintf('Average power loss: %.2f W\n', P_loss);
```

### Four-Quadrant Operation
Study regenerative operation with highly inductive loads (α > 90°)

---

## Example Results

### Firing Angle vs. Charging Time
The toolkit generates comparative plots showing the relationship between firing angle and charging time for all three configurations, helping you select the optimal topology for your application.

### Waveform Analysis
Simulink models provide detailed voltage and current waveforms, allowing you to:
- Verify conduction modes
- Identify discontinuous vs. continuous operation
- Analyze ripple characteristics
- Study transient behavior

---

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

### Development Guidelines
1. Follow MATLAB coding standards
2. Include comprehensive comments
3. Provide usage examples
4. Update documentation
5. Test with multiple parameter sets

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- Inspired by power electronics literature (Dr.-Ing. Moustafa Adly)
- Built with MATLAB and Simulink
- Uses Simscape Electrical for accurate component modeling

---

## Contact

**Mohammed Azab**
- Email: mohammed@azab.io
- GitHub: [@mohammedazab](https://github.com/mohammedazab)

---

## Related Resources

- [MATLAB Power Electronics Examples](https://www.mathworks.com/help/sps/power-electronics.html)
- [Thyristor Fundamentals](https://www.electronics-tutorials.ws/power/thyristor.html)
- [Simscape Electrical Documentation](https://www.mathworks.com/help/sps/)

---

## Roadmap

- [ ] Add three-phase rectifier support
- [ ] Implement PWM rectifier models
- [ ] Add GUI for parameter configuration
- [ ] Include harmonic analysis tools
- [ ] Expand to IGBT and MOSFET converters
- [ ] Add real-time control simulation
- [ ] Create video tutorials

---

**ThyristoLab** - Professional thyristor rectifier analysis made simple.

*Star ⭐ this repository if you find it useful!*
