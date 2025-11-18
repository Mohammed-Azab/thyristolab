function [alpha_deg, charging_time_hours, SoC_final, P_loss_avg] = full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, varargin)
% FULL_WAVE_CT_CHARGER Analyzes full-wave center-tapped controlled rectifier for battery charging
%
% Syntax:
%   [alpha_deg, charging_time_hours] = full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity)
%   [alpha_deg, charging_time_hours, SoC_final, P_loss_avg] = full_wave_ct_charger(..., t_charge, SoC_init)
%   [alpha_deg, charging_time_hours, P_loss_avg]= full_wave_ct_charger(..., t_charge, SoC_init, Vt, Ileak, t_rise, t_fall)
%
% Inputs:
%   Vrms        - Supply voltage RMS value (per half of CT transformer) [V]
%   f           - Supply frequency [Hz]
%   Vbat        - Battery nominal voltage [V]
%   Rbat        - Battery internal resistance [Ohm]
%   capacity    - Battery capacity value
%   capUnit     - Battery capacity [Ah or Wh] -string
%
% Optional Inputs (Name-Value pairs):
%   't_charge'  - Charging time [sec] (default: calculated for 20-80% SoC)
%   'SoC_init'  - Initial State of Charge [%] (default: 20)
%   'Vt'        - Thyristor forward voltage drop [V] (default: 0, ideal)
%   'Ileak'     - Thyristor reverse leakage current [A] (default: 0)
%   't_rise'    - Voltage/current rise time [sec] (default: 0)
%   't_fall'    - Voltage/current fall time [sec] (default: 0)
%
% Outputs:
%   alpha_deg           - Array of firing angles [degrees]
%   charging_time_hours - Corresponding charging times [hours]
%   SoC_final           - Final State of Charge [%] (if t_charge provided)
%   P_loss_avg          - Average power losses [W] (if thyristor params provided)
%
% Description:
%   This function analyzes a full-wave center-tapped controlled rectifier for
%   battery charging. It uses a center-tapped transformer with two thyristors,
%   one for each half-cycle. The output voltage is twice that of a half-wave
%   rectifier for the same firing angle.
%
% Circuit Configuration:
%   - Center-tapped transformer with two secondary windings
%   - Two thyristors (T1 and T2)
%   - T1 conducts during positive half-cycle
%   - T2 conducts during negative half-cycle
%   - Common cathode connection to load
%
% TODO: Implement the function according to project specifications
% TODO: Account for two pulses per cycle (double frequency compared to half-wave)
% TODO: Calculate average output voltage: Vdc = (Vm/pi) * (1 + cos(alpha))
% TODO: Compare results with half-wave configuration
% TODO: Generate comparative plots if needed

% Parse optional inputs
p = inputParser;
addParameter(p, 't_charge', [], @isnumeric);
addParameter(p, 'SoC_init', 20, @isnumeric);
addParameter(p, 'Vt', 0, @isnumeric);
addParameter(p, 'Ileak', 0, @isnumeric);
addParameter(p, 't_rise', 0, @isnumeric);
addParameter(p, 't_fall', 0, @isnumeric);
parse(p, varargin{:});

% Extract parameters
t_charge = p.Results.t_charge;
SoC_init = p.Results.SoC_init;
Vt = p.Results.Vt;
Ileak = p.Results.Ileak;
t_rise = p.Results.t_rise;
t_fall = p.Results.t_fall;

% Constants
Vm = sqrt(2) * Vrms;  % Peak voltage
omega = 2 * pi * f;   % Angular frequency

% TODO: Define range of firing angles to analyze
alpha_deg = 0:5:150;  % Placeholder
alpha_rad = deg2rad(alpha_deg);

% Initialize output arrays
charging_time_hours = zeros(size(alpha_deg));
SoC_final = [];
P_loss_avg = [];


























% Generate main output plot
figure('Name', 'Full-Wave CT Rectifier - Firing Angle vs Charging Time');
plot(alpha_deg, charging_time_hours, 'r-s', 'LineWidth', 2);
grid on;
xlabel('Firing Angle \alpha (degrees)');
ylabel('Charging Time (hours)');
title('Full-Wave Center-Tapped Controlled Rectifier for Battery Charging');
legend('Charging Time');

% Display summary
fprintf('\n=== Full-Wave Center-Tapped Rectifier Battery Charger ===\n');
fprintf('Supply: %.1f V RMS per half, %.1f Hz\n', Vrms, f);
fprintf('Battery: %.1f V, %.3f Ohm, %.1f Ah\n', Vbat, Rbat, capacity);
fprintf('Configuration: Center-Tapped Transformer with 2 SCRs\n');
fprintf('Initial SoC: %.1f%%\n', SoC_init);
if ~isempty(SoC_final)
    fprintf('Final SoC: %.1f%%\n', SoC_final);
end
fprintf('========================================================\n\n');

end
