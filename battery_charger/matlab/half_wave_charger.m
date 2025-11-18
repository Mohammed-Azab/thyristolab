function [alpha_deg, charging_time_hours, SoC_final, P_loss_avg] = half_wave_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, varargin)
%
% Syntax:
%   [alpha_deg, charging_time_hours] = half_wave_charger(Vrms, f, Vbat, Rbat, capacity)
%   [alpha_deg, charging_time_hours, SoC_final] =  half_wave_charger(..., t_charge, SoC_init)
%   [alpha_deg, charging_time_hours, P_loss_avg] = half_wave_charger(..., Vt, Ileak, t_rise, t_fall)
%
% Inputs:
%   Vrms        - Supply voltage RMS value [V]
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
%   This function analyzes a half-wave controlled rectifier for battery charging.
%   It calculates the relationship between firing angle and charging time, and
%   optionally tracks State of Charge and calculates power losses.
%
% TODO: Implement the function according to project specifications
% TODO: Calculate average charging current for each firing angle
% TODO: Determine charging time based on capacity and desired SoC change
% TODO: Generate plot of firing angle vs charging time
% TODO: Implement optional SoC tracking
% TODO: Implement optional power loss calculation

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

% TODO: Define range of firing angles to analyze (e.g., 0 to 150 degrees)
alpha_deg = 0:5:150;  % Placeholder
alpha_rad = deg2rad(alpha_deg);

% Initialize output arrays
charging_time_hours = zeros(size(alpha_deg));
SoC_final = [];
P_loss_avg = [];



























% Generate main output plot
figure('Name', 'Half-Wave Rectifier - Firing Angle vs Charging Time');
plot(alpha_deg, charging_time_hours, 'b-o', 'LineWidth', 2);
grid on;
xlabel('Firing Angle \alpha (degrees)');
ylabel('Charging Time (hours)');
title('Half-Wave Controlled Rectifier for Battery Charging');
legend('Charging Time');

% Display summary
fprintf('\n=== Half-Wave Rectifier Battery Charger ===\n');
fprintf('Supply: %.1f V RMS, %.1f Hz\n', Vrms, f);
fprintf('Battery: %.1f V, %.3f Ohm, %.1f Ah\n', Vbat, Rbat, capacity);
fprintf('Initial SoC: %.1f%%\n', SoC_init);
if ~isempty(SoC_final)
    fprintf('Final SoC: %.1f%%\n', SoC_final);
end
fprintf('==========================================\n\n');

end
