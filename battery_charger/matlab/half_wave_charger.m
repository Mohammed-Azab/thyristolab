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
I_charge = zeros(size(alpha_deg));

% Convert capacity from Wh to Coulomb Ah
if strcmpi(capUnit, 'Wh')
    capacity_Ah = capacity / Vbat;
else
    capacity_Ah = capacity;
end

% Initialized inputs for option 1 and 2
nbOfAlphas = numel(alpha_deg);
charging_time_hours = nan(1, nbOfAlphas);
SoC_final = nan(1, nbOfAlphas);
P_loss_avg = nan(nbOfAlphas);

% TASK 1 get charging currents
fprintf('Charging currents...\n');
disp("FUNCTION LOADED SUCCESSFULLY")

for i = 1 : nbOfAlphas
    alpha = alpha_rad(i);
   % Get Voltage first
    if alpha >= 0 && alpha <= pi
        V_avg = (Vm / (2*pi))  *(1 + cos(alpha)) - Vt/2;
    else
        V_avg = 0;
    end
   % since the voltage HAS to be positive we need to check

   V_avg = max(0,V_avg);

   if V_avg > Vbat % half wave we need to have the ouputvoltage bigger than battery voltage
       I_charge(i) = (V_avg - Vbat) / Rbat;
   else
       I_charge(i) = 0;
   end
   
   I_charge(i) = max(0,I_charge(i) - Ileak); % to ignore leakage current

end

% TASK 2 get charging times and SoC when t_charge is not given
fprintf('Charging times...\n');

% When the t_charge is not given
if isempty(t_charge)
    SoC_target = 80;
    charge_needed_Ah = capacity_Ah * (SoC_target - SoC_init) / 100;
else
    charge_needed_Ah = capacity_Ah * 0.6;
end

for i = 1 : length(alpha_deg)
    if I_charge(i) > 0.001
        charging_time_hours(i) = charge_needed_Ah / I_charge(i);
    else 
        charging_time_hours(i) = inf;
    end
end
 % when t_charge not give SoC
fprintf('Final State of Charge...\n');

if ~isempty(t_charge)

    valid_indices = find(I_charge > 0.001);
    if ~isempty(valid_indices)
        optimal_angle = valid_indices(1);
        I_optimal = I_charge(optimal_angle);

        charge_needed_Ah = I_optimal * (t_charge) / 3600;
        SoC_final = SoC_init + (charge_needed_Ah / capacity_Ah) * 100;
        SoC_final = min(SoC_final, 100);
    else
        SoC_final = SoC_init;
    end
end

% TASK 3 get SoC in case when the t_charge is given

if ~isempty(t_charge)
    for i = 1 : nbOfAlphas
        time_delivered = I_charge(i) * (t_charge / 3600); % Calculate time delivered based on charging current
        SoC_final  = SoC_init + (time_delivered / capacity_Ah) * 100;
        SoC_final = min(SoC_final(i),100);
    end
end

% TASK 4 get the power loss

fprintf('Power loss...\n');

for i = 1 : nbOfAlphas
    P_conduction = I_charge(i) * Vt;   
    P_leak = Vrms * Ileak;    
    P_loss_avg(i) = P_conduction + P_leak;
end





















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

% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
end
