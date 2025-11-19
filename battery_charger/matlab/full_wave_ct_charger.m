function [alpha_deg, charging_time_hours, SoC_final, P_loss_avg, metrics] = full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, varargin)
% Syntax:
%   [alpha_deg, charging_time_hours] = full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit)
%   [alpha_deg, charging_time_hours, SoC_final, P_loss_avg, metrics] = ...
%       full_wave_ct_charger(..., 't_charge', t_charge, 'SoC_init', SoC_init, 'Vt', Vt, 'Rth', Rth)
%
% Inputs:
%   Vrms        - Secondary RMS voltage per half of CT transformer [V]
%   f           - Supply frequency [Hz]
%   Vbat        - Battery nominal voltage [V]
%   Rbat        - Battery internal resistance [Ohm]
%   capacity    - Battery capacity value
%   capUnit     - Battery capacity unit 'Ah' or 'Wh'
%
% Name-Value Optional Inputs:
%   'alpha_deg' - Vector of firing angles to analyze [deg]
%   't_charge'  - Charging time [s]
%   'SoC_init'  - Initial state of charge [%]
%   'SoC_target'- Target state of charge [%]
%   'Vt'        - Thyristor forward drop per conducting device [V]
%   'Rth'       - Equivalent on-state resistance of thyristor [Ohm]
%   'Ileak'     - Reverse leakage Current [A]
%   't_rise'    - Rise time [s]
%   't_fall'    - Fall time[s]
%
% Outputs:
%   alpha_deg           - Analyzed firing angles [deg]
%   charging_time_hours - Charging time for each alpha [h]
%   SoC_final           - Final SoC if 't_charge' provided [%]
%   P_loss_avg          - Average SCR conduction loss for each alpha [W]
%   metrics             - Struct with fields per alpha: Vavg, Vrms, Iavg, Irms
%
% Notes:
% - Full-wave CT uses two SCRs; one drop Vt per half-cycle.
% - Ideal R-load average: Vdc = (Vm/pi)*(1+cos(alpha)).

% ---------- Load parameters from workspace ----------
% Try to load specific parameters from workspace, use defaults if not found
try
    SoC_target_ws = evalin('base', 'SoC_target');
catch
    SoC_target_ws = 80; 
end
try
    alpha_deg_ws = evalin('base', 'alpha_deg');
catch
    alpha_deg_ws = 0:5:175;
end
try
    enablePlots_ws = evalin('base', 'enablePlots');
catch
    enablePlots_ws = false; 
end
try
    Rth_ws = evalin('base', 'Rth');
catch
    Rth_ws = 0;
end
try
    alpha = evalin('base', 'alpha');
catch
    alpha = 90;
end


% ---------- Parse inputs ----------
p = inputParser;
addParameter(p, 't_charge', [], @isnumeric);
addParameter(p, 'SoC_init', 20, @isnumeric);
addParameter(p, 'SoC_target', SoC_target_ws, @isnumeric);
addParameter(p, 'Vt', 0, @isnumeric);
addParameter(p, 'alpha_deg', alpha_deg_ws, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'Rth', Rth_ws, @isnumeric);
addParameter(p, 'Ileak', 0, @isnumeric); 
addParameter(p, 't_rise', 0, @isnumeric); 
addParameter(p, 't_fall', 0, @isnumeric);
addParameter(p, 'enablePlots', enablePlots_ws, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

t_charge      = p.Results.t_charge;
SoC_init      = p.Results.SoC_init;
SoC_target    = p.Results.SoC_target;
Vt            = p.Results.Vt;
alpha_deg     = p.Results.alpha_deg;
Rth           = p.Results.Rth;
enablePlots   = logical(p.Results.enablePlots);

if isstring(capUnit) || ischar(capUnit)
    if strcmpi(string(capUnit), "Wh")
        capacity = capacity / Vbat; % Wh -> Ah
    end
end

% ---------- Constants ----------
Vm = sqrt(2) * Vrms;

% ---------- Sampling ----------
N = 4096;
theta = linspace(0, 2*pi, N+1); theta(end) = []; % exclude endpoint
theta_mod = mod(theta, pi);            
V_abs = Vm * abs(sin(theta));

Q_C  = capacity * 3600;     % Ah -> Coulombs

% ---------- Outputs ----------
na = numel(alpha_deg);
Vavg = zeros(1, na);
Vout_rms = zeros(1, na);
Iavg = zeros(1, na);
Irms = zeros(1, na);
P_loss_avg = zeros(1, na);
charging_time_hours = zeros(1, na);
SoC_final = [];

for k = 1:na
    a = deg2rad(alpha_deg(k));

    % Gate available only when theta_mod >= a (within each half-cycle)
    gate = theta_mod >= a & theta_mod <= pi;

    v_conv = (V_abs - Vt);

    % Conduction only when above battery clamp
    cond  = v_conv > Vbat;

    i_t = zeros(size(theta));
    on  = gate & cond;
    i_t(on) = (v_conv(on) - Vbat) ./ Rbat;

    v_out = zeros(size(theta));
    v_out(on) = v_conv(on);

    Vavg(k) = mean(v_out);
    Vout_rms(k) = sqrt(mean(v_out.^2));
    Iavg(k) = mean(i_t);
    Irms(k) = sqrt(mean(i_t.^2));

    P_loss_avg(k) = Vt*Iavg(k) + Rth*mean(i_t.^2);

    if isempty(t_charge)
        dSoC = max(SoC_target - SoC_init, 0)/100; 
        if Iavg(k) > 0
            t_sec = (Q_C * dSoC) / Iavg(k);
        else
            t_sec = inf;
        end
        charging_time_hours(k) = t_sec/3600;
    else
        t_sec = t_charge;
        SoC_final(k) = min(100, SoC_init + 100*(Iavg(k)*t_sec)/Q_C);
        charging_time_hours(k) = t_sec/3600;
    end
end

metrics = struct('Vavg', Vavg, 'Vrms', Vout_rms, 'Iavg', Iavg, 'Irms', Irms);

if enablePlots
    % Time vector
    t_ms = theta / (2*pi*f) * 1000;
    
    a_demo = deg2rad(alpha);

    gate_demo = theta_mod >= a_demo & theta_mod <= pi;
    v_conv_demo = (V_abs - Vt);
    cond_demo = v_conv_demo > Vbat;
    i_demo = zeros(size(theta));
    on_demo = gate_demo & cond_demo;
    i_demo(on_demo) = (v_conv_demo(on_demo) - Vbat) ./ Rbat;
    v_out_demo = zeros(size(theta));
    v_out_demo(on_demo) = v_conv_demo(on_demo);
    p_inst = i_demo.^2 * Rth + Vt * i_demo; % Instantaneous power loss
    
    % Thyristor voltage and current
    % When thyristor is ON: Vth = Vt + Rth*Ith, Ith = i_demo
    % When thyristor is OFF: Vth = blocking voltage, Ith = 0
    i_th = i_demo;
    v_th = zeros(size(theta));
    v_th(on_demo) = Vt + Rth * i_th(on_demo);
    v_th(~on_demo) = V_abs(~on_demo);
    
    % Figure 1: Voltages
    figure('Name', 'Full-Wave CT Rectifier - Voltages');
    t_layout1 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(t_layout1,sprintf('Full-Wave CT Rectifier Voltages (\\alpha = %.0f°, Vm=%.1fV)', alpha, Vm));
    
    nexttile; 
    plot(t_ms, V_abs, 'b--', t_ms, Vbat*ones(size(t_ms)), 'g:', t_ms, v_out_demo, 'b-','LineWidth',1.2); 
    grid on; xlabel('Time (ms)'); ylabel('Voltage (V)'); 
    legend('Source |V|', 'V_{bat}', 'V_{out}', 'Location', 'best');
    title('Source & Output Voltage');
    
    nexttile; 
    plot(t_ms, v_th, 'r-','LineWidth',1.2); 
    grid on; xlabel('Time (ms)'); ylabel('V_{th} (V)');
    title('Thyristor Voltage');
    
    nexttile; 
    plot(t_ms, v_out_demo, 'b-','LineWidth',1.2); 
    grid on; xlabel('Time (ms)'); ylabel('V_{out} (V)');
    title('Output Voltage');
    
    nexttile; axis off; 
    text(0.05, 0.6, sprintf('Battery:\n  Capacity = %.1f Ah\n  SoC: %.0f%% → %.0f%%\n  V_{bat} = %.1f V\n  R_{bat} = %.3f Ω', ...
        capacity, SoC_init, SoC_target, Vbat, Rbat), 'FontSize', 9);
    
    % Find index of alpha in alpha_deg array for metrics display
    [~, alpha_idx] = min(abs(alpha_deg - alpha));
    text(0.05, 0.15, sprintf('Averages (α=%.0f°):\n  V_{avg} = %.2f V\n  I_{avg} = %.2f A\n  P_{loss} = %.2f W', ...
        alpha, Vavg(alpha_idx), Iavg(alpha_idx), P_loss_avg(alpha_idx)), 'FontSize', 9);
    
    % Figure 2: Currents
    figure('Name', 'Full-Wave CT Rectifier - Currents');
    t_layout2 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(t_layout2,sprintf('Full-Wave CT Rectifier Currents (\\alpha = %.0f°)', alpha));
    
    nexttile; 
    plot(t_ms, i_demo, 'm-','LineWidth',1.2); 
    grid on; xlabel('Time (ms)'); ylabel('I_{batt} (A)');
    title('Battery Current');
    
    nexttile; 
    plot(t_ms, i_th, 'c-','LineWidth',1.2); 
    grid on; xlabel('Time (ms)'); ylabel('I_{th} (A)');
    title('Thyristor Current');
    
    nexttile; 
    plot(t_ms, p_inst, 'r-','LineWidth',1.2); 
    grid on; xlabel('Time (ms)'); ylabel('P_{loss} (W)');
    title('Instantaneous Power Loss');
    
    nexttile; axis off;
    text(0.05, 0.5, sprintf('Current Metrics:\n  I_{avg} = %.2f A\n  I_{rms} = %.2f A\n  P_{loss,avg} = %.2f W', ...
        Iavg(alpha_idx), Irms(alpha_idx), P_loss_avg(alpha_idx)), 'FontSize', 10);
end

figure('Name', 'Full-Wave CT Rectifier - Firing Angle vs Charging Time');
plot(alpha_deg, charging_time_hours, 'r-s', 'LineWidth', 2);
grid on;
xlabel('Firing Angle \\alpha (degrees)');
ylabel('Charging Time (hours)');
title('Full-Wave Center-Tapped Controlled Rectifier');
legend('Charging Time');

% Console summary
fprintf('\n=== Full-Wave Center-Tapped Rectifier Battery Charger ===\n');
fprintf('Supply : %.1f V RMS, %.1f Hz\n', Vrms, f);
fprintf('Battery: %.1f V, Rint=%.3f Ohm, Capacity=%.1f Ah\n', Vbat, Rbat, capacity);
if ~isempty(Vt)
    fprintf('Thyristor: Vt=%.2f V, Rth=%.4f Ohm\n', Vt, Rth);
end
fprintf('Analyzed alphas: [%d..%d] deg (n=%d)\n', round(alpha_deg(1)), round(alpha_deg(end)), numel(alpha_deg));
if isempty(t_charge)
    fprintf('Charging window: SoC %.1f%% -> %.1f%%\n', SoC_init, SoC_target);
else
    fprintf('t_charge=%.1f s; SoC_init=%.1f%%\n', t_charge, SoC_init);
end
fprintf('========================================================\n\n');

end
