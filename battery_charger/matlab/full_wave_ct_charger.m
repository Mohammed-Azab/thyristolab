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
%   't_charge'  - Charging time [s]
%   'SoC_init'  - Initial state of charge [%]
%   'Vt'        - Thyristor forward drop per conducting device [V]
%   'Rth'       - Equivalent on-state resistance of thyristor [Ohm]
%   'Ileak'     - Reverse leakage Current [A]
%   't_rise'    - Rise time [s]
%   't_fall'    - Fall time[s]
%   'alpha_deg' - Vector of firing angles to analyze [deg]
%   'SoC_target'- Target state of charge [%]
%
% Outputs:
%   alpha_deg           - Analyzed firing angles [deg]
%   charging_time_hours - Charging time for each alpha [h]
%   SoC_final           - Final SoC if 't_charge' provided [%]
%   P_loss_avg          - Average SCR conduction loss for each alpha [W]
%   metrics             - Struct with fields per alpha: Vavg, Vrms, Iavg, Irms,
%                         P_batt, P_thyristor, P_blocking, P_switching, P_total
%
% Notes:
% - Full-wave CT uses two SCRs; one drop Vt per half-cycle.
% - Ideal R-load average: Vdc = (Vm/(2*pi))*(1+cos(alpha)).
% - See README.md for detailed equations.

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
    alpha_ws = evalin('base', 'alpha');
catch
    alpha_ws = 90;
end
try
    SoC_init_ws = evalin('base', 'SoC_init');
catch
    SoC_init_ws = 90;
end


% ---------- Parse inputs ----------
p = inputParser;
addParameter(p, 't_charge', [], @isnumeric);
addParameter(p, 'SoC_init', SoC_init_ws, @isnumeric);
addParameter(p, 'SoC_target', SoC_target_ws, @isnumeric);
addParameter(p, 'Vt', 0, @isnumeric);
addParameter(p, 'alpha_deg', alpha_deg_ws, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'Rth', Rth_ws, @isnumeric);
addParameter(p, 'Ileak', 0, @isnumeric); 
addParameter(p, 't_rise', 0, @isnumeric); 
addParameter(p, 't_fall', 0, @isnumeric);
addParameter(p, 'alpha', alpha_ws, @isnumeric); 
addParameter(p, 'enablePlots', enablePlots_ws, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

t_charge      = p.Results.t_charge;
SoC_init      = p.Results.SoC_init;
SoC_target    = p.Results.SoC_target;
Vt            = p.Results.Vt;
alpha_deg     = p.Results.alpha_deg;
Rth           = p.Results.Rth;
Ileak         = p.Results.Ileak;
t_rise        = p.Results.t_rise;
t_fall        = p.Results.t_fall;
enablePlots   = logical(p.Results.enablePlots);
alpha         = p.Results.alpha;

if isstring(capUnit) || ischar(capUnit)
    if strcmpi(string(capUnit), "Wh")
        capacity = capacity / Vbat; % Wh -> Ah
    end
end

% ---------- Constants ----------
Vm = sqrt(2) * Vrms; % Peak voltage per half-winding

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
P_batt = zeros(1, na);        % Battery internal losses
P_thyristor = zeros(1, na);   % Thyristor conduction losses
P_blocking = zeros(1, na);    % Thyristor blocking/leakage losses
P_switching = zeros(1, na);   % Thyristor switching losses
P_total = zeros(1, na);       % Total power losses
charging_time_hours = zeros(1, na);
SoC_final = zeros(1, na);

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
    % Theoretical average for ideal resistive load (no battery clamp)
    %Vavg_theoretical(k) = (Vm/(2*pi)) * (1 + cos(a));
    Vout_rms(k) = sqrt(mean(v_out.^2));
    Iavg(k) = mean(i_t);
    Irms(k) = sqrt(mean(i_t.^2));

    % Power Loss Calculations
    % Battery internal resistance losses: P = I_rms^2 * Rbat
    P_batt(k) = Irms(k)^2 * Rbat;
    
    % Thyristor conduction losses: P = Vt*I_avg + Rth*I_rms^2
    P_thyristor(k) = Vt*Iavg(k) + Rth*Irms(k)^2;
    
    % Thyristor blocking/leakage losses: P = V_blocking * Ileak (average)
    if Ileak > 0
        v_blocking = V_abs;
        v_blocking(on) = 0;  % Zero when conducting
        P_blocking(k) = mean(v_blocking) * Ileak;
    else
        P_blocking(k) = 0;
    end
    
    % Thyristor switching losses: P = f * (E_on + E_off)
    if (t_rise > 0 || t_fall > 0) && Irms(k) > 0
        I_peak = max(i_t);
        V_block_avg = mean(V_abs(~on));
        if isnan(V_block_avg)
            V_block_avg = Vm;
        end
        E_on = (1/6) * V_block_avg * I_peak * t_rise;
        E_off = (1/6) * V_block_avg * I_peak * t_fall;
        P_switching(k) = f * (E_on + E_off);
    else
        P_switching(k) = 0;
    end
    
    % Total losses
    P_total(k) = P_batt(k) + P_thyristor(k) + P_blocking(k) + P_switching(k);
    
    % Legacy output for backward compatibility
    P_loss_avg(k) = P_thyristor(k);

    if isempty(t_charge)
        dSoC = max(SoC_target - SoC_init, 0)/100; 
        if Iavg(k) > 0
            t_sec = (Q_C * dSoC) / Iavg(k);
        else
            t_sec = inf;
        end
        charging_time_hours(k) = t_sec/3600;
        SoC_final(k) = SoC_target;
    else
        t_sec = t_charge;
        SoC_final(k) = min(100, SoC_init + 100*(Iavg(k)*t_sec)/Q_C);
        charging_time_hours(k) = t_sec/3600;
    end
end

metrics = struct('Vavg', Vavg, 'Vrms', Vout_rms, 'Iavg', Iavg, 'Irms', Irms, ...
                 'P_batt', P_batt, 'P_thyristor', P_thyristor, ...
                 'P_blocking', P_blocking, 'P_switching', P_switching, 'P_total', P_total);

[~, alpha_idx] = min(abs(alpha_deg - alpha));

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
    
    % Instantaneous power losses
    p_batt_inst = i_demo.^2 * Rbat;                    % Battery losses
    p_thyristor_inst = i_demo.^2 * Rth + Vt * i_demo;  % Thyristor losses
    
    % Blocking losses
    v_blocking_demo = V_abs;
    v_blocking_demo(on_demo) = 0;
    p_blocking_inst = v_blocking_demo * Ileak;
    
    p_total_inst = p_batt_inst + p_thyristor_inst + p_blocking_inst;     % Total losses
    
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
    if Ileak > 0 || Rth > 0 || Vt > 0
        plot(t_ms, p_batt_inst, 'b-', t_ms, p_thyristor_inst, 'r-', ...
             t_ms, p_blocking_inst, 'g-', t_ms, p_total_inst, 'k--', 'LineWidth', 1.2);
        legend('P_{batt}', 'P_{thyristor}', 'P_{blocking}', 'P_{total}', 'Location', 'best');
        title('Instantaneous Power Losses');
    else
        plot(t_ms, p_batt_inst, 'b-', 'LineWidth', 1.2);
        legend('P_{batt}', 'Location', 'best');
        title('Instantaneous Battery Power Loss');
    end
    grid on; xlabel('Time (ms)'); ylabel('Power Loss (W)');
    
end

figure('Name', 'Full-Wave CT Rectifier - Firing Angle vs Charging Time');
plot(alpha_deg, charging_time_hours, 'r-s', 'LineWidth', 2);
grid on;
xlabel('Firing Angle \alpha (degrees)');
ylabel('Charging Time (hours)');
title('Full-Wave Center-Tapped Controlled Rectifier');
legend('Charging Time');

% Power Losses vs Firing Angle
figure('Name', 'Full-Wave CT Rectifier - Power Losses vs Firing Angle');
if any(P_thyristor > 0) || any(P_blocking > 0) || any(P_switching > 0)
    hold on;
    plot(alpha_deg, P_batt, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    plot(alpha_deg, P_thyristor, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
    if any(P_blocking > 0)
        plot(alpha_deg, P_blocking, 'g-d', 'LineWidth', 2, 'MarkerSize', 6);
    end
    if any(P_switching > 0)
        plot(alpha_deg, P_switching, 'm-^', 'LineWidth', 2, 'MarkerSize', 6);
    end
    plot(alpha_deg, P_total, 'k-^', 'LineWidth', 2, 'MarkerSize', 6);
    hold off;
    
    leg_entries = {'Battery (I^2R_{bat})', 'Thyristor Conduction (Vt+Rth)'};
    if any(P_blocking > 0)
        leg_entries{end+1} = 'Blocking (V_{block}·I_{leak})';
    end
    if any(P_switching > 0)
        leg_entries{end+1} = 'Switching';
    end
    leg_entries{end+1} = 'Total Losses';
    legend(leg_entries, 'Location', 'best');
    title('Power Losses vs Firing Angle (Full-Wave CT)');
else
    % Plot only battery losses if thyristor losses are negligible
    plot(alpha_deg, P_batt, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    legend('Battery Losses (I^2R_{bat})', 'Location', 'best');
    title('Battery Power Losses vs Firing Angle (Full-Wave CT)');
end
grid on;
xlabel('Firing Angle \alpha (degrees)');
ylabel('Power Loss (W)');

if ~isempty(t_charge)
    t_charge_vec = linspace(0, t_charge/3600, 100); % hours
    SoC_vec = SoC_init + (SoC_final(alpha_idx) - SoC_init) * (t_charge_vec / (t_charge/3600));
    
    figure('Name', 'Battery State of Charge vs Time');
    plot(t_charge_vec, SoC_vec, 'b-', 'LineWidth', 2);
    hold on;
    plot(0, SoC_init, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(t_charge/3600, SoC_final(alpha_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('Time (hours)');
    ylabel('State of Charge (%)');
    title(sprintf('Battery Charging Profile (\\alpha = %.0f°)', alpha));
    legend('SoC', 'Initial State', 'Final State', 'Location', 'southeast');
    text(0, SoC_init-5, sprintf('  %.1f%%', SoC_init), 'FontSize', 10, 'Color', 'g');
    text(t_charge/3600, SoC_final(alpha_idx)+5, sprintf('  %.1f%%', SoC_final(alpha_idx)), 'FontSize', 10, 'Color', 'r');
    ylim([max(0, SoC_init-10), min(100, max(SoC_final(alpha_idx), SoC_target)+10)]);
    hold off;
else
    t_charge_hours = charging_time_hours(alpha_idx);
    t_charge_vec = linspace(0, t_charge_hours, 100);
    SoC_vec = SoC_init + (SoC_target - SoC_init) * (t_charge_vec / t_charge_hours);
    
    figure('Name', 'Battery State of Charge vs Time');
    plot(t_charge_vec, SoC_vec, 'b-', 'LineWidth', 2);
    hold on;
    plot(0, SoC_init, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(t_charge_hours, SoC_target, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('Time (hours)');
    ylabel('State of Charge (%)');
    title(sprintf('Battery Charging Profile (\\alpha = %.0f°)', alpha));
    legend('SoC', 'Initial State', 'Target State', 'Location', 'southeast');
    text(0, SoC_init-5, sprintf('  %.1f%%', SoC_init), 'FontSize', 10, 'Color', 'g');
    text(t_charge_hours, SoC_target+5, sprintf('  %.1f%%', SoC_target), 'FontSize', 10, 'Color', 'r');
    ylim([max(0, SoC_init-10), min(100, SoC_target+10)]);
    hold off;
end

% Display system parameters
fprintf('\n========== Full-Wave CT Rectifier Analysis ==========\n');
fprintf('Supply : %.1f V RMS, %.1f Hz\n', Vrms, f);
fprintf('Battery: %.1f V, Rint -> %.3f Ohm, Capacity -> %.1f Ah\n', Vbat, Rbat, capacity);
if ~isempty(Ileak)
    fprintf('Thyristor: Vt -> %.2f V, Rth -> %.4f Ohm, \n Ileak -> %.2f A, t_rise -> %.2e s, t_fall -> %.2e s\n', Vt, Rth, p.Results.Ileak, p.Results.t_rise, p.Results.t_fall);
end
fprintf('Alpha Range: [%d : %d] deg \n', min(alpha_deg), max(alpha_deg));
fprintf('======================================================\n');

fprintf('\n========== Battery CHARGER Params ==========\n');
fprintf('Firing Angle (α)    : %.0f°\n', alpha);
fprintf('Average Output (Vdc): %.2f V\n', Vavg(alpha_idx));
fprintf('Average Current     : %.2f A\n', Iavg(alpha_idx));
fprintf('RMS Current         : %.2f A\n', Irms(alpha_idx));
fprintf('\n--- Power Losses ---\n');
fprintf('Battery Loss        : %.3f W\n', P_batt(alpha_idx));
if P_thyristor(alpha_idx) > 0
    fprintf('Thyristor Conduction: %.3f W\n', P_thyristor(alpha_idx));
    fprintf('  - Forward drop    : %.3f W\n', Vt*Iavg(alpha_idx));
    fprintf('  - I²Rth           : %.3f W\n', Rth*Irms(alpha_idx)^2);
end
if P_blocking(alpha_idx) > 0
    fprintf('Blocking Loss       : %.3f W\n', P_blocking(alpha_idx));
end
if P_switching(alpha_idx) > 0
    fprintf('Switching Loss      : %.3f W\n', P_switching(alpha_idx));
end
fprintf('Total Power Loss    : %.3f W\n', P_total(alpha_idx));
if isempty(t_charge)
    fprintf('\nInitial SoC         : %.1f%%\n', SoC_init);
    fprintf('Target SoC          : %.1f%%\n', SoC_target);
    fprintf('Charging Time       : %.2f hours\n', charging_time_hours(alpha_idx));
else
    fprintf('\nInitial SoC         : %.1f%%\n', SoC_init);
    fprintf('Final SoC           : %.1f%%\n', SoC_final(alpha_idx));
    fprintf('Charging Time       : %.2f hours\n', t_charge/3600);
end
fprintf('=========================================\n\n');
end
