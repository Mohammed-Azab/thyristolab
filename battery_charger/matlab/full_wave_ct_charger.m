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
%   't_charge'  - Charging time [Hours].
%   'SoC_init'  - Initial state of charge [%]
%   'Vt'        - Thyristor forward drop per conducting device [V]
%   'Rth'       - Equivalent on-state resistance of thyristor [Ohm]
%   'Ileak'     - Reverse leakage Current [A]
%   't_rise'    - Rise time [s]
%   't_fall'    - Fall time[s]
%   'alpha_deg' - Vector of firing angles to analyze [deg] (for sweep analysis)
%   'alpha'     - Single firing angle [deg] (for detailed waveform plots)
%   'SoC_target'- Target state of charge [%]
%   'enablePlots'- Enable/disable all plots [boolean]
%
% Outputs:
%   alpha_deg           - Analyzed firing angles [deg] (vector)
%   charging_time_hours - Charging time for each alpha [h] (vector, same size as alpha_deg)
%   SoC_final           - Final SoC if 't_charge' provided [%] (vector)
%   P_loss_avg          - Average SCR conduction loss for each alpha [W] (vector)
%   metrics             - Struct with fields per alpha: Vavg, Vrms, Iavg, Irms,
%                         P_batt, P_thyristor, P_blocking, P_switching, P_total
%
% Notes:
% - Full-wave CT uses two SCRs; one drop Vt per half-cycle.
% - Ideal R-load average: Vdc = (Vm/(2*pi))*(1+cos(alpha)).
% - See README.md for detailed equations.

% ---------- Load parameters from workspace ----------
try
    SoC_target_ws = evalin('base', 'SoC_target');
catch
    SoC_target_ws = 80; 
end
try
    t_charge_ws = evalin('base', 't_charge');
    if isinf(t_charge_ws)
        t_charge_ws = [];
    end
catch
    t_charge_ws = []; 
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
    savePlots_ws = evalin('base', 'savePlots');
catch
    savePlots_ws = false;
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
addParameter(p, 't_charge', t_charge_ws, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
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
addParameter(p, 'savePlots', savePlots_ws, @(x) islogical(x) || isnumeric(x));
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
savePlots     = logical(p.Results.savePlots);
alpha         = p.Results.alpha;

if isstring(capUnit) || ischar(capUnit)
    if strcmpi(string(capUnit), "Wh")
        capacity = capacity / Vbat; % Wh -> Ah
    end
end

t_charge = t_charge*3600; % convert to Seconds

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
    
    if ~(Ileak == 0)
    % Thyristor conduction losses: P = Vt*I_avg + Rth*I_rms^2
        P_thyristor(k) = Vt*Iavg(k) + Rth*Irms(k)^2;
    end
    
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

    % Calculate charging time or final SoC
    % If t_charge is empty or inf, calculate time needed to reach SoC_target
    if isempty(t_charge) || isinf(t_charge)
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

% Create output directory if savePlots is enabled
if savePlots
    output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures', 'full_wave_ct');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
end

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
    fig1 = figure('Name', 'Full-Wave CT Rectifier - Voltages', 'Position', [100, 100, 1200, 800]);
    t_layout1 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(t_layout1,sprintf('Full-Wave CT Rectifier Voltages ($\\alpha = %.0f^\\circ$, $V_m=%.1f$ V)', alpha, Vm), 'Interpreter', 'latex', 'FontSize', 18);
    
    nexttile; 
    plot(t_ms, V_abs, 'b--', t_ms, Vbat*ones(size(t_ms)), 'g:', t_ms, v_out_demo, 'b-','LineWidth', 2); 
    grid on; 
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel('Voltage (V)', 'Interpreter', 'latex', 'FontSize', 14); 
    legend('$|V_{\mathrm{source}}|$', '$V_{\mathrm{bat}}$', '$V_{\mathrm{out}}$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
    title('Source \& Output Voltage', 'Interpreter', 'latex', 'FontSize', 14);
    
    nexttile; 
    plot(t_ms, v_th, 'r-','LineWidth', 2); 
    grid on; 
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel('$V_{\mathrm{th}}$ (V)', 'Interpreter', 'latex', 'FontSize', 14);
    title('Thyristor Voltage', 'Interpreter', 'latex', 'FontSize', 14);
    
    nexttile; 
    plot(t_ms, v_out_demo, 'b-','LineWidth', 2); 
    grid on; 
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel('$V_{\\mathrm{out}}$ (V)', 'Interpreter', 'latex', 'FontSize', 14);
    title('Output Voltage', 'Interpreter', 'latex', 'FontSize', 14);
    
    if savePlots
        saveas(fig1, fullfile(output_dir, sprintf('full_wave_ct_voltages_alpha_%d.png', round(alpha))));
        savefig(fig1, fullfile(output_dir, sprintf('full_wave_ct_voltages_alpha_%d.fig', round(alpha))));
    end
    
    % Figure 2: Currents
    fig2 = figure('Name', 'Full-Wave CT Rectifier - Currents', 'Position', [100, 100, 1200, 800]);
    t_layout2 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(t_layout2,sprintf('Full-Wave CT Rectifier Currents ($\\alpha = %.0f^\\circ$)', alpha), 'Interpreter', 'latex', 'FontSize', 18);
    
    nexttile; 
    plot(t_ms, i_demo, 'm-','LineWidth', 2); 
    grid on; 
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel('$I_{\mathrm{batt}}$ (A)', 'Interpreter', 'latex', 'FontSize', 14);
    title('Battery Current', 'Interpreter', 'latex', 'FontSize', 14);
    
    nexttile; 
    plot(t_ms, i_th, 'c-','LineWidth', 2); 
    grid on; 
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel('$I_{\mathrm{th}}$ (A)', 'Interpreter', 'latex', 'FontSize', 14);
    title('Thyristor Current', 'Interpreter', 'latex', 'FontSize', 14);
    
    nexttile;
    if Ileak > 0 || Rth > 0 || Vt > 0
        plot(t_ms, p_batt_inst, 'b-', t_ms, p_thyristor_inst, 'r-', ...
             t_ms, p_blocking_inst, 'g-', t_ms, p_total_inst, 'k--', 'LineWidth', 2);
        legend('$P_{\mathrm{batt}}$', '$P_{\mathrm{thyristor}}$', '$P_{\mathrm{blocking}}$', '$P_{\mathrm{total}}$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
        title('Instantaneous Power Losses', 'Interpreter', 'latex', 'FontSize', 14);
    else
        plot(t_ms, p_batt_inst, 'b-', 'LineWidth', 2);
        legend('$P_{\mathrm{batt}}$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
        title('Instantaneous Battery Power Loss', 'Interpreter', 'latex', 'FontSize', 14);
    end
    grid on; 
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel('Power Loss (W)', 'Interpreter', 'latex', 'FontSize', 14);
    
    if savePlots
        saveas(fig2, fullfile(output_dir, sprintf('full_wave_ct_currents_alpha_%d.png', round(alpha))));
        savefig(fig2, fullfile(output_dir, sprintf('full_wave_ct_currents_alpha_%d.fig', round(alpha))));
    end
    
end

% Plot charging time vs alpha only when t_charge is not provided (computing time to reach target SoC)
if enablePlots && (isempty(t_charge) || isinf(t_charge))
    figure('Name', 'Full-Wave CT Rectifier - Firing Angle vs Charging Time', 'Position', [100, 100, 900, 600]);
    
    % Use logarithmic scale for better visualization of large variations
    semilogy(alpha_deg, charging_time_hours, 'r-', 'LineWidth', 2.5, 'MarkerSize', 8, 'Marker', 's', 'MarkerFaceColor', 'r');
    grid on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);
    xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Charging Time (hours, log scale)', 'Interpreter', 'latex', 'FontSize', 16);
    title('Full-Wave Center-Tapped Controlled Rectifier', 'Interpreter', 'latex', 'FontSize', 18);
    legend('Charging Time', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
    
    % Add annotation for the current alpha value
    hold on;
    plot(alpha, charging_time_hours(alpha_idx), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
    text(alpha, charging_time_hours(alpha_idx)*1.3, sprintf('$\\alpha = %.0f^\\circ$\n%.2f h', alpha, charging_time_hours(alpha_idx)), ...
         'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', 'g', 'FontWeight', 'bold');
    hold off;
    
    % Set reasonable axis limits
    xlim([min(alpha_deg), max(alpha_deg)]);
    ylim([min(charging_time_hours(charging_time_hours>0))*0.5, max(charging_time_hours)*1.5]);
    
    if savePlots
        fig_temp = gcf;
        saveas(fig_temp, fullfile(output_dir, 'full_wave_ct_charging_time_vs_alpha.png'));
        savefig(fig_temp, fullfile(output_dir, 'full_wave_ct_charging_time_vs_alpha.fig'));
    end
end

% Power Losses vs Firing Angle
if enablePlots
    figure('Name', 'Full-Wave CT Rectifier - Power Losses vs Firing Angle', 'Position', [100, 100, 900, 600]);
    if any(P_thyristor > 0) || any(P_blocking > 0) || any(P_switching > 0)
        hold on;
        plot(alpha_deg, P_batt, 'b-', 'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        plot(alpha_deg, P_thyristor, 'r-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        if any(P_blocking > 0)
            plot(alpha_deg, P_blocking, 'g-', 'LineWidth', 2.5, 'Marker', 'd', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        end
        if any(P_switching > 0)
            plot(alpha_deg, P_switching, 'm-', 'LineWidth', 2.5, 'Marker', '^', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
        end
        plot(alpha_deg, P_total, 'k--', 'LineWidth', 3, 'Marker', 'pentagram', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
        hold off;
        
        leg_entries = {'$P_{\mathrm{batt}} = I_{\mathrm{rms}}^2 R_{\mathrm{bat}}$', '$P_{\mathrm{th,cond}} = V_t I_{\mathrm{avg}} + R_{\mathrm{th}} I_{\mathrm{rms}}^2$'};
        if any(P_blocking > 0)
            leg_entries{end+1} = '$P_{\mathrm{block}} = V_{\mathrm{block}} I_{\mathrm{leak}}$';
        end
        if any(P_switching > 0)
            leg_entries{end+1} = '$P_{\mathrm{switch}}$';
        end
        leg_entries{end+1} = '$P_{\mathrm{total}}$';
        legend(leg_entries, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
        title('Power Losses vs Firing Angle', 'Interpreter', 'latex', 'FontSize', 18);
    else
        % Plot only battery losses if thyristor losses are negligible
        plot(alpha_deg, P_batt, 'b-', 'LineWidth', 2.5, 'MarkerFaceColor', 'b');
        legend('$P_{\mathrm{batt}} = I_{\mathrm{rms}}^2 R_{\mathrm{bat}}$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
        title('Battery Power Losses vs Firing Angle', 'Interpreter', 'latex', 'FontSize', 18);
    end
    grid on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);
    xlabel('Firing Angle $\\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Power Loss (W)', 'Interpreter', 'latex', 'FontSize', 16);
    
    if savePlots
        fig_temp = gcf;
        saveas(fig_temp, fullfile(output_dir, 'full_wave_ct_power_losses_vs_alpha.png'));
        savefig(fig_temp, fullfile(output_dir, 'full_wave_ct_power_losses_vs_alpha.fig'));
    end
end

if enablePlots && ~isempty(t_charge) && ~isinf(t_charge)
    % Plot SoC vs time for ALL alpha values
    figure('Name', 'Battery State of Charge vs Time (All Alphas)', 'Position', [100, 100, 900, 600]);
    hold on;
    
    % Color map for different alpha values
    colors = jet(na);
    legend_entries = {};
    
    for k = 1:na
        % Calculate time to reach either 100% or the final SoC
        if SoC_final(k) >= 100
            % Calculate time to reach 100%
            time_to_100 = ((100 - SoC_init) / 100) * Q_C / Iavg(k) / 3600; % hours
            t_plot = linspace(0, time_to_100, 100);
            SoC_plot = SoC_init + (100 - SoC_init) * (t_plot / time_to_100);
        else
            % Use full charging time
            t_plot = linspace(0, t_charge/3600, 100);
            SoC_plot = SoC_init + (SoC_final(k) - SoC_init) * (t_plot / (t_charge/3600));
        end
        
        plot(t_plot, SoC_plot, 'LineWidth', 2, 'Color', colors(k,:));
        legend_entries{k} = sprintf('$\\alpha = %.0f^\\circ$', alpha_deg(k));
    end
    
    grid on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('State of Charge (\%)', 'Interpreter', 'latex', 'FontSize', 16);
    title('Battery Charging Profiles for Different Firing Angles', 'Interpreter', 'latex', 'FontSize', 18);
    
    % Limit legend entries to avoid warning
    if na <= 20
        legend(legend_entries, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'best', 'NumColumns', 2);
    else
        % Show legend for every other entry if there are many
        legend_subset = legend_entries(1:2:end);
        legend(legend_subset, 'Interpreter', 'latex', 'FontSize', 9, 'Location', 'best', 'NumColumns', 3);
    end
    
    ylim([max(0, SoC_init-5), 105]);
    hold off;
    
    if savePlots
        fig_temp = gcf;
        saveas(fig_temp, fullfile(output_dir, 'full_wave_ct_soc_vs_time_fixed_duration.png'));
        savefig(fig_temp, fullfile(output_dir, 'full_wave_ct_soc_vs_time_fixed_duration.fig'));
    end
elseif enablePlots
    % Plot SoC vs time for ALL alpha values (when t_charge not provided)
    figure('Name', 'Battery State of Charge vs Time (All Alphas)', 'Position', [100, 100, 900, 600]);
    hold on;
    
    % Color map for different alpha values
    colors = jet(na);
    legend_entries = {};
    
    for k = 1:na
        % Use the calculated charging time to reach target SoC for each alpha
        t_plot = linspace(0, charging_time_hours(k), 100);
        SoC_plot = SoC_init + (SoC_target - SoC_init) * (t_plot / charging_time_hours(k));
        
        plot(t_plot, SoC_plot, 'LineWidth', 2, 'Color', colors(k,:));
        legend_entries{k} = sprintf('$\\alpha = %.0f^\\circ$', alpha_deg(k));
    end
    
    grid on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('State of Charge (\%)', 'Interpreter', 'latex', 'FontSize', 16);
    title('Battery Charging Profiles for Different Firing Angles', 'Interpreter', 'latex', 'FontSize', 18);
    
    % Limit legend entries to avoid warning
    if na <= 20
        legend(legend_entries, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'best', 'NumColumns', 2);
    else
        % Show legend for every other entry if there are many
        legend_subset = legend_entries(1:2:end);
        legend(legend_subset, 'Interpreter', 'latex', 'FontSize', 9, 'Location', 'best', 'NumColumns', 3);
    end
    
    ylim([max(0, SoC_init-5), min(100, SoC_target+10)]);
    hold off;
    
    if savePlots
        fig_temp = gcf;
        saveas(fig_temp, fullfile(output_dir, 'full_wave_ct_soc_vs_time_target_soc.png'));
        savefig(fig_temp, fullfile(output_dir, 'full_wave_ct_soc_vs_time_target_soc.fig'));
    end
end

% Display system parameters
fprintf('\n========== Full-Wave CT Rectifier Analysis ==========\n');
fprintf('Supply : %.1f V RMS, %.1f Hz\n', Vrms, f);
fprintf('Battery: %.1f V, Rint -> %.3f Ohm, Capacity -> %.1f Ah\n', Vbat, Rbat, capacity);
if Vt > 0 || Ileak > 0 || t_rise > 0 || t_fall > 0
    fprintf('Thyristor: Vt -> %.2f V, Rth -> %.4f Ohm, \n Ileak -> %.2f A, t_rise -> %.2e s, t_fall -> %.2e s\n', Vt, Rth, Ileak, t_rise, t_fall);
end
fprintf('Alpha Range: [%d : %d] deg\n', min(alpha_deg), max(alpha_deg));
fprintf('======================================================\n');

fprintf('\n========== Battery CHARGER Params (α = %.0f°) ==========\n', alpha);
fprintf('Average Output (Vdc): %.2f V\n', Vavg(alpha_idx));
fprintf('Average Current     : %.2f A\n', Iavg(alpha_idx));
fprintf('RMS Current         : %.2f A\n', Irms(alpha_idx));
fprintf('\n--- Power Losses ---\n');
fprintf('Battery Loss        : %.3f KW\n', P_batt(alpha_idx)/1000);
if P_thyristor(alpha_idx) > 0
    fprintf('Thyristor Conduction: %.3f KW\n', P_thyristor(alpha_idx)/1000);
    fprintf('  - Forward drop    : %.3f KW\n', Vt*Iavg(alpha_idx)/1000);
end
if P_blocking(alpha_idx) > 0
    fprintf('  - Blocking Loss       : %.3f KW\n', P_blocking(alpha_idx)/1000);
end
if P_switching(alpha_idx) > 0
    fprintf('  - Switching Loss      : %.3f KW\n', P_switching(alpha_idx)/1000);
end
fprintf('Total Power Loss    : %.3f KW\n', P_total(alpha_idx)/1000);
if isempty(t_charge)
    fprintf('\nInitial SoC         : %.1f%%\n', SoC_init);
    fprintf('Target SoC          : %.1f%%\n', SoC_target);
    if charging_time_hours(alpha_idx) >= 0.25
        fprintf('Charging Time       : %.2f hours\n', charging_time_hours(alpha_idx));
    else 
        fprintf('Charging Time       : %.2f Minutes\n', charging_time_hours(alpha_idx)*60);
    end
else
    fprintf('\nInitial SoC         : %.1f%%\n', SoC_init);
    fprintf('Final SoC           : %.1f%%\n', SoC_final(alpha_idx));
    fprintf('Charging Time       : %.2f hours\n', t_charge/3600);
end
fprintf('=========================================\n\n');
end
