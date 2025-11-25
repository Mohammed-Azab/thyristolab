% COMPARE - Compares all three rectifier configurations
%
% Description:
%   This script loads parameters from params.m and calls all three 
%   rectifier configurations with three different scenarios, then 
%   generates comparative plots and tables to help analyze the 
%   differences in performance.
%

run params.m

fprintf('\n========================================\n');
fprintf('  Rectifier Configuration Comparison\n');
fprintf('========================================\n\n');
fprintf('Supply: %.1f V RMS, %.1f Hz\n', Vrms, f);
fprintf('Battery: %.1f V, %.3f Ohm, %.1f %s\n', Vbat, Rbat, capacity, capUnit);
fprintf('========================================\n\n');

%% ========== SCENARIO 0: Default (no optional parameters) ==========
fprintf('\n========== SCENARIO 0: Default Configuration ==========\n');
fprintf('No optional parameters (ideal thyristors)\n');
fprintf('=======================================================\n');

% Call all three chargers with enablePlots = false
[alpha_hw_s0, time_hw_s0, SoC_hw_s0, Ploss_hw_s0, metrics_hw_s0] = ...
    half_wave_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, 'enablePlots', false);

[alpha_ct_s0, time_ct_s0, SoC_ct_s0, Ploss_ct_s0, metrics_ct_s0] = ...
    full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, 'enablePlots', false);

[alpha_br_s0, time_br_s0, SoC_br_s0, Ploss_br_s0] = ...
    full_wave_bridge_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, 'enablePlots', false);

%% ========== SCENARIO 1: With t_charge and SoC_init ==========
fprintf('\n========== SCENARIO 1: Fixed Charging Time ==========\n');
fprintf('t_charge = %.2f hours, SoC_init = %.1f%%\n', t_charging_hours, SoC_init_ws);
fprintf('=====================================================\n');

[alpha_hw_s1, time_hw_s1, SoC_hw_s1, Ploss_hw_s1, metrics_hw_s1] = ...
    half_wave_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
    't_charge', t_charging_hours, 'SoC_init', SoC_init_ws, 'enablePlots', false);

[alpha_ct_s1, time_ct_s1, SoC_ct_s1, Ploss_ct_s1, metrics_ct_s1] = ...
    full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
    't_charge', t_charging_hours, 'SoC_init', SoC_init_ws, 'enablePlots', false);

[alpha_br_s1, time_br_s1, SoC_br_s1, Ploss_br_s1] = ...
    full_wave_bridge_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
    't_charge', t_charging_hours*3600, 'SoC_init', SoC_init_ws);

%% ========== SCENARIO 2: With non-ideal thyristor parameters ==========
fprintf('\n========== SCENARIO 2: Non-Ideal Thyristor Parameters ==========\n');
fprintf('Vt = %.2f V, Ileak = %.3f A, t_rise = %.2e s, t_fall = %.2e s\n', ...
    Vt, Ileak, t_rise, t_fall);
fprintf('================================================================\n');

[alpha_hw_s2, time_hw_s2, SoC_hw_s2, Ploss_hw_s2, metrics_hw_s2] = ...
    half_wave_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
    'Vt', Vt, 'Ileak', Ileak, 't_rise', t_rise, 't_fall', t_fall, 'enablePlots', false);

[alpha_ct_s2, time_ct_s2, SoC_ct_s2, Ploss_ct_s2, metrics_ct_s2] = ...
    full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
    'Vt', Vt, 'Ileak', Ileak, 't_rise', t_rise, 't_fall', t_fall, 'enablePlots', false);

[alpha_br_s2, time_br_s2, SoC_br_s2, Ploss_br_s2] = ...
    full_wave_bridge_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
    'Vt', Vt, 'Ileak', Ileak, 't_rise', t_rise, 't_fall', t_fall);

%% ========== GENERATE COMPARISON PLOTS ==========

% Plot 1: Scenario 0 - Charging Time vs Firing Angle
figure('Name', 'Scenario 0: Ideal Configuration', 'Position', [100, 100, 1000, 700]);
semilogy(alpha_hw_s0, time_hw_s0, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Half-Wave');
hold on;
semilogy(alpha_ct_s0, time_ct_s0, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave CT');
semilogy(alpha_br_s0, time_br_s0, 'g-d', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave Bridge');
hold off;
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Charging Time (hours, log scale)', 'Interpreter', 'latex', 'FontSize', 14);
title('Scenario 0: Ideal Thyristors - Charging Time Comparison', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 12);

% Plot 2: Scenario 1 - Final SoC Comparison
figure('Name', 'Scenario 1: Fixed Charging Time', 'Position', [150, 150, 1000, 700]);
plot(alpha_hw_s1, SoC_hw_s1, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Half-Wave');
hold on;
plot(alpha_ct_s1, SoC_ct_s1, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave CT');
plot(alpha_br_s1, SoC_br_s1, 'g-d', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave Bridge');
hold off;
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Final State of Charge (\%)', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('Scenario 1: Fixed Time (%.2f h) - Final SoC Comparison', t_charging_hours), 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 12);

% Plot 3: Scenario 2 - Charging Time with Non-Ideal Parameters
figure('Name', 'Scenario 2: Non-Ideal Thyristors', 'Position', [200, 200, 1000, 700]);
semilogy(alpha_hw_s2, time_hw_s2, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Half-Wave');
hold on;
semilogy(alpha_ct_s2, time_ct_s2, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave CT');
semilogy(alpha_br_s2, time_br_s2, 'g-d', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave Bridge');
hold off;
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Charging Time (hours, log scale)', 'Interpreter', 'latex', 'FontSize', 14);
title('Scenario 2: Non-Ideal Thyristors - Charging Time Comparison', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 12);

% Plot 4: Power Losses Comparison (Scenario 2)
figure('Name', 'Scenario 2: Power Losses', 'Position', [250, 250, 1000, 700]);
plot(alpha_hw_s2, metrics_hw_s2.P_total, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Half-Wave');
hold on;
plot(alpha_ct_s2, metrics_ct_s2.P_total, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave CT');
plot(alpha_br_s2, Ploss_br_s2, 'g-d', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Full-Wave Bridge');
hold off;
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Total Power Loss (W)', 'Interpreter', 'latex', 'FontSize', 14);
title('Scenario 2: Total Power Losses Comparison', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 12);

%% ========== COMPARISON TABLES ==========

% Select specific firing angles for comparison
alpha_compare = [0, 30, 60, 90, 120, 150];
fprintf('\n==================================================================\n');
fprintf('  SCENARIO 0: Performance Comparison at Selected Firing Angles\n');
fprintf('==================================================================\n');
fprintf('%-12s %-18s %-18s %-18s\n', 'Firing Angle', 'Half-Wave', 'Full-Wave CT', 'Full-Wave Bridge');
fprintf('%-12s %-18s %-18s %-18s\n', '(degrees)', '(hours)', '(hours)', '(hours)');
fprintf('------------------------------------------------------------------\n');

for alpha = alpha_compare
    [~, idx_hw] = min(abs(alpha_hw_s0 - alpha));
    [~, idx_ct] = min(abs(alpha_ct_s0 - alpha));
    [~, idx_br] = min(abs(alpha_br_s0 - alpha));
    
    fprintf('%-12d %-18.3f %-18.3f %-18.3f\n', alpha, ...
        time_hw_s0(idx_hw), time_ct_s0(idx_ct), time_br_s0(idx_br));
end

fprintf('\n==================================================================\n');
fprintf('  SCENARIO 1: Final SoC at Selected Firing Angles (%.2f h)\n', t_charging_hours);
fprintf('==================================================================\n');
fprintf('%-12s %-18s %-18s %-18s\n', 'Firing Angle', 'Half-Wave', 'Full-Wave CT', 'Full-Wave Bridge');
fprintf('%-12s %-18s %-18s %-18s\n', '(degrees)', '(%)', '(%)', '(%)');
fprintf('------------------------------------------------------------------\n');

for alpha = alpha_compare
    [~, idx_hw] = min(abs(alpha_hw_s1 - alpha));
    [~, idx_ct] = min(abs(alpha_ct_s1 - alpha));
    [~, idx_br] = min(abs(alpha_br_s1 - alpha));
    
    fprintf('%-12d %-18.2f %-18.2f %-18.2f\n', alpha, ...
        SoC_hw_s1(idx_hw), SoC_ct_s1(idx_ct), SoC_br_s1(idx_br));
end

fprintf('\n==================================================================\n');
fprintf('  SCENARIO 2: Performance with Non-Ideal Thyristors\n');
fprintf('==================================================================\n');
fprintf('%-12s %-18s %-18s %-18s\n', 'Firing Angle', 'Half-Wave', 'Full-Wave CT', 'Full-Wave Bridge');
fprintf('%-12s %-18s %-18s %-18s\n', '(degrees)', '(hours)', '(hours)', '(hours)');
fprintf('------------------------------------------------------------------\n');

for alpha = alpha_compare
    [~, idx_hw] = min(abs(alpha_hw_s2 - alpha));
    [~, idx_ct] = min(abs(alpha_ct_s2 - alpha));
    [~, idx_br] = min(abs(alpha_br_s2 - alpha));
    
    fprintf('%-12d %-18.3f %-18.3f %-18.3f\n', alpha, ...
        time_hw_s2(idx_hw), time_ct_s2(idx_ct), time_br_s2(idx_br));
end

%% ========== SUMMARY AND ANALYSIS ==========

fprintf('\n==================================================================\n');
fprintf('                         SUMMARY & ANALYSIS\n');
fprintf('==================================================================\n\n');

% Calculate average performance ratios at alpha = 30 degrees
[~, idx_30_hw] = min(abs(alpha_hw_s0 - 30));
[~, idx_30_ct] = min(abs(alpha_ct_s0 - 30));
[~, idx_30_br] = min(abs(alpha_br_s0 - 30));

ratio_ct_hw = time_hw_s0(idx_30_hw) / time_ct_s0(idx_30_ct);
ratio_br_hw = time_hw_s0(idx_30_hw) / time_br_s0(idx_30_br);
ratio_ct_br = time_ct_s0(idx_30_ct) / time_br_s0(idx_30_br);

fprintf('Key Observations (at α = 30°):\n\n');
fprintf('1. Charging Speed Comparison:\n');
fprintf('   - Full-wave configurations are %.2fx faster than half-wave\n', ratio_ct_hw);
fprintf('   - Center-tapped vs Half-wave: %.2fx speedup\n', ratio_ct_hw);
fprintf('   - Bridge vs Half-wave: %.2fx speedup\n', ratio_br_hw);
fprintf('   - Center-tapped vs Bridge: %.2fx (nearly identical)\n\n', ratio_ct_br);

fprintf('2. Configuration Trade-offs:\n');
fprintf('   - Half-Wave: Simplest (1 SCR), slowest charging\n');
fprintf('   - Full-Wave CT: 2 SCRs, requires center-tapped transformer\n');
fprintf('   - Full-Wave Bridge: 4 SCRs, standard transformer\n\n');

fprintf('3. Impact of Non-Ideal Parameters:\n');
avg_time_ideal = mean([time_hw_s0(idx_30_hw), time_ct_s0(idx_30_ct), time_br_s0(idx_30_br)]);
avg_time_nonideal = mean([time_hw_s2(idx_30_hw), time_ct_s2(idx_30_ct), time_br_s2(idx_30_br)]);
percent_increase = ((avg_time_nonideal - avg_time_ideal) / avg_time_ideal) * 100;
fprintf('   - Average charging time increase: %.2f%%\n', percent_increase);
fprintf('   - Thyristor voltage drop reduces effective charging voltage\n');
fprintf('   - Power losses increase with non-ideal parameters\n\n');

fprintf('4. Fixed Charging Time (Scenario 1):\n');
fprintf('   - At α = 30°, SoC after %.2f hours:\n', t_charging_hours);
fprintf('     * Half-Wave: %.1f%%\n', SoC_hw_s1(idx_30_hw));
fprintf('     * Full-Wave CT: %.1f%%\n', SoC_ct_s1(idx_30_ct));
fprintf('     * Full-Wave Bridge: %.1f%%\n\n', SoC_br_s1(idx_30_br));

fprintf('==================================================================\n\n');