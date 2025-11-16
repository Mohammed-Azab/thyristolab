function compare_rectifier_configs(Vrms, f, Vbat, Rbat, capacity)
% COMPARE_RECTIFIER_CONFIGS Compares all three rectifier configurations
%
% Syntax:
%   compare_rectifier_configs(Vrms, f, Vbat, Rbat, capacity)
%
% Inputs:
%   Vrms     - Supply voltage RMS value [V]
%   f        - Supply frequency [Hz]
%   Vbat     - Battery nominal voltage [V]
%   Rbat     - Battery internal resistance [Ohm]
%   capacity - Battery capacity [Ah]
%
% Description:
%   This utility function calls all three rectifier configurations and
%   generates comparative plots and tables to help analyze the differences
%   in performance.
%
% Example:
%   compare_rectifier_configs(230, 50, 12, 0.1, 50);
%
% TODO: Implement comparison function
% TODO: Call all three charger functions
% TODO: Generate side-by-side plots
% TODO: Create comparison table
% TODO: Calculate performance metrics (efficiency, charging time ratio, etc.)

fprintf('\n========================================\n');
fprintf('  Rectifier Configuration Comparison\n');
fprintf('========================================\n\n');

% TODO: Call half-wave charger
[alpha_hw, time_hw] = half_wave_charger(Vrms, f, Vbat, Rbat, capacity);

% TODO: Call full-wave center-tapped charger
[alpha_ct, time_ct] = full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity);

% TODO: Call full-wave bridge charger
[alpha_br, time_br] = full_wave_bridge_charger(Vrms, f, Vbat, Rbat, capacity);

% TODO: Generate comparison plot
figure('Name', 'Rectifier Configuration Comparison');
hold on;
plot(alpha_hw, time_hw, 'b-o', 'LineWidth', 2, 'DisplayName', 'Half-Wave');
plot(alpha_ct, time_ct, 'r-s', 'LineWidth', 2, 'DisplayName', 'Full-Wave CT');
plot(alpha_br, time_br, 'g-d', 'LineWidth', 2, 'DisplayName', 'Full-Wave Bridge');
hold off;
grid on;
xlabel('Firing Angle \alpha (degrees)');
ylabel('Charging Time (hours)');
title('Comparison of Rectifier Configurations');
legend('Location', 'best');

% TODO: Create comparison table for specific firing angles
fprintf('\nPerformance Comparison at Selected Firing Angles:\n');
fprintf('%-15s %-15s %-15s %-15s\n', 'Firing Angle', 'Half-Wave', 'Full-Wave CT', 'Full-Wave Bridge');
fprintf('%-15s %-15s %-15s %-15s\n', '(degrees)', '(hours)', '(hours)', '(hours)');
fprintf('---------------------------------------------------------------\n');
% TODO: Print comparison data

fprintf('\nKey Observations:\n');
fprintf('TODO: Add automatic analysis of results\n');
fprintf('- Full-wave configurations charge approximately 2x faster than half-wave\n');
fprintf('- Center-tapped and bridge have similar performance electrically\n');
fprintf('- Trade-off: CT requires special transformer, Bridge requires more SCRs\n\n');

end
