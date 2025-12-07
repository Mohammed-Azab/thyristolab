% default battery_charger params

% Supply
Vrms = 230; 
f = 50;
alphas_deg = 0:5:180;            
simTime = 0.1;
T = 1/f;

% Load scenarios
scenarios = struct();
scenarios(1).name = 'resistive_only';
scenarios(1).R = 10;           % Resistance only (ohms)
scenarios(1).L = 1e-6;         % Minimal inductance (H)

scenarios(2).name = 'R_L_load';
scenarios(2).R = 10;           % Resistance (ohms)
scenarios(2).L = 50e-3;        % Moderate inductance (H)

scenarios(3).name = 'highly_inductive';
scenarios(3).R = 5;            % Lower resistance (ohms)
scenarios(3).L = 200e-3;       % High inductance (H)

% Pulse Generator parameters
pulse_amplitude = 10;           % Gate signal amplitude (V)
pulse_width = 50;              % Pulse width (% of period)
pulse_period = 1/f;            % Period based on line frequency (s)

% Model selection: 'all', 'ct', 'bg', 'hf', or cell array {'ct', 'bg'}
% 'ct' = center_taped, 'bg' = bridge (Full_Wave_Bridge), 'hf' = half_wave
modelSelection = 'ct';        % Options: 'all', 'ct', 'bg', 'hf', or {'ct', 'bg'}
    

% Live plotting option
enableLivePlot = false;  % Set to false to disable live plotting

% Graph generation option
generateGraphs = true;   % Set to true to generate graphs after sweep
graphAlphas = [30, 90, 180];  % Alpha values to plot (degrees)

R = 30;
L = 0.1;

phase_delay = 0.005;
phase_delay2 = 0.015;