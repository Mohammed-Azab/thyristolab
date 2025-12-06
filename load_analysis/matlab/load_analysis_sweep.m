% Load parameters from params.m
run params.m

% Close all loaded Simulink models to avoid conflicts
bd = find_system('type', 'block_diagram');
for i = 1:length(bd)
    if ~strcmp(bd{i}, 'simulink')
        try
            close_system(bd{i}, 0);
        catch
        end
    end
end

saveFolder = fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

% locate models in ../simulink
scriptDir = fileparts(mfilename('fullpath'));
slxDir = fullfile(scriptDir,'..','simulink');
files = dir(fullfile(slxDir,'*.slx'));
if isempty(files)
    error('No .slx models found in %s. Place your models there and re-run.', slxDir);
end

models = {files.name};
fprintf('Found %d model(s) in %s:\n', numel(models), slxDir);
for k=1:numel(models)
    fprintf(' - %s\n', models{k});
end

fprintf('\nLoad scenarios defined:\n');
for s=1:numel(scenarios)
    fprintf(' %d. %s: R=%.2f Ω, L=%.1f mH\n', s, scenarios(s).name, scenarios(s).R, scenarios(s).L*1000);
end

if enableLivePlot
    figHandle = figure('Name', 'Live Waveforms', 'NumberTitle', 'off');
end

% Results container
results = struct();
% Main sweep loop - iterate over scenarios and models
for si = 1:numel(scenarios)
    scenario = scenarios(si);
    fprintf('\n========== Scenario %d: %s ==========\n', si, scenario.name);
    
    for mi = 1:numel(models)
        modelName = models{mi};
        modelPath = fullfile(slxDir, modelName);
        [~,modelBase] = fileparts(modelName);
        
        fprintf('\nRunning model: %s\n', modelName);
        
        % Close any existing model with the same name to avoid conflicts
        if bdIsLoaded(modelBase)
            fprintf('  Closing existing model instance...\n');
            try
                % Stop simulation if running
                set_param(modelBase, 'SimulationCommand', 'stop');
                pause(0.5); % Wait for simulation to stop
            catch
            end
            close_system(modelBase, 0);
        end
        
        % Load model from specified path
        load_system(modelPath);
        
        modelResults = cell(numel(alphas_deg),1);
        
        for ai = 1:numel(alphas_deg)
            alpha_deg = alphas_deg(ai);
            % convert alpha to phase delay time for Pulse Generator
            phase_delay = (alpha_deg/360) * (1/f);
            
            % push variables to base workspace so model blocks can read them
            assignin('base','alpha_deg',alpha_deg);
            assignin('base','phase_delay',phase_delay);
            assignin('base','f',f);
            assignin('base','Vrms',Vrms);
            
            % Load scenario parameters
            assignin('base','R',scenario.R);
            assignin('base','L',scenario.L);
            
            % Pulse generator parameters
            assignin('base','pulse_amplitude',pulse_amplitude);
            assignin('base','pulse_width',pulse_width);
            assignin('base','pulse_period',pulse_period);
            
            fprintf('  alpha = %3d deg -> phase_delay = %.6f s ... ', alpha_deg, phase_delay);
            
            % run simulation and collect SimulationOutput
            try
                simOut = sim(modelBase, ...
                    'StopTime', num2str(simTime), ...
                    'SaveOutput','on', ...
                    'SaveTime','on', ...
                    'SaveState','off', ...
                    'ReturnWorkspaceOutputs','on');
            catch ME
                fprintf('FAILED\n');
                warning('Simulation failed for %s alpha=%d: %s', modelName, alpha_deg, ME.message);
                modelResults{ai} = struct('alpha_deg',alpha_deg,'error',ME);
                continue
            end
            
            fprintf('OK\n');
            
            % package useful outputs
            entry.alpha_deg = alpha_deg;
            entry.phase_delay = phase_delay;
            entry.R = scenario.R;
            entry.L = scenario.L;
            entry.simOut = simOut; %#ok<STRNU>
            
            % Post-process voltage and current from To Workspace blocks
            V_data = [];
            I_data = [];
            t_data = [];
            
            try
                % Check what's available
                baseVars = evalin('base', 'who');
                hasVout = ismember('Vout', baseVars);
                hasIout = ismember('Iout', baseVars);
                
                if hasVout && hasIout
                    % Get from base workspace
                    Vout = evalin('base', 'Vout');
                    Iout = evalin('base', 'Iout');
                    
                    % Extract data and time
                    if isstruct(Vout) && isfield(Vout, 'Data')
                        V_data = Vout.Data;
                        t_data = Vout.Time;
                    elseif isstruct(Vout) && isfield(Vout, 'signals')
                        V_data = Vout.signals.values;
                        t_data = Vout.time;
                    else
                        V_data = Vout;
                        t_data = simOut.tout;
                    end
                    
                    if isstruct(Iout) && isfield(Iout, 'Data')
                        I_data = Iout.Data;
                    elseif isstruct(Iout) && isfield(Iout, 'signals')
                        I_data = Iout.signals.values;
                    else
                        I_data = Iout;
                    end
                    
                    % Calculate average and RMS values
                    entry.Vavg = mean(V_data);
                    entry.Vrms = sqrt(mean(V_data.^2));
                    entry.Iavg = mean(I_data);
                    entry.Irms = sqrt(mean(I_data.^2));
                    
                    % Live plotting
                    if enableLivePlot && mod(ai, 6) == 1  % Plot every 6th alpha (0, 30, 60, 90...)
                        figure(figHandle);
                        subplot(2,1,1);
                        plot(t_data, V_data, 'b-', 'LineWidth', 1.5);
                        xlabel('Time (s)'); ylabel('Voltage (V)');
                        title(sprintf('%s - %s: Voltage at α=%d°', modelBase, scenario.name, alpha_deg));
                        grid on;
                        
                        subplot(2,1,2);
                        plot(t_data, I_data, 'r-', 'LineWidth', 1.5);
                        xlabel('Time (s)'); ylabel('Current (A)');
                        title(sprintf('Current at α=%d°', alpha_deg));
                        grid on;
                        drawnow;
                    end
                else
                    if alpha_deg == 0
                        fprintf('\n  WARNING: Vout/Iout not in base workspace. Available vars: %s\n', strjoin(baseVars, ', '));
                        fprintf('  Make sure To Workspace blocks are configured properly!\n');
                    end
                    entry.Vavg = NaN;
                    entry.Vrms = NaN;
                    entry.Iavg = NaN;
                    entry.Irms = NaN;
                end
            catch ME
                if alpha_deg == 0
                    fprintf('\n  ERROR extracting data: %s\n', ME.message);
                end
                entry.Vavg = NaN;
                entry.Vrms = NaN;
                entry.Iavg = NaN;
                entry.Irms = NaN;
            end
            
            % If model has logged signals in logsout or To Workspace variables,
            % they will appear inside simOut. Save the full simOut for postprocessing.
            modelResults{ai} = entry;
        end
        
        % store and save per-model result file with scenario name
        fieldName = sprintf('%s_%s', modelBase, scenario.name);
        results.(fieldName) = modelResults; %#ok<STRNU>
        saveFile = fullfile(saveFolder, sprintf('%s_%s_results.mat', modelBase, scenario.name));
        save(saveFile, 'modelResults','scenario','-v7.3');
        fprintf('Saved results to %s\n', saveFile);
        
        % Display summary statistics
        fprintf('\nSummary for %s - %s:\n', modelName, scenario.name);
        fprintf('  Alpha(deg)  Vavg(V)   Vrms(V)   Iavg(A)   Irms(A)\n');
        fprintf('  -----------------------------------------------\n');
        for ai = 1:min(5, numel(modelResults))  % Show first 5 alpha values
            if ~isempty(modelResults{ai}) && isfield(modelResults{ai}, 'Vavg')
                e = modelResults{ai};
                fprintf('  %6d     %7.2f   %7.2f   %7.3f   %7.3f\n', ...
                    e.alpha_deg, e.Vavg, e.Vrms, e.Iavg, e.Irms);
            end
        end
        if numel(modelResults) > 5
            fprintf('  ... (%d more alpha values)\n', numel(modelResults)-5);
        end
        
        % close model to keep environment clean
        try
            close_system(modelPath,0);
        catch
        end
    end
end
% Save aggregated results
save(fullfile(saveFolder,'all_results.mat'),'results','-v7.3');
fprintf('\nAll sweeps finished. Results in: %s\n', saveFolder);
fprintf('Load the MAT files and inspect `simOut` fields or logged variables.\n');
