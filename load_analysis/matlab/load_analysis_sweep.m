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

% Filter models based on modelSelection parameter
if ~strcmp(modelSelection, 'all')
    modelMap = struct('ct', 'center_taped', 'bg', 'Full_Wave_Bridge', 'hf', 'half_wave');
    
    if ischar(modelSelection)
        % Single model selection
        if isfield(modelMap, modelSelection)
            targetName = modelMap.(modelSelection);
            models = models(contains(lower(models), lower(targetName)));
        else
            error('Invalid modelSelection: %s. Use ''all'', ''ct'', ''bg'', ''hf'', or cell array.', modelSelection);
        end
    elseif iscell(modelSelection)
        % Multiple model selection
        selectedModels = {};
        for i = 1:length(modelSelection)
            sel = modelSelection{i};
            if isfield(modelMap, sel)
                targetName = modelMap.(sel);
                selectedModels = [selectedModels; models(contains(lower(models), lower(targetName)))];
            end
        end
        models = unique(selectedModels);
    end
end

if isempty(models)
    error('No models match the selection criteria: %s', mat2str(modelSelection));
end

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
                % Check simulation status
                simStatus = get_param(modelBase, 'SimulationStatus');
                fprintf('    Current simulation status: %s\n', simStatus);
                
                % If compiling or running, stop it
                if ~strcmp(simStatus, 'stopped')
                    fprintf('    Stopping simulation...\n');
                    set_param(modelBase, 'SimulationCommand', 'stop');
                    % Wait until simulation actually stops
                    timeout = 10; % seconds
                    startTime = tic;
                    while ~strcmp(get_param(modelBase, 'SimulationStatus'), 'stopped') && toc(startTime) < timeout
                        pause(0.2);
                    end
                    simStatus = get_param(modelBase, 'SimulationStatus');
                    if ~strcmp(simStatus, 'stopped')
                        warning('Simulation did not stop within timeout (status: %s)', simStatus);
                    else
                        fprintf('    Simulation stopped.\n');
                    end
                end
            catch ME
                warning('Error checking/stopping simulation: %s', ME.message);
            end
            
            % Give it a moment before closing
            pause(0.5);
            
            % Now try to close
            try
                close_system(modelBase, 0);
                fprintf('    Model closed successfully.\n');
            catch ME
                warning('Could not close model %s: %s', modelBase, ME.message);
                fprintf('    Skipping model close - will try to load anyway...\n');
            end
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
                
                % Try common variable names for To Workspace blocks
                vNames = {'Vout', 'vout', 'V_out', 'voltage', 'out'};
                iNames = {'Iout', 'iout', 'I_out', 'current', 'out'};
                
                hasVout = false;
                hasIout = false;
                VoutVar = '';
                IoutVar = '';
                
                % Find voltage variable
                for v = 1:length(vNames)
                    if ismember(vNames{v}, baseVars)
                        VoutVar = vNames{v};
                        hasVout = true;
                        break;
                    end
                end
                
                % Find current variable
                for v = 1:length(iNames)
                    if ismember(iNames{v}, baseVars)
                        IoutVar = iNames{v};
                        hasIout = true;
                        break;
                    end
                end
                
                % Check if 'out' is a structure with multiple signals
                if ismember('out', baseVars) && ~hasVout && ~hasIout
                    outData = evalin('base', 'out');
                    if isstruct(outData)
                        % Check for common field names in 'out' structure
                        fnames = fieldnames(outData);
                        if alpha_deg == 0
                            fprintf('\n  DEBUG: ''out'' structure has fields: %s\n', strjoin(fnames, ', '));
                        end
                        % You may need to adapt field names based on your model
                        if ismember('Vout', fnames) || ismember('voltage', fnames)
                            VoutVar = 'out';
                            hasVout = true;
                        end
                        if ismember('Iout', fnames) || ismember('current', fnames)
                            IoutVar = 'out';
                            hasIout = true;
                        end
                    end
                end
                
                if hasVout && hasIout
                    % Get from base workspace
                    Vout = evalin('base', VoutVar);
                    Iout = evalin('base', IoutVar);
                    
                    % Debug output for first alpha
                    if alpha_deg == 0
                        fprintf('\n  DEBUG: Found %s (class: %s)', VoutVar, class(Vout));
                        if isstruct(Vout)
                            fprintf(', fields: %s', strjoin(fieldnames(Vout), ', '));
                        end
                        fprintf('\n  DEBUG: Found %s (class: %s)', IoutVar, class(Iout));
                        if isstruct(Iout)
                            fprintf(', fields: %s', strjoin(fieldnames(Iout), ', '));
                        end
                        fprintf('\n');
                    end
                    
                    % Extract data and time
                    if isstruct(Vout) && isfield(Vout, 'Data')
                        V_data = Vout.Data;
                        t_data = Vout.Time;
                    elseif isstruct(Vout) && isfield(Vout, 'signals')
                        V_data = Vout.signals.values;
                        t_data = Vout.time;
                    elseif isstruct(Vout) && isfield(Vout, 'time') && isfield(Vout, 'signals')
                        V_data = Vout.signals(1).values;
                        t_data = Vout.time;
                    elseif isnumeric(Vout) || islogical(Vout)
                        V_data = Vout;
                        t_data = simOut.tout;
                    else
                        error('Unrecognized Vout structure');
                    end
                    
                    if isstruct(Iout) && isfield(Iout, 'Data')
                        I_data = Iout.Data;
                    elseif isstruct(Iout) && isfield(Iout, 'signals')
                        I_data = Iout.signals.values;
                    elseif isstruct(Iout) && isfield(Iout, 'time') && isfield(Iout, 'signals')
                        I_data = Iout.signals(1).values;
                    elseif isnumeric(Iout) || islogical(Iout)
                        I_data = Iout;
                    else
                        error('Unrecognized Iout structure');
                    end
                    
                    % Ensure data is numeric
                    if ~isnumeric(V_data) || ~isnumeric(I_data)
                        error('V_data or I_data is not numeric');
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
