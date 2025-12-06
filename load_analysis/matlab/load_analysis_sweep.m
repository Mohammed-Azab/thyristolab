run params.m

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

scriptDir = fileparts(mfilename('fullpath'));
slxDir = fullfile(scriptDir,'..','simulink');
files = dir(fullfile(slxDir,'*.slx'));
if isempty(files)
    error('No .slx models found in %s. Place your models there and re-run.', slxDir);
end

models = {files.name};

if ~strcmp(modelSelection, 'all')
    modelMap = struct('ct', 'center_taped', 'bg', 'Full_Wave_Bridge', 'hf', 'half_wave');
    
    if ischar(modelSelection)
        if isfield(modelMap, modelSelection)
            targetName = modelMap.(modelSelection);
            models = models(contains(lower(models), lower(targetName)));
        else
            error('Invalid modelSelection: %s. Use ''all'', ''ct'', ''bg'', ''hf'', or cell array.', modelSelection);
        end
    elseif iscell(modelSelection)
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

results = struct();

for si = 1:numel(scenarios)
    scenario = scenarios(si);
    fprintf('\n========== Scenario %d: %s ==========\n', si, scenario.name);
    
    for mi = 1:numel(models)
        modelName = models{mi};
        modelPath = fullfile(slxDir, modelName);
        [~,modelBase] = fileparts(modelName);
        
        fprintf('\nRunning model: %s\n', modelName);
        
        if bdIsLoaded(modelBase)
            fprintf('  Closing existing model instance...\n');
            try
                simStatus = get_param(modelBase, 'SimulationStatus');
                fprintf('    Current simulation status: %s\n', simStatus);
                
                if ~strcmp(simStatus, 'stopped')
                    fprintf('    Stopping simulation...\n');
                    set_param(modelBase, 'SimulationCommand', 'stop');
                    timeout = 10;
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
            
            pause(0.5);
            
            try
                close_system(modelBase, 0);
                fprintf('    Model closed successfully.\n');
            catch ME
                warning('Could not close model %s: %s', modelBase, ME.message);
                fprintf('    Skipping model close - will try to load anyway...\n');
            end
        end
        
        load_system(modelPath);
        
        modelResults = cell(numel(alphas_deg),1);
        
        for ai = 1:numel(alphas_deg)
            alpha_deg = alphas_deg(ai);
            phase_delay = (alpha_deg/360) * (1/f);
            phase_delay2 = phase_delay + (1/(2*f));
            
            assignin('base','alpha_deg',alpha_deg);
            assignin('base','phase_delay',phase_delay);
            assignin('base','f',f);
            assignin('base','Vrms',Vrms);
            assignin('base','R',scenario.R);
            assignin('base','L',scenario.L);
            assignin('base','pulse_amplitude',pulse_amplitude);
            assignin('base','pulse_width',pulse_width);
            assignin('base','pulse_period',pulse_period);
            
            evalin('base', 'clear Vout Iout V_data I_data vout iout');
            
            fprintf('  alpha = %3d deg -> phase_delay = %.6f s ... ', alpha_deg, phase_delay);
            
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
            
            entry.alpha_deg = alpha_deg;
            entry.phase_delay = phase_delay;
            entry.R = scenario.R;
            entry.L = scenario.L;
            entry.simOut = simOut;
            
            V_data = [];
            I_data = [];
            t_data = simOut.tout;
            
            try
                hasVout = false;
                hasIout = false;
                
                if isprop(simOut, 'Vout') && isprop(simOut, 'Iout')
                    Vout = simOut.Vout;
                    Iout = simOut.Iout;
                    
                    if ~isempty(Vout) && ~isempty(Iout)
                        hasVout = true;
                        hasIout = true;
                    end
                end
                
                if ~hasVout || ~hasIout
                    baseVars = evalin('base', 'who');
                    
                    VoutVar = '';
                    IoutVar = '';
                    
                    vNames = {'V_data', 'Vout', 'vout', 'V_out', 'voltage'};
                    iNames = {'I_data', 'Iout', 'iout', 'I_out', 'current'};
                    
                    for v = 1:length(vNames)
                        if ismember(vNames{v}, baseVars)
                            VoutVar = vNames{v};
                            Vout = evalin('base', VoutVar);
                            hasVout = true;
                            break;
                        end
                    end
                    
                    for v = 1:length(iNames)
                        if ismember(iNames{v}, baseVars)
                            IoutVar = iNames{v};
                            Iout = evalin('base', IoutVar);
                            hasIout = true;
                            break;
                        end
                    end
                end
                
                if hasVout && hasIout
                    if isa(Vout, 'timeseries')
                        V_data = Vout.Data;
                        t_data = Vout.Time;
                    elseif isnumeric(Vout) || islogical(Vout)
                        if ~isvector(Vout) && size(Vout, 2) > 1
                            V_data = Vout(:, end);
                        else
                            V_data = Vout(:);
                        end
                        t_data = simOut.tout;
                    elseif isa(Vout, 'Simulink.SimulationOutput')
                        props = properties(Vout);
                        if ismember('logsout', props) && ~isempty(Vout.logsout)
                            V_data = Vout.logsout.getElement(1).Values.Data;
                            t_data = Vout.logsout.getElement(1).Values.Time;
                        elseif ismember('yout', props) && ~isempty(Vout.yout)
                            V_data = Vout.yout.Data;
                            t_data = Vout.yout.Time;
                        else
                            varNames = who(Vout);
                            if ~isempty(varNames)
                                V_data = Vout.(varNames{1});
                                if isstruct(V_data) && isfield(V_data, 'Data')
                                    t_data = V_data.Time;
                                    V_data = V_data.Data;
                                elseif isstruct(V_data) && isfield(V_data, 'signals')
                                    t_data = V_data.time;
                                    V_data = V_data.signals.values;
                                else
                                    t_data = simOut.tout;
                                end
                            else
                                error('Cannot find data in Vout SimulationOutput');
                            end
                        end
                    elseif isstruct(Vout) && isfield(Vout, 'Data')
                        V_data = Vout.Data;
                        t_data = Vout.Time;
                    elseif isstruct(Vout) && isfield(Vout, 'signals')
                        V_data = Vout.signals.values;
                        t_data = Vout.time;
                    elseif isstruct(Vout) && isfield(Vout, 'time') && isfield(Vout, 'signals')
                        V_data = Vout.signals(1).values;
                        t_data = Vout.time;
                    else
                        error('Unrecognized Vout structure');
                    end
                    
                    if isa(Iout, 'timeseries')
                        I_data = Iout.Data;
                    elseif isnumeric(Iout) || islogical(Iout)
                        if ~isvector(Iout) && size(Iout, 2) > 1
                            I_data = Iout(:, end);
                        else
                            I_data = Iout(:);
                        end
                    elseif isa(Iout, 'Simulink.SimulationOutput')
                        props = properties(Iout);
                        if ismember('logsout', props) && ~isempty(Iout.logsout)
                            I_data = Iout.logsout.getElement(1).Values.Data;
                        elseif ismember('yout', props) && ~isempty(Iout.yout)
                            I_data = Iout.yout.Data;
                        else
                            varNames = who(Iout);
                            if ~isempty(varNames)
                                I_data = Iout.(varNames{1});
                                if isstruct(I_data) && isfield(I_data, 'Data')
                                    I_data = I_data.Data;
                                elseif isstruct(I_data) && isfield(I_data, 'signals')
                                    I_data = I_data.signals.values;
                                end
                            else
                                error('Cannot find data in Iout SimulationOutput');
                            end
                        end
                    elseif isstruct(Iout) && isfield(Iout, 'Data')
                        I_data = Iout.Data;
                    elseif isstruct(Iout) && isfield(Iout, 'signals')
                        I_data = Iout.signals.values;
                    elseif isstruct(Iout) && isfield(Iout, 'time') && isfield(Iout, 'signals')
                        I_data = Iout.signals(1).values;
                    else
                        error('Unrecognized Iout structure');
                    end
                    
                    if ~isnumeric(V_data) || ~isnumeric(I_data)
                        error('V_data or I_data is not numeric');
                    end
                    if ~isnumeric(t_data) || isempty(t_data)
                        error('t_data is not a valid numeric vector');
                    end
                    
                    V_data = V_data(:);
                    I_data = I_data(:);
                    t_data = t_data(:);
                    
                    if length(t_data) ~= length(V_data) || length(t_data) ~= length(I_data)
                        error('Data length mismatch: t_data(%d), V_data(%d), I_data(%d)', ...
                            length(t_data), length(V_data), length(I_data));
                    end
                    
                    entry.Vavg = mean(V_data);
                    entry.Vrms = sqrt(mean(V_data.^2));
                    entry.Iavg = mean(I_data);
                    entry.Irms = sqrt(mean(I_data.^2));
                    
                    if enableLivePlot && mod(ai, 6) == 1
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
                    warning('Vout/Iout not found in workspace for %s at alpha=%d', modelBase, alpha_deg);
                    entry.Vavg = NaN;
                    entry.Vrms = NaN;
                    entry.Iavg = NaN;
                    entry.Irms = NaN;
                end
            catch ME
                warning('Error extracting data for %s at alpha=%d: %s', modelBase, alpha_deg, ME.message);
                entry.Vavg = NaN;
                entry.Vrms = NaN;
                entry.Iavg = NaN;
                entry.Irms = NaN;
            end
            
            modelResults{ai} = entry;
        end
        
        fieldName = sprintf('%s_%s', modelBase, scenario.name);
        results.(fieldName) = modelResults;
        saveFile = fullfile(saveFolder, sprintf('%s_%s_results.mat', modelBase, scenario.name));
        save(saveFile, 'modelResults','scenario','-v7.3');
        fprintf('Saved results to %s\n', saveFile);
        
        fprintf('\nSummary for %s - %s:\n', modelName, scenario.name);
        fprintf('  Alpha(deg)  Vavg(V)   Vrms(V)   Iavg(A)   Irms(A)\n');
        fprintf('  -----------------------------------------------\n');
        for ai = 1:min(5, numel(modelResults))
            if ~isempty(modelResults{ai}) && isfield(modelResults{ai}, 'Vavg')
                e = modelResults{ai};
                fprintf('  %6d     %7.2f   %7.2f   %7.3f   %7.3f\n', ...
                    e.alpha_deg, e.Vavg, e.Vrms, e.Iavg, e.Irms);
            end
        end
        if numel(modelResults) > 5
            fprintf('  ... (%d more alpha values)\n', numel(modelResults)-5);
        end
        
        if generateGraphs
            fprintf('\nGenerating graphs for selected alpha values...\n');
            figFolder = fullfile(saveFolder, '..', '..', 'figures', modelBase, scenario.name);
            if ~exist(figFolder, 'dir')
                mkdir(figFolder);
            end
            
            for alpha_val = graphAlphas
                resultIdx = -1;
                for ai = 1:numel(modelResults)
                    if ~isempty(modelResults{ai}) && isfield(modelResults{ai}, 'alpha_deg')
                        if modelResults{ai}.alpha_deg == alpha_val
                            resultIdx = ai;
                            break;
                        end
                    end
                end
                
                if resultIdx > 0 && isfield(modelResults{resultIdx}, 'simOut')
                    result = modelResults{resultIdx};
                    simOut = result.simOut;
                    
                    V_plot = [];
                    I_plot = [];
                    t_plot = simOut.tout;
                    
                    if isprop(simOut, 'Vout') && isprop(simOut, 'Iout')
                        VoutData = simOut.Vout;
                        IoutData = simOut.Iout;
                        if isa(VoutData, 'timeseries')
                            V_plot = VoutData.Data;
                            t_plot = VoutData.Time;
                        elseif isnumeric(VoutData) && ~isempty(VoutData)
                            if size(VoutData, 2) > 1
                                V_plot = VoutData(:, end);
                            else
                                V_plot = VoutData(:);
                            end
                        end
                        
                        if isa(IoutData, 'timeseries')
                            I_plot = IoutData.Data;
                        elseif isnumeric(IoutData) && ~isempty(IoutData)
                            if size(IoutData, 2) > 1
                                I_plot = IoutData(:, end);
                            else
                                I_plot = IoutData(:);
                            end
                        end
                    end
                    
                    if ~isempty(V_plot) && ~isempty(I_plot)
                        fig = figure('Name', sprintf('%s - %s - Alpha=%d', modelBase, scenario.name, alpha_val), ...
                                     'Position', [100, 100, 1200, 600]);
                        
                        subplot(2, 1, 1);
                        ax1 = gca;
                        ax1.Toolbar.Visible = 'off';
                        plot(t_plot, V_plot, 'b-', 'LineWidth', 1.5);
                        grid on;
                        xlabel('Time (s)');
                        ylabel('Voltage (V)');
                        title(sprintf('%s - %s: Voltage Waveform (\\alpha = %d°)', ...
                            strrep(modelBase, '_', ' '), strrep(scenario.name, '_', ' '), alpha_val));
                        
                        subplot(2, 1, 2);
                        ax2 = gca;
                        ax2.Toolbar.Visible = 'off';
                        plot(t_plot, I_plot, 'r-', 'LineWidth', 1.5);
                        grid on;
                        xlabel('Time (s)');
                        ylabel('Current (A)');
                        title(sprintf('Current Waveform (\\alpha = %d°)', alpha_val));
                        
                        figName = sprintf('%s_%s_alpha_%d', modelBase, scenario.name, alpha_val);
                        saveas(fig, fullfile(figFolder, [figName '.png']));
                        saveas(fig, fullfile(figFolder, [figName '.fig']));
                        fprintf('  Saved graph for alpha=%d to %s\n', alpha_val, figFolder);
                        close(fig);
                    else
                        fprintf('  Warning: No data available for alpha=%d\n', alpha_val);
                    end
                else
                    fprintf('  Warning: Result not found for alpha=%d\n', alpha_val);
                end
            end
        end
        
        try
            close_system(modelPath,0);
        catch
        end
    end
end

save(fullfile(saveFolder,'all_results.mat'),'results','-v7.3');
fprintf('\nAll sweeps finished. Results in: %s\n', saveFolder);
fprintf('Load the MAT files and inspect `simOut` fields or logged variables.\n');
