function make_F7_startup_cross_maps() %#ok<*UNRCH>
%MAKE_F7_STARTUP_CROSS_MAPS Compute and plot all pairwise startup maps.

% ---- Manual choices: edit these lines in MATLAB ----
parameterFile = @startup_reference_parameters;
formalMode = true;            % false: 13x13 quick maps; true: 51x51 maps
observationCycles = 40;       % use 40 or more for a formal manuscript run
reuseExistingResults = false; % true only after confirming inputs are unchanged
resumeFromCheckpoint = true;  % resume after interruption, one completed y-row at a time

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot,fullfile(repoRoot,'configs'));
p = parameterFile();
outDir = fullfile(repoRoot,'results','F7_startup_cross_maps');
figDir = fullfile(repoRoot,'figures');
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end

if formalMode
    mode = "formal";
    n = 51;
    observationCycles = max(observationCycles,40);
    outputDt = min(p.solver.output_dt,0.005);
    suffix = "formal";
else
    mode = "quick";
    n = 13;
    outputDt = 0.01; % saved sampling only; the solver remains adaptive ODE45
    suffix = "quick";
end

omega0 = sqrt(p.g/p.ell);
opts = struct( ...
    'A',linspace(0.001,0.008,n), ...
    'L',linspace(0.14,0.24,n), ...
    'theta_seed_deg',linspace(0.1,2.0,n), ...
    'omega',linspace(1.4*omega0,2.6*omega0,n), ...
    'A_fixed',p.drive.amplitude, ...
    'L_fixed',p.ell, ...
    'theta_seed_fixed_deg',abs(rad2deg(p.initial_state.theta(2))), ...
    'observation_cycles',observationCycles, ...
    'mode',mode, ...
    'output_dt',outputDt, ...
    'checkpoint_dir',string(outDir), ...
    'resume',resumeFromCheckpoint, ...
    'continue_on_error',true);
dataFile = fullfile(outDir,"F7_startup_cross_maps_"+suffix+".mat");

if reuseExistingResults && isfile(dataFile)
    load(dataFile,'maps');
else
    maps = lato.scan.startup_cross_maps(p,opts);
    parameterSnapshot = p;
    save(dataFile,'maps','parameterSnapshot','-v7.3');
    writetable(maps.A_omega.table, ...
        fullfile(outDir,"F7_A_omega_"+suffix+".csv"));
    writetable(maps.L_omega.table, ...
        fullfile(outDir,"F7_L_omega_"+suffix+".csv"));
    writetable(maps.theta_seed_omega.table, ...
        fullfile(outDir,"F7_theta_seed_omega_"+suffix+".csv"));
end

[heatmapFig,surfaceFig] = lato.plot.startup_cross_maps( ...
    maps,struct('visible',"on"));
exportgraphics(heatmapFig,fullfile(figDir,'F7_startup_cross_heatmaps.png'), ...
    'Resolution',300);
exportgraphics(heatmapFig,fullfile(figDir,'F7_startup_cross_heatmaps.pdf'), ...
    'ContentType','vector');
exportgraphics(surfaceFig,fullfile(figDir,'F7_startup_cross_surfaces.png'), ...
    'Resolution',300);
exportgraphics(surfaceFig,fullfile(figDir,'F7_startup_cross_surfaces.pdf'), ...
    'ContentType','vector');

fprintf(['F7 complete: 3 x %d x %d ODE45 startup runs (%s).\n', ...
    'Finite-time responses: %d out-of-scope, %d solver failures.\n'], ...
    n,n,mode,count_field(maps,'out_of_scope'),count_field(maps,'failed'));
end

function n = count_field(maps,name)
pairs = {'A_omega','L_omega','theta_seed_omega'};
n = 0;
for k = 1:numel(pairs)
    n = n + sum(maps.(pairs{k}).(name),'all');
end
end
