function maps = startup_cross_maps(p, opts)
%STARTUP_CROSS_MAPS Pairwise finite-time startup maps using the ODE45 core.
%   MAPS = lato.scan.startup_cross_maps(P, OPTS) computes A-Omega,
%   L-Omega and theta_seed-Omega maps.  These are collision-suppressed
%   startup-reference diagnostics: the plotted amplitude is the largest
%   response reached in a fixed number of drive cycles, not an assumed
%   steady amplitude or a geometrically realized two-ball pre-impact path.
%
%   Formal mode requires at least 51 samples on every varied axis.  The
%   returned structures retain explicit steady/transient, out-of-scope and
%   solver-failure masks so downstream plots cannot silently treat a
%   transient or a trajectory beyond 90 degrees as a steady result.

arguments
    p (1,1) struct
    opts (1,1) struct = struct()
end

omega0 = sqrt(p.g/p.ell);
defaults = struct( ...
    'A', linspace(0.001,0.008,51), ...
    'L', linspace(0.14,0.24,51), ...
    'theta_seed_deg', linspace(0.1,2.0,51), ...
    'omega', linspace(1.4*omega0,2.6*omega0,51), ...
    'omega_A', [], ...
    'omega_L', [], ...
    'omega_seed', [], ...
    'A_fixed', p.drive.amplitude, ...
    'L_fixed', p.ell, ...
    'theta_seed_fixed_deg', abs(rad2deg(p.initial_state.theta(2))), ...
    'observation_cycles', 40, ...
    'observation_duration_s', NaN, ...
    'terminal_steady_window_cycles', 10, ...
    'terminal_steady_relative_mean_deviation_tolerance', 0.10, ...
    'mode', "quick", ...
    'output_dt', p.solver.output_dt, ...
    'checkpoint_dir', "", ...
    'resume', true, ...
    'continue_on_error', true, ...
    'pairs', ["A","L","theta_seed_deg"], ...
    'reuse_maps', struct());
opts = lato.scan.options(opts, defaults);
opts.mode = string(opts.mode);
opts.checkpoint_dir = string(opts.checkpoint_dir);
mustBeMember(opts.mode,["quick","formal"]);
opts.pairs = reshape(string(opts.pairs),1,[]);
mustBeMember(opts.pairs,["A","L","theta_seed_deg"]);
if numel(unique(opts.pairs)) ~= numel(opts.pairs)
    error('lato:scan:DuplicatePair','opts.pairs must not contain duplicates.');
end

if isempty(opts.omega_A), opts.omega_A = opts.omega; end
if isempty(opts.omega_L), opts.omega_L = opts.omega; end
if isempty(opts.omega_seed), opts.omega_seed = opts.omega; end

validate_scan_axis(opts.A,'A',0);
validate_scan_axis(opts.L,'L',realmin);
validate_scan_axis(opts.theta_seed_deg,'theta_seed_deg',realmin);
validate_scan_axis(opts.omega,'omega',realmin);
validate_scan_axis(opts.omega_A,'omega_A',realmin);
validate_scan_axis(opts.omega_L,'omega_L',realmin);
validate_scan_axis(opts.omega_seed,'omega_seed',realmin);

if opts.mode == "formal"
    formalCounts = zeros(0,2);
    if any(opts.pairs=="A")
        formalCounts(end+1,:) = [numel(opts.A),numel(opts.omega_A)]; %#ok<AGROW>
    end
    if any(opts.pairs=="L")
        formalCounts(end+1,:) = [numel(opts.L),numel(opts.omega_L)]; %#ok<AGROW>
    end
    if any(opts.pairs=="theta_seed_deg")
        formalCounts(end+1,:) = ...
            [numel(opts.theta_seed_deg),numel(opts.omega_seed)]; %#ok<AGROW>
    end
    if any(formalCounts < 51,'all')
        error('lato:scan:FormalResolution', ...
            'Formal startup cross maps require >=51 points on each selected axis.');
    end
end
if opts.observation_cycles < 1 || fix(opts.observation_cycles) ~= opts.observation_cycles
    error('lato:scan:InvalidCycles', ...
        'observation_cycles must be a positive integer.');
end
if opts.terminal_steady_window_cycles < 1 || ...
        fix(opts.terminal_steady_window_cycles) ~= ...
        opts.terminal_steady_window_cycles
    error('lato:scan:InvalidSteadyWindow', ...
        'terminal_steady_window_cycles must be a positive integer.');
end
if ~isscalar(opts.terminal_steady_relative_mean_deviation_tolerance) || ...
        ~isfinite(opts.terminal_steady_relative_mean_deviation_tolerance) || ...
        opts.terminal_steady_relative_mean_deviation_tolerance < 0
    error('lato:scan:InvalidSteadyTolerance', ...
        'terminal steady relative mean deviation tolerance must be finite and nonnegative.');
end
if ~(isnan(opts.observation_duration_s) || ...
        (isscalar(opts.observation_duration_s) && ...
         isfinite(opts.observation_duration_s) && opts.observation_duration_s>0))
    error('lato:scan:InvalidDuration', ...
        'observation_duration_s must be NaN or a positive finite scalar.');
end
if opts.output_dt <= 0
    error('lato:scan:InvalidOutputStep','output_dt must be positive.');
end
if strlength(opts.checkpoint_dir) > 0 && ~isfolder(opts.checkpoint_dir)
    mkdir(opts.checkpoint_dir);
end

maps = struct();
maps.kind = "startup_cross_maps";
maps.mode = opts.mode;
maps.g = p.g;
maps.beta = p.beta;
maps.omega0_reference = omega0;
maps.options = rmfield(opts,'reuse_maps');
maps.resume_enabled = logical(opts.resume);

if any(opts.pairs=="A")
    maps.A_omega = run_pair("A", opts.A, opts.omega_A);
    save_checkpoint('F7_A_omega_checkpoint.mat','A_omega',maps.A_omega);
end
if any(opts.pairs=="L")
    maps.L_omega = run_pair("L", opts.L, opts.omega_L);
    save_checkpoint('F7_L_omega_checkpoint.mat','L_omega',maps.L_omega);
end
if any(opts.pairs=="theta_seed_deg")
    maps.theta_seed_omega = run_pair( ...
        "theta_seed_deg", opts.theta_seed_deg, opts.omega_seed);
    save_checkpoint('F7_theta_seed_omega_checkpoint.mat', ...
        'theta_seed_omega', maps.theta_seed_omega);
end

    function pair = run_pair(parameter, yValues, omegaValues)
        ny = numel(yValues);
        nx = numel(omegaValues);
        maximum = nan(ny,nx);
        mu = nan(ny,nx);
        time45 = nan(ny,nx);
        cycles = zeros(ny,nx);
        steady = false(ny,nx);
        terminalSteady = false(ny,nx);
        terminalSteadyAmplitude = nan(ny,nx);
        terminalSteadyRelativeDeviation = nan(ny,nx);
        large = false(ny,nx);
        outOfScope = false(ny,nx);
        completedInterval = false(ny,nx);
        incomplete = false(ny,nx);
        failed = false(ny,nx);
        status = strings(ny,nx);
        records = repmat(record_template(),ny*nx,1);
        reuseMask = false(ny,nx);
        reusedPoints = 0;
        [maximum,mu,time45,cycles,steady,terminalSteady, ...
            terminalSteadyAmplitude,terminalSteadyRelativeDeviation, ...
            large,outOfScope, ...
            completedInterval,incomplete,failed,status,records, ...
            reuseMask,reusedPoints] = reuse_existing_pair(parameter, ...
            yValues,omegaValues,maximum,mu,time45,cycles,steady, ...
            terminalSteady,terminalSteadyAmplitude, ...
            terminalSteadyRelativeDeviation,large, ...
            outOfScope,completedInterval,incomplete,failed,status,records, ...
            reuseMask);
        completedRows = false(ny,1);
        resumedRows = 0;
        progressFile = row_checkpoint_file(parameter);
        signature = checkpoint_signature(parameter,yValues,omegaValues);
        if opts.resume && strlength(progressFile)>0 && isfile(progressFile)
            try
                loaded = load(progressFile,'progress');
                if isfield(loaded,'progress') && ...
                        isequaln(loaded.progress.signature,signature)
                    progress = loaded.progress;
                    maximum = progress.maximum;
                    mu = progress.mu;
                    time45 = progress.time45;
                    cycles = progress.cycles;
                    steady = progress.steady;
                    terminalSteady = progress.terminal_steady;
                    terminalSteadyAmplitude = ...
                        progress.terminal_steady_amplitude_deg;
                    terminalSteadyRelativeDeviation = ...
                        progress.terminal_steady_relative_mean_deviation;
                    large = progress.large;
                    outOfScope = progress.out_of_scope;
                    completedInterval = progress.completed_requested_interval;
                    incomplete = progress.incomplete;
                    failed = progress.failed;
                    status = progress.status;
                    records = progress.records;
                    completedRows = progress.completed_rows;
                    resumedRows = sum(completedRows);
                    fprintf('F7 %s-Omega: resumed %d/%d completed rows.\n', ...
                        parameter,resumedRows,ny);
                else
                    warning('lato:scan:CheckpointMismatch', ...
                        ['Ignoring %s because axes, fixed parameters or ', ...
                         'solver settings changed.'],progressFile);
                end
            catch ME
                warning('lato:scan:CheckpointUnreadable', ...
                    'Ignoring unreadable checkpoint %s: %s',progressFile,ME.message);
            end
        end

        q = 0;
        for iy = 1:ny
            if completedRows(iy)
                q = q + nx;
                continue
            end
            for ix = 1:nx
                q = q + 1;
                if reuseMask(iy,ix)
                    continue
                end
                omega = omegaValues(ix);
                A = opts.A_fixed;
                L = opts.L_fixed;
                seedDeg = opts.theta_seed_fixed_deg;
                switch parameter
                    case "A"
                        A = yValues(iy);
                    case "L"
                        L = yValues(iy);
                    case "theta_seed_deg"
                        seedDeg = yValues(iy);
                end
                duration = opts.observation_cycles * 2*pi / omega;
                if isfinite(opts.observation_duration_s)
                    duration = opts.observation_duration_s;
                end
                pk = lato.scan.configure_point(p,omega,A,L,duration);
                pk.initial_state.theta = deg2rad([-seedDeg,seedDeg]);
                pk.initial_state.omega = [0,0];
                pk.initial_state.contact = false;
                pk.initial_state.slack = [false,false];
                pk.solver.output_dt = opts.output_dt;

                rec = record_template();
                rec.parameter = parameter;
                rec.parameter_value = yValues(iy);
                rec.omega = omega;
                rec.A = A;
                rec.L = L;
                rec.theta_seed_deg = seedDeg;
                rec.omega_ratio_local = omega/sqrt(pk.g/L);
                rec.h = A*omega^2/pk.g;
                try
                    result = lato.simulate_hybrid(pk);
                    summary = lato.scan.summarize_result(result, ...
                        'Omega',omega,'A',A,'L',L);
                    rec.max_amplitude_deg = summary.max_amplitude_deg;
                    rec.envelope_growth_rate_per_s = envelope_growth_rate(result);
                    rec.time_to_45_deg_s = first_crossing_time(result,45);
                    rec.completed_cycles = numel(result.cycle.index);
                    rec.steady = summary.steady;
                    [rec.terminal_steady, ...
                        rec.terminal_steady_amplitude_deg, ...
                        rec.terminal_steady_relative_mean_deviation, ...
                        terminalCyclePeaksDeg] = terminal_steady_amplitude( ...
                        result,opts.terminal_steady_window_cycles, ...
                        opts.terminal_steady_relative_mean_deviation_tolerance);
                    rec.terminal_cycle_peaks_deg = {terminalCyclePeaksDeg};
                    rec.large_amplitude = summary.large_amplitude;
                    rec.out_of_scope = summary.out_of_scope;
                    rec.status = summary.status;
                    rec.completed_requested_interval = logical( ...
                        result.termination.completed_requested_interval);
                    rec.incomplete = ~rec.completed_requested_interval && ...
                        result.termination.reason ~= "angle_limit";
                    if rec.incomplete
                        rec.steady = false;
                        rec.status = "incomplete";
                    end
                    if ~rec.completed_requested_interval || rec.out_of_scope
                        rec.terminal_steady = false;
                        rec.terminal_steady_amplitude_deg = NaN;
                    end
                catch ME
                    rec.status = "solver_error";
                    rec.failed = true;
                    if ~opts.continue_on_error
                        rethrow(ME);
                    end
                    warning('lato:scan:StartupPointFailed', ...
                        '%s=%g, omega=%g failed: %s', ...
                        parameter,yValues(iy),omega,ME.message);
                end

                records(q) = rec;
                maximum(iy,ix) = rec.max_amplitude_deg;
                mu(iy,ix) = rec.envelope_growth_rate_per_s;
                time45(iy,ix) = rec.time_to_45_deg_s;
                cycles(iy,ix) = rec.completed_cycles;
                steady(iy,ix) = rec.steady;
                terminalSteady(iy,ix) = rec.terminal_steady;
                terminalSteadyAmplitude(iy,ix) = ...
                    rec.terminal_steady_amplitude_deg;
                terminalSteadyRelativeDeviation(iy,ix) = ...
                    rec.terminal_steady_relative_mean_deviation;
                large(iy,ix) = rec.large_amplitude;
                outOfScope(iy,ix) = rec.out_of_scope;
                completedInterval(iy,ix) = rec.completed_requested_interval;
                incomplete(iy,ix) = rec.incomplete;
                failed(iy,ix) = rec.failed;
                status(iy,ix) = rec.status;
            end
            completedRows(iy) = true;
            save_row_progress(progressFile,signature,maximum,mu,time45, ...
                cycles,steady,terminalSteady,terminalSteadyAmplitude, ...
                terminalSteadyRelativeDeviation,large,outOfScope, ...
                completedInterval,incomplete,failed,status,records,completedRows);
        end

        pair = struct();
        pair.kind = "finite_time_startup_pair";
        pair.parameter = parameter;
        pair.parameter_values = yValues;
        pair.omega = omegaValues;
        pair.max_amplitude_deg = maximum;
        pair.envelope_growth_rate_per_s = mu;
        pair.time_to_45_deg_s = time45;
        pair.completed_cycles = cycles;
        pair.steady = steady;
        pair.terminal_steady = terminalSteady;
        pair.terminal_steady_amplitude_deg = terminalSteadyAmplitude;
        pair.terminal_steady_relative_mean_deviation = ...
            terminalSteadyRelativeDeviation;
        pair.transient = ~steady & completedInterval & ~outOfScope & ~failed;
        pair.large_amplitude = large;
        pair.out_of_scope = outOfScope;
        pair.completed_requested_interval = completedInterval;
        pair.incomplete = incomplete;
        pair.failed = failed;
        pair.valid = completedInterval & ~outOfScope & ~failed & isfinite(maximum);
        pair.status = status;
        pair.table = struct2table(records);
        pair.completed_rows = completedRows;
        pair.resumed_rows = resumedRows;
        pair.reused_points = reusedPoints;
        pair.progress_file = progressFile;
    end

    function [maximum,mu,time45,cycles,steady,terminalSteady, ...
            terminalSteadyAmplitude,terminalSteadyRelativeDeviation, ...
            large,outOfScope, ...
            completedInterval,incomplete,failed,status,records,reuseMask, ...
            reusedPoints] = reuse_existing_pair(parameter,yValues,omegaValues, ...
            maximum,mu,time45,cycles,steady,terminalSteady, ...
            terminalSteadyAmplitude,terminalSteadyRelativeDeviation, ...
            large,outOfScope, ...
            completedInterval,incomplete,failed,status,records,reuseMask)
        reusedPoints = 0;
        if ~isstruct(opts.reuse_maps) || isempty(fieldnames(opts.reuse_maps))
            return
        end
        fieldName = pair_field_name(parameter);
        if ~isfield(opts.reuse_maps,fieldName)
            return
        end
        old = opts.reuse_maps.(fieldName);
        requiredFields = {'parameter_values','omega','max_amplitude_deg', ...
            'envelope_growth_rate_per_s','time_to_45_deg_s', ...
            'completed_cycles','steady','terminal_steady', ...
            'terminal_steady_amplitude_deg', ...
            'terminal_steady_relative_mean_deviation', ...
            'large_amplitude','out_of_scope', ...
            'completed_requested_interval','incomplete','failed','status','table'};
        if ~all(isfield(old,requiredFields))
            warning('lato:scan:ReuseMapIncomplete', ...
                'Ignoring incomplete reusable data for %s.',fieldName);
            return
        end
        oldNx = numel(old.omega);
        for iy = 1:numel(yValues)
            oldIy = matching_index(old.parameter_values,yValues(iy));
            if isempty(oldIy), continue; end
            for ix = 1:numel(omegaValues)
                oldIx = matching_index(old.omega,omegaValues(ix));
                if isempty(oldIx), continue; end
                q = (iy-1)*numel(omegaValues)+ix;
                oldQ = (oldIy-1)*oldNx+oldIx;
                maximum(iy,ix) = old.max_amplitude_deg(oldIy,oldIx);
                mu(iy,ix) = old.envelope_growth_rate_per_s(oldIy,oldIx);
                time45(iy,ix) = old.time_to_45_deg_s(oldIy,oldIx);
                cycles(iy,ix) = old.completed_cycles(oldIy,oldIx);
                steady(iy,ix) = old.steady(oldIy,oldIx);
                terminalSteady(iy,ix) = old.terminal_steady(oldIy,oldIx);
                terminalSteadyAmplitude(iy,ix) = ...
                    old.terminal_steady_amplitude_deg(oldIy,oldIx);
                terminalSteadyRelativeDeviation(iy,ix) = ...
                    old.terminal_steady_relative_mean_deviation(oldIy,oldIx);
                large(iy,ix) = old.large_amplitude(oldIy,oldIx);
                outOfScope(iy,ix) = old.out_of_scope(oldIy,oldIx);
                completedInterval(iy,ix) = ...
                    old.completed_requested_interval(oldIy,oldIx);
                incomplete(iy,ix) = old.incomplete(oldIy,oldIx);
                failed(iy,ix) = old.failed(oldIy,oldIx);
                status(iy,ix) = old.status(oldIy,oldIx);
                records(q) = table2struct(old.table(oldQ,:));
                reuseMask(iy,ix) = true;
                reusedPoints = reusedPoints + 1;
            end
        end
        fprintf('F7 %s-Omega: reused %d existing points.\n', ...
            parameter,reusedPoints);
    end

    function name = pair_field_name(parameter)
        switch string(parameter)
            case "A"
                name = 'A_omega';
            case "L"
                name = 'L_omega';
            otherwise
                name = 'theta_seed_omega';
        end
    end

    function index = matching_index(values,target)
        tolerance = 5e-11*max(1,abs(target));
        index = find(abs(values-target)<=tolerance,1,'first');
    end

    function file = row_checkpoint_file(parameter)
        file = "";
        if strlength(opts.checkpoint_dir)==0
            return
        end
        safeName = replace(string(parameter),"theta_seed_deg","theta_seed");
        file = string(fullfile(opts.checkpoint_dir, ...
            "F7_"+safeName+"_omega_progress.mat"));
    end

    function signature = checkpoint_signature(parameter,yValues,omegaValues)
        signature = struct();
        signature.version = 6;
        signature.parameter = string(parameter);
        signature.parameter_values = reshape(yValues,1,[]);
        signature.omega = reshape(omegaValues,1,[]);
        signature.A_fixed = opts.A_fixed;
        signature.L_fixed = opts.L_fixed;
        signature.theta_seed_fixed_deg = opts.theta_seed_fixed_deg;
        signature.observation_cycles = opts.observation_cycles;
        signature.observation_duration_s = opts.observation_duration_s;
        signature.terminal_steady_window_cycles = ...
            opts.terminal_steady_window_cycles;
        signature.terminal_steady_relative_mean_deviation_tolerance = ...
            opts.terminal_steady_relative_mean_deviation_tolerance;
        signature.output_dt = opts.output_dt;
        signature.g = p.g;
        signature.mass = p.mass;
        signature.radius = p.radius;
        signature.ell_reference = p.ell;
        signature.beta = p.beta;
        signature.drive_function = string(func2str(p.drive.function));
        signature.drive_phase0 = p.drive.phase0;
        signature.air_enabled = logical(p.air.enabled);
        signature.air_model = string(p.air.model);
        signature.air_density = p.air.density;
        signature.air_viscosity = p.air.dynamic_viscosity;
        signature.air_constant_cd = p.air.constant_cd;
        signature.collision_enabled = logical(p.events.enable_collision);
        signature.slack_enabled = logical(p.events.enable_slack);
        signature.t_start = p.solver.t_start;
        signature.relative_tolerance = p.solver.relative_tolerance;
        signature.absolute_tolerance = p.solver.absolute_tolerance;
        signature.maximum_step = p.solver.maximum_step;
        signature.stop_at_angle_limit = logical(p.solver.stop_at_angle_limit);
        signature.angle_limit = p.solver.angle_limit;
        signature.large_angle = p.analysis.large_angle;
        signature.maximum_interpretable_angle = ...
            p.analysis.maximum_interpretable_angle;
        signature.stability = p.stability;
    end

    function save_row_progress(progressFile,signature,maximum,mu,time45, ...
            cycles,steady,terminalSteady,terminalSteadyAmplitude, ...
            terminalSteadyRelativeDeviation,large,outOfScope, ...
            completedInterval,incomplete,failed,status,records,completedRows)
        if strlength(progressFile)==0
            return
        end
        progress = struct();
        progress.signature = signature;
        progress.maximum = maximum;
        progress.mu = mu;
        progress.time45 = time45;
        progress.cycles = cycles;
        progress.steady = steady;
        progress.terminal_steady = terminalSteady;
        progress.terminal_steady_amplitude_deg = terminalSteadyAmplitude;
        progress.terminal_steady_relative_mean_deviation = ...
            terminalSteadyRelativeDeviation;
        progress.large = large;
        progress.out_of_scope = outOfScope;
        progress.completed_requested_interval = completedInterval;
        progress.incomplete = incomplete;
        progress.failed = failed;
        progress.status = status;
        progress.records = records;
        progress.completed_rows = completedRows;
        temporary = progressFile + ".tmp";
        save(temporary,'progress','-v7');
        [ok,message] = movefile(temporary,progressFile,'f');
        if ~ok
            error('lato:scan:CheckpointWriteFailed', ...
                'Could not replace checkpoint %s: %s',progressFile,message);
        end
    end

    function save_checkpoint(fileName, variableName, value)
        if strlength(opts.checkpoint_dir) == 0
            return
        end
        payload = struct();
        payload.(variableName) = value;
        save(fullfile(opts.checkpoint_dir,fileName),'-struct','payload','-v7.3');
    end
end

function mu = envelope_growth_rate(result)
% A finite-time peak-envelope slope, not a periodic-orbit Floquet exponent.
peak = result.cycle.peak_half_opening(:);
time = result.cycle.mid_time(:);
valid = isfinite(peak) & peak > 1e-10 & isfinite(time);
peak = peak(valid);
time = time(valid);
mu = NaN;
if numel(peak) < 6
    return
end
nWindow = min(5,floor(numel(peak)/2));
first = 1:nWindow;
last = (numel(peak)-nWindow+1):numel(peak);
dt = mean(time(last))-mean(time(first));
if dt <= 0
    return
end
mu = (mean(log(peak(last)))-mean(log(peak(first))))/dt;
end

function value = first_crossing_time(result, thresholdDeg)
amplitudeDeg = rad2deg(result.response.half_opening_angle(:));
index = find(amplitudeDeg >= thresholdDeg,1,'first');
if isempty(index)
    value = NaN;
else
    value = result.time(index)-result.time(1);
end
end

function [isSteady,amplitudeDeg,relativeMeanDeviation,recentPeakDeg] = ...
        terminal_steady_amplitude(result,windowCycles,tolerance)
% At the cycle cap, average the final window only when its mean absolute
% relative deviation satisfies the requested tolerance. Preserve that
% window so a later plotting pass can audit or change the criterion.
peak = result.cycle.peak_half_opening(:);
peak = peak(isfinite(peak));
isSteady = false;
amplitudeDeg = NaN;
relativeMeanDeviation = NaN;
recentPeakDeg = nan(1,windowCycles);
if numel(peak) < windowCycles
    return
end
recent = peak(end-windowCycles+1:end);
recentPeakDeg = reshape(rad2deg(recent),1,[]);
meanPeak = mean(recent);
scale = max(abs(meanPeak),sqrt(eps));
relativeMeanDeviation = mean(abs(recent-meanPeak))/scale;
isSteady = relativeMeanDeviation <= tolerance;
if isSteady
    amplitudeDeg = rad2deg(meanPeak);
end
end

function rec = record_template()
rec = struct( ...
    'parameter',"", ...
    'parameter_value',NaN, ...
    'omega',NaN, ...
    'omega_ratio_local',NaN, ...
    'A',NaN, ...
    'L',NaN, ...
    'theta_seed_deg',NaN, ...
    'h',NaN, ...
    'max_amplitude_deg',NaN, ...
    'envelope_growth_rate_per_s',NaN, ...
    'time_to_45_deg_s',NaN, ...
    'completed_cycles',0, ...
    'steady',false, ...
    'terminal_steady',false, ...
    'terminal_steady_amplitude_deg',NaN, ...
    'terminal_steady_relative_mean_deviation',NaN, ...
    'terminal_cycle_peaks_deg',{{nan(1,10)}}, ...
    'large_amplitude',false, ...
    'out_of_scope',false, ...
    'completed_requested_interval',false, ...
    'incomplete',false, ...
    'failed',false, ...
    'status',"not_run");
end

function validate_scan_axis(values,name,minimum)
if ~isnumeric(values) || ~isvector(values) || isempty(values) || ...
        any(~isfinite(values))
    error('lato:scan:InvalidAxis', ...
        '%s must be a nonempty finite numeric vector.',name);
end
values = values(:);
if any(values < minimum) || any(diff(values) <= 0) || ...
        numel(unique(values)) ~= numel(values)
    error('lato:scan:InvalidAxis', ...
        '%s must be physical, strictly increasing and unique.',name);
end
end
