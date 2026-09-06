function app = lato_gui(varargin)
%LATO_GUI Virtual Lato experiment using the shared ODE45 hybrid core.
%   lato_gui opens a MATLAB-first interface for fixed-frequency,
%   continuous-frequency and stepped-frequency runs. The interface follows
%   the predecessor new_logic.m workflow, but contains no independent
%   dynamics, PID controller or artificial kick.
%
%   APP = lato_gui('Visible','off') is useful for tests. APP.readParameters
%   converts the current controls to the canonical parameter struct and
%   APP.simulate evaluates it with lato.simulate_hybrid without animation.
%
%   The solver completes before playback starts. Pause and Stop therefore
%   control playback only. The A/Omega +/- controls edit the next run.

ip = inputParser;
ip.FunctionName = 'lato_gui';
addParameter(ip, 'Visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(ip, 'ParameterFile', 'paper_parameters', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(ip, varargin{:});

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'configs'));

base = load_parameter_file(string(ip.Results.ParameterFile));
base = lato.validate_parameters(base);

fig = figure('Name', 'Lato virtual experiment platform', ...
    'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure', ...
    'Color', [0.96 0.96 0.96], 'Position', [80 80 1180 720], ...
    'Visible', char(ip.Results.Visible), 'CloseRequestFcn', @close_app);

panel = uipanel('Parent', fig, 'Title', 'Run configuration', ...
    'Units', 'pixels', 'Position', [15 15 290 690], ...
    'FontWeight', 'bold');
animation_axes = axes('Parent', fig, 'Units', 'pixels', ...
    'Position', [335 315 400 365]);
trace_axes = axes('Parent', fig, 'Units', 'pixels', ...
    'Position', [780 365 370 315]);
drive_axes = axes('Parent', fig, 'Units', 'pixels', ...
    'Position', [780 60 370 220]);

controls = struct();
y = 625;
controls.parameter_file = labelled_edit(panel, 'Parameter function', ...
    char(ip.Results.ParameterFile), y); y = y - 45;

uicontrol(panel, 'Style', 'text', 'String', 'Drive mode', ...
    'HorizontalAlignment', 'left', 'Position', [15 y+22 120 18]);
controls.mode = uicontrol(panel, 'Style', 'popupmenu', ...
    'String', {'Fixed frequency','Continuous sweep','Step sweep'}, ...
    'Value', 1, 'Position', [15 y 250 24], 'Callback', @mode_changed);
y = y - 48;

controls.amplitude = labelled_edit(panel, 'Amplitude A (m)', ...
    num2str(base.drive.amplitude, '%.6g'), y); y = y - 45;
controls.omega_start = labelled_edit(panel, 'Start Omega (rad/s)', ...
    num2str(base.drive.omega, '%.6g'), y); y = y - 45;
controls.omega_end = labelled_edit(panel, 'End Omega (rad/s)', ...
    num2str(base.drive.omega, '%.6g'), y); y = y - 45;
controls.step_duration = labelled_edit(panel, 'Step interval (s)', ...
    '2', y); y = y - 45;
controls.duration = labelled_edit(panel, 'Simulation duration (s)', ...
    num2str(base.solver.t_final, '%.6g'), y); y = y - 45;
controls.theta1 = labelled_edit(panel, 'Initial theta_1 (deg)', ...
    num2str(rad2deg(base.initial_state.theta(1)), '%.6g'), y); y = y - 45;
controls.theta2 = labelled_edit(panel, 'Initial theta_2 (deg)', ...
    num2str(rad2deg(base.initial_state.theta(2)), '%.6g'), y); y = y - 45;

uicontrol(panel, 'Style', 'text', 'String', 'Manual changes (next run)', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold', ...
    'Position', [15 y+8 250 20]);
y = y - 27;
controls.omega_step = labelled_edit(panel, 'Omega increment', '0.1', y);
y = y - 45;
controls.amplitude_step = labelled_edit(panel, 'A increment', '0.002', y);
y = y - 42;

controls.omega_minus = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'Omega -', 'Position', [15 y 57 28], ...
    'Callback', @(~,~) adjust_edit(controls.omega_start, ...
    -read_scalar(controls.omega_step, 'Omega increment')));
controls.omega_plus = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'Omega +', 'Position', [78 y 57 28], ...
    'Callback', @(~,~) adjust_edit(controls.omega_start, ...
    read_scalar(controls.omega_step, 'Omega increment')));
controls.amplitude_minus = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'A -', 'Position', [145 y 55 28], ...
    'Callback', @(~,~) adjust_edit(controls.amplitude, ...
    -read_scalar(controls.amplitude_step, 'A increment')));
controls.amplitude_plus = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'A +', 'Position', [206 y 55 28], ...
    'Callback', @(~,~) adjust_edit(controls.amplitude, ...
    read_scalar(controls.amplitude_step, 'A increment')));

controls.run = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Run / rerun', ...
    'FontSize', 11, 'FontWeight', 'bold', 'Position', [350 245 115 36], ...
    'Callback', @run_callback);
controls.pause = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Pause', ...
    'FontSize', 10, 'Enable', 'off', 'Position', [480 245 90 36], ...
    'Callback', @pause_callback);
controls.stop = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop playback', ...
    'FontSize', 10, 'Enable', 'off', 'Position', [585 245 125 36], ...
    'Callback', @stop_callback);
controls.status = uicontrol(fig, 'Style', 'text', ...
    'String', ['Ready. Simulation uses lato.simulate_hybrid (ODE45). ' ...
    'Pause/stop act on playback.'], 'HorizontalAlignment', 'left', ...
    'BackgroundColor', get(fig, 'Color'), 'Position', [335 185 405 50]);
controls.note = uicontrol(fig, 'Style', 'text', ...
    'String', ['Angles above 45 deg are labelled large-amplitude. ' ...
    'The shared core stops interpretation at 90 deg.'], ...
    'HorizontalAlignment', 'left', 'ForegroundColor', [0.35 0.20 0.05], ...
    'BackgroundColor', get(fig, 'Color'), 'Position', [335 125 405 45]);

state = struct('paused', false, 'stopped', false, 'animating', false);
setappdata(fig, 'lato_gui_state', state);
setappdata(fig, 'lato_gui_result', []);

title(animation_axes, 'Computed motion playback');
xlabel(animation_axes, 'x (m)'); ylabel(animation_axes, 'y (m)');
axis(animation_axes, 'equal'); grid(animation_axes, 'on');
title(trace_axes, 'Angular response');
xlabel(trace_axes, 'Time (s)'); ylabel(trace_axes, 'Angle (deg)');
grid(trace_axes, 'on');
title(drive_axes, 'Drive history');
xlabel(drive_axes, 'Time (s)'); ylabel(drive_axes, 'Omega (rad/s)');
grid(drive_axes, 'on');
mode_changed();

app.figure = fig;
app.controls = controls;
app.axes = struct('animation', animation_axes, 'trace', trace_axes, ...
    'drive', drive_axes);
app.readParameters = @read_parameters;
app.simulate = @simulate_current;
app.getResult = @() getappdata(fig, 'lato_gui_result');

    function h = labelled_edit(parent, label, value, y_position)
        uicontrol(parent, 'Style', 'text', 'String', label, ...
            'HorizontalAlignment', 'left', 'Position', [15 y_position+22 180 18]);
        h = uicontrol(parent, 'Style', 'edit', 'String', value, ...
            'BackgroundColor', 'white', 'Position', [15 y_position 250 24]);
    end

    function mode_changed(~,~)
        mode = selected_mode();
        variable_mode = ~strcmp(mode, 'fixed');
        set(controls.omega_end, 'Enable', on_off(variable_mode));
        set(controls.step_duration, 'Enable', on_off(strcmp(mode, 'step')));
    end

    function mode = selected_mode()
        value = get(controls.mode, 'Value');
        labels = {'fixed','continuous','step'};
        mode = labels{value};
    end

    function p = read_parameters()
        parameter_name = string(strtrim(get(controls.parameter_file, 'String')));
        p = load_parameter_file(parameter_name);

        amplitude = read_scalar(controls.amplitude, 'Amplitude A');
        omega_start = read_scalar(controls.omega_start, 'Start Omega');
        omega_end = read_scalar(controls.omega_end, 'End Omega');
        duration = read_scalar(controls.duration, 'Simulation duration');
        theta = deg2rad([read_scalar(controls.theta1, 'Initial theta_1'), ...
            read_scalar(controls.theta2, 'Initial theta_2')]);

        assert(amplitude >= 0, 'lato:gui:InvalidInput', ...
            'Amplitude A must be non-negative.');
        assert(omega_start > 0 && omega_end > 0, 'lato:gui:InvalidInput', ...
            'Drive frequencies must be positive.');
        assert(duration > 0, 'lato:gui:InvalidInput', ...
            'Simulation duration must be positive.');
        assert(all(abs(theta) < deg2rad(90)), 'lato:gui:AngleLimit', ...
            'Initial angles must lie strictly between -90 and 90 degrees.');

        t0 = p.solver.t_start;
        p.solver.t_final = duration;
        p.solver.t_end = t0 + duration;
        p.initial_state.theta = theta;
        p.initial_state.omega = [0 0];
        p.initial_state.slack = [false false];
        p.initial_state.x = [NaN NaN];
        p.initial_state.y = [NaN NaN];
        p.initial_state.vx = [NaN NaN];
        p.initial_state.vy = [NaN NaN];

        p.drive.amplitude = amplitude;
        p.drive.A = amplitude;
        p.drive.omega = omega_start;
        p.drive.phase0 = 0;
        p.drive.time_origin = t0;
        p.drive.breakpoints = [];

        switch selected_mode()
            case 'fixed'
                p.drive.function = @lato.drive.harmonic;
            case 'continuous'
                p.drive.function = @lato.drive.linear_chirp;
                p.drive.omega_start = omega_start;
                p.drive.sweep_rate = (omega_end - omega_start) / duration;
                p.drive.omega_min = min(omega_start, omega_end);
                p.drive.omega_max = max(omega_start, omega_end);
            case 'step'
                interval = read_scalar(controls.step_duration, 'Step interval');
                assert(interval > 0, 'lato:gui:InvalidInput', ...
                    'Step interval must be positive.');
                times = t0:interval:(t0 + duration - 10*eps(t0 + duration));
                if isempty(times), times = t0; end
                p.drive.function = @lato.drive.step_schedule;
                p.drive.times = times;
                p.drive.omega_values = linspace(omega_start, omega_end, numel(times));
                p.drive.amplitude_values = amplitude;
                p.drive.breakpoints = times(2:end);
        end
        p = lato.validate_parameters(p);
    end

    function [result, p] = simulate_current()
        assert(isgraphics(fig), 'lato:gui:Closed', 'The GUI has been closed.');
        p = read_parameters();
        set(controls.run, 'Enable', 'off');
        set(controls.status, 'String', 'Computing with ODE45; playback follows ...');
        drawnow;
        restore = onCleanup(@() restore_run_button());
        result = lato.simulate_hybrid(p);
        setappdata(fig, 'lato_gui_result', result);
        plot_result_summary(result);
        set(controls.status, 'String', result_status(result));
        clear restore
    end

    function run_callback(~,~)
        try
            result = simulate_current();
            play_result(result);
        catch exception
            if isgraphics(fig)
                set(controls.status, 'String', ['Error: ' exception.message]);
                errordlg(exception.message, 'Lato simulation error', 'modal');
            end
        end
    end

    function plot_result_summary(result)
        t = result.time(:);
        theta_deg = rad2deg(result.states.theta);
        cla(trace_axes);
        plot(trace_axes, t, theta_deg(:,1), 'r-', 'DisplayName', '\theta_1');
        hold(trace_axes, 'on');
        plot(trace_axes, t, theta_deg(:,2), 'b-', 'DisplayName', '\theta_2');
        yline(trace_axes, 45, ':', '45 deg', 'Color', [0.65 0.45 0.05], ...
            'HandleVisibility', 'off');
        yline(trace_axes, -45, ':', 'Color', [0.65 0.45 0.05], ...
            'HandleVisibility', 'off');
        hold(trace_axes, 'off');
        grid(trace_axes, 'on'); legend(trace_axes, 'Location', 'best');
        title(trace_axes, 'Angular response');
        xlabel(trace_axes, 'Time (s)'); ylabel(trace_axes, 'Angle (deg)');

        cla(drive_axes);
        yyaxis(drive_axes, 'left');
        omegaColor = [0.10 0.40 0.75];
        amplitudeColor = [0.85 0.33 0.10];
        plot(drive_axes, t, result.drive.omega(:), '-', ...
            'Color', omegaColor, 'LineWidth', 1.1);
        drive_axes.YAxis(1).Color = omegaColor;
        ylabel(drive_axes, 'Omega (rad/s)');
        yyaxis(drive_axes, 'right');
        plot(drive_axes, t, result.drive.A(:) * 1000, '-', ...
            'Color', amplitudeColor, 'LineWidth', 1.1);
        drive_axes.YAxis(2).Color = amplitudeColor;
        ylabel(drive_axes, 'A (mm)'); xlabel(drive_axes, 'Time (s)');
        title(drive_axes, 'Drive history'); grid(drive_axes, 'on');
    end

    function play_result(result)
        state = getappdata(fig, 'lato_gui_state');
        state.paused = false; state.stopped = false; state.animating = true;
        setappdata(fig, 'lato_gui_state', state);
        set(controls.pause, 'Enable', 'on', 'String', 'Pause');
        set(controls.stop, 'Enable', 'on');

        t = result.time(:);
        x = result.states.x;
        y_ball = result.states.y;
        y_pivot = result.drive.y(:);
        n = numel(t);
        stride = max(1, ceil(n / 900));
        indices = unique([1:stride:n, n]);
        radius = result.parameters.radius;
        extent = max(result.parameters.ell + abs(result.parameters.drive.amplitude), ...
            max(abs(x(:))) + 2*radius);

        cla(animation_axes);
        hold(animation_axes, 'on'); grid(animation_axes, 'on');
        axis(animation_axes, 'equal');
        xlim(animation_axes, 1.15 * [-extent extent]);
        ymin = min([y_ball(:); y_pivot(:)]) - 3*radius;
        ymax = max([y_ball(:); y_pivot(:)]) + 3*radius;
        ylim(animation_axes, [ymin ymax]);
        pivot_marker = plot(animation_axes, 0, y_pivot(1), 'ko', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 5);
        strings = plot(animation_axes, [0 x(1,1); 0 x(1,2)], ...
            [y_pivot(1) y_ball(1,1); y_pivot(1) y_ball(1,2)], 'k-');
        ball1 = rectangle(animation_axes, 'Position', ...
            [x(1,1)-radius, y_ball(1,1)-radius, 2*radius, 2*radius], ...
            'Curvature', [1 1], 'FaceColor', [0.85 0.18 0.14]);
        ball2 = rectangle(animation_axes, 'Position', ...
            [x(1,2)-radius, y_ball(1,2)-radius, 2*radius, 2*radius], ...
            'Curvature', [1 1], 'FaceColor', [0.12 0.35 0.85]);
        hold(animation_axes, 'off');

        for k = indices
            if ~isgraphics(fig), return; end
            state = getappdata(fig, 'lato_gui_state');
            if state.stopped, break; end
            while state.paused && ~state.stopped && isgraphics(fig)
                drawnow;
                pause(0.03);
                state = getappdata(fig, 'lato_gui_state');
            end
            if ~isgraphics(fig) || state.stopped, break; end
            set(pivot_marker, 'YData', y_pivot(k));
            set(strings(1), 'XData', [0 x(k,1)], ...
                'YData', [y_pivot(k) y_ball(k,1)]);
            set(strings(2), 'XData', [0 x(k,2)], ...
                'YData', [y_pivot(k) y_ball(k,2)]);
            set(ball1, 'Position', [x(k,1)-radius, y_ball(k,1)-radius, ...
                2*radius, 2*radius]);
            set(ball2, 'Position', [x(k,2)-radius, y_ball(k,2)-radius, ...
                2*radius, 2*radius]);
            title(animation_axes, sprintf('Computed playback, t = %.2f s', t(k)));
            drawnow limitrate;
            pause(0.004);
        end

        if isgraphics(fig)
            state = getappdata(fig, 'lato_gui_state');
            state.animating = false; state.paused = false;
            setappdata(fig, 'lato_gui_state', state);
            set(controls.pause, 'Enable', 'off', 'String', 'Pause');
            set(controls.stop, 'Enable', 'off');
            if state.stopped
                set(controls.status, 'String', ...
                    ['Playback stopped. Computed result is retained. ' ...
                    result_status(result)]);
            else
                set(controls.status, 'String', result_status(result));
            end
        end
    end

    function pause_callback(~,~)
        state = getappdata(fig, 'lato_gui_state');
        if ~state.animating, return; end
        state.paused = ~state.paused;
        setappdata(fig, 'lato_gui_state', state);
        if state.paused
            set(controls.pause, 'String', 'Resume');
        else
            set(controls.pause, 'String', 'Pause');
        end
    end

    function stop_callback(~,~)
        state = getappdata(fig, 'lato_gui_state');
        state.stopped = true; state.paused = false;
        setappdata(fig, 'lato_gui_state', state);
    end

    function close_app(~,~)
        if ~isgraphics(fig), return; end
        state = getappdata(fig, 'lato_gui_state');
        state.stopped = true; state.paused = false;
        setappdata(fig, 'lato_gui_state', state);
        delete(fig);
    end

    function restore_run_button()
        if isgraphics(fig), set(controls.run, 'Enable', 'on'); end
    end
end

function p = load_parameter_file(name)
name = erase(strtrim(string(name)), ".m");
assert(strlength(name) > 0 && isvarname(char(name)), ...
    'lato:gui:InvalidParameterFile', ...
    'Parameter function must be a MATLAB function name, for example paper_parameters.');
f = str2func(char(name));
assert(~isempty(which(char(name))), 'lato:gui:MissingParameterFile', ...
    'Parameter function "%s" was not found on the MATLAB path.', name);
p = f();
assert(isstruct(p) && isscalar(p), 'lato:gui:InvalidParameterFile', ...
    'Parameter function "%s" must return one scalar struct.', name);
end

function value = read_scalar(control, label)
value = str2double(get(control, 'String'));
assert(isscalar(value) && isfinite(value), 'lato:gui:InvalidInput', ...
    '%s must be a finite scalar.', label);
end

function adjust_edit(control, increment)
value = read_scalar(control, 'Current value') + increment;
set(control, 'String', num2str(value, '%.8g'));
end

function value = on_off(condition)
if condition, value = 'on'; else, value = 'off'; end
end

function message = result_status(result)
theta_max = max(abs(rad2deg(result.states.theta)), [], 'all');
if theta_max > 45
    amplitude_label = 'large-amplitude (>45 deg)';
else
    amplitude_label = 'small/moderate-amplitude (<=45 deg)';
end
message = sprintf('Computed to %.3f s; %s; termination: %s; events: %d.', ...
    result.time(end), amplitude_label, char(result.termination.reason), ...
    result.events.count);
end
