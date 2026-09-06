function cycle = cycle_trace(result)
%CYCLE_TRACE Extract per-drive-cycle process quantities from a run.

t = result.time(:);
theta = get_theta(result);
[phase, omega, A] = get_drive_history(result, t);

% Prefer the core's completed-cycle bookkeeping when available.
if isfield(result,'cycle') && isfield(result.cycle,'index') && ...
        ~isempty(result.cycle.index)
    rc=result.cycle;
    amplitudeDeg=rad2deg(rc.peak_half_opening(:));
    n=numel(amplitudeDeg);
    maxAngleDeg=amplitudeDeg;
    cycle=table(rc.index(:),rc.start_time(:),rc.end_time(:), ...
        rc.mean_omega(:),rc.mean_A(:),amplitudeDeg,maxAngleDeg, ...
        'VariableNames',{'cycle_index','start_time','end_time','omega_mean', ...
        'A_mean','amplitude_deg','max_angle_deg'});
    cycle.amplitude_class=lato.scan.classify_amplitude(cycle.amplitude_deg);
    cycle.out_of_scope=cycle.amplitude_deg>90;
    cycle.large_amplitude=cycle.amplitude_deg>45 & cycle.amplitude_deg<=90;
    return
end

index = floor((phase-phase(1))/(2*pi));
ids = unique(index, 'stable');
ids = ids(ids >= 0);
n = numel(ids);

startTime = nan(n,1); endTime = nan(n,1); meanOmega = nan(n,1);
meanA = nan(n,1); amplitudeDeg = nan(n,1); maxAngleDeg = nan(n,1);
for k = 1:n
    mask = index == ids(k);
    tk = t(mask);
    if numel(tk) < 2, continue; end
    startTime(k) = tk(1); endTime(k) = tk(end);
    meanOmega(k) = trapz(tk, omega(mask)) / (tk(end)-tk(1));
    meanA(k) = trapz(tk, A(mask)) / (tk(end)-tk(1));
    th = theta(mask,:);
    if size(th,2)>=2
        halfOpening=0.5*abs(th(:,1)-th(:,2));
        amplitudeDeg(k)=max(rad2deg(halfOpening),[],'omitnan');
    else
        amplitudeDeg(k)=max(abs(rad2deg(th)),[],'all','omitnan');
    end
    maxAngleDeg(k) = max(abs(rad2deg(th)),[],'all','omitnan');
end

complete = isfinite(startTime) & isfinite(endTime);
cycle = table(ids(complete),startTime(complete),endTime(complete), ...
    meanOmega(complete),meanA(complete),amplitudeDeg(complete), ...
    maxAngleDeg(complete), ...
    'VariableNames',{'cycle_index','start_time','end_time','omega_mean', ...
    'A_mean','amplitude_deg','max_angle_deg'});
cycle.amplitude_class = lato.scan.classify_amplitude(cycle.amplitude_deg);
cycle.out_of_scope = cycle.amplitude_deg > 90;
cycle.large_amplitude = cycle.amplitude_deg > 45 & cycle.amplitude_deg <= 90;
end

function theta = get_theta(result)
if isfield(result,'theta1') && isfield(result,'theta2')
    theta = [result.theta1(:),result.theta2(:)];
elseif isfield(result,'states') && size(result.states,2) >= 3
    theta = result.states(:,[1 3]);
elseif isfield(result,'states')
    theta = result.states(:,1:min(2,size(result.states,2)));
else
    theta = nan(numel(result.time),1);
end
end

function [phase,omega,A] = get_drive_history(result,t)
phase = []; omega = []; A = [];
if isfield(result,'drive') && isstruct(result.drive)
    if isfield(result.drive,'phase'), phase = result.drive.phase(:); end
    if isfield(result.drive,'omega'), omega = result.drive.omega(:); end
    if isfield(result.drive,'A'), A = result.drive.A(:); end
end
if isempty(omega)
    if isempty(phase)
        omega = ones(size(t));
    else
        omega = gradient(phase,t);
    end
end
if isempty(phase)
    phase = cumtrapz(t,omega);
end
if isempty(A), A = nan(size(t)); end
if isscalar(omega), omega = repmat(omega,size(t)); end
if isscalar(A), A = repmat(A,size(t)); end
end
