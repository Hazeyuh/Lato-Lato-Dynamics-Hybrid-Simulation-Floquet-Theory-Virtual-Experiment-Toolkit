function s = startup_parameters(p)
%STARTUP_PARAMETERS Resolve parameters used by the linear onset model.
%   S = lato.theory.startup_parameters(P) accepts the public parameter
%   structure and returns a small canonical structure.  Several legacy
%   field names are accepted so the onset code can also inspect older
%   parameter files without copying their dynamics.

arguments
    p (1, 1) struct
end

s.g_m_s2 = require_scalar(p, ...
    {'g_m_s2', 'g', 'gravity_m_s2', 'gravity', ...
     'physics.g_m_s2', 'physics.g', 'constants.g_m_s2'}, ...
    'gravitational acceleration');
s.ell_m = require_scalar(p, ...
    {'ell_m', 'ell', 'length_m', 'pendulum_length_m', ...
     'geometry.ell_m', 'geometry.length_m'}, ...
    'pivot-to-center length');
s.beta_s_inv = require_scalar(p, ...
    {'beta_s_inv', 'mechanical_damping_beta', 'beta', ...
     'mechanical_damping_beta_s_inv', 'damping.beta_s_inv', ...
     'dissipation.beta_s_inv'}, ...
    'linear damping coefficient beta');

s.drive_amplitude_m = optional_scalar(p, ...
    {'drive_amplitude_m', 'drive_amplitude', 'A_m', 'A', ...
     'drive.amplitude', 'drive.amplitude_m', 'drive.A_m', ...
     'excitation.amplitude_m'}, NaN);
s.drive_omega_rad_s = optional_scalar(p, ...
    {'drive_omega_rad_s', 'drive_omega', 'Omega_rad_s', 'Omega', ...
     'omega_rad_s', 'drive.omega', 'drive.omega_rad_s', ...
     'drive.angular_frequency_rad_s', 'excitation.omega_rad_s'}, NaN);

validateattributes(s.g_m_s2, {'numeric'}, {'positive', 'finite'});
validateattributes(s.ell_m, {'numeric'}, {'positive', 'finite'});
validateattributes(s.beta_s_inv, {'numeric'}, {'nonnegative', 'finite'});
if isfinite(s.drive_amplitude_m)
    validateattributes(s.drive_amplitude_m, {'numeric'}, {'nonnegative'});
end
if isfinite(s.drive_omega_rad_s)
    validateattributes(s.drive_omega_rad_s, {'numeric'}, {'positive'});
end

s.omega0_rad_s = sqrt(s.g_m_s2 / s.ell_m);
s.damping_ratio_beta_over_omega0 = s.beta_s_inv / s.omega0_rad_s;

if isfinite(s.drive_amplitude_m) && isfinite(s.drive_omega_rad_s)
    s.current_x = s.drive_omega_rad_s / s.omega0_rad_s;
    s.current_h = s.drive_amplitude_m * s.drive_omega_rad_s^2 / s.g_m_s2;
else
    s.current_x = NaN;
    s.current_h = NaN;
end

end

function value = require_scalar(p, paths, label)
value = optional_scalar(p, paths, NaN);
if ~isfinite(value)
    error('lato:theory:MissingParameter', ...
        'Could not resolve %s from the parameter structure.', label);
end
end

function value = optional_scalar(p, paths, default_value)
value = default_value;
for k = 1:numel(paths)
    [found, candidate] = value_at_path(p, paths{k});
    if found && isnumeric(candidate) && isscalar(candidate) && isreal(candidate)
        value = double(candidate);
        return;
    end
end
end

function [found, value] = value_at_path(s, path)
parts = strsplit(path, '.');
value = s;
found = true;
for k = 1:numel(parts)
    if ~isstruct(value) || ~isfield(value, parts{k})
        found = false;
        value = [];
        return;
    end
    value = value.(parts{k});
end
end
