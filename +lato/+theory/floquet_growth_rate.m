function result = floquet_growth_rate(varargin)
%FLOQUET_GROWTH_RATE Floquet growth of the taut non-impact reference model.
%   RESULT = lato.theory.floquet_growth_rate(OMEGA,A,L,BETA) integrates
%
%     theta'' + 2*beta*theta' +
%       [g/L - A*Omega^2/L*cos(Omega*t)]*theta = 0
%
%   over one drive period. A is the one-sided pivot displacement amplitude,
%   L is the pivot-to-ball-center length, and beta is the measured amplitude
%   decay rate.
%
%   A legacy call RESULT = floquet_growth_rate(P,OMEGA,H,...) is accepted.
%   In that form H=A*Omega^2/g and P.beta is interpreted as the decay rate.

[omegaRadS,amplitudeM,lengthM,betaSInv,options,legacy] = ...
    parse_inputs(varargin{:});

g = options.Gravity;
periodS = 2*pi/omegaRadS;
z0 = reshape(eye(2),4,1);
odeOptions = odeset( ...
    'RelTol',options.RelTol, ...
    'AbsTol',options.AbsTol, ...
    'MaxStep',periodS*options.MaxStepFraction);
rhs = @(t,z) fundamental_rhs(t,z,omegaRadS,amplitudeM,lengthM, ...
    betaSInv,g,options.PhaseRad);

success = true;
message = "ok";
try
    [tSolution,zSolution] = ode45(rhs,[0,periodS],z0,odeOptions);
    monodromy = reshape(zSolution(end,:),2,2);
    multipliers = eig(monodromy);
    spectralRadius = max(abs(multipliers));
    muMax = log(spectralRadius)/periodS;
catch exception
    success = false;
    message = string(exception.message);
    tSolution = zeros(0,1);
    monodromy = nan(2);
    multipliers = nan(2,1);
    spectralRadius = NaN;
    muMax = NaN;
end

omega0 = sqrt(g/lengthM);
muNormalized = muMax/omega0;
detM = det(monodromy);
expectedDetM = exp(-2*betaSInv*periodS);
detRelativeError = abs(detM-expectedDetM)/max(abs(expectedDetM),eps);

tol = options.StabilityTolerance;
if ~success || ~isfinite(muMax)
    stabilityCode = int8(0);
    stabilityLabel = "solver_failure";
elseif muMax > tol
    stabilityCode = int8(1);
    stabilityLabel = "unstable";
elseif muMax < -tol
    stabilityCode = int8(-1);
    stabilityLabel = "stable";
else
    stabilityCode = int8(0);
    stabilityLabel = "neutral/boundary";
end

signature = sprintf( ...
    'floquet-v2|Omega=%.12g|A=%.12g|L=%.12g|beta=%.12g|rtol=%.3g|atol=%.3g|maxstep=%.6g', ...
    omegaRadS,amplitudeM,lengthM,betaSInv,options.RelTol, ...
    options.AbsTol,options.MaxStepFraction);

result.mu_max = muMax;
result.multipliers = multipliers;
result.detM = detM;
result.period = periodS;
result.success = success;
result.signature = string(signature);
result.message = message;
result.monodromy = monodromy;
result.expected_detM = expectedDetM;
result.detM_relative_error = detRelativeError;
result.Omega_rad_s = omegaRadS;
result.A_m = amplitudeM;
result.L_m = lengthM;
result.beta_s_inv = betaSInv;
result.omega0_rad_s = omega0;
result.Omega_over_omega0 = omegaRadS/omega0;
result.A_over_L = amplitudeM/lengthM;
result.h = amplitudeM*omegaRadS^2/g;
result.stability_code = stabilityCode;
result.stability_label = stabilityLabel;
result.solver = struct( ...
    'name',"ode45", ...
    'RelTol',options.RelTol, ...
    'AbsTol',options.AbsTol, ...
    'MaxStepFraction',options.MaxStepFraction, ...
    'accepted_time_points',numel(tSolution));

% Compatibility aliases used by the historical onset scripts.
result.period_s = periodS;
result.mu_max_s_inv = muMax;
result.mu_max_normalized = muNormalized;
result.spectral_radius = spectralRadius;
result.expected_determinant = expectedDetM;
result.determinant_relative_error = detRelativeError;
result.omega_rad_s = omegaRadS;
result.x_omega_over_omega0 = omegaRadS/omega0;
result.h_a0_over_g = result.h;
result.a0_m_s2 = amplitudeM*omegaRadS^2;
result.equivalent_amplitude_m = amplitudeM;
result.legacy_input = legacy;

end

function dz = fundamental_rhs(t,z,omegaRadS,amplitudeM,lengthM,betaSInv,g,phaseRad)
phi = reshape(z,2,2);
coefficient = g/lengthM - ...
    amplitudeM*omegaRadS^2/lengthM*cos(omegaRadS*t+phaseRad);
B = [0,1;-coefficient,-2*betaSInv];
dz = reshape(B*phi,4,1);
end

function [omegaRadS,amplitudeM,lengthM,betaSInv,options,legacy] = ...
        parse_inputs(varargin)
if nargin < 3
    error('lato:theory:InvalidFloquetInput', ...
        'Expected (Omega,A,L,beta) or the legacy (p,Omega,h) call.');
end

options = default_options();
legacy = isstruct(varargin{1});
if legacy
    p = varargin{1};
    omegaRadS = double(varargin{2});
    h = double(varargin{3});
    s = lato.theory.startup_parameters(p);
    amplitudeM = h*s.g_m_s2/omegaRadS^2;
    lengthM = s.ell_m;
    betaSInv = s.beta_s_inv;
    optionArgs = varargin(4:end);
    options.Gravity = s.g_m_s2;
else
    if nargin < 4
        error('lato:theory:InvalidFloquetInput', ...
            'The physical interface requires Omega, A, L, and beta.');
    end
    omegaRadS = double(varargin{1});
    amplitudeM = double(varargin{2});
    lengthM = double(varargin{3});
    betaSInv = double(varargin{4});
    optionArgs = varargin(5:end);
end

validateattributes(omegaRadS,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(amplitudeM,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(lengthM,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(betaSInv,{'numeric'},{'scalar','real','finite','nonnegative'});

if numel(optionArgs) == 1 && isstruct(optionArgs{1})
    supplied = optionArgs{1};
    names = fieldnames(supplied);
    for k = 1:numel(names)
        options = assign_option(options,names{k},supplied.(names{k}));
    end
elseif mod(numel(optionArgs),2) == 0
    for k = 1:2:numel(optionArgs)
        options = assign_option(options,string(optionArgs{k}),optionArgs{k+1});
    end
else
    error('lato:theory:InvalidFloquetOptions', ...
        'Options must be a structure or name-value pairs.');
end
end

function options = default_options()
options.RelTol = 1e-8;
options.AbsTol = 1e-10;
options.MaxStepFraction = 1/80;
options.PhaseRad = 0;
options.StabilityTolerance = 1e-10;
options.Gravity = 9.81;
end

function options = assign_option(options,name,value)
name = lower(strrep(char(name),'_',''));
switch name
    case 'reltol'
        options.RelTol = double(value);
    case 'abstol'
        options.AbsTol = double(value);
    case 'maxstepfraction'
        options.MaxStepFraction = double(value);
    case 'phaserad'
        options.PhaseRad = double(value);
    case {'stabilitytolerance','stabilitytolerancenormalized'}
        options.StabilityTolerance = double(value);
    case {'gravity','g'}
        options.Gravity = double(value);
    otherwise
        error('lato:theory:UnknownFloquetOption', ...
            'Unknown Floquet option: %s',name);
end
end
