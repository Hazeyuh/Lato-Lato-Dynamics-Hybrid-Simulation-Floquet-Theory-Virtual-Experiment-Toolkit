function theory = averaged_startup_boundary(p,xValues,hValues)
%AVERAGED_STARTUP_BOUNDARY Weak-modulation principal-onset approximation.
%   The damping parameter p.beta is the measured amplitude decay rate in
%   theta''+2*beta*theta'+...=0.

arguments
    p (1,1) struct
    xValues (1,:) double {mustBeFinite,mustBePositive} = linspace(1.4,2.6,121)
    hValues (1,:) double {mustBeFinite,mustBeNonnegative} = double.empty(1,0)
end

s = lato.theory.startup_parameters(p);
xValues = double(xValues(:).');
b = s.beta_s_inv/s.omega0_rad_s;

% With h=A*Omega^2/g, first-order averaging gives
% h_boundary = 2*x*sqrt((x-2)^2+4*b^2).
theory.x = xValues;
theory.h_boundary = 2*xValues.*sqrt((xValues-2).^2+4*b^2);
theory.omega0_rad_s = s.omega0_rad_s;
theory.beta_s_inv = s.beta_s_inv;
theory.beta_normalized = b;
theory.h_tip = 8*b;
theory.x_tip = 2;
theory.assumptions = ["small angle","taut string", ...
    "linear damping","weak modulation","near principal resonance"];

if isempty(hValues)
    theory.h = double.empty(1,0);
    theory.mu_normalized = double.empty(0,numel(xValues));
    theory.mu_s_inv = double.empty(0,numel(xValues));
    return
end

hValues = double(hValues(:).');
[X,H] = meshgrid(xValues,hValues);
coupling = H./(4*X);
detuning = (X-2)/2;
muNormalized = -b+real(sqrt(complex(coupling.^2-detuning.^2,0)));

theory.h = hValues;
theory.mu_normalized = muNormalized;
theory.mu_s_inv = muNormalized*s.omega0_rad_s;
end
