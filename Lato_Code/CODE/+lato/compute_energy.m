function energy = compute_energy(result, p)
%COMPUTE_ENERGY Mechanical-energy diagnostics in the laboratory frame.
%   Potential energy has an arbitrary zero at y=0. The relative-pivot form
%   is also returned because it equals m*g*ell*(1-cos(theta)) for a taut ball.

arguments
    result (1,1) struct
    p (1,1) struct
end

s = result.states;
speed_squared = s.vx.^2 + s.vy.^2;
energy.kinetic = 0.5*p.mass*speed_squared;
energy.gravitational = p.mass*p.g*s.y;
energy.mechanical = energy.kinetic + energy.gravitational;
energy.total_mechanical = sum(energy.mechanical, 2);

pivot_y = result.drive.y(:);
energy.potential_relative_to_pivot = p.mass*p.g .* ...
    (s.y - pivot_y + p.ell);
energy.relative_pendulum = energy.kinetic + ...
    energy.potential_relative_to_pivot;
energy.total_relative_pendulum = sum(energy.relative_pendulum, 2);

% These powers are diagnostics, not a complete collision/work balance.
energy.linear_damping_power = zeros(size(s.theta));
taut = ~s.slack;
energy.linear_damping_power(taut) = ...
    -2*p.mass*p.ell^2*p.beta .* s.omega(taut).^2;
energy.total_linear_damping_power = sum(energy.linear_damping_power, 2);

energy.reference_note = [ ...
    "Laboratory mechanical energy uses U=m*g*y; its zero is arbitrary. ", ...
    "Instantaneous collision losses and pivot constraint work are not ", ...
    "included in the damping-power diagnostic."];

end
