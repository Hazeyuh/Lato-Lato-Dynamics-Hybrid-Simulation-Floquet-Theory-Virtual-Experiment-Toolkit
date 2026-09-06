function d = evaluate(function_handle, t, cfg)
%EVALUATE Call and validate a Lato pivot-drive function.

assert(isa(function_handle, 'function_handle'), 'lato:InvalidDrive', ...
    'Drive must be a function handle.');
d = function_handle(t, cfg);
assert(isstruct(d), 'lato:InvalidDriveOutput', ...
    'Drive function must return a struct.');
names = {'y', 'v', 'a', 'phase', 'omega'};
for k = 1:numel(names)
    name = names{k};
    assert(isfield(d, name), 'lato:InvalidDriveOutput', ...
        'Drive output is missing field "%s".', name);
    if isscalar(d.(name)) && ~isscalar(t)
        d.(name) = d.(name) + zeros(size(t));
    end
    assert(isequal(size(d.(name)), size(t)) && ...
        (all(isfinite(d.(name)), 'all') || ...
        any(strcmp(name, {'phase', 'omega'}))), ...
        'lato:InvalidDriveOutput', ...
        'Drive field "%s" must match t and be finite.', name);
end
if ~isfield(d, 'A')
    if isfield(cfg, 'amplitude')
        d.A = cfg.amplitude + zeros(size(t));
    elseif isfield(cfg, 'A')
        d.A = cfg.A + zeros(size(t));
    else
        d.A = NaN(size(t));
    end
end

end
