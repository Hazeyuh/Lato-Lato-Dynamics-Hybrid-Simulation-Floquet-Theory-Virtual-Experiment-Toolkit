function opts = options(opts, defaults)
%OPTIONS Fill missing fields of a user-editable scan options struct.

if nargin < 1 || isempty(opts), opts = struct(); end
if ~isstruct(opts) || ~isscalar(opts)
    error('lato:scan:InvalidOptions','Options must be a scalar struct.');
end
names = fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k})=defaults.(names{k});
    end
end
end
