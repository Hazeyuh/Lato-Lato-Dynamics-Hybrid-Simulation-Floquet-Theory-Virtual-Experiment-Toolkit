function T = summary_table(items)
%SUMMARY_TABLE Convert a struct array returned by scan functions to table.

if isempty(items)
    T = table();
    return
end

fields = {'omega','A','L','h','steady','steady_amplitude_deg', ...
    'max_amplitude_deg','amplitude_class','out_of_scope', ...
    'large_amplitude','collision_rate_hz','slack_fraction', ...
    'midline_offset_deg','mean_energy_J','status'};
T = struct2table(items);
present = fields(ismember(fields, T.Properties.VariableNames));
T = T(:, present);
end
