function make_F1_mechanism_timeline()
%MAKE_F1_MECHANISM_TIMELINE Parameter modulation and one 20 s hybrid run.
parameter_file = @paper_parameters; % Change this line to select another file.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'configs'));
p = parameter_file();
p.solver.t_end = 20;
result = lato.simulate_hybrid(p);
run_id = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
results_dir = fullfile(root, 'results');
figures_dir = fullfile(root, 'figures');
if ~isfolder(results_dir), mkdir(results_dir); end
if ~isfolder(figures_dir), mkdir(figures_dir); end
result_file = fullfile(results_dir, "F1_single_run_" + run_id + ".mat");
save(result_file, 'result', '-v7.3');
output_base = fullfile(figures_dir, "F1_mechanism_timeline");
lato.plot.mechanism_timeline(result, output_base);
fprintf('Saved %s and matching figures.\n', result_file);

end
