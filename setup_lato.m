function root = setup_lato()
%SETUP_LATO Add the user-facing folders of the Lato toolkit to the path.
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'configs'));
addpath(fullfile(root, 'gui'));
addpath(fullfile(root, 'scripts'));
fprintf('Lato simulation toolkit ready: %s\n', root);
end
