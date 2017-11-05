% Set number of subjects and result matrix
GROUP_SIZE = 12;
results = zeros(6,4,GROUP_SIZE);

% Get parameters for automaticityModel
loaded = load('fmri/particleswarm_target_17_10_31.mat');
arg_vector = loaded.x;
CONFIG = 'FMRI';

% Create parallel pool
gcp;
parfor i=1:GROUP_SIZE
    arg_struct = argvectortostruct(arg_vector, CONFIG);
    [r_x_vals, r_y_vals] = createFMRIInput(11520);
    visualinput = [r_x_vals, r_y_vals];
    optional_parms = struct('FMRI_META_GROUP_RUN', 1, ...
                            'VIS_INPUT_FROM_PARM', 1, ...
                            'visualinput', visualinput);
    [~,results(:,:,i)] = automaticityModelFast_mex(arg_struct, optional_parms);
end
% Delete parallel pool
delete(gcp('nocreate'));

% Get correlation information from results
permuted_results = permute(results(:,:,:), [3 2 1]);
corr_mat = zeros(5,4);
for i=1:5
    corr_mat(i,:) = corrneuron(permuted_results(:,:,i), permuted_results(:,:,6));
end