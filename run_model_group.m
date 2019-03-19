% Add createFMRIInput function to path
addpath('visinputgen', 'functions', 'classes');

% Set number of subjects and result matrix
GROUP_SIZE = 12;
results = zeros(4,4,GROUP_SIZE);

% Get parameters for automaticityModel
loaded = load('fmri/particleswarm_target_17_12_14.mat');
arg_vector = loaded.x;
CONFIG = AutomaticityConfiguration.FMRI;

% Create parallel pool
tic;
poolobj = gcp;
parfor i=1:GROUP_SIZE
    arg_struct = argvectortostruct(arg_vector, CONFIG);
    [r_x_vals, r_y_vals] = createFMRIInput(11520);
    visualinput = [r_x_vals, r_y_vals];
    optional_parms = struct('FMRI_META_GROUP_RUN', 1, ...
                            'VIS_INPUT_FROM_PARM', 1, ...
                            'visualinput', visualinput);
    [~,results(:,:,i)] = automaticityModel_mex(arg_struct, optional_parms);
end
% Delete parallel pool
delete(gcp('nocreate'));
toc;

% Get correlation information from results
permuted_results = permute(results(:,:,:), [3 2 1]);
corr_mat = zeros(3,4);
for i=1:4
    corr_mat(i,:) = corrneuron(permuted_results(:,:,i), permuted_results(:,:,4));
end

% Plot correlation data
figure;
rows = 2; columns = 3;

subplot(rows,columns,1);
plot(corr_mat(1,:));
title('CN');
subplot(rows,columns,2);
plot(corr_mat(2,:));
title('MDN');
subplot(rows,columns,3);
plot(corr_mat(3,:));
title('PMC');

% Plot target data
figure;
target = load('fmri/targetFMRICorrelations.mat');

subplot(rows,columns,1);
plot(target.CN);
title('CN');
subplot(rows,columns,2);
plot(target.MDN);
title('MDN');
subplot(rows,columns,3);
plot(target.PMC);
title('PMC');