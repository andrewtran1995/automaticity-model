RUNS = 8;
PMC_BOLD_S1  = zeros(RUNS,1);
PMC_BOLD_S4  = zeros(RUNS,1);
PMC_BOLD_S10 = zeros(RUNS,1);
PMC_BOLD_S20 = zeros(RUNS,1);
Accuracy_S1  = zeros(RUNS,1);
Accuracy_S4  = zeros(RUNS,1);
Accuracy_S10 = zeros(RUNS,1);
Accuracy_S20 = zeros(RUNS,1);

loaded_value = load('fmri/particleSwarmX_17_6_2.mat');
arg_vector = loaded_value.x;

gcp;
parfor i=1:RUNS
    [r_x_vals, r_y_vals] = createFMRIInput(11520);
    visualinput = [r_x_vals, r_y_vals];
    optional_parms = struct('FMRI_META_GROUP_RUN', 1, ...
                            'VIS_INPUT_FROM_PARM', 1, ...
                            'visualinput', visualinput);
    [unused_val,results] = automaticityModelFast_mex(arg_vector, optional_parms);
    PMC_BOLD_S1(i)  = results(1);
    PMC_BOLD_S4(i)  = results(2);
    PMC_BOLD_S10(i) = results(3);
    PMC_BOLD_S20(i) = results(4);
    Accuracy_S1(i)  = results(5);
    Accuracy_S4(i)  = results(6);
    Accuracy_S10(i) = results(7);
    Accuracy_S20(i) = results(8);
end
delete(gcp('nocreate'));

PMC_corr = [corr(PMC_BOLD_S1,Accuracy_S1), corr(PMC_BOLD_S4,Accuracy_S4), corr(PMC_BOLD_S10,Accuracy_S10), corr(PMC_BOLD_S20,Accuracy_S20)];

figure;
plot(PMC_corr);
title('PMC Correlations (Average Bold vs. Accuracy)');
axis([0 4 -1.1 1.1]);