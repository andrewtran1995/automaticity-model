RUNS = 8;
PMC_BOLD_S1  = zeros(1,RUNS);
PMC_BOLD_S4  = zeros(1,RUNS);
PMC_BOLD_S10 = zeros(1,RUNS);
PMC_BOLD_S20 = zeros(1,RUNS);
Accuracy_S1  = zeros(1,RUNS);
Accuracy_S4  = zeros(1,RUNS);
Accuracy_S10 = zeros(1,RUNS);
Accuracy_S20 = zeros(1,RUNS);

loaded_value = load('fmri/particleSwarmX_17_6_2.mat');
arg_vector = loaded_value.x;

for i=1:RUNS
    [r_x_vals, r_y_vals] = createFMRIInput(11520);
    visualinput = [r_x_vals, r_y_vals];
    optional_parms = struct('FMRI_META_GROUP_RUN', 1, ...
                            'VIS_INPUT_FROM_PARM', 1, ...
                            'visualinput', visualinput);
    [unused_val,results] = automaticityModelFast(arg_vector, optional_parms);
    PMC_BOLD_S1(i)  = results(1);
    PMC_BOLD_S4(i)  = results(2);
    PMC_BOLD_S10(i) = results(3);
    PMC_BOLD_S20(i) = results(4);
    Accuracy_S1(i)  = results(5);
    Accuracy_S4(i)  = results(6);
    Accuracy_S10(i) = results(7);
    Accuracy_S20(i) = results(8);
end