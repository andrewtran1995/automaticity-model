function [params] = getconfigparams(configuration)
%GETPARAMS Get parameters that are dependent on configuration.
    param_names          = {'PRE_LEARNING_TRIALS'; 'LEARNING_TRIALS'; 'POST_LEARNING_TRIALS'; 'PFC_DECISION_PT'; 'PMC_DECISION_PT'; 'MC_DECISION_PT'};
    IMAGE_CORR_CONFIG    = {                  100;               200;                    100;               400;               400;              400};
    BUTTON_SWITCH_CONFIG = {                    0;             11520;                      0;               700;               700;              700};
    DUAL_TASK_CONFIG     = {                    0;             11520;                      0;               700;               700;              700};
    switch class(configuration)
        case 'ModelConfigImageCorr'
            param_vals = IMAGE_CORR_CONFIG;
        case 'ModelConfigButtonSwitch'
            param_vals = BUTTON_SWITCH_CONFIG;
        case 'ModelConfigDualTask'
            param_vals = DUAL_TASK_CONFIG;
        otherwise
            error('Improper configuration requested in get_parameters(configuration): %s!', configuration);
    end
    params = cell2struct(vertcat( param_vals), ...  % Parameter values
                         vertcat(param_names), ... % Parameter names
                         1);
end

