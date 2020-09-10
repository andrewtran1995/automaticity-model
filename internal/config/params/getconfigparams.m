function [params] = getconfigparams(configuration)
%GETPARAMS Get parameters that are dependent on configuration.
    param_names          = {'PRE_LEARNING_TRIALS'; 'LEARNING_TRIALS'; 'POST_LEARNING_TRIALS'; 'W_MAX'};
    ELECTRO_CONFIG       = {                  5;              500;                    5;       5};
    BUTTON_SWITCH_CONFIG = {                    0;             11520;                      0;       5};
    DUAL_TASK_CONFIG     = {                    0;             11520;                      0;       5};
    switch class(configuration)
        case 'ModelConfigElectro'
            param_vals = ELECTRO_CONFIG;
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

