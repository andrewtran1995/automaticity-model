function [param_struct] = getAutomaticityParams(configuration)
%GETAUTOMATICITYPARAMS Get parameters for automaticityModel
%   Get parameters based on configuration
    % Necessary for codegen
    % coder.extrinsic('cell2struct');
    % Preinitialize param_struct to allow codegen to infer type
    % param_struct = struct('PRE_LEARNING_TRIALS',0, 'LEARNING_TRIALS',0, 'POST_LEARNING_TRIALS',0, 'NOISE',0, 'PFC_DECISION_PT',0, 'PMC_DECISION_PT',0,'HEB_CONSTS',0,'ANTI_HEB_CONSTS',0,'NMDA',0,'AMPA',0,'W_MAX',0);
    
    % Initialize parameters that depend on configuration
    param_names   = {'PRE_LEARNING_TRIALS'; 'LEARNING_TRIALS'; 'POST_LEARNING_TRIALS'; 'NOISE'; 'PFC_DECISION_PT'; 'PMC_DECISION_PT'};
    MADDOX_CONFIG = {                    0;               500;                      0;       0;                 4;                 4};
    WALLIS_CONFIG = {                  100;               200;                    100;       2;               400;               400};
    FMRI_CONFIG   = {                    0;             11520;                      0;       2;               400;               400};
    if strcmp(configuration,'MADDOX')
        params = MADDOX_CONFIG;
    elseif strcmp(configuration,'WALLIS')
        params = WALLIS_CONFIG;
    elseif strcmp(configuration,'FMRI')
        params = FMRI_CONFIG;
    else
        error('Improper configuration requested in get_parameters(configuration)!');
    end
    % Initialize parameters used in optimization
    % Note that PMC_DECISION_PT & NOISE are used in optimization as well, but their default value is dependent on configuration
    % Therefore, we do not re-initialize it here
    optim_param_names    = {'HEB_CONSTS';'ANTI_HEB_CONSTS';'NMDA';'AMPA';'W_MAX'};
    optim_param_defaults = {        1e-8;             1e-8;  600;      0;     10};
    param_struct = cell2struct(vertcat(params, optim_param_defaults), ...   % Parameter values
                               vertcat(param_names, optim_param_names), ... % Parameter names
                               1);

end

