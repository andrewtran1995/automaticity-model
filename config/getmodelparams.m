function [param_struct] = getmodelparams(configuration)
%GETAUTOPARAMS Get parameters for automaticityModel based on configuration
%   Get parameters for automaticityModel. Parameters initialized here are
%   either dependent on the chosen configuration or are exposed in the
%   param_struct for the purpose of optimizing their value.
    % Load dependencies
    addpath('classes');
    % Initialize parameters that depend on configuration
    param_names   = {'PRE_LEARNING_TRIALS'; 'LEARNING_TRIALS'; 'POST_LEARNING_TRIALS'; 'PFC_DECISION_PT'; 'PMC_DECISION_PT'; 'MC_DECISION_PT'};
    MADDOX_CONFIG = {                    0;               500;                      0;                 4;                 4;                4};
    WALLIS_CONFIG = {                  100;               200;                    100;               400;               400;              400};
    FMRI_CONFIG   = {                    0;             11520;                      0;               700;               700;              700};
    IMAGE_CONFIG  = {                    0;             11520;                      0;               700;               700;              700};
    switch configuration
        case AutomaticityConfiguration.MADDOX
            param_vals = MADDOX_CONFIG;
        case AutomaticityConfiguration.WALLIS
            param_vals = WALLIS_CONFIG;
        case AutomaticityConfiguration.FMRI
            param_vals = FMRI_CONFIG;
        case AutomaticityConfiguration.IMAGE
            param_vals = IMAGE_CONFIG;
        otherwise
            error('Improper configuration requested in get_parameters(configuration): %s!', configuration);
    end
    % Initialize configuration-agnostic parameters used in optimization
    agn_names = {'HEB_CONSTS';'NMDA';'AMPA';'W_MAX';'NOISE_PFC';'NOISE_PMC';'NOISE_MC';'PMC_A_W_OUT';'PMC_B_W_OUT'};
    agn_vals  = {       1e-12;  400;      0;     10;        2.4;          2;         5;            1;            1};
    % FROST Params
    frost_names = {'PFC_A_W_OUT_MDN';'PFC_B_W_OUT_MDN';'DRIV_PFC_W_OUT';'MDN_A_W_OUT';'MDN_B_W_OUT'};
    frost_vals  = {                5;                5;               5;            5;            5};
    % COVIS Params
    covis_names = {'COVIS_DELTA_C';'COVIS_DELTA_E';'COVIS_PERSEV';'COVIS_LAMBDA'};
    covis_vals  = {              1;            3.5;           3.5;           3.5};
    param_struct = cell2struct(vertcat( param_vals,  agn_vals,  frost_vals, covis_vals), ...  % Parameter values
                               vertcat(param_names, agn_names, frost_names, covis_names), ... % Parameter names
                               1);
end

