function [ opt_val ] = automaticityModelOpt( arg_vector )
%AUTOMATICITYMODELOPT Cost function for global optimization
%   Wrapper for automaticityModel with a function signature that adheres to
%   global optimization function requirements
    GROUP_SIZE = 12;
    PFC = zeros(GROUP_SIZE, 4);
    CN  = zeros(GROUP_SIZE, 4);
    GP  = zeros(GROUP_SIZE, 4);
    MDN = zeros(GROUP_SIZE, 4);
    PMC = zeros(GROUP_SIZE, 4);
    acc = zeros(GROUP_SIZE, 4);
    
    parfor i = 1:GROUP_SIZE
        % Get default parameters for configuration (hardcoded to FMRI)
        params = getAutomaticityParams('FMRI');

        % Set static values
        params.PMC_DECISION_PT   = arg_vector(1);
        params.NOISE             = arg_vector(2);
        params.NMDA              = arg_vector(3);
        params.AMPA              = arg_vector(4);
        params.W_MAX             = arg_vector(5);

        % Get normal distribution parameters for certain model parms
        mu_1                     = arg_vector(6);
        sigma_1                  = arg_vector(7);
        mu_2                     = arg_vector(8);
        sigma_2                  = arg_vector(9);
        params.HEB_CONSTS        = normrnd(mu_1, sigma_1);
        params.PFC_A_W_OUT_MDN   = normrnd(mu_2, sigma_2);
        params.PFC_B_W_OUT_MDN   = normrnd(mu_2, sigma_2);
        params.DRIV_PFC_W_OUT    = normrnd(mu_2, sigma_2);
        params.MDN_A_W_OUT       = normrnd(mu_2, sigma_2);
        params.MDN_B_W_OUT       = normrnd(mu_2, sigma_2);

        % Populate second function argument
        optional_parms = struct('FMRI_META_GROUP_RUN', 1, ...
                                'VIS_INPUT_FROM_PARM', 0, ...
                                'visualinput', zeros(2));

        % Call automaticity model function
        [~, retval] = automaticityModelFast_mex(params, optional_parms);
        
        % Assign values for correlation
        PFC(i) = retval(1);
        CN(i)  = retval(2);
        GP(i)  = retval(3);
        MDN(i) = retval(4);
        PMC(i) = retval(5);
        acc(i) = retval(6);
    end
    
    % Calculate SSE
    corr_neuron = @(x) [corr(x(1), acc(1)), corr(x(2), acc(2)), corr(x(3), acc(3)), corr(x(4), acc(4))];
    actual_corr = [arrayfun(corr_neuron, PFC); ...
                   arrayfun(corr_neuron, CN); ...
                   arrayfun(corr_neuron, GP); ...
                   arrayfun(corr_neuron, MDN); ...
                   arrayfun(corr_neuron, PMC)];
    target = load('fmri/targetFMRICorrelations.mat');
    target_diff = actual_corr - [target.PFC; target.CN; target.GP; target.MDN; target.PMC];
    opt_val = sum(sum(target_diff.^2));
    keyboard;
end