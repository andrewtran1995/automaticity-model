function [ arg_struct ] = argvectortostruct( arg_vector, configuration )
%ARGVECTORTOSTRUCT Convert arg_vector into arg_struct that is compatible
%with automaticityModel
    % Get default parameters for configuration (hardcoded to FMRI)
    params = getAutomaticityParams(configuration);

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
    params.HEB_CONSTS        = max(eps, normrnd(mu_1, sigma_1));
    params.PFC_A_W_OUT_MDN   = max(eps, normrnd(mu_2, sigma_2));
    params.PFC_B_W_OUT_MDN   = max(eps, normrnd(mu_2, sigma_2));
    params.DRIV_PFC_W_OUT    = max(eps, normrnd(mu_2, sigma_2));
    params.MDN_A_W_OUT       = max(eps, normrnd(mu_2, sigma_2));
    params.MDN_B_W_OUT       = max(eps, normrnd(mu_2, sigma_2));
    
    arg_struct = params;
end

