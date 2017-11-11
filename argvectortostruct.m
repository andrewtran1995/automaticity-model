function [ arg_struct ] = argvectortostruct( arg_vector, configuration )
%ARGVECTORTOSTRUCT Convert arg_vector into arg_struct that is compatible
%with automaticityModel
    % Get default parameters for configuration (hardcoded to FMRI)
    params = getAutomaticityParams(configuration);

    % Set static values for base parameters
    params.PMC_DECISION_PT = arg_vector(1);
    params.NOISE_PFC       = arg_vector(2);
    params.NOISE_PMC       = arg_vector(3);
    params.NMDA            = arg_vector(4);
    params.AMPA            = arg_vector(5);
    params.W_MAX           = arg_vector(6);

    % Get normal distribution parameters for certain model parms
    mu_1                   = arg_vector(7);
    sigma_1                = arg_vector(8);
    mu_2                   = arg_vector(9);
    sigma_2                = arg_vector(10);
    params.HEB_CONSTS      = max(eps, normrnd(mu_1, sigma_1));
    params.PFC_A_W_OUT_MDN = max(eps, normrnd(mu_2, sigma_2));
    params.PFC_B_W_OUT_MDN = max(eps, normrnd(mu_2, sigma_2));
    params.DRIV_PFC_W_OUT  = max(eps, normrnd(mu_2, sigma_2));
    params.MDN_A_W_OUT     = max(eps, normrnd(mu_2, sigma_2));
    params.MDN_B_W_OUT     = max(eps, normrnd(mu_2, sigma_2));
    
    % Set static values for COVIS
    params.COVIS_DELTA_C   = arg_vector(11);
    params.COVIS_DELTA_E   = arg_vector(12);
    params.COVIS_PERSEV    = arg_vector(13);
    params.COVIS_LAMBDA    = arg_vector(14);
    
    arg_struct = params;
end

