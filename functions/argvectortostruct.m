function [ arg_struct ] = argvectortostruct( arg_vector, configuration )
%ARGVECTORTOSTRUCT Convert arg_vector into arg_struct compatible w/ automaticityModel
    % Get default parameters for configuration
    params = getmodelparams(configuration);

    % Set static values for base parameters
    params.PMC_DECISION_PT = arg_vector(1);
    params.MC_DECISION_PT  = arg_vector(2);
    params.NOISE_PFC       = arg_vector(3);
    params.NOISE_PMC       = arg_vector(4);
    params.NOISE_MC        = arg_vector(5);
    params.NMDA            = arg_vector(6);
    params.AMPA            = arg_vector(7);
    params.W_MAX           = arg_vector(8);
    params.PMC_A_W_OUT     = arg_vector(9);
    params.PMC_B_W_OUT     = arg_vector(10);

    % Get normal distribution parameters for certain model parms
    rnd_1                  = arg_vector(10);
    rnd_2                  = arg_vector(11);
    rnd_3                  = arg_vector(12);
    rnd_4                  = arg_vector(13);
    params.HEB_CONSTS      = autohebrnd(rnd_1, rnd_2);
    params.PFC_A_W_OUT_MDN = autowmrnd(rnd_3, rnd_4);
    params.PFC_B_W_OUT_MDN = autowmrnd(rnd_3, rnd_4);
    params.DRIV_PFC_W_OUT  = autowmrnd(rnd_3, rnd_4);
    params.MDN_A_W_OUT     = autowmrnd(rnd_3, rnd_4);
    params.MDN_B_W_OUT     = autowmrnd(rnd_3, rnd_4);
    
    % Set static values for COVIS
    params.COVIS_DELTA_C   = arg_vector(14);
    params.COVIS_DELTA_E   = arg_vector(15);
    params.COVIS_PERSEV    = arg_vector(16);
    params.COVIS_LAMBDA    = arg_vector(17);
    
    arg_struct = params;
end

