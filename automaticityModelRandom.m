function [ opt_val_1, opt_val_2, new_vals ] = automaticityModelRandom( arg_vector, optional_parms )
%AUTOMATICITYMODELRANDOM Automaticity model with pre-defined variability
%   Automaticity model with some pre-defined variability, given the
%   arg_vector, to represent an individual within a group run
    SIGMA_HEBB = 5e-9;
    SIGMA_ANTI_HEBB = 5e-8;
    SIGMA_PERS = 1;
    
    new_vals = zeros(3,1);
    new_vals(1) = normrnd(arg_vector(1), SIGMA_HEBB);
    new_vals(2) = normrnd(arg_vector(2), SIGMA_ANTI_HEBB);
    new_vals(3) = normrnd(5, SIGMA_PERS);
    
    arg_vector(1) = new_vals(1);
    arg_vector(2) = new_vals(2);
    optional_parms.COVIS_PERSEV_PARAM = new_vals(3);
    
    [opt_val_1, opt_val_2] = automaticityModelFast_mex(arg_vector, optional_parms);
end

