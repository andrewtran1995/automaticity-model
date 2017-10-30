function [ opt_val_1, opt_val_2, new_vals ] = automaticityModelRandom( arg_struct, optional_parms )
%AUTOMATICITYMODELRANDOM Automaticity model with pre-defined variability
%   Automaticity model with some pre-defined variability, given the
%   arg_vector, to represent an individual within a group run
    SIGMA_HEBB = 5e-9;
    SIGMA_PERS = 1;
    
    params = getAutomaticityParams('FMRI');
    
    new_vals = zeros(2,1);
    new_vals(1) = normrnd(arg_struct.HEB_CONSTS, SIGMA_HEBB);
    new_vals(2) = normrnd(5, SIGMA_PERS);
    
    params.HEB_CONSTS   = new_vals(1);
    params.COVIS_PERSEV = new_vals(2);
    
    [opt_val_1, opt_val_2] = automaticityModelFast_mex(params, optional_parms);
end

