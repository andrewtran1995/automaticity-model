function [ opt_val ] = automaticityModelOpt( arg_vector )
%AUTOMATICITYMODELOPT Cost function for global optimization
%   Wrapper for automaticityModel with a function signature that adheres to
%   global optimization function requirements
    optional_parms = struct('FMRI_META_GROUP_RUN', 0, ...
                            'VIS_INPUT_FROM_PARM', 0, ...
                            'visualinput', zeros(2));
    [opt_val, ~] = automaticityModelFast_mex(arg_vector, optional_parms);
end

