function [params] = getmodelparams(configuration)
%GETAUTOPARAMS Get parameters for automaticityModel based on configuration
%   Get parameters for automaticityModel. Parameters initialized here are
%   either dependent on the chosen configuration or are exposed in the
%   param_struct for the purpose of optimizing their value.
    addpath(genpath('.'));
    params = mergestructs(getconfigparams(configuration), getconstants());
end

