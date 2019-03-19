function [ params_valid ] = areparamsvalid( param_struct )
%VALIDATEPARAMS Validates that params are valid biologically/functionally
    params_valid = all(structfun(@(x) x >= 0, param_struct));
    return;
end