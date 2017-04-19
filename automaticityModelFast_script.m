% AUTOMATICITYMODELFAST_SCRIPT   Generate MEX-function automaticityModelFast_mex
%  from automaticityModelFast.
% 
% Script generated from project 'automaticityModelFast.prj' on 11-Apr-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.
%
% Usage: run this script in order to compile automaticityModelFast.m. This
% should result in a file named "automaticityModelFast_mex.*", where the
% file extension is dependent on the type of system this is run on.
% Note that ARGS is a structure that defines the input type of
% automaticityModelFast, and must be of the same length as the input
% vector.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;

%% Define argument types for entry-point 'automaticityModelFast'.
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARGS{1}{1} = coder.typeof(0,[1 6]);

%% Invoke MATLAB Coder.
codegen -config cfg automaticityModelFast -args ARGS{1}
