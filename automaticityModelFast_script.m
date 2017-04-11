% AUTOMATICITYMODELFAST_SCRIPT   Generate MEX-function automaticityModelFast_mex
%  from automaticityModelFast.
% 
% Script generated from project 'automaticityModelFast.prj' on 11-Apr-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;

%% Define argument types for entry-point 'automaticityModelFast'.
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARGS{1}{1} = coder.typeof(0,[1 2]);

%% Invoke MATLAB Coder.
codegen -config cfg automaticityModelFast -args ARGS{1}
