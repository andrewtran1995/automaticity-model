% AUTOMATICITYMODELFAST_SCRIPT   Generate MEX-function automaticityModelFast_mex
%  from automaticityModelFast.
% 
% Script generated from project 'automaticityModelFast.prj' on 06-Oct-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.
% Generated by command: coder -tocode automaticityModelFast -script
% automaticityModelFast_script.m

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.EnableJIT = true;

%% Define argument types for entry-point 'automaticityModelFast'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[1 7]);
ARGS{1}{2} = struct;
ARGS{1}{2}.FMRI_META_GROUP_RUN = coder.typeof(0);
ARGS{1}{2}.VIS_INPUT_FROM_PARM = coder.typeof(0);
ARGS{1}{2}.visualinput = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2});

%% Invoke MATLAB Coder.
codegen -config cfg automaticityModelFast -args ARGS{1}
