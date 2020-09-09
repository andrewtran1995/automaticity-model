% AUTOMATICITYMODEL_SCRIPT   Generate MEX-function automaticityModel_mex from
%  automaticityModel.
% 5
% Script generated from project 'automaticityModel.prj' on 07-Jan-2018.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.
% Generated by: coder -tocode automaticityModel -script automaticityModel_script.m

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.EnableJIT = true;

%% Define argument types for entry-point 'automaticityModel'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = ModelConfig;
ARGS{1}{2} = struct;
ARGS{1}{2}.PRE_LEARNING_TRIALS = coder.typeof(0);
ARGS{1}{2}.LEARNING_TRIALS = coder.typeof(0);
ARGS{1}{2}.POST_LEARNING_TRIALS = coder.typeof(0);
ARGS{1}{2}.PFC_DECISION_PT = coder.typeof(0);
ARGS{1}{2}.PMC_DECISION_PT = coder.typeof(0);
ARGS{1}{2}.MC_DECISION_PT = coder.typeof(0);
ARGS{1}{2}.HEB_CONSTS = coder.typeof(0);
ARGS{1}{2}.NMDA = coder.typeof(0);
ARGS{1}{2}.AMPA = coder.typeof(0);
ARGS{1}{2}.W_MAX = coder.typeof(0);
ARGS{1}{2}.NOISE_PFC = coder.typeof(0);
ARGS{1}{2}.NOISE_PMC = coder.typeof(0);
ARGS{1}{2}.NOISE_MC = coder.typeof(0);
ARGS{1}{2}.PMC_A_W_OUT = coder.typeof(0);
ARGS{1}{2}.PMC_B_W_OUT = coder.typeof(0);
ARGS{1}{2}.PFC_A_W_OUT_MDN = coder.typeof(0);
ARGS{1}{2}.PFC_B_W_OUT_MDN = coder.typeof(0);
ARGS{1}{2}.DRIV_PFC_W_OUT = coder.typeof(0);
ARGS{1}{2}.MDN_A_W_OUT = coder.typeof(0);
ARGS{1}{2}.MDN_B_W_OUT = coder.typeof(0);
ARGS{1}{2}.COVIS_DELTA_C = coder.typeof(0);
ARGS{1}{2}.COVIS_DELTA_E = coder.typeof(0);
ARGS{1}{2}.COVIS_PERSEV = coder.typeof(0);
ARGS{1}{2}.COVIS_LAMBDA = coder.typeof(0);
ARGS{1}{2} = coder.typeof(ARGS{1}{1});
ARGS{1}{3} = struct;
ARGS{1}{3}.FMRI_META_GROUP_RUN = coder.typeof(0);
ARGS{1}{3}.VIS_INPUT_FROM_PARM = coder.typeof(0);
ARGS{1}{3}.visualinput = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{2});

%% Add necessary folders to path
addpath(genpath('.'));

%% Invoke MATLAB Coder.
codegen -config cfg automaticityModel -args ARGS{1}
