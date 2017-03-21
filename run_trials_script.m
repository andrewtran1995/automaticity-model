% RUN_TRIALS_SCRIPT   Generate MEX-function run_trials_mex from run_trials.
% 
% Script generated from project 'run_trials.prj' on 19-Mar-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.MATLABSourceComments = true;
cfg.GenerateReport = true;
cfg.EnableDebugging = true;
cfg.EnableJIT = true;

%% Define argument types for entry-point 'run_trials'.
ARGS = cell(1,1);
ARGS{1} = cell(11,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.V_SCALE = coder.typeof(0);
ARGS{1}{1}.W_LI = coder.typeof(0);
ARGS{1}{1}.DECISION_PT = coder.typeof(0);
ARGS{1}{1}.rx_matrix = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1});
ARGS{1}{1} = coder.cstructname(ARGS{1}{1},'PFC');
ARGS{1}{2} = struct;
ARGS{1}{2}.V_SCALE = coder.typeof(0);
ARGS{1}{2}.W_LI = coder.typeof(0);
ARGS{1}{2}.DECISION_PT = coder.typeof(0);
ARGS{1}{2}.rx_matrix = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2});
ARGS{1}{2} = coder.cstructname(ARGS{1}{2},'PMC');
ARGS{1}{3} = struct;
ARGS{1}{3}.W_OUT = coder.typeof(0);
ARGS{1}{3}.out = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3}.out_all = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{3}.spikes = coder.typeof(0);
ARGS{1}{3}.v = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3}.u = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3}.pos_volt = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3}.v_stim = coder.typeof(0);
ARGS{1}{3} = coder.typeof(ARGS{1}{3});
ARGS{1}{3} = coder.cstructname(ARGS{1}{3},'PFC_A');
ARGS{1}{4} = struct;
ARGS{1}{4}.W_OUT = coder.typeof(0);
ARGS{1}{4}.out = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.out_all = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{4}.spikes = coder.typeof(0);
ARGS{1}{4}.v = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.u = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.pos_volt = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.v_stim = coder.typeof(0);
ARGS{1}{4} = coder.typeof(ARGS{1}{4});
ARGS{1}{4} = coder.cstructname(ARGS{1}{4},'PFC_B');
ARGS{1}{5} = struct;
ARGS{1}{5}.W_OUT = coder.typeof(0);
ARGS{1}{5}.out = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{5}.out_all = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{5}.spikes = coder.typeof(0);
ARGS{1}{5}.v = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{5}.u = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{5}.pos_volt = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{5}.v_stim = coder.typeof(0);
ARGS{1}{5}.weights = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{5}.weights_avg = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{5} = coder.typeof(ARGS{1}{5});
ARGS{1}{5} = coder.cstructname(ARGS{1}{5},'PMC_A');
ARGS{1}{6} = struct;
ARGS{1}{6}.W_OUT = coder.typeof(0);
ARGS{1}{6}.out = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{6}.out_all = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{6}.spikes = coder.typeof(0);
ARGS{1}{6}.v = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{6}.u = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{6}.pos_volt = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{6}.v_stim = coder.typeof(0);
ARGS{1}{6}.weights = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{6}.weights_avg = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{6} = coder.typeof(ARGS{1}{6});
ARGS{1}{6} = coder.cstructname(ARGS{1}{6},'PMC_B');
ARGS{1}{7} = struct;
ARGS{1}{7}.C = coder.typeof(0);
ARGS{1}{7}.rv = coder.typeof(0);
ARGS{1}{7}.vt = coder.typeof(0);
ARGS{1}{7}.k = coder.typeof(0);
ARGS{1}{7}.a = coder.typeof(0);
ARGS{1}{7}.b = coder.typeof(0);
ARGS{1}{7}.c = coder.typeof(0);
ARGS{1}{7}.d = coder.typeof(0);
ARGS{1}{7}.vpeak = coder.typeof(0);
ARGS{1}{7}.E = coder.typeof(0);
ARGS{1}{7} = coder.typeof(ARGS{1}{7});
ARGS{1}{7} = coder.cstructname(ARGS{1}{7},'RSN');
ARGS{1}{8} = coder.typeof(0);
ARGS{1}{9} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10} = coder.typeof(0);
ARGS{1}{11} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg run_trials -args ARGS{1}

