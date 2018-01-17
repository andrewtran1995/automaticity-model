%% Particle swarm script
%{
Calls automaticityModelOpt, a wrapper to automaticityModel_mex, to
globally optimize the function for a set of parameters.
BEFORE RUNNING: automaticityModel_mex must be successfully compiled. In order
to compile it, run automaticityModel_script first, and ensure that the code
is successfully generated (no error report).
If adding an additional parameter in the parameter vector, update argvectortostruct.
%}
%% Script body
% Declare particleswarm arguments, starting with function and num variables
fun = @automaticityModelOpt;
nvars = 18;
% Lower and upper bounds of
% PMC_DECISION_PT, MC_DECISION_PT, NOISE_PFC, NOISE_PMC, NOISE_MC, NDMA, AMPA, W_MAX, PMC_A_W_OUT, PMC_B_W_OUT, rnd_1, rnd_2, rnd_3, rnd_4, DELTA_C, DELTA_E, PERSEV, LAMBDA
lb = [  0,   0, eps, eps, eps,   0,   0,   0, -15, 0, 1, -inf,      0,      0,   0,   0,   0,   0];
ub = [inf, inf, inf, inf, inf, inf, inf, inf,  -5, 1, 9,  inf, 100000, 100000, inf, inf, inf, inf];
% Declare optimization options
options = optimoptions(@particleswarm, 'UseParallel', true, ...
                                       'Display', 'iter', ...
                                       'MaxStallIterations', 6, ...
                                       'FunctionTolerance', 1e-13, ...
                                       'SwarmSize', 12);
% Call particleswarm function
[x,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
