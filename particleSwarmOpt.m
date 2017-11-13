%% Particle swarm script
%{
Calls automaticityModelOpt, a wrapper to automaticityModelFast_mex, to
globally optimize the function for a set of parameters.
BEFORE RUNNING: automaticityModelFast_mex must be successfully compiled. In order
to compile it, run automaticityModelFast_script first, and ensure that the code
is successfully generated (no error report).
%}
%% Script body
% Declare particleswarm arguments, starting with function and num variables
fun = @automaticityModelOpt;
nvars = 14;
% Lower and upper bounds of
% PMC_DECISION_PT, NOISE_PFC, NOISE_PMC, NDMA, AMPA, W_MAX, mu_1, sigma_1, mu_2, sigma_2, DELTA_C, DELTA_E, PERSEV, LAMBDA
% lb = [700, eps, eps, 600, 0, 10, 1e-10, 1e-10,  1, 1e-10, 10, 1, 5, 1];
% ub = [700,   5,   5, 600, 0, 10,  1e-7,  1e-7, 10,  1e-7, 10, 1, 5, 1];
lb = [700, eps, eps, 600, 0,  0, 1e-10, 1e-10,  1,  1,  1,  1,  1,  1];
ub = [700,   5,   5, 600, 0, 10,  1e-7,  1e-7, 10, 10, 10, 10, 10, 10];
% Declare optimization options
options = optimoptions(@particleswarm, 'UseParallel', true, ...
                                       'Display', 'iter', ...
                                       'MaxStallIterations', 6, ...
                                       'FunctionTolerance', 1e-13, ...
                                       'SwarmSize', 12);
% Call particleswarm function
[x,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
