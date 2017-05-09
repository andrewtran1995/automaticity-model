% Particle swarm stuff goes here!
% Declare particleswarm arguments, starting with function and num variables
fun = @automaticityModelFast_mex;
nvars = 7;
% Lower and upper bounds of
% HEB_CONSTS, ANTI_HEB_CONSTS, PMC_DECISION_PT, NOISE, NMDA, AMPA, W_MAX
lb = [ 1e-16, 1e-16,   10, 0 + eps, 0 + eps,    0,  0];
ub = [  1e-3,  1e-3, 2000,       5,    2000, 2000, 50];
% Declare optimization options
% options = optimoptions(@particleswarm, 'HybridFcn', @fmincon, 'UseParallel', true, 'Display', 'iter', 'MaxStallIterations', 4, 'FunctionTolerance', 1e-8, 'SwarmSize', 12);
options = optimoptions(@particleswarm, 'UseParallel', true, 'Display', 'iter', 'MaxStallIterations', 6, 'FunctionTolerance', 1e-13, 'SwarmSize', 12);
% Call particleswarm function (requires Matlab 2014 or higher)
[x,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);