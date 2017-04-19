% Particle swarm stuff goes here!
% Declare particleswarm arguments, starting with function and num variables
fun = @automaticityModelFast_mex;
nvars = 6;
% Lower and upper bounds of
% heb_consts, anti_heb_consts, decision_pt,
% noise, nmda, ampa
lb = [ 1e-16, 1e-16,   10, 0 + eps, 0 + eps,    0];
ub = [  1e-6,  1e-6, 2000,       5,    2000, 2000];
% Declare optimization options
% options = optimoptions(@particleswarm, 'HybridFcn', @fmincon, 'UseParallel', true, 'Display', 'iter', 'MaxStallIterations', 4, 'FunctionTolerance', 1e-8, 'SwarmSize', 12);
options = optimoptions(@particleswarm, 'UseParallel', true, 'Display', 'iter', 'MaxStallIterations', 6, 'FunctionTolerance', 1e-13, 'SwarmSize', 12);
% Call particleswarm function (requires Matlab 2014 or higher)
[x,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
