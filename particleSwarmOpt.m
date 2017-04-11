% Particle swarm stuff goes here!
% Declare particleswarm arguments, starting with function and num variables
fun = @automaticityModelFast_mex;
nvars = 2;
% Lower and upper bounds of Hebbian coefficient(s) and PMC decision point
lb = [0 + eps,  200];
ub = [   1e-6, 1000];
% Declare optimization options
options = optimoptions(@particleswarm, 'HybridFcn', @fmincon, 'UseParallel', true, 'Display', 'iter', 'MaxStallIterations', 4, 'FunctionTolerance', 1e-8, 'SwarmSize', 12);
% Call particleswarm function (requires Matlab 2014 or higher)
[x,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
