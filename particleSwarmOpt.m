% Declare particleswarm arguments, starting with function and num variables
fun = @automaticityModelOpt;
nvars = 9;
% Lower and upper bounds of
% PMC_DECISION_PT, NOISE, NDMA, AMPA, W_MAX, mu_1, sigma_1, mu_2, sigma_2
lb = [700, 2, 600, 0, 10,  1, 1e-10,  1, 1e-10];
ub = [700, 2, 600, 0, 10, 10,  1e-7, 10,  1e-7];
% Declare optimization options
options = optimoptions(@particleswarm, 'UseParallel', true, ...
                                       'Display', 'iter', ...
                                       'MaxStallIterations', 6, ...
                                       'FunctionTolerance', 1e-13, ...
                                       'SwarmSize', 12);
% Call particleswarm function (requires Matlab 2014 or higher)
[x,fval,exitflag,output] = particleswarm(fun, nvars, lb, ub, options);
