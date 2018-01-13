function [ autoparm ] = autohebrnd( rndval_1, rndval_2 )
%AUTOPARMRND Return randomly generated value for automaticity parameter
%   Uses input arguments to randomly generate a value for use in the
%   automaticity model. Random values are distributed uniformly across
%   log-space in the range defined implicitly by the exponent.
    coeff = rndval_1;
    exponent = rndval_2;
    lb = exponent - 1;
    ub = exponent + 1;
    autoparm = coeff * 10^(lb + (ub-lb)*rand);
    
    % Ensure parm is positive
    autoparm = max(eps,autoparm);
    
    %% Alternative random value
    % Backwards compatibility for target vectors before 1/5/2018
    % autoparm = max(eps, normrnd(rndval_1, rndval_2));
end

