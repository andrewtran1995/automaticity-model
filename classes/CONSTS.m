classdef CONSTS
    properties (Constant)
        n = 1000
        TAU = 1
        LAMBDA = 20
        % Calculate lambda beforehand for performance reasons
        % LAMDA_PRECALC = (t/LAMBDA).*exp((LAMBDA-t)/LAMBDA), where t = (0:n)'
        LAMBDA_PRECALC = ((0:1000)'/20).*exp((20-(0:1000)')/20)
    end
end

