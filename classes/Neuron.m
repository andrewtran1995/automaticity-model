classdef Neuron
	%Represents a neuron, carrying commonly used properties

	properties
        out % output vector
		v % voltage matrix (positive)
		u % voltage matrix (negative)
        n % number of iterations (of time) during a trial
		i = 1% time iteration during a trial
        TAU
        LAMBDA
        LAMBDA_PRECALC
		spikes = 0 % spiking rate per trial
    end
    
    methods
        function obj = Neuron(n, TAU, LAMBDA)
           obj.n = n;
           obj.TAU = TAU;
           obj.LAMBDA = LAMBDA;

           % Pre-calcate the lambda vector for performance reasons
           t = (0:n)';
           obj.LAMBDA_PRECALC = (t/LAMBDA).*exp((LAMBDA-t)/LAMBDA);
           
           obj.out = zeros(n,1);
           obj.u = zeros(n,1);
        end
        function obj = reset(obj)
            obj.i = 1;
            obj.out(:) = 0;
            obj.u(:) = 0;
            obj.spikes = 0;
        end
        function arr = pos_volt(obj)
            arr = obj.v;
            arr(obj.v > 0) = obj.v(obj.v > 0);
        end
    end
end