classdef (Abstract) Neuron
	%Represents a neuron, carrying commonly used properties

	properties
        out % output vector
        v % voltage matrix (positive)
        u % voltage matrix (negative)
        i = 1% time iteration during a trial
        LAMBDA_PRECALC
        spikes = 0 % spiking rate per trial
    end
    
    properties (Constant)
        n      = 1000  % Number of iterations (of time) during a trial.
        TAU    = 1
        LAMBDA = 20
    end
    
    properties (Abstract, Constant)
        rv
    end
    
    properties (Dependent)
        integralPosVolt
    end
    
    methods
        function obj = Neuron()
            % Pre-calcate the lambda vector for performance reasons
            t = (0:obj.n)';
            obj.LAMBDA_PRECALC = coder.const((t/obj.LAMBDA).*exp((obj.LAMBDA-t)/obj.LAMBDA));

            obj.out = zeros(obj.n,1);
            obj.u = zeros(obj.n,1);
        end
        function obj = reset(obj)
            obj.i = 1;
            obj.out(:) = 0;
            obj.v(:) = obj.rv;
            obj.u(:) = 0;
            obj.spikes = 0;
        end
        function arr = posVolt(obj)
            arr = zeros(obj.n,1);
            arr(obj.v > 0) = obj.v(obj.v > 0);
        end
        function val = get.integralPosVolt(obj)
            val = trapz(obj.posVolt());
        end
        function dispVoltage(obj, neuron_name)
            plot(obj.TAU * (1:obj.n), obj.v);
            axis([0 obj.n -100 100]);
            title(sprintf('%s Neuron Voltage', neuron_name));
        end
        function dispOutput(obj, neuron_name)
            plot(obj.TAU * (1:obj.n), obj.out);
            axis([0 obj.n -1 10]);
            title(sprintf('%s Neuron Output', neuron_name));
        end
    end
end