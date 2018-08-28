classdef MSN < Neuron
    %Medium Spiny Neuron
    %   Medium spiny neuron in caudate nucleus.
    properties (Constant)
        C = 50
        rv = -80
        vt = -25
        k = 1
        a = 0.01
        b = -20
        c = -55
        d = 150
        vpeak = 40
        E = 100
    end

    methods
        function obj = MSN(n, TAU, LAMBDA)
            obj@Neuron(n, TAU, LAMBDA);
        end
        function obj = reset(obj)
            obj.v(:) = obj.rv;
            obj.u(:) = 0;
            obj.spikes = 0;
            obj.out(:) = 0;
            obj.restartTime();
        end
    end
end

