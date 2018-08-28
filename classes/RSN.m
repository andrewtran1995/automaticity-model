classdef RSN < Neuron
    %Regular Spiking Neuron
    %   Cortical regular spiking neuron.
    properties (Constant)
        C = 100
        rv = -60
        vt = -40
        k = 0.7
        a = 0.03
        b = -2
        c = -50
        d = 100
        vpeak = 35
        E = 60
    end

    methods
        function obj = RSN(n, TAU, LAMBDA)
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