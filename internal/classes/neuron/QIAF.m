classdef QIAF < Neuron
    %Quadratic Integrate and Fire Neuron
    %   Stimulates neurons in the Globus Pallidus. For FROST.
    properties (Constant)
        beta = 11.83
        gamma = 0.117
        vt = -40
        rv = -60
        vpeak = 35
        vreset = -50
    end

    methods
        function obj = QIAF()
        end
    end
end

