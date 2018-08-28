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
end