classdef MCNeuron
    %Motor Cortex Neuron
    %   This class represents a neuron from the motor cortex (MC).
    
    properties
        W_OUT = 0
        out
        spikes = 0
        v
        u
        pos_volt
        weights
        v_stim = 0
    end
    
    properties (Constant)
        V_SCALE = 1
        W_LI = 2
        INIT_WEIGHT = 1
        W_MAX = 100;
    end
    
    methods
        function obj = MC(n, trials)
            obj.out = zeros(n,1);
            obj.v = repmat(RSN.rv,n,1);
            obj.u = zeros(n,1);
            obj.pos_volt = zeros(n,1);
            obj.weights = MC.INIT_WEIGHT*ones(2,trials);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.pos_volt(:) = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
        end
    end
    
end

