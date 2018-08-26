classdef PMCNeuron
    %Primary Motor Cortex Neuron
    %   This class represents a neuron from the Primary Motor Cortex (PMC).
    
    properties
        W_OUT
        out
        spikes = 0
        v
        u
        pos_volt
        v_stim = 0
        weights
        weights_avg
    end
    
    properties (Constant)
        V_SCALE = 1 % can use to scale PMC visual input value if it comes out too high
        W_LI = 2 % lateral inhibition between PMC A / PMC B
        INIT_WEIGHT = 0.08
    end
    
    methods
        function obj = PMC(W_OUT, n, trials, SAVE_MEM, COVIS_ENABLED, GRID_SIZE)
            obj.W_OUT = W_OUT;
            obj.out = zeros(n,1);
            obj.v = repmat(RSN.rv,n,1);
            obj.u = zeros(n,1);
            obj.pos_volt = zeros(n,1);
            
            % Create weights matrix condiitonally
            if SAVE_MEM
                weight_length = 1;
            else
                weight_length = trials;
            end
            if COVIS_ENABLED
                obj.weights = PMC.INIT_WEIGHT*ones(GRID_SIZE, GRID_SIZE, weight_length, 4);
            else
                obj.weights = PMC.INIT_WEIGHT*ones(GRID_SIZE, GRID_SIZE, weight_length);
            end
            obj.weights_avg = zeros(trials,1);
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

