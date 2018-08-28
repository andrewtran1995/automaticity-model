classdef PMCNeuron < Neuron
    %Primary Motor Cortex Neuron
    %   This class represents a neuron from the Primary Motor Cortex (PMC).
    
    properties
        W_OUT
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
        function obj = PMCNeuron(n, TAU, LAMBDA, trials, W_OUT, SAVE_MEM, COVIS_ENABLED, GRID_SIZE)
            obj@Neuron(n, TAU, LAMBDA);
            obj.W_OUT = W_OUT;
            obj.v = repmat(RSN.rv,n,1);
            
            % Create weights matrix condiitonally
            if SAVE_MEM
                weight_length = 1;
            else
                weight_length = trials;
            end
            if COVIS_ENABLED
                obj.weights = obj.INIT_WEIGHT*ones(GRID_SIZE, GRID_SIZE, weight_length, 4);
            else
                obj.weights = obj.INIT_WEIGHT*ones(GRID_SIZE, GRID_SIZE, weight_length);
            end
            obj.weights_avg = zeros(trials,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
            obj.restartTime();
        end

        function obj = iterate(obj, NOISE_PMC, PMC_OTHER, PFC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=(obj.v(i) + TAU*(RSN.k*(obj.v(i)-RSN.rv)*(obj.v(i)-RSN.vt)-obj.u(i)+ RSN.E + obj.v_stim + PFC.W_OUT*PFC.out(i) - obj.W_LI*PMC_OTHER.out(i) )/RSN.C) + normrnd(0,NOISE_PMC);
            obj.u(i+1)=obj.u(i)+TAU*RSN.a*(RSN.b*(obj.v(i)-RSN.rv)-obj.u(i));
            if obj.v(i+1)>=RSN.vpeak
                obj.v(i)= RSN.vpeak;
                obj.v(i+1)= RSN.c;
                obj.u(i+1)= obj.u(i+1)+ RSN.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

