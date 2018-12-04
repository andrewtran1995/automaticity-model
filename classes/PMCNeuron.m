classdef PMCNeuron < RSN
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
            obj@RSN(n, TAU, LAMBDA);
            obj.W_OUT = W_OUT;
            obj.v = repmat(obj.rv,n,1);
            
            % Create weights matrix conditionally
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
            obj = reset@RSN(obj);
        end

        function obj = iterate(obj, NOISE_PMC, PMC_OTHER, PFC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=(obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i)+ obj.E + obj.v_stim + PFC.W_OUT*PFC.out(i) - obj.W_LI*PMC_OTHER.out(i) )/obj.C) + normrnd(0,NOISE_PMC);
            obj.u(i+1)=obj.u(i)+TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1)>=obj.vpeak
                obj.v(i)= obj.vpeak;
                obj.v(i+1)= obj.c;
                obj.u(i+1)= obj.u(i+1)+ obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

