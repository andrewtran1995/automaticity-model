classdef MCNeuron < RSN
    %Motor Cortex Neuron
    %   This class represents a neuron from the motor cortex (MC).
    
    properties
        W_OUT = 0
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
        function obj = MCNeuron(n, TAU, LAMBDA, trials)
            obj@RSN(n, TAU, LAMBDA);
            obj.v = repmat(obj.rv,n,1);
            obj.weights = obj.INIT_WEIGHT*ones(2,trials);
        end
        
        function obj = reset(obj)
            obj = reset@RSN(obj);
        end

        function obj = iterate(obj, TRIAL, NOISE_MC, MC_OTHER, PMC, PMC_OTHER, MC_PRIMARY_WEIGHT, MC_SECONDARY_WEIGHT)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;

            obj.v(i+1)=(obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i)+ obj.E + (PMC.W_OUT*MC_PRIMARY_WEIGHT*obj.weights(1,TRIAL)*PMC.out(i) + PMC_OTHER.W_OUT*MC_SECONDARY_WEIGHT*MC_OTHER.weights(2,TRIAL)*PMC_OTHER.out(i)) - obj.W_LI*MC_OTHER.out(i) )/obj.C) + normrnd(0,NOISE_MC);
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

