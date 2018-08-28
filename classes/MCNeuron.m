classdef MCNeuron < Neuron
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
            obj@Neuron(n, TAU, LAMBDA);
            obj.v = repmat(RSN.rv,n,1);
            obj.weights = obj.INIT_WEIGHT*ones(2,trials);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
            obj.restartTime();
        end

        function obj = iterate(obj, TRIAL, NOISE_MC, MC_OTHER, PMC, PMC_OTHER, MC_PRIMARY_WEIGHT, MC_SECONDARY_WEIGHT)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;

            obj.v(i+1)=(obj.v(i) + TAU*(RSN.k*(obj.v(i)-RSN.rv)*(obj.v(i)-RSN.vt)-obj.u(i)+ RSN.E + (PMC.W_OUT*MC_PRIMARY_WEIGHT*obj.weights(1,TRIAL)*PMC.out(i) + PMC_OTHER.W_OUT*MC_SECONDARY_WEIGHT*MC_OTHER.weights(2,TRIAL)*PMC_OTHER.out(i)) - obj.W_LI*MC_OTHER.out(i) )/RSN.C) + normrnd(0,NOISE_MC);
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

