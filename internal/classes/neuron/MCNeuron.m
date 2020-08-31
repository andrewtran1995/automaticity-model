classdef MCNeuron < HebbianLearningNeuron
    %Motor Cortex Neuron
    %   This class represents a neuron from the motor cortex (MC).
    
    properties
        W_OUT = 0
        v_stim = 0
    end
    
    properties (Constant)
        V_SCALE = 1
        W_LI = 2
        INIT_WEIGHT = 1
        W_MAX = 100;
        NOISE = getconstants().NOISE_MC
        WEIGHT = struct('PRIMARY', 0.9, 'SECONDARY', 0.1)
    end
    
    methods
        function obj = MCNeuron(trials, max_weight, hebbianConsts)
            obj@HebbianLearningNeuron(hebbianConsts, max_weight);
            obj.v = repmat(obj.rv,obj.n,1);
            obj.weights = obj.INIT_WEIGHT*ones(2,trials);
        end

        function obj = iterate(obj, TRIAL, MC_OTHER, PMC, PMC_OTHER)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            if TRIAL == 1
                OBJ_WEIGHT = obj.INIT_WEIGHT;
                OTHER_WEIGHT = obj.INIT_WEIGHT;
            else
                OBJ_WEIGHT = obj.weights(1, TRIAL-1);
                OTHER_WEIGHT = MC_OTHER.weights(2, TRIAL-1);
            end

            obj.v(i+1) = obj.v(i) ...
                       + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + (PMC.W_OUT*obj.WEIGHT.PRIMARY*OBJ_WEIGHT*PMC.out(i) + PMC_OTHER.W_OUT*obj.WEIGHT.SECONDARY*OTHER_WEIGHT*PMC_OTHER.out(i)) - obj.W_LI*MC_OTHER.out(i) )/obj.C ...
                       + normrnd(0, obj.NOISE);
            obj.u(i+1) = obj.u(i)+TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1) >= obj.vpeak
                obj.v(i:i+1) = [obj.vpeak, obj.c];
                obj.u(i+1) = obj.u(i+1)+ obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end

        function weights = previousweights(obj, trial)
            if trial == 1
                previous_trial = 1;
            else
                previous_trial = trial - 1;
            end
            weights = obj.weights(:, previous_trial);
        end
        
        function obj = doHebbianLearning(obj, config, preSynapticOutput)
            % Default the scalar matrix to the constant weight values.
            obj.weights(:,config.trial) = calcHebbianWeights( ...
                obj, ...
                config, ...
                [obj.WEIGHT.PRIMARY; obj.WEIGHT.SECONDARY], ...
                preSynapticOutput ...
            );
        end
        
        function w = weightsForTrial(obj, config)
            w = obj.weights(:,config.trial);
        end
    end
    
end

