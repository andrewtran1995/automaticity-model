classdef PMCNeuron < HebbianLearningNeuron
    %Primary Motor Cortex Neuron
    %   This class represents a neuron from the Primary Motor Cortex (PMC).
    
    properties
        W_OUT
        v_stim = 0
        weights_avg
    end
    
    properties (Constant)
        V_SCALE = 1 % can use to scale PMC visual input value if it comes out too high
        W_LI = 5 % lateral inhibition between PMC A / PMC B
        INIT_WEIGHT = 0.08
        NOISE = getconstants().NOISE_PMC
    end
    
    methods
        function obj = PMCNeuron(trials, W_OUT, max_weight, COVIS_ENABLED, GRID_SIZE, hebbianConsts)
            obj@HebbianLearningNeuron(hebbianConsts, max_weight);
            obj.W_OUT = W_OUT;
            obj.v = repmat(obj.rv,obj.n,1);
            
            % Create weights matrix conditionally
            if trials > HebbianLearningNeuron.LARGE_TRIAL_BOUNDARY
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

        function obj = iterate(obj, PMC_OTHER, PFC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1) = obj.v(i) ...
                       + TAU * ( ...
                                 obj.k * (obj.v(i) - obj.rv) * (obj.v(i) - obj.vt) ...
                                 - obj.u(i) ...
                                 + obj.E ...
                                 + obj.v_stim ...
                                 + PFC.W_OUT*PFC.out(i) ...
                                 - obj.W_LI*PMC_OTHER.out(i) ...
                        ) / obj.C ...
                       + normrnd(0, obj.NOISE);
            obj.u(i+1)=obj.u(i)+TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1)>=obj.vpeak
                obj.v(i:i+1) = [obj.vpeak, obj.c];
                obj.u(i+1)= obj.u(i+1)+ obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
        
        % The pre-synaptic output is directly related to the visual input
        % (represented by the radial basis vector). More precisely, the
        % pre-synaptic output to do Hebbian learning on are the existing
        % PMC weights multiplied by Neuron.n (the amount of time where the
        % visual stimulus is exposed), which is then multiplied by the
        % scalar (the radial basis vector).
        function obj = doHebbianLearning(obj, config, scalar, preSynapticOutput)
            w = obj.calcHebbianWeights(config, scalar, preSynapticOutput);
            if config.isCOVISEnabled
                obj.weights(:,:,config.weightIdx,config.COVISRules.chosen) = w;
            else
                obj.weights(:,:,config.weightIdx,1) = w;
            end
        end
        
        function w = weightsForTrial(obj, config)
            % If first trial, or weights have no time dimension, set to
            % initial weights.
            if config.trials > HebbianLearningNeuron.LARGE_TRIAL_BOUNDARY || config.trial==1
                if config.isCOVISEnabled
                    w = obj.weights(:,:,1,config.COVISRules.chosen);
                else
                    w = obj.weights(:,:,1);
                end
            % Else, set weights to results of previous trial.    
            else
                if config.isCOVISEnabled
                    w = obj.weights(:,:,config.trial-1,config.COVISRules.chosen);
                else
                    w = obj.weights(:,:,config.trial-1);
                end
            end
        end
    end
    
end

