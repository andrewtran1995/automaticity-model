classdef (Abstract) HebbianLearningNeuron < Neuron
    %HEBBIANLEARNINGNEURON A neuron that engages in learning per Hebbian
    %theory.
    
    properties
        weights
    end
    
    properties (SetAccess = immutable)
        HEBBIAN (1,1) HebbianConst
        MAX_WEIGHT
    end
    
    properties (Constant)
        % Maximum number of trials where per-trial weight information can
        % still be stored until memory constraints are reached.
        LARGE_TRIAL_BOUNDARY = 2000 
    end
    
    properties (Dependent)
        g_t_1
        g_t_2
    end
    
    methods
        function obj = HebbianLearningNeuron(hebbianConsts, max_weight)
            obj.HEBBIAN = hebbianConsts;
            obj.MAX_WEIGHT = max_weight;
        end
        
        function g_t_1 = get.g_t_1(obj)
            g_t_1 = max(0, obj.integralPosVolt - obj.HEBBIAN.NMDA);
        end
        
        function g_t_2 = get.g_t_2(obj)
            g_t_2 = max(0, obj.HEBBIAN.NMDA - obj.integralPosVolt - obj.HEBBIAN.AMPA);
        end
    end
    
    methods     
        function w = calcHebbianWeights(obj, config, scalar, inputNeuron)
            arguments
                obj
                config (1,1) ModelConfig
                scalar (:,:) {mustBeNumeric}
                inputNeuron (1,1) Neuron
            end
            w = obj.weightsForTrial(config);
            w = w + scalar.*( ...
                obj.HEBBIAN.COEF * inputNeuron.integralPosVolt * obj.g_t_1.*(obj.MAX_WEIGHT - w) ...
                - obj.HEBBIAN.ANTI * inputNeuron.integralPosVolt * obj.g_t_2.*w ...
            );
            w = bound_array(w, 0, obj.MAX_WEIGHT);
        end
    end
    
    methods (Abstract)
        doHebbianLearning(obj, config, scalar, inputNeuron)
        w = weightsForTrial(obj, config)
    end
end

