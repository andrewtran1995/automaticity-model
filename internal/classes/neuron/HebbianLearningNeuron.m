classdef (Abstract) HebbianLearningNeuron < RSN
    %HEBBIANLEARNINGNEURON A neuron that engages in learning per Hebbian
    %theory.
    
    properties
        weights
    end
    
    properties (SetAccess = immutable)
        HEBBIAN (1,1) HebbianConst
        MAX_WEIGHT (1,1) {mustBeNonnegative}
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
            g_t_2 = max(0, obj.HEBBIAN.NMDA - obj.integralPosVolt);
        end
    end

    methods     
        % "inputNeuron" should be the integral of the alpha function of the pre-synaptic unit.
        % This can be simplified to the integral of the positive voltage.
        % For the visual units, since they're not spiking, they have a
        % constant output which is the same for all Neuron.n time points
        % which is essentially the radial basis function vector multiplied
        % by the weights of the PMC.
        % Multiply radial basis values by weights, then 1000.
        function w = calcHebbianWeights(obj, config, scalar, preSynapticOutput)
            w = obj.weightsForTrial(config);
            w = w + scalar.*( ...
                obj.HEBBIAN.COEF * preSynapticOutput * obj.g_t_1.*(obj.MAX_WEIGHT - w) ...
                - obj.HEBBIAN.COEF * preSynapticOutput * obj.g_t_2.*w ...
            );
            w = bound_array(w, 0, obj.MAX_WEIGHT);
        end
    end
    
    methods (Abstract)
        doHebbianLearning(obj, config, scalar, inputNeuron)
        w = weightsForTrial(obj, config)
    end
end

