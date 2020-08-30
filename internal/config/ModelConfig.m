classdef ModelConfig
    properties (Abstract, Constant)
        isFROSTEnabled
        isCOVISEnabled
        isMCLearningEnabled
        hasCriterialNoise
    end
    
    properties (Constant)
        STIMULUS_GRID_SIZE = 100;      % Length of side of square grid used for visual input; should be an even number
        BORDER_SIZE        = 20;       % Width of border used to pad the grid such that visual stimulus on the edge still has an appropriate effect
        GRID_SIZE          = 100+2*20; % Total length of grid, i.e., the stimulus grid size and the border
    end
    
    properties
        trials (1,1) {mustBeNumeric} % Number of trials for the configuration.3
        trial (1,1) {mustBeNonzero} = 1
        accuracy (:,1) {mustBeNumeric} % Boolean array recording if the reacting neuron is correct or not.
        
        visual (1,1) VisualStimulus % Struct containing information on visual stimulus and rules.
        RBF (1,1) RadialBasisFunction
        COVISRules (1,1) COVISRuleSet % Struct for information related to COVIS.
        meta % Struct for metadata related to a specific configuration.
    end
    
    properties (Dependent)
        weightIdx
        correctRuleIdx
        chosenRuleIdx
    end
    
    methods
        function c = ModelConfig()
            c.visual = VisualStimulus();
            c.RBF = RadialBasisFunction(ModelConfig.GRID_SIZE, VisualStimulus.STIM);
            c.COVISRules = COVISRuleSet();
        end
    end
    
    methods (Abstract)
        [x_coords, y_coords, coord_groups] = loadCoords(obj)
        config = doPreprocessing(obj)
    end
    
    methods
        function obj = setTrials(obj, trials)
            obj.trials = trials;
            obj.accuracy = zeros(trials,1);
        end
        
        function [config, RULE] = chooseCOVISRule(obj)
            if obj.trial <= obj.COVISRules.GUESSES
                config.COVISRules.chosen = randi(obj.COVISRules.NUM);
            elseif obj.accuracy(obj.trial-1) == 1 || obj.trial > PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS
                config.COVISRules.chosen = obj.COVISRules.log(obj.trial-1);
            else
                config.COVISRules.chosen = rand_discrete(obj.COVISRules.prob);
            end
            RULE = VISUAL.RULES(config.COVISRules.chosen);
            config.COVISRules.log(config.trial) = config.COVISRules.chosen;
        end
        
        function tf = shouldButtonSwitch(obj)
            tf = isa(obj, 'ModelConfigButtonSwitch') && ...
                 obj.trial == (obj.trials - obj.meta.trialsAfterSwitch + 1);
        end
        
        function obj = doButtonSwitch(obj, pmcA, pmcB)
            % Switch the "correct" rule to its inverse.
            obj.COVISRules.correct = obj.visual.RULES(obj.COVISRules.correct).INVERSE;
            
            % Record the PMC weights at this instant.
            obj.meta.PMC_A_weights(:,:,1,:) = pmcA.weights(:,:,obj.weightIdx,:);
            obj.meta.PMC_B_weights(:,:,1,:) = pmcB.weights(:,:,obj.weightIdx,:);
        end
        
        function idx = get.weightIdx(obj)
            if obj.trials > HebbianLearningNeuron.LARGE_TRIAL_BOUNDARY
                idx = 1;
            else
                idx = obj.trial;
            end
        end
        
        function idx = get.correctRuleIdx(obj)
            if obj.isCOVISEnabled
                idx = obj.COVISRules.correct;
            else
                idx = 1;
            end
        end
        
        function idx = get.chosenRuleIdx(obj)
            if obj.isCOVISEnabled
                idx = obj.COVISRules.chosen;
            else
                idx = 1;
            end
        end
        
    end
    
    methods
        function obj = initCOVISRules(obj, params, chosen_rule, correct_rule, trials)
            obj.COVISRules = COVISRuleSet(params, 4, chosen_rule, correct_rule, trials);
        end
    end
end