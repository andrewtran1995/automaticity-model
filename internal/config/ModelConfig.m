classdef ModelConfig
    properties (Abstract, Constant)
        name
        isFROSTEnabled
        isCOVISEnabled
        isMCLearningEnabled
        hasCriterialNoise
        hasStroopInterference
    end
    
    properties (Constant)
        STIMULUS_GRID_SIZE = 100;      % Length of side of square grid used for visual input; should be an even number
        BORDER_SIZE        = 20;       % Width of border used to pad the grid such that visual stimulus on the edge still has an appropriate effect
        GRID_SIZE          = 100+2*20; % Total length of grid, i.e., the stimulus grid size and the border
    end
    
    properties
        trials (1,1) {mustBeNumeric} % Number of trials for the configuration.3
        preLearningTrials (1,1) {mustBeNonnegative}
        learningTrials (1,1) {mustBeNonnegative}
        postLearningTrials (1,1) {mustBeNonnegative}
        trial (1,1) {mustBeNonzero} = 1
        accuracy (:,1) {mustBeNumeric} % Boolean array recording if the reacting neuron is correct or not.
        
        visual (1,1) VisualStimulus % Struct containing information on visual stimulus and rules.
        RBF (1,1) RadialBasisFunction
        COVISRules (1,1) COVISRuleSet % Struct for information related to COVIS.
        meta % Struct for metadata related to a specific configuration.
    end
    
    properties (Dependent)
        isLearningTrial
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
    
    methods (Static)
        function m = byName(config_name)
            switch config_name
                case ModelConfigElectro.name()
                    m = ModelConfigElectro();
                case ModelConfigDualTask.name()
                    m = ModelConfigDualTask();
                case ModelConfigButtonSwitch.name()
                    m = ModelConfigButtonSwitch();
            end
        end
    end
    
    methods (Abstract)
        [x_coords, y_coords, coord_groups] = loadCoords(obj)
        config = doPreprocessing(obj)
        dispResults(obj)
    end
    
    methods
        function config = setTrials(config, preLearningTrials, learningTrials, postLearningTrials)
            config.preLearningTrials = preLearningTrials;
            config.learningTrials = learningTrials;
            config.postLearningTrials = postLearningTrials;
            config.trials = preLearningTrials + learningTrials + postLearningTrials;
            config.accuracy = zeros(config.trials,1);
        end
        
        function obj = initCOVISRules(obj, params, chosen_rule, correct_rule, trials)
            obj.COVISRules = COVISRuleSet(params, 2, chosen_rule, correct_rule, trials);
        end
        
        function [config, RULE] = chooseCOVISRule(config)
            if config.trial <= config.COVISRules.GUESSES
                config.COVISRules.chosen = randi(config.COVISRules.NUM);
            elseif config.accuracy(config.trial-1) == 1 || config.trial > config.preLearningTrials + config.learningTrials + config.postLearningTrials
                config.COVISRules.chosen = config.COVISRules.log(config.trial-1);
            else
                config.COVISRules.chosen = rand_discrete(config.COVISRules.prob);
            end
            RULE = config.visual.RULES(config.COVISRules.chosen);
            config.COVISRules.log(config.trial) = config.COVISRules.chosen;
        end
        
        function dispCOVISLog(config)
            figure;
            
            for i=1:config.COVISRules.NUM
                plot(smooth(config.COVISRules.log == i, 500));
                hold on;
            end
            colormap(lines);
            xlim([0, config.trials]);
            legend(arrayfun(@(i) sprintf('Rule %d', i), config.COVISRules.idxs, 'UniformOutput', false));
            title('COVIS Rule Log Frequency');
            
            if isa(config, 'ModelConfigButtonSwitch')
                config.dispButtonSwitchLine();
            end
        end
        
        function dispAccuracy(config)
            figure;
            plot(smooth(config.accuracy, 200), 'b');
            xlim([0, config.trials]);
            ylim([0,1]);
            if isa(config, 'ModelConfigButtonSwitch')
                config.dispButtonSwitchLine();
            end
            title('Accuracy');
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
        
        function tf = get.isLearningTrial(obj)
            tf = obj.trial <= obj.trials && ...
                obj.preLearningTrials < obj.trial && ...
                obj.trial <= obj.preLearningTrials + obj.learningTrials;
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
end