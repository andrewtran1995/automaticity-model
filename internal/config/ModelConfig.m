classdef ModelConfig < handle
    properties (SetAccess = immutable)
        isFROSTEnabled
        isCOVISEnabled
    end
    
    properties (Constant)
        STIMULUS_GRID_SIZE = 100;      % Length of side of square grid used for visual input; should be an even number
        BORDER_SIZE        = 20;       % Width of border used to pad the grid such that visual stimulus on the edge still has an appropriate effect
        GRID_SIZE          = 100+2*20; % Total length of grid, i.e., the stimulus grid size and the border
    end
    
    properties (SetAccess = private)
        trials (1,1) {mustBeNumeric} % Number of trials for the configuration.
        accuracy (:,1) {mustBeNumeric} % Boolean array recording if the reacting neuron is correct or not.
        visual = createvisualstimulusrules() % Struct containing information on visual stimulus and rules.
        RBF    = RadialBasisFunction(ModelConfig.GRID_SIZE, createvisualstimulusrules().STIM)
        COVISRules (1,1) COVISRuleSet % Struct for information related to COVIS.
        meta % Struct for metadata related to a specific configuration.
    end
    
    properties
        trial (1,1) {mustBeNumeric,mustBeNonzero} = 1
    end
    
    properties (Dependent)
        weightIdx
        correctRuleIdx
        chosenRuleIdx
    end
    
    methods
        function c = ModelConfig(isFROSTEnabled, isCOVISEnabled, metadata)
            c.isFROSTEnabled = isFROSTEnabled;
            c.isCOVISEnabled = isCOVISEnabled;
            c.meta = metadata;
        end
    end
    
    enumeration
        %{
        Corresponds to the Helie, S. paper.
        Observes a button switch effect with some late-stage dual-task
        constraints.
        %}
        BUTTON_SWITCH(false, true, struct( ...
            'trialsAfterSwitch', 600, ...
            'PMC_A_weights', ones(ModelConfig.GRID_SIZE,ModelConfig.GRID_SIZE,1,4), ...
            'PMC_B_weights', ones(ModelConfig.GRID_SIZE,ModelConfig.GRID_SIZE, 1,4), ...
            'optimization', struct( ...
                'NUM_TRIALS', 11520, 'GROUP_RUN', 0, ...
                'SES_1',      1:480, 'SES_4',    1681:2160, ...
                'SES_10', 5161:5640, 'SES_20', 11041:11520 ...
            ) ...
        )),
        %{
        Corresponds to the Wallis, Jonathan D. paper.
        Observes image correlation effects.
        %}
        IMAGE_CORR   (false, false, struct()),
        %{
        Corresponds to the Zeithamova, D. paper.
        Observes dual task interference in category learning.
        %}
        DUAL_TASK    (true, true, struct())
    end
    
    methods
        function setTrials(obj, trials)
            obj.trials = trials;
            obj.accuracy = zeros(trials,1);
        end
        
        function tf = shouldButtonSwitch(obj)
            arguments
                obj
            end
            tf = obj == ModelConfig.BUTTON_SWITCH && ...
                 obj.trial == obj.meta.trialsAfterSwitch + 1;
        end
        
        function doButtonSwitch(obj, pmcA, pmcB)
            arguments
                obj
                pmcA (1,1) PMCNeuron
                pmcB (1,1) PMCNeuron
            end
            % Switch the "correct" rule to its inverse.
            ruleProps = obj.visual.RULES(obj.COVIS.rule.correct);
            obj.COVIS.rule.correct = obj.visual.RULES(ruleProps.INVERSE);
            
            % Record the PMC weights at this instant.
            obj.meta.PMC_A_weights(:,:,1,:) = pmcA.weights(:,:,obj.weightIndex(obj.trial),:);
            obj.meta.PMC_B_weights(:,:,1,:) = pmcB.weights(:,:,obj.weightIndex(obj.trial),:);
        end
        
        function idx = get.weightIdx(obj)
            if obj.trials > PMCNeuron.LARGE_TRIAL_BOUNDARY
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
        function initCOVISRules(obj, params, chosen_rule, correct_rule, trials)
            obj.COVISRules = COVISRuleSet(params, 4, chosen_rule, correct_rule, trials);
        end
    end
end