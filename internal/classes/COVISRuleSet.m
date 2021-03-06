classdef COVISRuleSet
    %COVISRULESET Contains stored data and other information regarding
    %COVIS rulesets and neuron-specific information per rule.
    
    properties (SetAccess = immutable)
        DELTA_C (1,1) {mustBeNumeric} = 1
        DELTA_E (1,1) {mustBeNumeric} = 1
        PERSEV (1,1) {mustBeNumeric} = 1
        LAMBDA (1,1) {mustBeNumeric} = 1
        GUESSES (1,1) {mustBeNumeric} = 5
        
        NUM (1,1) {mustBeNumeric}
    end
    
    properties
        chosen (1,1) {mustBeNumeric} % The chosen rule from the rule set.
        correct (1,1) {mustBeNumeric} % The correct rule from the rule set.
        
        saliences (:,1)
        weights (:,1)
        log (:,1)
    end
    
    properties (Dependent)
        idxs
        prob
    end
    
    methods        
        function obj = COVISRuleSet(params, num_rules, chosen_rule, correct_rule, trials)
            if nargin > 0
                % Set number of rules.
                obj.NUM = num_rules;

                % Set chosen and correct rules.
                obj.chosen = chosen_rule;
                obj.correct = correct_rule;

                % Set logged/saved information per rule.
                obj.saliences = ones(num_rules,1);
                obj.weights = ones(num_rules,1);
                obj.log = ones(trials,1);
            end
        end
        
        function idxs = get.idxs(obj)
            idxs = 1:obj.NUM;
        end
        
        function prob = get.prob(obj)
            prob = obj.weights./sum(obj.weights);
        end
        
        function obj = processRuleAttempt(obj, is_accurate)
            % Step 1: Re-adjust saliences.
            if is_accurate
                obj.saliences(obj.chosen) = obj.saliences(obj.chosen) + obj.DELTA_C;
            else
                obj.saliences(obj.chosen) = max(obj.saliences(obj.chosen) - obj.DELTA_E,1);
            end
            
            % Step 2: Produce weight according to salience.
            obj.weights(obj.chosen) = obj.saliences(obj.chosen) + obj.PERSEV;
            
            % Step 3: Update randomly chosen rule.
            random_rule = randi(obj.NUM);
            obj.weights(random_rule) = obj.weights(random_rule) + poissrnd(obj.LAMBDA);
            
            % Step 4: Carry over all other saliences into weights.
            other_rules = ~ismember(obj.idxs, [obj.chosen, random_rule]);
            obj.weights(other_rules) = obj.saliences(other_rules);
        end
    end
end

