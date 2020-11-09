classdef (Abstract) ModelConfig
    % Ideally, these would be abstract properties. However, it is
    % necessary to make ModelConfig non-abstract for code generation.
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
    
    methods (Abstract)
        [x_coords, y_coords, coord_groups] = loadCoords(obj)
        config = doPreprocessing(obj)
    end

    methods (Static)
        function vals = hebbianValues()
            vals = Hebbian();
        end
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
        
        function dispResults(config, neurons)
            %% Unpack neurons from struct.
            PFC = neurons.PFC;
            PFC_A = neurons.PFC_A;
            PFC_B = neurons.PFC_B;
            PFC_Stroop_A = neurons.PFC_Stroop_A;
            PFC_Stroop_B = neurons.PFC_Stroop_B;
            PMC = neurons.PMC;
            PMC_A = neurons.PMC_A;
            PMC_B = neurons.PMC_B;
            MC = neurons.MC;
            MC_A = neurons.MC_A;
            MC_B = neurons.MC_B;
            Driv_PFC = neurons.Driv_PFC;
            CN = neurons.CN;
            GP = neurons.GP;
            MDN_A = neurons.MDN_A;
            MDN_B = neurons.MDN_B;
            AC_A = neurons.AC_A;
            AC_B = neurons.AC_B;
            %% Neuron information from last trial or throughout trials
            figure;
            rows = 3; columns = 4;

            subplot(rows,columns,1); PFC_A.dispVoltage('PFC_A');
            subplot(rows,columns,2); PFC_B.dispVoltage('PFC_B');
            subplot(rows,columns,3); PMC_A.dispVoltage('PMC_A');
            subplot(rows,columns,4); PMC_B.dispVoltage('PMC_B');
            subplot(rows,columns,5); PFC_A.dispOutput('PFC_A');
            subplot(rows,columns,6); PFC_B.dispOutput('PFC_B');
            subplot(rows,columns,7); PMC_A.dispOutput('PMC_A');
            subplot(rows,columns,8); PMC_B.dispOutput('PMC_B');
            
            subplot(rows,columns,9);
            config.RBF.dispStimulus(config.visual.coord);

            subplot(rows,columns,10);
            x_axis = linspace(1, config.trials, config.trials);
            plot(x_axis, PMC_A.weights_avg, 'r', x_axis, PMC_B.weights_avg, 'b');
            legend('PMC_A', 'PMC_B', 'Location', 'southeast');
            title('PMC_A & PMC_B Weight Average');

            subplot(rows,columns,11);
            linspace(1, config.trials, config.trials);
            PMC_A_Rx = PMC.reactions(:,1) == 1;
            PMC_B_Rx = ~PMC_A_Rx;
            scatter(find(PMC_A_Rx), PMC.reactions(PMC_A_Rx,2), 10, 'r', 'filled'); hold on;
            scatter(find(PMC_B_Rx), PMC.reactions(PMC_B_Rx,2), 10, 'b', 'filled');
            xlim([0, config.trials]);
            legend('PMC_A', 'PMC_B');
            suplabel('PMC_A & PMC_B Reaction Time', 't');
            
            if config.isFROSTEnabled
                figure;
                rows = 7; columns = 2;

                subplot(rows,columns,1); Driv_PFC.dispVoltage('Driv_PFC');
                subplot(rows,columns,3); CN.dispVoltage('CN');
                subplot(rows,columns,5); GP.dispVoltage('GP');
                subplot(rows,columns,7); MDN_A.dispVoltage('MDN_A');
                subplot(rows,columns,9); MDN_B.dispVoltage('MDN_B');
                subplot(rows,columns,11); AC_A.dispVoltage('AC_A');
                subplot(rows,columns,13); AC_B.dispVoltage('AC_B');

                subplot(rows,columns,2); Driv_PFC.dispOutput('Driv_PFC');
                subplot(rows,columns,4); CN.dispOutput('CN');
                subplot(rows,columns,6); GP.dispOutput('GP');
                subplot(rows,columns,8); MDN_A.dispOutput('MDN_A');
                subplot(rows,columns,10); MDN_B.dispOutput('MDN_B');
                subplot(rows,columns,12); AC_A.dispOutput('AC_A');
                subplot(rows,columns,14); AC_B.dispOutput('AC_B');

                suplabel('Neuron Information from Last Trial, Rx Times, Etc.', 't');
            end
            
            %% COVIS Figures
            if config.isCOVISEnabled
                config.dispCOVISLog();
            end
            
            % If not BUTTONS_SWITCH, assume there is a record for the weights for each trial
            if not(isa(config, 'ModelConfigButtonSwitch'))
                figure;
                rows = 1; columns = 2;
                PMC_A.dispWeightsWithSlider('PMC_A', rows, columns, 1, config);
                PMC_B.dispWeightsWithSlider('PMC_B', rows, columns, 2, config);
            % Create figures for button switch
            elseif isa(config, 'ModelConfigButtonSwitch')
                rows = 2; columns = 4;
                figure;
                for i=1:4
                    subplot(rows,columns,i);
                    colormap('hot');
                    imagesc(config.meta.PMC_A_weights(config.BORDER_SIZE:end-config.BORDER_SIZE, config.BORDER_SIZE:end-config.BORDER_SIZE, 1, i));
                    title(sprintf('PMC_A Rule %d', i));
                end
                for i=1:4
                    subplot(rows,columns,i + 4);
                    colormap('hot');
                    imagesc(config.meta.PMC_B_weights(config.BORDER_SIZE:end-config.BORDER_SIZE, config.BORDER_SIZE:end-config.BORDER_SIZE, 1, i));
                    title(sprintf('PMC_B Rule %d', i));
                end
                suplabel('Initial Heatmaps (Upon Button Switch)', 't');
                figure;
                subplot(rows,columns,1);
                for i=1:4
                    subplot(rows,columns,i);
                    colormap('hot');
                    imagesc(PMC_A.weights(config.BORDER_SIZE:end-config.BORDER_SIZE, config.BORDER_SIZE:end-config.BORDER_SIZE, 1, i));
                    title(sprintf('PMC_A Rule %d', i));
                end
                for i=1:4
                    subplot(rows,columns,i + 4);
                    colormap('hot');
                    imagesc(PMC_B.weights(config.BORDER_SIZE:end-config.BORDER_SIZE, config.BORDER_SIZE:end-config.BORDER_SIZE, 1, i));
                    title(sprintf('PMC_B Rule %d', i));
                end
                suplabel('Final Heatmaps', 't');
            end
            
            %% MC Weights
            figure;
            rows = 2; columns = 1;
            subplot(rows,columns,1);
            MC_A.dispWeights(config, 'MC_A');
            
            subplot(rows,columns,2);
            MC_B.dispWeights(config, 'MC_B');
            suplabel('MC Weights', 't');
            
            %% Reaction Latency
            % Histograms of reaction latencies by neuron and trial subsets
            figure;
            LEARNING_IDX = (config.preLearningTrials+1):(config.preLearningTrials+config.learningTrials);
            latencies = {PFC.reactions(1:config.preLearningTrials,2), 'PFC Latencies (Pre-Learning)'; ...
                         PMC.reactions(1:config.preLearningTrials,2), 'PMC Latencies (Pre-Learning)'; ...
                         PFC.reactions(LEARNING_IDX,2), 'PFC Latencies (Learning)'; ...
                         PMC.reactions(LEARNING_IDX,2), 'PMC Latencies (Learning)'; ...
                         PFC.reactions(end-config.postLearningTrials+1:end, 2), 'PFC Latencies (Post-Learning)'; ...
                         PMC.reactions(end-config.postLearningTrials+1:end, 2), 'PMC Latencies (Post-Learning)'};
            latencies(cellfun(@isempty, latencies(:,1)), :) = [];
            rows = length(latencies(:,1))/2; columns = 2;
            numBins = 20;
            for i=1:rows*columns
                subplot(rows,columns,i);
                histogram(latencies{i,1}, numBins);
                xlabel('Latency of selectivity for the behavioral response (ms)');
                ylabel('Number of neurons');
                title(latencies{i,2});
            end
            suplabel('Reaction Latency Histograms', 't');

            %% Reaction Plots
            figure;
            rows = 3; columns = 1;
            latencies = {PFC.reactions(:,2), PMC.reactions(:,2), MC.reactions(:,2)};
            latencyTitles = {'PFC Reaction Time', 'PMC Reaction Time', 'MC Reaction Time'};
            for i=1:3
                subplot(rows,columns,i);
                plot(smooth(latencies{i},30));
                xlim([0, config.trials]);
                if isa(config, 'ModelConfigButtonSwitch')
                    config.dispButtonSwitchLine();
                end
                title(latencyTitles{i});
            end
            suplabel('Reaction Latencies Over Time', 't');
            
            %% Accuracy
            config.dispAccuracy()
        end
    end
end