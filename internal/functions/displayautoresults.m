function displayautoresults(config, RBF, PFC, PMC, MC, PFC_A, PFC_B, PMC_A, PMC_B, MC_A, MC_B, Driv_PFC, CN, GP, MDN_A, MDN_B, AC_A, AC_B )
%DISPLAYAUTORESULTS Display results an Automaticity Model run
%   Display results from an Automaticity Model run. Requires *all* (relevant)
%   variables from the Automaticity Model workspace to be passed in.
%   Separated for code clarity and ease of code-generation.
    BORDER_SIZE = config.BORDER_SIZE;
    TRIALS = config.trials;
    PRE_LEARNING_TRIALS = config.preLearningTrials;
    LEARNING_TRIALS = config.learningTrials;
    POST_LEARNING_TRIALS = config.postLearningTrials;
    LEARNING_IDX = (PRE_LEARNING_TRIALS+1):(PRE_LEARNING_TRIALS+LEARNING_TRIALS);

    config.dispResults(RBF, PFC, PMC, MC, PFC_A, PFC_B, PMC_A, PMC_B, MC_A, MC_B, Driv_PFC, CN, GP, MDN_A, MDN_B, AC_A, AC_B);

    %% Figure 2 - Synaptic Weight Heatmaps
    % Only relevant if any learning trials were conducted
    if LEARNING_TRIALS > 0
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
                imagesc(config.meta.PMC_A_weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, i));
                title(sprintf('PMC_A Rule %d', i));
            end
            for i=1:4
                subplot(rows,columns,i + 4);
                colormap('hot');
                imagesc(config.meta.PMC_B_weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, i));
                title(sprintf('PMC_B Rule %d', i));
            end
            suplabel('Initial Heatmaps (Upon Button Switch)', 't');
            figure;
            subplot(rows,columns,1);
            for i=1:4
                subplot(rows,columns,i);
                colormap('hot');
                imagesc(PMC_A.weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, i));
                title(sprintf('PMC_A Rule %d', i));
            end
            for i=1:4
                subplot(rows,columns,i + 4);
                colormap('hot');
                imagesc(PMC_B.weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, i));
                title(sprintf('PMC_B Rule %d', i));
            end
            suplabel('Final Heatmaps', 't');
        end
        % MC Weights
        figure;
        rows = 2; columns = 1;
        subplot(rows,columns,1);
        plot(smooth(MC_A.weights(1,:),50), 'r'); hold on;
        plot(smooth(MC_A.weights(2,:),50), 'b');
        if isa(config, 'ModelConfigButtonSwitch')
            config.dispButtonSwitchLine();
        end
        xlim([0, TRIALS]);
        legend('Weight to PMC_A', 'Weight to PMC_B');
        title('MC_A');
        subplot(rows,columns,2);
        plot(smooth(MC_B.weights(1,:),50), 'r'); hold on;
        plot(smooth(MC_B.weights(2,:),50), 'b');
        if isa(config, 'ModelConfigButtonSwitch')
            config.dispButtonSwitchLine();
        end
        xlim([0, TRIALS]);
        legend('Weight to PMC_B', 'Weight to PMC_A');
        title('MC_B');
        suplabel('MC Weights', 't');
    end

    %% Figure 5 - Reaction Latency
    % Histograms of reaction latencies by neuron and trial subsets
    figure;
    latencies = {PFC.reactions(1:PRE_LEARNING_TRIALS,2), 'PFC Latencies (Pre-Learning)'; ...
                 PMC.reactions(1:PRE_LEARNING_TRIALS,2), 'PMC Latencies (Pre-Learning)'; ...
                 PFC.reactions(LEARNING_IDX,2), 'PFC Latencies (Learning)'; ...
                 PMC.reactions(LEARNING_IDX,2), 'PMC Latencies (Learning)'; ...
                 PFC.reactions(end-POST_LEARNING_TRIALS+1:end, 2), 'PFC Latencies (Post-Learning)'; ...
                 PMC.reactions(end-POST_LEARNING_TRIALS+1:end, 2), 'PMC Latencies (Post-Learning)'};
    latencies(cellfun(@isempty, latencies(:,1)), :) = [];
    rows = length(latencies(:,1))/2; columns = 2;
    numBins = 20;
    for i=1:rows*columns
        subplot(rows,columns,i);
        hist(latencies{i,1}, numBins);
        xlabel('Latency of selectivity for the behavioral response (ms)');
        ylabel('Number of neurons');
        title(latencies{i,2});
    end
    suplabel('Reaction Latency Histograms', 't');
    
    %% Figure 6 - Reaction Plots
    figure;
    rows = 3; columns = 1;
    latencies = {PFC.reactions(:,2), PMC.reactions(:,2), MC.reactions(:,2)};
    latencyTitles = {'PFC Reaction Time', 'PMC Reaction Time', 'MC Reaction Time'};
    for i=1:3
        subplot(rows,columns,i);
        plot(smooth(latencies{i},30));
        xlim([0, TRIALS]);
        if isa(config, 'ModelConfigButtonSwitch')
            config.dispButtonSwitchLine();
        end
        title(latencyTitles{i});
    end
    suplabel('Reaction Latencies Over Time', 't');
    
    %% Figure 7 - Accuracy
    config.dispAccuracy()

    %% Starts debug mode, allowing variables to be observed before the function ends
    keyboard;

end

% Find the hazard function as defined by Hazard = f(t)/S(t),
% where f(t) is the PDF and S(t) is the survivor function
function [f] = get_hazard_estimate(x, pts)
    [f_pdf, ~] = ksdensity(x, pts, 'function', 'pdf');
    [f_sur, ~] = ksdensity(x, pts, 'function', 'survivor');
    f = f_pdf./f_sur;
end