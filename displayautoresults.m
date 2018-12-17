function displayautoresults( FROST_ENABLED, COVIS_ENABLED, BUTTON_SWITCH_ENABLED, BUTTON_SWITCH, COVIS_VARS, FMRI_META, configuration, TAU, n, RBF, BORDER_SIZE, VISUAL, TRIALS, PRE_LEARNING_TRIALS, LEARNING_TRIALS, POST_LEARNING_TRIALS, accuracy, PFC, PMC, MC, PFC_A, PFC_B, PMC_A, PMC_B, MC_A, MC_B, Driv_PFC, CN, GP, MDN_A, MDN_B, AC_A, AC_B, PERF_OUTPUT, start_time, loop_times, trial_times, rt_calc_times, chosen_rule )
%DISPLAYAUTORESULTS Display results an Automaticity Model run
%   Display results from an Automaticity Model run. Requires *all* (relevant)
%   variables from the Automaticity Model workspace to be passed in.
%   Separated for code clarity and ease of code-generation.
    addpath('classes','libraries');
    LEARNING_IDX = (PRE_LEARNING_TRIALS+1):(PRE_LEARNING_TRIALS+LEARNING_TRIALS);

    %% Starts debug mode if you want to view the variables before executing code
    % keyboard;

    %% Figure 1 - neuron information from last trial or throughout trials
    figure;
    rows = 3; columns = 4;

    subplot(rows,columns,1); plot(TAU*(1:n),PFC_A.v);
    axis([0 n -100 100]); title('PFC_A Neuron Voltage');

    subplot(rows,columns,2); plot(TAU*(1:n),PFC_B.v);
    axis([0 n -100 100]); title('PFC_B Neuron Voltage');

    subplot(rows,columns,3); plot(TAU*(1:n),PMC_A.v);
    axis([0 n -100 100]); title('PMC_A Neuron Voltage');

    subplot(rows,columns,4); plot(TAU*(1:n),PMC_B.v);
    axis([0 n -100 100]); title('PMC_B Neuron Voltage');

    subplot(rows,columns,5); plot(TAU*(1:n),PFC_A.out);
    axis([0 n -1 10]); title('PFC_A Neuron Output');

    subplot(rows,columns,6); plot(TAU*(1:n),PFC_B.out);
    axis([0 n -1 10]); title('PFC_B Neuron Output');

    subplot(rows,columns,7); plot(TAU*(1:n),PMC_A.out);
    axis([0 n -1 10]); title('PMC_A Neuron Output');

    subplot(rows,columns,8); plot(TAU*(1:n),PMC_B.out);
    axis([0 n -1 10]); title('PMC_B Neuron Output');

    subplot(rows,columns,9);
    colormap('hot');
    imagesc(RBF.rbv(BORDER_SIZE:end-BORDER_SIZE-1,BORDER_SIZE:end-BORDER_SIZE-1,:));
    title(sprintf('Stimulus: (%d,%d)', VISUAL.y_coord, VISUAL.x_coord));

    subplot(rows,columns,10);
    x_axis = linspace(1, TRIALS, TRIALS);
    plot(x_axis, PMC_A.weights_avg, 'r', x_axis, PMC_B.weights_avg, 'b');
    legend('PMC_A', 'PMC_B', 'Location', 'southeast');
    title('PMC_A & PMC_B Weight Average');

    subplot(rows,columns,11);
    x_axis = linspace(1, TRIALS, TRIALS);
    PMC_A_Rx = PMC.reactions(:,1) == 1;
    PMC_B_Rx = ~PMC_A_Rx;
    scatter(find(PMC_A_Rx), PMC.reactions(PMC_A_Rx,2), 10, 'r', 'filled'); hold on;
    scatter(find(PMC_B_Rx), PMC.reactions(PMC_B_Rx,2), 10, 'b', 'filled');
    xlim([0, TRIALS]);
    legend('PMC_A', 'PMC_B');
    suplabel('PMC_A & PMC_B Reaction Time', 't');

    %% Figure 1B - neuron information from last trial or throughout trials
    if FROST_ENABLED
        figure;
        rows = 7; columns = 2;

        subplot(rows,columns,1); plot(TAU*(1:n),Driv_PFC.v);
        axis([0 n -100 100]); title('Driv_PFC Voltage');

        subplot(rows,columns,3); plot(TAU*(1:n),CN.v);
        axis([0 n -100 100]); title('CN Voltage');

        subplot(rows,columns,5); plot(TAU*(1:n),GP.v);
        axis([0 n -100 100]); title('GP Voltage');

        subplot(rows,columns,7); plot(TAU*(1:n),MDN_A.v);
        axis([0 n -100 100]); title('MDN_A Voltage');

        subplot(rows,columns,9); plot(TAU*(1:n),MDN_B.v);
        axis([0 n -100 100]); title('MDN_B Voltage');

        subplot(rows,columns,11); plot(TAU*(1:n),AC_A.v);
        axis([0 n -100 100]); title('AC_A Voltage');

        subplot(rows,columns,13); plot(TAU*(1:n),AC_B.v);
        axis([0 n -100 100]); title('AC_B Voltage');

        subplot(rows,columns,2); plot(TAU*(1:n),Driv_PFC.out);
        axis([0 n 0 30]); title('Driv PFC Output');

        subplot(rows,columns,4); plot(TAU*(1:n),CN.out);
        axis([0 n 0 30]); title('CN Output');

        subplot(rows,columns,6); plot(TAU*(1:n),GP.out);
        axis([0 n 0 30]); title('GP Output');

        subplot(rows,columns,8); plot(TAU*(1:n),MDN_A.out);
        axis([0 n 0 30]); title('MDN_A Output');

        subplot(rows,columns,10); plot(TAU*(1:n),MDN_B.out);
        axis([0 n 0 30]); title('MDN_B Output');

        subplot(rows,columns,12); plot(TAU*(1:n),AC_A.out);
        axis([0 n 0 30]); title('AC_A Output');

        subplot(rows,columns,14); plot(TAU*(1:n),AC_B.out);
        axis([0 n 0 30]); title('AC_B Output');

        suplabel('Neuron Information from Last Trial, Rx Times, Etc.', 't');
    end

    %% COVIS Figures
    if COVIS_ENABLED
        figure;
        rule_legend = { 'r', 'Rule 1'; ...
                        'b', 'Rule 2'; ...
                        'g', 'Rule 3'; ...
                        'm', 'Rule 4' ...
        };
        rules = { COVIS_VARS.rule_log == 1; ...
                  COVIS_VARS.rule_log == 2; ...
                  COVIS_VARS.rule_log == 3; ...
                  COVIS_VARS.rule_log == 4 ...
        };
        for i=1:4
            plot(smooth(rules{i}, 500), rule_legend{i,1}); hold on;
        end
        xlim([0, TRIALS]);
        legend(rule_legend{:,2});
        title('COVIS Rule Log Frequency');
    end

    %% Figure 2 - Synaptic Weight Heatmaps
    % Only relevant if any learning trials were conducted
    if LEARNING_TRIALS > 0
        % If not FMRI, assume there is a record for the weights for each trial
        if configuration ~= AutomaticityConfiguration.FMRI
            figure;
            rows = 1; columns = 2;
            % Force slider to integer/discrete value:
            % https://www.mathworks.com/matlabcentral/answers/45769-forcing-slider-values-to-round-to-a-valid-number

            PMC_A_trial_num = 1;
            PMC_A_no_border = PMC_A.weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, LEARNING_IDX);
            subplot(rows,columns,1);
            data3 = PMC_A_no_border(:,:,PMC_A_trial_num);
            colormap('hot');
            imagesc(data3);
            colorbar;
            title(sprintf('PMC_A Synaptic Heatmap, Trial %d\n', PMC_A_trial_num));
            slider_PMC_A = uicontrol('Style', 'slider', ...
                                     'Min', 1, 'Max', LEARNING_TRIALS, ...
                                     'Value', 1, ...
                                     'Position', [100 50 300 20]);
            set(slider_PMC_A, 'Callback', {@synaptic_slider_callback, 1, PMC_A_no_border, 'PMC_A'});

            PMC_B_trial_num = 1;
            PMC_B_no_border = PMC_B.weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, LEARNING_IDX);
            subplot(rows,columns,2);
            data4 = PMC_B_no_border(:,:,PMC_B_trial_num);
            colormap('hot');
            imagesc(data4);
            colorbar;
            title(sprintf('PMC_B Synaptic Heatmap, Trial %d\n', PMC_B_trial_num));
            slider_PMC_B = uicontrol('Style', 'slider', ...
                                     'Min', 1, 'Max', LEARNING_TRIALS, ...
                                     'Value', 1, ...
                                     'Position', [500 50 300 20]);
            set(slider_PMC_B, 'Callback', {@synaptic_slider_callback, 2, PMC_B_no_border, 'PMC_B'});
            
            suplabel('Synaptic Heatmaps', 't');
        % Create figures for button switch
        elseif configuration == AutomaticityConfiguration.FMRI && BUTTON_SWITCH_ENABLED
            rows = 2; columns = 4;
            figure;
            for i=1:4
                subplot(rows,columns,i);
                colormap('hot');
                imagesc(BUTTON_SWITCH.PMC_A_weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, i));
                title(sprintf('PMC_A Rule %d', i));
            end
            for i=1:4
                subplot(rows,columns,i + 4);
                colormap('hot');
                imagesc(BUTTON_SWITCH.PMC_B_weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, i));
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
        if BUTTON_SWITCH_ENABLED
            hold on;
            plot([TRIALS - BUTTON_SWITCH.TRIALS; TRIALS - BUTTON_SWITCH.TRIALS], get(gca,'ylim'), 'r');
        end
        xlim([0, TRIALS]);
        legend('Weight to PMC_A', 'Weight to PMC_B');
        title('MC_A');
        subplot(rows,columns,2);
        plot(smooth(MC_B.weights(1,:),50), 'r'); hold on;
        plot(smooth(MC_B.weights(2,:),50), 'b');
        if BUTTON_SWITCH_ENABLED
            hold on;
            plot([TRIALS - BUTTON_SWITCH.TRIALS; TRIALS - BUTTON_SWITCH.TRIALS], get(gca,'ylim'), 'r');
        end
        xlim([0, TRIALS]);
        legend('Weight to PMC_B', 'Weight to PMC_A');
        title('MC_B');
        suplabel('MC Weights', 't');
    end

    if configuration == AutomaticityConfiguration.MADDOX
        %% Figure 3
        % CDFs of RTs (reaction times) dependent on stimulus type - Short, Medium, or Long
        % CDF = P(RT <= t), for each specific value t
        PMC_S = PMC.reactions(LEARNING_IDX,3) == 'S';
        PMC_M = PMC.reactions(LEARNING_IDX,3) == 'M';
        PMC_L = PMC.reactions(LEARNING_IDX,3) == 'L';

        figure;
        p1 = cdfplot(PMC.reactions(PMC_S, 2)); set(p1, 'Color', 'r'); hold on;
        p2 = cdfplot(PMC.reactions(PMC_M, 2)); set(p2, 'Color', 'b'); hold on;
        p3 = cdfplot(PMC.reactions(PMC_L, 2)); set(p3, 'Color', 'g');
        legend('S', 'M', 'L', 'Location', 'southeast');
        suplabel('CDFs of PMC Rx Times (Grouped by Distance)');

        %% Figure 4 - Hazard Functions
        % Hazard Function = f(t)/[1-F(t)], where f(t) = PDF, F(t) = CDF
        % https://www.mathworks.com/help/stats/survival-analysis.html#btnxirj-1
        figure; title('PMC Rx Times Hazard Functions');

        % Reuse vars from CDF plot
        pts = (min(PMC.reactions(LEARNING_IDX, 2)):0.25:max(PMC.reactions(LEARNING_IDX, 2)));
        plot(pts, get_hazard_estimate(PMC.reactions(PMC_S, 2), pts), 'Color', 'r'); hold on;
        plot(pts, get_hazard_estimate(PMC.reactions(PMC_M, 2), pts), 'Color', 'b'); hold on;
        plot(pts, get_hazard_estimate(PMC.reactions(PMC_L, 2), pts), 'Color', 'g');
        legend('S', 'M', 'L', 'Location', 'southeast');
        title('Hazard Functions');
    end

    %% Figure 5 - Reaction Latency
    % Histograms of reaction latencies by neuron and trial subsets
    figure;
    latencies = {PFC.reactions(1:PRE_LEARNING_TRIALS,2), 'PFC Latencies (Pre-Learning)'; ...
                 PMC.reactions(1:PRE_LEARNING_TRIALS,2), 'PMC Latencies (Pre-Learning)'; ...
                 PFC.reactions(LEARNING_IDX,2), 'PFC Latencies (Learning)'; ...
                 PMC.reactions(LEARNING_IDX,2), 'PMC Latencies (Learning)'; ...
                 PFC.reactions(end-POST_LEARNING_TRIALS+1:end, 2), 'PFC Latencies (Post-Learning)'; ...
                 PMC.reactions(end-POST_LEARNING_TRIALS+1:end, 2), 'PMC Latencies (Post-Learning'};
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
        if BUTTON_SWITCH_ENABLED
            hold on;
            plot([TRIALS - BUTTON_SWITCH.TRIALS; TRIALS - BUTTON_SWITCH.TRIALS], get(gca,'ylim'), 'r');
        end
        title(latencyTitles{i});
    end
    suplabel('Reaction Latencies Over Time', 't');
    
    %% Figure 7 - Accuracy
    figure;
    plot(smooth(accuracy, 200), 'b');
    xlim([0, TRIALS]); ylim([0, 1]);
    if BUTTON_SWITCH_ENABLED
       hold on;
       plot([TRIALS - BUTTON_SWITCH.TRIALS; TRIALS - BUTTON_SWITCH.TRIALS], get(gca,'ylim'), 'r');
    end
    title('Accuracy');

    %% Figure 8 - Performance Tests
    % Information regarding the run-time of this program
    if PERF_OUTPUT
        elapsedTime = toc(start_time);
        figure; title(sprintf('TOTAL: %d, MEAN(LOOP): %d', elapsedTime, mean(loop_times))); hold on;
        plot(loop_times, 'b'); hold on;
        plot(trial_times, 'r'); hold on;
        plot(rt_calc_times, 'g');
    end

    %% Starts debug mode, allowing variables to be observed before the function ends
    keyboard;

end

% Handles the slider functionality for the synaptic weight heatmaps
function synaptic_slider_callback(src, ~, position, data, neuron_name)
    subplot(1, 2, position, 'replace');
    trial_num = round(get(src, 'value'));
    set(src, 'value', trial_num);
    colormap('hot');
    imagesc(data(:,:,trial_num));
    colorbar;
    title(sprintf('%s Synaptic Heatmap, Trial %d\n', neuron_name, trial_num'));
end

% Find the hazard function as defined by Hazard = f(t)/S(t),
% where f(t) is the PDF and S(t) is the survivor function
function [f] = get_hazard_estimate(x, pts)
    [f_pdf, ~] = ksdensity(x, pts, 'function', 'pdf');
    [f_sur, ~] = ksdensity(x, pts, 'function', 'survivor');
    f = f_pdf./f_sur;
end