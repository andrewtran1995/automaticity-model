function displayautoresults( FROST_ENABLED, COVIS_ENABLED, FMRI_META, CONFIGURATION, MADDOX, WALLIS, FMRI, TAU, n, RBF, BORDER_SIZE, VISUAL, TRIALS, LEARNING_TRIALS, LEARNING_IDX, PFC_A, PFC_B, PMC_A, PMC_B, Driv_PFC, CN, GP, MDN_A, MDN_B, AC_A, AC_B, PERF_TEST, start_time, loop_times, trial_times, rt_calc_times, chosen_rule )
%DISPLAYAUTORESULTS Display results an Automaticity Model run
%   Display results from an Automaticity Model run. Requires *all*
%   variables from the Automaticity Model workspace to be passed in.
%   Separated for code clarity and ease of code-generation.
    %% Figure 1 - neuron information from last trial or throughout trials
    figure; title('Neuron Information from Last Trial, Rx Times, Etc.');
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
    title(sprintf('Stimulus: (%d,%d); Weight: %d', r_y, r_x, VISUAL.STIM));

    subplot(rows,columns,10);
    x_axis = linspace(1, TRIALS, TRIALS);
    plot(x_axis, PMC_A.weights_avg, 'r', x_axis, PMC_B.weights_avg, 'b');
    legend('PMC_A', 'PMC_B', 'Location', 'southeast');
    title('PMC_A & PMC_B Weight Average');

    subplot(rows,columns,11);
    x_axis = linspace(1, TRIALS, TRIALS);
    PMC_A_Rx = PMC.rx_matrix(:,1) == 1;
    PMC_B_Rx = ~PMC_A_Rx;
    scatter(find(PMC_A_Rx), PMC.rx_matrix(PMC_A_Rx,2), 10, 'r', 'filled');
    hold on;
    scatter(find(PMC_B_Rx), PMC.rx_matrix(PMC_B_Rx,2), 10, 'b', 'filled');
    legend('PMC_A', 'PMC_B');
    title('PMC_A & PMC_B Reaction Time');

    %% Figure 1B - neuron information from last trial or throughout trials
    if FROST_ENABLED
        figure; title('Neuron Information from Last Trial, Rx Times, Etc.');
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
        axis([0 n 0 30]); title('Driv_PFC Output');

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
    end

    %% Figure 1C - COVIS
    if COVIS_ENABLED
        figure; title('COVIS Information');
        scatter(1:max(size(COVIS_VARS.rule_log)), COVIS_VARS.rule_log);
    end

    %% Figure 2
    % Synaptic weight heatmaps with sliders to allow the observation of the heatmap at different intervals in time
    % Only relevant if any learning trials were conducted
    if CONFIGURATION ~= FMRI && LEARNING_TRIALS > 0
        figure; title('Synaptic Heatmaps');
        rows = 1; columns = 2;
        % Force slider to integer/discrete value:
        % https://www.mathworks.com/matlabcentral/answers/45769-forcing-slider-values-to-round-to-a-valid-number
        PMC_A_trial_num = 1;
        PMC_A_no_border = PMC_A.weights(BORDER_SIZE:end-BORDER_SIZE, ...
                                        BORDER_SIZE:end-BORDER_SIZE, ...
                                        LEARNING_IDX);
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
        PMC_B_no_border = PMC_B.weights(BORDER_SIZE:end-BORDER_SIZE, ...
                                        BORDER_SIZE:end-BORDER_SIZE, ...
                                        LEARNING_IDX);
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
    end

    if CONFIGURATION == MADDOX
        %% Figure 3
        % CDFs of RTs (reaction times) dependent on stimulus type -- Short, Medium, or Long
        % CDF = P(RT <= t), for each specific value t
        % Set-up
        PMC_S = PMC.rx_matrix(LEARNING_IDX,3) == 'S';
        PMC_M = PMC.rx_matrix(LEARNING_IDX,3) == 'M';
        PMC_L = PMC.rx_matrix(LEARNING_IDX,3) == 'L';        
        figure; title('CDFs of PMC Rx Times (Grouped by Distance)');

        p1 = cdfplot(PMC.rx_matrix(PMC_S, 2));
        set(p1, 'Color', 'r');
        hold on;
        p2 = cdfplot(PMC.rx_matrix(PMC_M, 2));
        set(p2, 'Color', 'b');
        hold on;
        p3 = cdfplot(PMC.rx_matrix(PMC_L, 2));
        set(p3, 'Color', 'g');
        legend('S', 'M', 'L', 'Location', 'southeast');
        title('CDFs of RTs by Grouping');

        %% Figure 4 - Hazard Functions
        % Hazard Function = f(t)/[1-F(t)], where f(t) = PDF, F(t) = CDF
        % https://www.mathworks.com/help/stats/survival-analysis.html#btnxirj-1
        figure; title('PMC Rx Times Hazard Functions');

        % Reuse vars from CDF plot
        pts = (min(PMC.rx_matrix(LEARNING_IDX, 2)):0.25:max(PMC.rx_matrix(LEARNING_IDX, 2)));
        plot(pts, get_hazard_estimate(PMC.rx_matrix(PMC_S, 2), pts), 'Color', 'r');
        hold on;
        plot(pts, get_hazard_estimate(PMC.rx_matrix(PMC_M, 2), pts), 'Color', 'b');
        hold on;
        plot(pts, get_hazard_estimate(PMC.rx_matrix(PMC_L, 2), pts), 'Color', 'g');
        legend('S', 'M', 'L', 'Location', 'southeast');
        title('Hazard Functions');
    end

    %% Figure 5 - Reaction Latency
    % Compare the latency of the PFC versus the PMC
    % TODO: factor out x/y labeling
    f = figure;
    title('Reaction Latency Histograms');
    rows = 3; columns = 2;
    numBins = 20;

    % Pre-Learning
    subplot(rows,columns,1);
    hist(PFC.rx_matrix(1:PRE_LEARNING_TRIALS,2), numBins);
    xlabel('Latency of selectivity for the behavioral response (ms)');
    ylabel('Number of neurons');
    title('PFC Latencies (Pre-Learning)');
    subplot(rows,columns,2);
    hist(PMC.rx_matrix(1:PRE_LEARNING_TRIALS,2), numBins);
    xlabel('Latency of selectivity for the behavioral response (ms)');
    ylabel('Number of neurons');
    title('PMC Latencies (Pre-Learning)');
    % Learning
    subplot(rows,columns,3);
    hist(PFC.rx_matrix(LEARNING_IDX,2), numBins);
    xlabel('Latency of selectivity for the behavioral response (ms)');
    ylabel('Number of neurons');
    title('PFC Latencies (Learning)');
    subplot(rows,columns,4);
    hist(PMC.rx_matrix(LEARNING_IDX,2), numBins);
    xlabel('Latency of selectivity for the behavioral response (ms)');
    ylabel('Number of neurons');
    title('PMC Latencies (Learning)');
    % Post-Learning
    subplot(rows,columns,5);
    hist(PFC.rx_matrix(end-POST_LEARNING_TRIALS+1:end, 2), numBins);
    xlabel('Latency of selectivity for the behavioral response (ms)');
    ylabel('Number of neurons');
    title('PFC Latencies (No Learning)');
    subplot(rows,columns,6);
    hist(PMC.rx_matrix(end-POST_LEARNING_TRIALS+1:end, 2), numBins);
    xlabel('Latency of selectivity for the behavioral response (ms)');
    ylabel('Number of neurons');
    title('PMC Latencies (No Learning)');

    %% Figure 6 - Accuracy
    figure;
    plot(smooth(accuracy, 11), 'b');
    title('Accuracy');

    %% Figure 7 - Performance Tests
    % Information regarding the performance, or run-time, of this program
    if PERF_TEST
        elapsedTime = toc(start_time);
        figure; title(sprintf('TOTAL: %d, MEAN(LOOP): %d', elapsedTime, mean(loop_times))); hold on;
        plot(loop_times, 'b'); hold on;
        plot(trial_times, 'r'); hold on;
        plot(rt_calc_times, 'g');
    end

    %% Starts debug mode, allowing variables to be observed before the function ends
    keyboard;
    %% Following code can be run (copy-paste in terminal should work) to generate heat maps
    % Without border, COVIS_ENABLED
    figure;
    title('Synaptic Heatmaps');
    rows = 1; columns = 2;
    subplot(rows,columns,1);
    colormap('hot');
    imagesc(PMC_A.weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, chosen_rule));
    colorbar;
    subplot(rows,columns,2);
    colormap('hot');
    imagesc(PMC_B.weights(BORDER_SIZE:end-BORDER_SIZE, BORDER_SIZE:end-BORDER_SIZE, 1, chosen_rule));
    colorbar;
    % With border, COVIS_ENABLED
    figure;
    title('Synaptic Heatmaps');
    rows = 1; columns = 2;
    subplot(rows,columns,1);
    colormap('hot');
    imagesc(PMC_A.weights(:,:,1,chosen_rule));
    colorbar;
    subplot(rows,columns,2);
    colormap('hot');
    imagesc(PMC_B.weights(:,:,1,chosen_rule));
    colorbar;

end

% Handles the slider functionality for the synaptic weight heatmaps
% REMOVE FOR CODEGEN
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
% REMOVE FOR CODEGEN
function [f] = get_hazard_estimate(x, pts)
    [f_pdf, ~] = ksdensity(x, pts, 'function', 'pdf');
    [f_sur, ~] = ksdensity(x, pts, 'function', 'survivor');
    f = f_pdf./f_sur;
end