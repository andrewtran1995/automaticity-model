%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Model of Automaticity in Rule Based Learning %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% NOTES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In general, variables that are written in all capital letters are meant
% to be constant values -- set once in the beginning of the program
% (variable initialization) and nowhere else

% Structures are used in the program by calling struct(...).
% These structures are meant to "bundle" or "group" together variables that
% have some kind of commonality, e.g., belonging to the same neuron, behavior,
% model, etc.
% The goal is to enforce readability by standardizing the names of grouped variables
% and making relationships between variables more apparent

% If debugging, one can observe the workspace of the function by issuing the following
% command before execution: "dbstop if error"

function automaticityModel()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Programming Parameters
    PERF_TEST = 1;      % Enable/disable performance output
    SANDBOX = 0;        % Controls whether "sandbox" area executes, or main func
    INPUT_GROUPING = 0; % Set to 1 (true) if loaded visual input groups stimulus
    if PERF_TEST
        startTime = tic;
    end
    
    %% Load visual stimulus matrix
    % %% Random Visual Input, 100 x 100 %%
    if 1
        load('randomVisualInput.mat');
        r_x_vals = r_x_mat;
        r_y_vals = r_y_mat;
    elseif 0
    % %% Random Visual Input to Maddox Grid, 100 X 100 %%
        load('maddoxVisualInput.mat');
        INPUT_GROUPING = 1;
        r_x_vals = maddoxVisualInput(:, 1);
        r_y_vals = maddoxVisualInput(:, 2);
        r_groups = maddoxVisualInput(:, 3);
    end
    
    % If input grouping not enabled, set r_groups to zeros
    if not(INPUT_GROUPING)
        r_groups = zeros(1, length(r_x_vals));
    end

    %% Initialize/configure constants (though some data structure specific constants are initialized below)
    % Set behavior and number of trials
    PRE_LEARNING_TRIALS = 100;  % Number of control trials run before learning trials
    LEARNING_TRIALS = 100;        % Number of learning trials in automaticity experiment
    POST_LEARNING_TRIALS = 100; % Number of trials where no learning is involved after learning trials
    TRIALS = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS; % Total number of trials
    % Create matrix to store information on when learning should occur
    LEARNING = [zeros(1, PRE_LEARNING_TRIALS), ...
                ones(1,LEARNING_TRIALS), ...
                zeros(1,POST_LEARNING_TRIALS)]; 
    % Convenience variables
    LEARNING_IDX = PRE_LEARNING_TRIALS+1:PRE_LEARNING_TRIALS+LEARNING_TRIALS;

    % Set properties of grid
    STIM_GRID_SIZE = 100;      % Length of side of square grid used for visual input; shoudl be an even number
    BORDER_SIZE = 20;          % Width of border used to pad the grid such that visual stimulus on the edge still has an appropriate effect
    GRID_SIZE = STIM_GRID_SIZE + 2*BORDER_SIZE; % Total length of grid, i.e., the stimulus grid size and the border
    
    % Other parameters
    n = 1000;                  % Time period for one trial (in milliseconds)
    TAU = 1;
    LAMBDA = 20;               % Lambda Value
    W_MAX = 10;                % maximum possible weight for Hebbian Synapses
    INIT_PMC_WEIGHT = 0.08;    % Initial weight for PMC neurons
    NOISE = 2;                 % Std. dev. of noise given to PFC/PMC v; set to 0 for no noise
    
    loop_times = zeros(1, TRIALS); % Records how much time was needed for each loop

    % Quantity of Visual Stimulus
    Visual = struct( ...
        'stim', 50 ...
    );

    % Radial Basis Function
    [X, Y] = meshgrid(1:GRID_SIZE, 1:GRID_SIZE);
    RBF = struct( ...
        'RADIUS', 2, ...
        'rbv', zeros(GRID_SIZE), ...
        'X', X, ...
        'Y', Y, ...
        'HALF_NUM_WEIGHTS', GRID_SIZE/2 * GRID_SIZE, ...
        'NUM_WEIGHTS', GRID_SIZE * GRID_SIZE ...
    );

    %% General settings for PFC, PMC neurons
    % Note that rx_matrix is big enough for both learning trials and no-learning trials to allow for comparisons
    % PFC scaling information
    PFC = struct( ...
        'V_SCALE', 1, ...                                   % scaling factor for visual input into PFC neurons
        'W_LI', 2, ...                                      % lateral inhibition between PFC A / PFC B
        'DECISION_PT', 4, ...                               % Integral value which determines which PFC neuron acts on a visual input
        'rx_matrix', zeros(TRIALS,3) ... % Stores information about PFC neuron reacting during trial
    );

    % PMC scaling information
    PMC = struct( ...
        'V_SCALE', 1, ...                                   % can use to scale PMC visual input value if it comes out way too high
        'W_LI', 2, ...                                      % lateral inhibition between PMC A / PMC B
        'DECISION_PT', 4, ...                               % Integral value which determines which PMC neuron acts on a visual input
        'rx_matrix', zeros(TRIALS,3) ... % Stores information about PMC neuron reacting during trial
    );

    %% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses)
    % NMDA - upper threshold
    % AMPA - lower threshold
    % Strengthening occurs if integral_PMCAvoltage > Hebbian.NMDA
    % Weakening occurs if Hebbian.NMDA - integral_PMCAvoltage - Hebbian.AMPA > 0, i.e., only if integral_PMCAvoltage < Hebbian.NMDA + Hebbian.AMPA
    Hebbian = struct( ...
        'heb_coef', 0.00000001, ...
        'anti_heb', 0.00000001, ...
        'NMDA', 1500, ...
        'AMPA', 0 ...
    );

    %% Initialize data structures (some constants contained in data structures)
    % Neuron constants (RSN: Regular Spiking Neuron)
    % Set for a cortical regular spiking neuron
    % Usage: RSN.C to get C value
    RSN = struct( ...
        'C', 100, ...
        'rv', -60, ...
        'vt', -40, ...
        'k', 0.7, ...
        'a', 0.03, ...
        'b', -2, ...
        'c', -50, ...
        'd', 100, ...
        'vpeak', 35, ...
        'E', 60 ...
    );

    % Neuron-related variables and matrices contained in
    % "structures", such as output matrix, output weight, etc.
    % Certain variables are also initialized here and nowhere
    % else in the program (W_OUT)
    % Neuron.W_OUT: weight of output from one neuron to another (PFC to PMC or PMC to PFC)
    % Neuron.out: output array
    % Neuron.spikes: variables tracking spiking rate per trial
    % Neuron.v: voltage matrix (positive)
    % Neuron.u: voltage matrix (negative)
    PFC_A = struct( ...
        'W_OUT', 9, ...
        'out', zeros(1,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0 ...
    );

    PFC_B = struct( ...
        'W_OUT', 9, ...
        'out', zeros(1,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0 ...
    );

    % Please note that PMC_A.weights and PMC_B.weights has TRIALS+1 amount of entries
    % where weights(1) is the intialization matrix, and weights(i+1) is the true
    % synaptic weight values for trial i
    % The initialization matrix is removed at the end of the calculation section so that
    % the index matches the synaptic weight values for that trial
    PMC_A = struct( ...
        'W_OUT', 0, ...
        'out', zeros(1,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0, ...
        'weights', INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,TRIALS), ...
        'weights_avg', zeros(1,TRIALS) ...
    );

    PMC_B = struct( ...
        'W_OUT', 0, ...
        'out', zeros(1,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0, ...
        'weights', INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,TRIALS), ...
        'weights_avg', zeros(1,TRIALS) ...
    );

    %% Sandbox area
    % Placed after all values are initalized, and serves as an area where code can be prototyped and tested
    % (for validity or performance reasons) before being implemented into the main body of the function
    if PERF_TEST && SANDBOX
        return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Pre-calculations
    % Precompute RBF.rbv matrices
%     parfor j=1:TRIALS
%         
%     end
    
    %% Learning trials
    trial_number = 0;

    for j=1:TRIALS
        if PERF_TEST
            tic;
        end
        %% Initialize appropriate variables for each loop
        trial_number = trial_number + 1;    % track number of current trial

        % variables tracking spiking rate in each neuron
        PFC_A.spikes = 0;       PMC_A.spikes = 0;
        PFC_B.spikes = 0;       PMC_B.spikes = 0;

        % variables keep track of positive voltage values for calculation of
        % integral (for Hebbian learning equation)
        PFC_A.pos_volt(:) = 0;  PMC_A.pos_volt(:) = 0;
        PFC_B.pos_volt(:) = 0;  PMC_B.pos_volt(:) = 0;

        % Re-initialize Neuron Voltage Matrices (set all v to RSN.rv; all u to 0)
        PFC_A.v(:) = RSN.rv;     PFC_A.u(:) = 0;
        PFC_B.v(:) = RSN.rv;     PFC_B.u(:) = 0;
        PMC_A.v(:) = RSN.rv;     PMC_A.u(:) = 0;
        PMC_B.v(:) = RSN.rv;     PMC_B.u(:) = 0;

        % Re-initialize Neuron Output Matrices
        PFC_A.out(:) = 0;        PMC_A.out(:) = 0;
        PFC_B.out(:) = 0;        PMC_B.out(:) = 0;
        
        % Set PMC weights; if first trial, set to initial weights
        if j==1
            PMC_A_weights = PMC_A.weights(:,:,1);
            PMC_B_weights = PMC_B.weights(:,:,1);
        % Else, set weights to results of previous trial
        else
            PMC_A_weights = PMC_A.weights(:,:,j-1);
            PMC_B_weights = PMC_B.weights(:,:,j-1);
        end

        % Determine visual stimulus in range [1, GRID_SIZE] to pick random gabor for each trial, padded with
        % the BORDER_SIZE such that the visual stimulus is accounted for properly
        r_y = r_y_vals(j) + BORDER_SIZE;
        r_x = r_x_vals(j) + BORDER_SIZE;
        r_group = r_groups(j);

        %% Radial Basis Function (RBF) Implementation
        % Calculate RBF grid
        RBF.rbv(:, :) = exp( -(sqrt((r_y-RBF.Y).^2 + (r_x-RBF.X).^2))/RBF.RADIUS ) * Visual.stim;
        % Sum appropriate RBF values to find PFC_A and PFC_B v_stim values
        PFC_A.v_stim = sum(reshape(    RBF.rbv(:, 1:GRID_SIZE/2), [1 RBF.HALF_NUM_WEIGHTS]));
        PFC_B.v_stim = sum(reshape(RBF.rbv(:, GRID_SIZE/2+1:end), [1 RBF.HALF_NUM_WEIGHTS]));
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_A_weights,      [1 RBF.NUM_WEIGHTS]));
        PMC_B.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_B_weights,      [1 RBF.NUM_WEIGHTS]));
        % Scale v_stim values to prevent them from becoming too large
        PFC_A.v_stim = PFC_A.v_stim * PFC.V_SCALE;
        PFC_B.v_stim = PFC_B.v_stim * PFC.V_SCALE;
        PMC_A.v_stim = PMC_A.v_stim * PMC.V_SCALE;
        PMC_B.v_stim = PMC_B.v_stim * PMC.V_SCALE;

        %% Individual Time Trial
        for i=1:n-1
            % Neuron Equations
            % PFC A Neuron
            PFC_A.v(i+1)=(PFC_A.v(i) + TAU*(RSN.k*(PFC_A.v(i)-RSN.rv)*(PFC_A.v(i)-RSN.vt)-PFC_A.u(i)+ RSN.E + PFC_A.v_stim + (PMC_A.W_OUT*PMC_A.out(i)) - PFC.W_LI*PFC_B.out(i))/RSN.C) + normrnd(0,NOISE);
            PFC_A.u(i+1)=PFC_A.u(i)+TAU*RSN.a*(RSN.b*(PFC_A.v(i)-RSN.rv)-PFC_A.u(i));
            if PFC_A.v(i+1)>=RSN.vpeak;
                PFC_A.v(i)= RSN.vpeak;
                PFC_A.v(i+1)= RSN.c;
                PFC_A.u(i+1)= PFC_A.u(i+1)+ RSN.d;
            end

            if (PFC_A.v(i) >= RSN.vpeak)
                PFC_A.spikes = PFC_A.spikes + 1;
                for k=i:n
                   t= k-i;
                   PFC_A.out(k)= PFC_A.out(k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
                end
            end

            % PFC B Neuron
            PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out(i)) - PFC.W_LI*PFC_A.out(i))/RSN.C) + normrnd(0,NOISE);
            PFC_B.u(i+1)=PFC_B.u(i)+TAU*RSN.a*(RSN.b*(PFC_B.v(i)-RSN.rv)-PFC_B.u(i));
            if PFC_B.v(i+1)>=RSN.vpeak;
                PFC_B.v(i)= RSN.vpeak;
                PFC_B.v(i+1)= RSN.c;
                PFC_B.u(i+1)= PFC_B.u(i+1)+ RSN.d;
            end

            if (PFC_B.v(i) >= RSN.vpeak)
                PFC_B.spikes = PFC_B.spikes + 1;
                for k=i:n
                   t= k-i;
                   PFC_B.out(k)= PFC_B.out(k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
                end
            end

            % PMC_A Neuron
            PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + PMC_A.v_stim + (PFC_A.W_OUT*PFC_A.out(i)) - PMC.W_LI*PMC_B.out(i) )/RSN.C) + normrnd(0,NOISE);
            PMC_A.u(i+1)=PMC_A.u(i)+TAU*RSN.a*(RSN.b*(PMC_A.v(i)-RSN.rv)-PMC_A.u(i));
            if PMC_A.v(i+1)>=RSN.vpeak;
                PMC_A.v(i)= RSN.vpeak;
                PMC_A.v(i+1)= RSN.c;
                PMC_A.u(i+1)= PMC_A.u(i+1)+ RSN.d;
            end

            if (PMC_A.v(i) >= RSN.vpeak)
                PMC_A.spikes = PMC_A.spikes + 1;
                for k=i:n
                   t= k-i;
                   PMC_A.out(k)= PMC_A.out(k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
                end
            end

            % PMC_B Neuron
            PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + PMC_B.v_stim + (PFC_B.W_OUT*PFC_B.out(i)) - PMC.W_LI*PMC_A.out(i) )/RSN.C) + normrnd(0,NOISE);
            PMC_B.u(i+1)=PMC_B.u(i)+TAU*RSN.a*(RSN.b*(PMC_B.v(i)-RSN.rv)-PMC_B.u(i));
            if PMC_B.v(i+1)>=RSN.vpeak;
                PMC_B.v(i)= RSN.vpeak;
                PMC_B.v(i+1)= RSN.c;
                PMC_B.u(i+1)= PMC_B.u(i+1)+ RSN.d;
            end

            if (PMC_B.v(i) >= RSN.vpeak)
                PMC_B.spikes = PMC_B.spikes + 1;
                for k=i:n
                   t= k-i;
                   PMC_B.out(k)= PMC_B.out(k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
                end
            end

            % Record voltage value if positive. Else, do nothing.
            % For computation of integral
            if PFC_A.v(i) > 0; PFC_A.pos_volt(i) = PFC_A.v(i); end
            if PFC_B.v(i) > 0; PFC_B.pos_volt(i) = PFC_B.v(i); end
            if PMC_A.v(i) > 0; PMC_A.pos_volt(i) = PMC_A.v(i); end
            if PMC_B.v(i) > 0; PMC_B.pos_volt(i) = PMC_B.v(i); end

        end

        %% Determine decision neuron and reaction time
%         PMC.rx_matrix(j,1:2) = determine_reacting_neuron(PMC_A.out, PMC_B.out, PMC.DECISION_PT);
%         PMC.rx_matrix(j,3) = r_group;
        for i=1:n
            if trapz(PFC_A.out(1:i)) >= PFC.DECISION_PT
                PFC.rx_matrix(j,:) = [1, i, r_group];
                break;
            elseif trapz(PFC_B.out(1:i)) >= PFC.DECISION_PT
                PFC.rx_matrix(j,:) = [2, i, r_group];
                break;
            else
                continue;
            end
        end
        for i=1:n
            % If PMC_A meets the decision point sooner, indicate it in the
            % first column with a '1'
            if trapz(PMC_A.out(1:i)) >= PMC.DECISION_PT
                PMC.rx_matrix(j,:) = [1, i, r_group];
                break;
            % Else, indicate PMC_B with '2'
            elseif trapz(PMC_B.out(1:i)) >= PMC.DECISION_PT
                PMC.rx_matrix(j,:) = [2, i, r_group];
                break;
            else
                continue;
            end
        end

        %% Weight change calculations
        if LEARNING(j)
            %% Calculation of Hebbian Weight for PMC_A
            % Visual input to PMC_A neuron (presynaptic)
            integral_visinputA   = trapz(PFC_A.pos_volt);
            % Activation of PMC_A neuron   (post-synaptic)
            integral_PMCAvoltage = trapz(PMC_A.pos_volt);

            % Ensure g(t)-1 and g(2)-2 are never less than zero
            g_t_1_A = max(0, integral_PMCAvoltage - Hebbian.NMDA);
            g_t_2_A = max(0, Hebbian.NMDA - integral_PMCAvoltage - Hebbian.AMPA);

            % Determine new weights of visual PMC_A synapses
            PMC_A.weights(:,:,j) = PMC_A_weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputA)*g_t_1_A.*(W_MAX - PMC_A_weights) - (Hebbian.anti_heb)*(integral_visinputA)*g_t_2_A.*PMC_A_weights);

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_A.weights(:,:,j) = max(PMC_A.weights(:,:,j), 0);
            PMC_A.weights(:,:,j) = min(PMC_A.weights(:,:,j), W_MAX);

            %% Calculation of Hebbian Weight for PMC_B
            % Visual input to PMC_B neuron (presynaptic)
            integral_visinputB   = trapz(PFC_B.pos_volt);
            % Activation of PMC_B neuron   (post-synaptic)
            integral_PMCBvoltage = trapz(PMC_B.pos_volt);

            % Ensures g(t)-1 and g(2)-2 are never less than zero
            g_t_1_B = max(0, integral_PMCBvoltage - Hebbian.NMDA);
            g_t_2_B = max(0, Hebbian.NMDA - integral_PMCBvoltage - Hebbian.AMPA);

            % Determine new weights of visual PMC_B synapses
            PMC_B.weights(:,:,j) = PMC_B_weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputB)*g_t_1_B.*(W_MAX - PMC_B_weights) - (Hebbian.anti_heb)*(integral_visinputB)*g_t_2_B.*PMC_B_weights);

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_B.weights(:,:,j) = max(PMC_B.weights(:,:,j), 0);
            PMC_B.weights(:,:,j) = min(PMC_B.weights(:,:,j), W_MAX); 
        % Else, if not learning, set new weights to previous weights
        else
            PMC_A.weights(:,:,j) = PMC_A_weights;
            PMC_B.weights(:,:,j) = PMC_B_weights;
        end
        
        % Record average weight for PMC_A and PMC_B
        PMC_A.weights_avg(j) = mean(mean(PMC_A.weights(:,:,j)));
        PMC_B.weights_avg(j) = mean(mean(PMC_B.weights(:,:,j)));

        %% Print data to console
        fprintf('~~~ TRIAL #: %d ~~~\n', trial_number);
        fprintf('r_y: %d\n', r_y);
        fprintf('r_x: %d\n', r_x);
        fprintf('PFC_A.v_stim: %d\n', PFC_A.v_stim);
        fprintf('PFC_B.v_stim: %d\n', PFC_B.v_stim);
        fprintf('PMC_A.v_stim: %d\n', PMC_A.v_stim);
        fprintf('PMC_B.v_stim: %d\n', PMC_B.v_stim);
        if PERF_TEST
            loop_times(j) = toc;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Figure 1 - neuron information from last trial or throughout trials
    figure;
    title('Neuron Information from Last Trial, Rx Times, Etc.');

    % Plot items
    rows = 3;
    columns = 4;

    subplot(rows,columns,1);
    plot(TAU*(1:n),PFC_A.v);
    axis([0 n -100 100]);
    title('PFC_A Neuron Voltage');

    subplot(rows,columns,2);
    plot(TAU*(1:n),PFC_B.v);
    axis([0 n -100 100]);
    title('PFC_B Neuron Voltage');

    subplot(rows,columns,3);
    plot(TAU*(1:n),PMC_A.v);
    axis([0 n -100 100]);
    title('PMC_A Neuron Voltage');

    subplot(rows,columns,4);
    plot(TAU*(1:n),PMC_B.v);
    axis([0 n -100 100]);
    title('PMC_B Neuron Voltage');

    subplot(rows,columns,5);
    plot(TAU*(1:n),PFC_A.out);
    axis([0 n -1 10]);
    title('PFC_A Neuron Output');    % PMC_B.out

    subplot(rows,columns,6);
    plot(TAU*(1:n),PFC_B.out);
    axis([0 n -1 10]);
    title('PFC_B Neuron Output');

    subplot(rows,columns,7);
    plot(TAU*(1:n),PMC_A.out);
    axis([0 n -1 10]);
    title('PMC_A Neuron Output');

    subplot(rows,columns,8);
    plot(TAU*(1:n),PMC_B.out);
    axis([0 n -1 10]);
    title('PMC_B Neuron Output');

    subplot(rows,columns,9);
    colormap('hot');
    imagesc(RBF.rbv(BORDER_SIZE:end-BORDER_SIZE-1,BORDER_SIZE:end-BORDER_SIZE-1,:));
    title(sprintf('Stimulus: (%d,%d); Weight: %d', r_y, r_x, Visual.stim));

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
    

    
    %% Figure 2
    % Synaptic weight heatmaps with sliders to allow the observation of the heatmap at different intervals in time
    % Only relevant if any learning trials were conducted
    if LEARNING_TRIALS > 0
        figure;
        title('Synaptic Heatmaps');
        rows = 1;
        columns = 2;
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
    
    %% Figure 3
    % CDFs of RTs (reaction times) dependent on stimulus type -- Short, Medium, or Long
    % CDF = P(RT <= t), for each specific value t
    if INPUT_GROUPING
        % Set-up
        PMC_S = PMC.rx_matrix(LEARNING_IDX,3) == 'S';
        PMC_M = PMC.rx_matrix(LEARNING_IDX,3) == 'M';
        PMC_L = PMC.rx_matrix(LEARNING_IDX,3) == 'L';        
        figure;
        title('CDFs of PMC Rx Times (Grouped by Distance)');

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
        figure;
        title('PMC Rx Times Hazard Functions');

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
    rows = 3;
    columns = 2;
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
    
    %% Figure 6 - Performance Tests
    % Information regarding the performance, or run-time, of this program
    if PERF_TEST
        elapsedTime = toc(startTime);
        figure;
        plot(loop_times);
        title(sprintf('TOTAL: %d, MEAN(LOOP): %d', elapsedTime, mean(loop_times)));
    end
    
    %% Starts debug mode, allowing variables to be observed before the function ends
    keyboard;
end

%% Return what neuron reacts to the stimuli, and the latency
% Returns neuron_id = 1 for n1, neuron_id = 2 for n2
% Not currently used -- potentially slower, and inaccurate results
function [neuron_id, latency] = determine_reacting_neuron(n1, n2, decision_pt)
    n1_latency = find(cumtrapz(n1) >= decision_pt, 1);
    n2_latency = find(cumtrapz(n2) >= decision_pt, 1);
    % If n1 or n2 only contains zeroes, the alternative must be the reacting neuron
    if any(n1) == 0
        neuron_id = 2;
        latency = n2_latency;
    elseif any(n2) == 0
        neuron_id = 1;
        latency = n1_latency;
    % Else, both n1 and n2 have valid values -- compare latencies
    elseif n2_latency < n1_latency
        neuron_id = 2;
        latency = n2_latency;
    else
        neuron_id = 1;
        latency = n1_latency;
    end
end

%% Handles the slider functionality for the synaptic weight heatmaps
function synaptic_slider_callback(src, ~, position, data, neuron_name)
    subplot(1, 2, position, 'replace');
    trial_num = round(get(src, 'value'));
    set(src, 'value', trial_num);
    colormap('hot');
    imagesc(data(:,:,trial_num));
    colorbar;
    title(sprintf('%s Synaptic Heatmap, Trial %d\n', neuron_name, trial_num'));
end

%% Find the hazard function as defined by Hazard = f(t)/S(t),
% where f(t) is the PDF and S(t) is the survivor function
function [f] = get_hazard_estimate(x, pts)
    [f_pdf, ~] = ksdensity(x, pts, 'function', 'pdf');
    [f_sur, ~] = ksdensity(x, pts, 'function', 'survivor');
    f = f_pdf./f_sur;
end