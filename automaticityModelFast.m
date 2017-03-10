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

function automaticityModelFast()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load parameters
    MADDOX = 'MADDOX';
    WALLIS = 'WALLIS';
    CONFIGURATIONS = {MADDOX, WALLIS};
    CONFIGURATION = CONFIGURATIONS{2};
    PARAM_CONFS = get_parameter_configurations();
    PARAMS = PARAM_CONFS(CONFIGURATION);
    
    % Programming Parameters
    PERF_TEST = 1;      % Enable/disable performance output
    SANDBOX = 0;        % Controls whether "sandbox" area executes, or main func
    PARALLEL = 0;       % Enable for parallel computing -- NOT YET SUPPORTED
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
    elseif 0
    % %% Wallis Visual Input, 100 X 100 %%
        load('wallisVisualInput.mat');
        r_x_vals = wallisVisualInput5(:,1);
        r_y_vals = wallisVisualInput5(:,2);
    end
    
    % If input grouping not enabled, set r_groups to zeros
    if not(INPUT_GROUPING)
        r_groups = zeros(1, length(r_x_vals));
    end

    %% Initialize/configure constants (though some data structure specific constants are initialized below)
    % Set behavior and number of trials
    PRE_LEARNING_TRIALS = PARAMS.PRE_LEARNING_TRIALS;   % Number of control trials run before learning trials
    LEARNING_TRIALS = PARAMS.LEARNING_TRIALS;           % Number of learning trials in automaticity experiment
    POST_LEARNING_TRIALS = PARAMS.POST_LEARNING_TRIALS; % Number of trials where no learning is involved after learning trials
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
    NOISE = PARAMS.NOISE;      % Std. dev. of noise given to PFC/PMC v; set to 0 for no noise
    
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
        'V_SCALE', 1, ...                            % scaling factor for visual input into PFC neurons
        'W_LI', 2, ...                               % lateral inhibition between PFC A / PFC B
        'DECISION_PT', PARAMS.PFC_DECISION_PT, ...   % Integral value which determines which PFC neuron acts on a visual input
        'rx_matrix', zeros(TRIALS,3) ...             % Stores information about PFC neuron reacting during trial
    );

    % PMC scaling information
    PMC = struct( ...
        'V_SCALE', 1, ...                            % can use to scale PMC visual input value if it comes out way too high
        'W_LI', 2, ...                               % lateral inhibition between PMC A / PMC B
        'DECISION_PT', PARAMS.PMC_DECISION_PT, ...   % Integral value which determines which PMC neuron acts on a visual input
        'rx_matrix', zeros(TRIALS,3) ...             % Stores information about PMC neuron reacting during trial
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

    %% Neuron constants (RSN: Regular Spiking Neuron), set for a cortical regular spiking neuron
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

    %% Neuron-related variables and matrices contained in structures, such as output matrix, output weight, etc.
    % Certain variables are also initialized here and nowhere else in the program (W_OUT)
    % Neuron.W_OUT: weight of output from one neuron to another (PFC to PMC or PMC to PFC)
    % Neuron.out: output array
    % Neuron.spikes: variables tracking spiking rate per trial
    % Neuron.v: voltage matrix (positive)
    % Neuron.u: voltage matrix (negative)
    PFC_A = struct( ...
        'W_OUT', 9, ...
        'out', zeros(1,n), ...
        'out_all', zeros(TRIALS,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0 ...
    );

    PFC_B = struct( ...
        'W_OUT', 9, ...
        'out', zeros(1,n), ...
        'out_all', zeros(TRIALS,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0 ...
    );

    PMC_A = struct( ...
        'W_OUT', 0, ...
        'out', zeros(1,n), ...
        'out_all', zeros(TRIALS,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0, ...
        'weights', INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE), ...
        'weights_avg', zeros(1,TRIALS) ...
    );

    PMC_B = struct( ...
        'W_OUT', 0, ...
        'out', zeros(1,n), ...
        'out_all', zeros(TRIALS,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0, ...
        'weights', INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE), ...
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
        PMC_A.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_A.weights,      [1 RBF.NUM_WEIGHTS]));
        PMC_B.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_B.weights,      [1 RBF.NUM_WEIGHTS]));
        % Scale v_stim values to prevent them from becoming too large
        PFC_A.v_stim = PFC_A.v_stim * PFC.V_SCALE;
        PFC_B.v_stim = PFC_B.v_stim * PFC.V_SCALE;
        PMC_A.v_stim = PMC_A.v_stim * PMC.V_SCALE;
        PMC_B.v_stim = PMC_B.v_stim * PMC.V_SCALE;

        %% Individual Time Trial
        for i=1:n-1
            % Neuron Equations
            % PFC A Neuron
            PFC_A.v(i+1)=(PFC_A.v(i) + TAU*(RSN.k*(PFC_A.v(i)-RSN.rv)*(PFC_A.v(i)-RSN.vt)-PFC_A.u(i)+ RSN.E + PFC_A.v_stim + (PMC_A.W_OUT*PMC_A.out_all(j,i)) - PFC.W_LI*PFC_B.out_all(j,i))/RSN.C) + normrnd(0,NOISE);
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
                   PFC_A.out_all(j,k)= PFC_A.out_all(j,k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
                end
            end

            % PFC B Neuron
            PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out_all(j,i)) - PFC.W_LI*PFC_A.out_all(j,i))/RSN.C) + normrnd(0,NOISE);
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
                   PFC_B.out_all(j,k)= PFC_B.out_all(j,k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
                end
            end

            % PMC_A Neuron
            PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + PMC_A.v_stim + (PFC_A.W_OUT*PFC_A.out_all(j,i)) - PMC.W_LI*PMC_B.out_all(j,i) )/RSN.C) + normrnd(0,NOISE);
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
                   PMC_A.out_all(j,k)= PMC_A.out_all(j,k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
                end
            end

            % PMC_B Neuron
            PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + PMC_B.v_stim + (PFC_B.W_OUT*PFC_B.out_all(j,i)) - PMC.W_LI*PMC_A.out_all(j,i) )/RSN.C) + normrnd(0,NOISE);
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
                   PMC_B.out_all(j,k)= PMC_B.out_all(j,k)+((t/LAMBDA)*exp((LAMBDA-t)/LAMBDA));
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
        if PARALLEL
            PFC_A.out_all(j,:) = PFC_A.out(:);
            PFC_B.out_all(j,:) = PFC_B.out(:);
            PMC_A.out_all(j,:) = PMC_A.out(:);
            PMC_B.out_all(j,:) = PMC_B.out(:);
        else
            [neuron_id, latency] = determine_reacting_neuron(PFC_A.out_all(j), PFC_B.out_all(j), PFC.DECISION_PT);
            PFC.rx_matrix(j,1:2) = [neuron_id, latency];
            PFC.rx_matrix(j,3) = r_group;
            [neuron_id, latency] = determine_reacting_neuron(PMC_A.out_all(j), PMC_B.out_all(j), PMC.DECISION_PT);
            PMC.rx_matrix(j,1:2) = [neuron_id, latency];
            PMC.rx_matrix(j,3) = r_group;
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
            PMC_A.weights = PMC_A.weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputA)*g_t_1_A.*(W_MAX - PMC_A.weights) - (Hebbian.anti_heb)*(integral_visinputA)*g_t_2_A.*PMC_A.weights);

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_A.weights = max(PMC_A.weights, 0);
            PMC_A.weights = min(PMC_A.weights, W_MAX);

            %% Calculation of Hebbian Weight for PMC_B
            % Visual input to PMC_B neuron (presynaptic)
            integral_visinputB   = trapz(PFC_B.pos_volt);
            % Activation of PMC_B neuron   (post-synaptic)
            integral_PMCBvoltage = trapz(PMC_B.pos_volt);

            % Ensures g(t)-1 and g(2)-2 are never less than zero
            g_t_1_B = max(0, integral_PMCBvoltage - Hebbian.NMDA);
            g_t_2_B = max(0, Hebbian.NMDA - integral_PMCBvoltage - Hebbian.AMPA);

            % Determine new weights of visual PMC_B synapses
            PMC_B.weights = PMC_B.weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputB)*g_t_1_B.*(W_MAX - PMC_B.weights) - (Hebbian.anti_heb)*(integral_visinputB)*g_t_2_B.*PMC_B.weights);

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_B.weights = max(PMC_B.weights, 0);
            PMC_B.weights = min(PMC_B.weights, W_MAX); 
        % Else, if not learning, do nothing
        end
        
        % Record average weight for PMC_A and PMC_B
        PMC_A.weights_avg(j) = mean(mean(PMC_A.weights));
        PMC_B.weights_avg(j) = mean(mean(PMC_B.weights));

        %% Print data to console
        fprintf('~~~ TRIAL #: %d ~~~\n', trial_number);
%         fprintf('r_y: %d\n', r_y);
%         fprintf('r_x: %d\n', r_x);
%         fprintf('PFC_A.v_stim: %d\n', PFC_A.v_stim);
%         fprintf('PFC_B.v_stim: %d\n', PFC_B.v_stim);
%         fprintf('PMC_A.v_stim: %d\n', PMC_A.v_stim);
%         fprintf('PMC_B.v_stim: %d\n', PMC_B.v_stim);
        if PERF_TEST
            loop_times(j) = toc;
        end
    end
    
    %% Post-calculations
    if PARALLEL
        % Copy variables over to temporary vars for parallel computation
        PFC_rx_matrix = PFC.rx_matrix(:,1:2);
        PMC_rx_matrix = PMC.rx_matrix(:,1:2);
        PFC_A_out_all = PFC_A.out_all(:,:);
        PFC_B_out_all = PFC_B.out_all(:,:);
        PMC_A_out_all = PMC_A.out_all(:,:);
        PMC_B_out_all = PMC_B.out_all(:,:);
        PFC_DECISION_PT = PFC.DECISION_PT;
        PMC_DECISION_PT = PMC.DECISION_PT;
        neuron_id = 0;
        latency = 0;
        keyboard;
        % Run parallel for loop to determine reaction times
        parfor j = 1:TRIALS
            [neuron_id, latency] = determine_reacting_neuron(PFC_A_out_all(j), PFC_B_out_all(j), PFC_DECISION_PT);
            PFC_rx_matrix(j,:) = [neuron_id, latency];
            [neuron_id, latency] = determine_reacting_neuron(PMC_A_out_all(j), PMC_B_out_all(j), PMC_DECISION_PT);
            PMC_rx_matrix(j,:) = [neuron_id, latency];
        end
        % Copy results into rx_matrix
        PFC.rx_matrix(:,1:2) = PFC_rx_matrix;
        PMC.rx_matrix(:,1:2) = PMC_rx_matrix;
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
    plot(TAU*(1:n),PFC_A.out_all(end));
    axis([0 n -1 10]);
    title('PFC_A Neuron Output');    % PMC_B.out

    subplot(rows,columns,6);
    plot(TAU*(1:n),PFC_B.out_all(end));
    axis([0 n -1 10]);
    title('PFC_B Neuron Output');

    subplot(rows,columns,7);
    plot(TAU*(1:n),PMC_A.out_all(end));
    axis([0 n -1 10]);
    title('PMC_A Neuron Output');

    subplot(rows,columns,8);
    plot(TAU*(1:n),PMC_B.out_all(end));
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

    if strcmp(CONFIGURATION, MADDOX)
        %% Figure 3
        % CDFs of RTs (reaction times) dependent on stimulus type -- Short, Medium, or Long
        % CDF = P(RT <= t), for each specific value t
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Return set of parameters based on argument?
% https://www.mathworks.com/matlabcentral/answers/59686-strategies-to-store-load-configuration-data
function [param_map] = get_parameter_configurations()
    % Field names of parameters (used in returned structure)
    param_names = {'CONF_NAME', 'PRE_LEARNING_TRIALS', 'LEARNING_TRIALS', 'POST_LEARNING_TRIALS', ...
                   'NOISE', 'PFC_DECISION_PT', 'PMC_DECISION_PT'};
    % Different configurations of parameters
    configurations = {'MADDOX',   0, 100,   0, 0,   4,   4; ...
                      'WALLIS',   0, 200, 100, 2, 400, 400; ...
                     };
    % Join parameter names and specified parameter configuration as structure
    param_struct = cell2struct(configurations(:,2:end), param_names(:,2:end), 2);
    param_struct_in_cells = arrayfun(@(x) x, param_struct', 'UniformOutput', false);
    % Create mapping from strings to parameter configurations
    param_map = containers.Map(configurations(:,1), param_struct_in_cells);
end

%% Return what neuron reacts to the stimuli, and the latency
% Returns neuron_id = 1 for n1, neuron_id = 2 for n2
% Not currently used -- potentially slower, and inaccurate results
function [neuron_id, latency] = determine_reacting_neuron(n1, n2, decision_pt)
    n1_latency = find(cumtrapz(n1) >= decision_pt, 1);
    n2_latency = find(cumtrapz(n2) >= decision_pt, 1);
    % n1_latency or n2_latency could be empty if the decision_pt was never reached
    % If so, set it to the maximum allowed value
    if isempty(n1_latency)
        n1_latency = length(n1);
    end
    if isempty(n2_latency)
        n2_latency = length(n2);
    end
    % Else, both n1 and n2 have valid values -- compare latencies
    if n2_latency < n1_latency
        neuron_id = 2;
        latency = n2_latency;
    else
        neuron_id = 1;
        latency = n1_latency;
    end
end

%% Find the hazard function as defined by Hazard = f(t)/S(t),
% where f(t) is the PDF and S(t) is the survivor function
function [f] = get_hazard_estimate(x, pts)
    [f_pdf, ~] = ksdensity(x, pts, 'function', 'pdf');
    [f_sur, ~] = ksdensity(x, pts, 'function', 'survivor');
    f = f_pdf./f_sur;
end