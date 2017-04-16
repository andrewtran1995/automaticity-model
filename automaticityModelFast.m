%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Model of Automaticity in Rule Based Learning %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
In general, variables that are written in all capital letters are meant
to be constant values -- set once in the beginning of the program
(variable initialization) and nowhere else

Structures are used in the program by calling struct(...).
These structures are meant to "bundle" or "group" together variables that
have some kind of commonality, e.g., belonging to the same neuron, behavior,
model, etc.
The goal is to enforce readability by standardizing the names of grouped variables
and making relationships between variables more apparent

If debugging, one can observe the workspace of the function by issuing the following
command before execution: "dbstop if error"
%}

function [sse_val] = automaticityModelFast(arg_vector) %#codegen
    %% ============================= %%
    %%%%%%%%%% INPUT PARSING %%%%%%%%%%
    %  =============================  %
    % Use default arguments if arg_vector empty
    if isempty(arg_vector)
        heb_consts = 1e-6;
        pmc_dec_pt = 400;
        noise_param = 4;
    else
        heb_consts = arg_vector(1);
        pmc_dec_pt = arg_vector(2);
        noise_param = arg_vector(3);
    end

    %% ======================================= %%
    %%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
    %  =======================================  %
    
    % Load configuration and config parameters
    MADDOX = 1;
    WALLIS = 2;
    FMRI = 3;
    CONFIGURATION = FMRI;
    PARAMS = struct( ...
        'PRE_LEARNING_TRIALS', 0, ...
        'LEARNING_TRIALS', 11520, ...
        'POST_LEARNING_TRIALS', 0, ...
        'NOISE', 2, ...
        'PFC_DECISION_PT', 400, ...
        'PMC_DECISION_PT', 400 ...
    );
    
    % Override parameter values if they were specified as inputs
    if nargin ~= 0
        PARAMS.HEB_CONSTS = heb_consts;
        PARAMS.PMC_DECISION_PT = pmc_dec_pt;
        PARAMS.NOISE = noise_param;
    end
    
    % Struct to contain meta-data of FMRI configuration
    FMRI_META = struct('NUM_TRIALS', 11520, ...
                       'SES_1',      1:480, 'SES_4',    1681:2160, ...
                       'SES_10', 5161:5640, 'SES_20', 11041:11520);
    
    % Programming Parameters
    OPTIMIZATION_RUN = 1;
    
    %% Load visual stimulus matrix
    % %% Random Visual Input, 100 x 100 %%
    if 0
        loaded_input = load('datasets/randomVisualInput.mat');
        r_x_vals = loaded_input.r_x_mat;
        r_y_vals = loaded_input.r_y_mat;
    elseif CONFIGURATION == MADDOX
    % %% Random Visual Input to Maddox Grid, 100 X 100 %%
        loaded_input = load('datasets/maddoxVisualInput.mat');
        r_x_vals = loaded_input.maddoxVisualInput(:, 1);
        r_y_vals = loaded_input.maddoxVisualInput(:, 2);
        r_groups = loaded_input.maddoxVisualInput(:, 3);
    elseif CONFIGURATION == WALLIS
    % %% Wallis Visual Input, 100 X 100 %%
        loaded_input = load('datasets/wallisVisualInput.mat');
        r_x_vals = loaded_input.wallisVisualInput5(:,1);
        r_y_vals = loaded_input.wallisVisualInput5(:,2);
    elseif CONFIGURATION == FMRI
        loaded_input = load('datasets/fMRI_data.mat');
        r_x_vals = loaded_input.r_x_mat;
        r_y_vals = loaded_input.r_y_mat;
    end
    
    % If input grouping not enabled, set r_groups to zeros
    if CONFIGURATION ~= MADDOX
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

    % Accuracy matrix, where first dimension has two rows (1 = PFC; 2 = PMC)
    % and second dimension is trial number
    % Each element is a boolean indicating if the correct neuron reacted that trial
    accuracy = zeros(TRIALS, 1);

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
        'heb_coef', PARAMS.HEB_CONSTS, ...
        'anti_heb', PARAMS.HEB_CONSTS, ...
        'NMDA',     1500, ...
        'AMPA',     750 ...
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
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0 ...
    );
    PFC_B = struct( ...
        'W_OUT', 9, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0 ...
    );
    % Use a simplified weights matrix for FMRI (due to large amount of trials)
    if CONFIGURATION == FMRI
        WEIGHTS_MATRIX = INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE);
    else
        WEIGHTS_MATRIX = INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,TRIALS);
    end
    PMC_A = struct( ...
        'W_OUT', 0, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0, ...
        'weights', WEIGHTS_MATRIX, ...
        'weights_avg', zeros(TRIALS,1) ...
    );
    PMC_B = struct( ...
        'W_OUT', 0, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0, ...
        'weights', WEIGHTS_MATRIX, ...
        'weights_avg', zeros(TRIALS,1) ...
    );

    %% ============================ %%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %  ============================  %

    %% Pre-calculations (for performance reasons)
    % Calculate lambda values for individual trials
    t = (0:n)';
    LAMBDA_PRECALC = (t/LAMBDA).*exp((LAMBDA-t)/LAMBDA);
    
    %% Learning trials
    for j=1:TRIALS
        %% Initialize appropriate variables for each loop
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
        if CONFIGURATION == FMRI || j==1
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
        PFC_A.v_stim = sum(sum(    RBF.rbv(:, 1:GRID_SIZE/2)));
        PFC_B.v_stim = sum(sum(RBF.rbv(:, GRID_SIZE/2+1:end)));
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(sum(RBF.rbv(:,:).*PMC_A_weights));
        PMC_B.v_stim = sum(sum(RBF.rbv(:,:).*PMC_B_weights));
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
            if PFC_A.v(i+1)>=RSN.vpeak
                PFC_A.v(i)= RSN.vpeak;
                PFC_A.v(i+1)= RSN.c;
                PFC_A.u(i+1)= PFC_A.u(i+1)+ RSN.d;
            end
            if PFC_A.v(i) >= RSN.vpeak
                PFC_A.out(i:n) = PFC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end

            % PFC B Neuron
            PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out(i)) - PFC.W_LI*PFC_A.out(i))/RSN.C) + normrnd(0,NOISE);
            PFC_B.u(i+1)=PFC_B.u(i)+TAU*RSN.a*(RSN.b*(PFC_B.v(i)-RSN.rv)-PFC_B.u(i));
            if PFC_B.v(i+1)>=RSN.vpeak
                PFC_B.v(i)= RSN.vpeak;
                PFC_B.v(i+1)= RSN.c;
                PFC_B.u(i+1)= PFC_B.u(i+1)+ RSN.d;
            end
            if PFC_B.v(i) >= RSN.vpeak
                PFC_B.out(i:n) = PFC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end

            % PMC_A Neuron
            PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + PMC_A.v_stim + (PFC_A.W_OUT*PFC_A.out(i)) - PMC.W_LI*PMC_B.out(i) )/RSN.C) + normrnd(0,NOISE);
            PMC_A.u(i+1)=PMC_A.u(i)+TAU*RSN.a*(RSN.b*(PMC_A.v(i)-RSN.rv)-PMC_A.u(i));
            if PMC_A.v(i+1)>=RSN.vpeak
                PMC_A.v(i)= RSN.vpeak;
                PMC_A.v(i+1)= RSN.c;
                PMC_A.u(i+1)= PMC_A.u(i+1)+ RSN.d;
            end
            if PMC_A.v(i) >= RSN.vpeak
                PMC_A.out(i:n) = PMC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end

            % PMC_B Neuron
            PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + PMC_B.v_stim + (PFC_B.W_OUT*PFC_B.out(i)) - PMC.W_LI*PMC_A.out(i) )/RSN.C) + normrnd(0,NOISE);
            PMC_B.u(i+1)=PMC_B.u(i)+TAU*RSN.a*(RSN.b*(PMC_B.v(i)-RSN.rv)-PMC_B.u(i));
            if PMC_B.v(i+1)>=RSN.vpeak
                PMC_B.v(i)= RSN.vpeak;
                PMC_B.v(i+1)= RSN.c;
                PMC_B.u(i+1)= PMC_B.u(i+1)+ RSN.d;
            end
            if PMC_B.v(i) >= RSN.vpeak
                PMC_B.out(i:n) = PMC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
        end
        % Count number of spikes
        PFC_A.spikes = nnz(PFC_A.v >= RSN.vpeak);
        PFC_B.spikes = nnz(PFC_B.v >= RSN.vpeak);
        PMC_A.spikes = nnz(PMC_A.v >= RSN.vpeak);
        PMC_B.spikes = nnz(PMC_B.v >= RSN.vpeak);
        % Record voltage value if positive. Else, do nothing.
        % For computation of integral
        PFC_A.pos_volt(PFC_A.v > 0) = PFC_A.v(PFC_A.v > 0);
        PFC_B.pos_volt(PFC_B.v > 0) = PFC_B.v(PFC_B.v > 0);
        PMC_A.pos_volt(PMC_A.v > 0) = PMC_A.v(PMC_A.v > 0);
        PMC_B.pos_volt(PMC_B.v > 0) = PMC_B.v(PMC_B.v > 0);

        %% Determine decision neuron and reaction time, and record accuracy
        [neuron_id_PFC, latency] = determine_reacting_neuron(PFC_A.out, PFC_B.out, PFC.DECISION_PT);
        PFC.rx_matrix(j,1:2) = [neuron_id_PFC, latency];
        PFC.rx_matrix(j,3) = r_group;
        [neuron_id_PMC, latency] = determine_reacting_neuron(PMC_A.out, PMC_B.out, PMC.DECISION_PT);
        PMC.rx_matrix(j,1:2) = [neuron_id_PMC, latency];
        PMC.rx_matrix(j,3) = r_group;
        accuracy(j) = neuron_id_PFC == neuron_id_PMC;

        %% Weight change calculations
        if CONFIGURATION == FMRI
            k = 1;
        else
            k = j;
        end
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
            PMC_A.weights(:,:,k) = PMC_A_weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputA)*g_t_1_A.*(W_MAX - PMC_A_weights) - (Hebbian.anti_heb)*(integral_visinputA)*g_t_2_A.*PMC_A_weights);

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_A.weights(:,:,k) = max(PMC_A.weights(:,:,k), 0);
            PMC_A.weights(:,:,k) = min(PMC_A.weights(:,:,k), W_MAX);

            %% Calculation of Hebbian Weight for PMC_B
            % Visual input to PMC_B neuron (presynaptic)
            integral_visinputB   = trapz(PFC_B.pos_volt);
            % Activation of PMC_B neuron   (post-synaptic)
            integral_PMCBvoltage = trapz(PMC_B.pos_volt);

            % Ensures g(t)-1 and g(2)-2 are never less than zero
            g_t_1_B = max(0, integral_PMCBvoltage - Hebbian.NMDA);
            g_t_2_B = max(0, Hebbian.NMDA - integral_PMCBvoltage - Hebbian.AMPA);

            % Determine new weights of visual PMC_B synapses
            PMC_B.weights(:,:,k) = PMC_B_weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputB)*g_t_1_B.*(W_MAX - PMC_B_weights) - (Hebbian.anti_heb)*(integral_visinputB)*g_t_2_B.*PMC_B_weights);

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_B.weights(:,:,k) = max(PMC_B.weights(:,:,k), 0);
            PMC_B.weights(:,:,k) = min(PMC_B.weights(:,:,k), W_MAX); 
        % Else, if not learning, set new weights to previous weights
        else
            PMC_A.weights(:,:,k) = PMC_A_weights;
            PMC_B.weights(:,:,k) = PMC_B_weights;
        end
        
        % Record average weight for PMC_A and PMC_B
        PMC_A.weights_avg(j) = mean(mean(PMC_A.weights(:,:,k)));
        PMC_B.weights_avg(j) = mean(mean(PMC_B.weights(:,:,k)));

        %% Print data to console
%         if mod(j,1) == 0
%             fprintf('~~~ TRIAL #: %d ~~~\n', int32(j));
%         end
    end
    
    %% ========================================= %%
    %%%%%%%%%% OPTIMIZATION CALCULATIONS %%%%%%%%%%
    %  =========================================  %
    % Return prematurely if we are optimizing (e.g., particle swarm optimization)
    % Calculate Sum of Squared Errors of Prediction (SSE)
    if CONFIGURATION == FMRI && OPTIMIZATION_RUN
        target = load('fmri_data/means1dCondition.mat');
        % Calculate Mean Accuracy for trials from Session 4, 10, and 20
        output_acc = [mean(accuracy(FMRI_META.SES_1)), ...
                      mean(accuracy(FMRI_META.SES_4)), ...
                      mean(accuracy(FMRI_META.SES_10)), ...
                      mean(accuracy(FMRI_META.SES_20))];
        % Calculate Mean Median RT for trials from Session 4, 10, and 20
        % Reaction times must be converted from ms to seconds
        norm_output_rt = [median(PMC.rx_matrix(FMRI_META.SES_1,2)), ...
                          median(PMC.rx_matrix(FMRI_META.SES_4,2)), ...
                          median(PMC.rx_matrix(FMRI_META.SES_10,2)), ...
                          median(PMC.rx_matrix(FMRI_META.SES_20,2))]./1000;
        disp([output_acc; norm_output_rt]);
        sse_val = sum(sum((target.means1dCondition - [output_acc;norm_output_rt]).^2));
        return
    end
end

%% =============================== %%
%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%
%  ===============================  %

% Return what neuron reacts to the stimuli, and the latency
% Returns neuron_id = 1 for n1, neuron_id = 2 for n2
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