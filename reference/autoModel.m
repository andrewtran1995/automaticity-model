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



function [Results] = autoModel(arg_vector) %#codegen
    

    %% ======================================= %%
    %%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
    %  =======================================  %
    
    
    
    % Load configuration and config parameters
    MADDOX = 1; WALLIS = 2; FMRI = 3;
    CONFIGURATIONS = {'MADDOX', 'WALLIS', 'FMRI'};
    CONFIGURATION = FMRI;
    PARAMS = get_parameters(CONFIGURATIONS{CONFIGURATION});
    
    % Override parameter values if they were specified as inputs
    if nargin ~= 0 && ~isempty(arg_vector)
        PARAMS.HEB_CONSTS      = arg_vector(1);
        PARAMS.ANTI_HEB_CONSTS = arg_vector(2);
        PARAMS.PMC_DECISION_PT = arg_vector(3);
        PARAMS.NOISE           = arg_vector(4);
        PARAMS.NMDA            = arg_vector(5);
        PARAMS.AMPA            = arg_vector(6);
        PARAMS.W_MAX           = arg_vector(7);
    end
    
    % Programming Parameters
    PERF_TEST = 1;      % Enable/disable performance output
    SANDBOX = 0;        % Controls whether "sandbox" area executes, or main func
    OPTIMIZATION_RUN = 0;
    if PERF_TEST; startTime = tic; end;
    
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
        r_groups = zeros(1, length(r_x_vals));
    elseif CONFIGURATION == FMRI
        loaded_input = load('/Users/gregashby/Desktop/Correlation_Analysis_Automaticity_Project/fMRI_data.mat');
        r_x_vals = loaded_input.r_x_mat;
        r_y_vals = loaded_input.r_y_mat;
        r_groups = zeros(1, length(r_x_vals));
    end
    
    % Struct to contain meta-data of FMRI configuration
    FMRI_META = struct('NUM_TRIALS', 11520, ...
                       'SES_1',      1:480, 'SES_4',    1681:2160, ...
                       'SES_10', 5161:5640, 'SES_20', 11041:11520);

    %% Initialize/configure constants (though some data structure specific constants are initialized below)
    % Set behavior and number of trials
    PRE_LEARNING_TRIALS  = PARAMS.PRE_LEARNING_TRIALS;                                   % Number of control trials run before learning trials
    LEARNING_TRIALS      = PARAMS.LEARNING_TRIALS;                                       % Number of learning trials in automaticity experiment
    POST_LEARNING_TRIALS = PARAMS.POST_LEARNING_TRIALS;                                  % Number of trials where no learning is involved after learning trials
    TRIALS               = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS; % Total number of trials
    % Create matrix to store information on when learning should occur
    LEARNING = [zeros(1, PRE_LEARNING_TRIALS), ...
                ones( 1, LEARNING_TRIALS), ...
                zeros(1, POST_LEARNING_TRIALS)]; 
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
    W_MAX = PARAMS.W_MAX;      % maximum possible weight for Hebbian Synapses
    INIT_PMC_WEIGHT = 0.08;    % Initial weight for PMC neurons
    NOISE = PARAMS.NOISE;      % Std. dev. of noise given to PFC/PMC v; set to 0 for no noise
    
    % Performance parameters
    loop_times = zeros(1, TRIALS); % Records how much time was needed for each loop
    trial_times = zeros(1, TRIALS);
    rt_calc_times = zeros(1, TRIALS);

    % Quantity of Visual Stimulus
    Visual = struct( ...
        'stim', 50 ...
    );

    % Using COVIS to select Rule
    correct_rule = 1;
    Rule_Matrix = [1 2 3 4];
    Salience_Matrix = ones(1,4);
    Rule_Weight_Matrix = ones(1,4);
    Rule_Probabilities_Matrix = ones(1,4);
    Probability_Space_Matrix = ones(1,3);        % This matrix saves the boundaries between the rule categories in the Rule Probability Matrix (used to select rule)
    
    % COVIS parameters
    
    Delta_C = 10;
    Delta_E = 1;
    Persev_Const = 5;
    lambda_COVIS = 1;
    n_guess = 5;               % number of random rule selections before COVIS model "kicks in" and begins selecting rule on each trial. 
    
    % Rule Log - Keeps track of which rule is used on each trial
    
    Rule_Log = ones(1,TRIALS);
    
    
    % Intialize Variables for Keeping Track of Total Activation in PFC, PMC, Head of Caudate, GPi, and Thalamus
    
    PFC_activations  = zeros(1, LEARNING_TRIALS);
    PMC_activations  = zeros(1, LEARNING_TRIALS);
    Caud_activations = zeros(1, LEARNING_TRIALS);
    GPI_activations  = zeros(1, LEARNING_TRIALS);
    Thal_activations = zeros(1, LEARNING_TRIALS);
    

    % Stimulus Rules
    RULE_1D_1   = struct('A', 1:GRID_SIZE/2, 'B', GRID_SIZE/2+1:GRID_SIZE, ...
                       'A_NUM_WEIGHTS', GRID_SIZE/2 * GRID_SIZE, 'B_NUM_WEIGHTS', GRID_SIZE/2 * GRID_SIZE);
    RULE_1D_2   = struct('B', 1:GRID_SIZE/2, 'A', GRID_SIZE/2+1:GRID_SIZE, ...
                       'A_NUM_WEIGHTS', GRID_SIZE/2 * GRID_SIZE, 'B_NUM_WEIGHTS', GRID_SIZE/2 * GRID_SIZE);

    % Radial Basis Function
    [X, Y] = meshgrid(1:GRID_SIZE, 1:GRID_SIZE);
    RBF = struct( ...
        'RADIUS', 0.8, ...
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
        'rx_matrix', zeros(TRIALS,3), ...            % Stores information about PMC neuron reacting during trial
        'alpha', zeros(TRIALS,n) ...
    );

    %% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses)
    % NMDA - upper threshold
    % AMPA - lower threshold
    % Strengthening occurs if integral_PMCAvoltage > Hebbian.NMDA
    % Weakening occurs if Hebbian.NMDA - integral_PMCAvoltage - Hebbian.AMPA > 0, i.e., only if integral_PMCAvoltage < Hebbian.NMDA + Hebbian.AMPA
    Hebbian = struct( ...
        'heb_coef', PARAMS.HEB_CONSTS, ...
        'anti_heb', PARAMS.ANTI_HEB_CONSTS, ...
        'NMDA',     PARAMS.NMDA, ...
        'AMPA',     PARAMS.AMPA ...
    );

    %% Neuron constants 
    %%(RSN: Regular Spiking Neuron)                set for a cortical regular spiking neuron
    %%(MSN: Medium Spiny Neuron)                   set for a medium spiny neuron in the caudate nucleus
    %%(QIAF: Quadratic Integrate and Fire Neuron)  set to simulate neurons in the Globus Pallidus
    
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

    MSN = struct( ...
        'C', 50, ...
        'rv', -80, ...
        'vt', -25, ...
        'k', 1, ...
        'a', 0.01, ...
        'b', -20, ...
        'c', -55, ...
        'd', 150, ...
        'vpeak', 40, ...
        'E', 100 ...
    );

    QIAF = struct( ...
        'beta', 11.83, ...
        'gamma', 0.117, ...
        'vt', -40, ...
        'rv', -60, ...
        'vpeak', 35, ...
        'vreset', -50 ...
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
        'W_OUT2MDN', 1, ...
        'W_OUT2AC', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0 ...
    );
    PFC_B = struct( ...
        'W_OUT', 9, ...
        'W_OUT2MDN', 1, ...
        'W_OUT2AC', 1, ...
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



%% Neurons that are part of the FROST component of the Model

    Driv_PFC = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...    % RSN.rv*ones(n,1)
        'u', zeros(n,1), ...
        'rule_stim', 0 ...
        );

    CN = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', MSN.rv*ones(n,1), ...
        'u', zeros(n,1) ...
        );
    
    GP = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', QIAF.rv*ones(n,1)...
        );
    
    MDN_A = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1) ...
        );
    
    MDN_B = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1) ...
        );
    
    AC_A = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1), ...
        'rule_stim', 0.1 ...
        );
    
    AC_B = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(n,1), ...
        'u', zeros(n,1), ...
        'rule_stim', 0.1 ...
        );


    %% Sandbox area
    % Placed after all values are initalized, and serves as an area where code can be prototyped and tested
    % (for validity or performance reasons) before being implemented into the main body of the function
    if PERF_TEST && SANDBOX; return; end;

    %% ============================ %%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %  ============================  %

    %% Pre-calculations (for performance reasons)
    % Calculate lambda values for individual trials
    t = (0:n)';
    LAMBDA_PRECALC = (t/LAMBDA).*exp((LAMBDA-t)/LAMBDA);
    
    %% Learning trials
    
    
    for j=1:TRIALS
        if PERF_TEST; loopStart = tic; end;
        
        % Select Rule for Trial

        random_number = rand;
        
        if j <= n_guess
            rule = randi(4);
        end   
       
        if j > n_guess && accuracy(j-1) == 1
            rule = Rule_Log(j-1);
        else
            if random_number <= Probability_Space_Matrix(1)
                rule = 1;
            elseif random_number <= Probability_Space_Matrix(2) && random_number > Probability_Space_Matrix(1)
                rule = 2;
            elseif random_number <= Probability_Space_Matrix(3) && random_number > Probability_Space_Matrix(2)
                rule = 3;
            else
                rule = 4;
            end
        end
             
        
    %  if rule == 1 || 3
    %      RULE = RULE_1D_1;
    %  else
    %      RULE = RULE_1D_2;
    %  end
        
        
        
        
        %% Initialize appropriate variables for each loop
        % variables tracking spiking rate in each neuron
        PFC_A.spikes = 0;       PMC_A.spikes = 0;
        PFC_B.spikes = 0;       PMC_B.spikes = 0;
        Driv_PFC.spikes = 0;
        CN.spikes = 0;
        GP.spikes = 0;
        MDN_A.spikes = 0;
        MDN_B.spikes = 0;
        AC_A.spikes = 0;
        AC_B.spikes = 0;

        % variables keep track of positive voltage values for calculation of
        % integral (for Hebbian learning equation)
        PFC_A.pos_volt(:) = 0;  PMC_A.pos_volt(:) = 0;
        PFC_B.pos_volt(:) = 0;  PMC_B.pos_volt(:) = 0;

        % Re-initialize Neuron Voltage Matrices (set all v to RSN.rv; all u to 0)
        PFC_A.v(:) = RSN.rv;     PFC_A.u(:) = 0;
        PFC_B.v(:) = RSN.rv;     PFC_B.u(:) = 0;
        PMC_A.v(:) = RSN.rv;     PMC_A.u(:) = 0;
        PMC_B.v(:) = RSN.rv;     PMC_B.u(:) = 0;
        Driv_PFC.v(:) = RSN.rv;  Driv_PFC.u(:) = 0;      %RSN.rv
        CN.v(:) = RSN.rv;        CN.u(:) = 0;
        GP.v(:) = RSN.rv;
        MDN_A.v(:) = RSN.rv;     MDN_A.u(:) = 0;
        MDN_B.v(:) = RSN.rv;     MDN_B.u(:) = 0;
        AC_A.v(:) = RSN.rv;      AC_A.u(:) = 0;
        AC_B.v(:) = RSN.rv;      AC_B.u(:) = 0;

        % Re-initialize Neuron Output Matrices
        PFC_A.out(:) = 0;        PMC_A.out(:) = 0;
        PFC_B.out(:) = 0;        PMC_B.out(:) = 0;
        Driv_PFC.out(:) = 0;
        CN.out(:) = 0;
        GP.out(:) = 0;
        MDN_A.out(:) = 0;
        MDN_B.out(:) = 0;
        AC_A.out(:) = 0;
        AC_B.out(:) = 0;
        
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
        
        %% Calculate which Rule to Use
        
        % First Trial, randomly select rule
        
        % Every trials after that
            % If previous trial was correct   (stick with same rule)
            % If previous trial was incorrect (select new rule)
        
        
        
        

        %% Radial Basis Function (RBF) Implementation
        
        % Calculate RBF grid
        RBF.rbv(:, :) = exp( -(sqrt((r_y-RBF.Y).^2 + (r_x-RBF.X).^2))/RBF.RADIUS ) * Visual.stim;

        if rule == 1
        
        % Sum appropriate RBF values to find PFC_A and PFC_B v_stim values
        PFC_A.v_stim = sum(reshape(RBF.rbv(:, RULE_1D_1.A), [RULE_1D_1.A_NUM_WEIGHTS  1]));                     
        PFC_B.v_stim = sum(reshape(RBF.rbv(:, RULE_1D_1.B), [RULE_1D_1.B_NUM_WEIGHTS  1]));
        
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_A_weights, [RBF.NUM_WEIGHTS 1]));
        PMC_B.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_B_weights, [RBF.NUM_WEIGHTS 1]));
        
        end
        
        
        if rule == 2
        
        % Sum appropriate RBF values to find PFC_A and PFC_B v_stim values
        PFC_A.v_stim = sum(reshape(RBF.rbv(:, RULE_1D_2.A), [RULE_1D_2.A_NUM_WEIGHTS  1]));                     
        PFC_B.v_stim = sum(reshape(RBF.rbv(:, RULE_1D_2.B), [RULE_1D_2.B_NUM_WEIGHTS  1]));
        
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_A_weights, [RBF.NUM_WEIGHTS 1]));
        PMC_B.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_B_weights, [RBF.NUM_WEIGHTS 1]));
        end
        
        if rule == 3
        
        % Sum appropriate RBF values to find PFC_A and PFC_B v_stim values
        PFC_A.v_stim = sum(reshape(RBF.rbv(RULE_1D_1.A, :), [RULE_1D_1.A_NUM_WEIGHTS  1]));
        PFC_B.v_stim = sum(reshape(RBF.rbv(RULE_1D_1.B, :), [RULE_1D_1.B_NUM_WEIGHTS  1]));
        
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_A_weights, [RBF.NUM_WEIGHTS 1]));
        PMC_B.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_B_weights, [RBF.NUM_WEIGHTS 1]));
        end
        
        if rule == 4
        
        % Sum appropriate RBF values to find PFC_A and PFC_B v_stim values
        PFC_A.v_stim = sum(reshape(RBF.rbv(RULE_1D_2.A, :), [RULE_1D_2.A_NUM_WEIGHTS  1]));
        PFC_B.v_stim = sum(reshape(RBF.rbv(RULE_1D_2.B, :), [RULE_1D_2.B_NUM_WEIGHTS  1]));
        
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_A_weights, [RBF.NUM_WEIGHTS 1]));
        PMC_B.v_stim = sum(reshape(RBF.rbv(:,:).*PMC_B_weights, [RBF.NUM_WEIGHTS 1]));
        end
        
        % Scale v_stim values to prevent them from becoming too large
        PFC_A.v_stim = PFC_A.v_stim * PFC.V_SCALE;
        PFC_B.v_stim = PFC_B.v_stim * PFC.V_SCALE;
        PMC_A.v_stim = PMC_A.v_stim * PMC.V_SCALE;
        PMC_B.v_stim = PMC_B.v_stim * PMC.V_SCALE;

        %% Individual Time Trial
        timeTrialStart = tic;
        for i=1:n-1
            % Neuron Equations
            % PFC A Neuron
            PFC_A.v(i+1)=(PFC_A.v(i) + TAU*(RSN.k*(PFC_A.v(i)-RSN.rv)*(PFC_A.v(i)-RSN.vt)-PFC_A.u(i)+ RSN.E + MDN_A.W_OUT*MDN_A.out(i) + AC_A.W_OUT*AC_A.out(i) + PFC_A.v_stim + (PMC_A.W_OUT*PMC_A.out(i)) - PFC.W_LI*PFC_B.out(i))/RSN.C) + normrnd(0,NOISE);
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
            PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + MDN_B.W_OUT*MDN_A.out(i) + AC_B.W_OUT*AC_A.out(i) + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out(i)) - PFC.W_LI*PFC_A.out(i))/RSN.C) + normrnd(0,NOISE);
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
            
            %% Add Neurons from FROST Model
            
            % Driv_PFC Neuron
                    % Input from Rule Stimulus (arbitrary value - constant)
                    % Output to CN Neuron
          
            Driv_PFC.v(i+1)=((Driv_PFC.v(i) + TAU*(RSN.k*(Driv_PFC.v(i)-RSN.rv)*(Driv_PFC.v(i)-RSN.vt)-Driv_PFC.u(i) + RSN.E + Driv_PFC.rule_stim))/RSN.C);
            Driv_PFC.u(i+1)=Driv_PFC.u(i)+TAU*RSN.a*(RSN.b*(Driv_PFC.v(i)-RSN.rv)-Driv_PFC.u(i));
            if Driv_PFC.v(i+1)>=RSN.vpeak
                Driv_PFC.v(i)= RSN.vpeak;
                Driv_PFC.v(i+1)= RSN.c;
                Driv_PFC.u(i+1)= Driv_PFC.u(i+1)+ RSN.d;
            end
            if Driv_PFC.v(i) >= RSN.vpeak
                Driv_PFC.out(i:n) = Driv_PFC.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
            
            % CN Neuron
                    % Input from Driv_PFC Neuron
                    % Output to GP Neuron
            
            CN.v(i+1)=((CN.v(i) + TAU*(MSN.k*(CN.v(i)-MSN.rv)*(CN.v(i)-MSN.vt)- CN.u(i) + Driv_PFC.W_OUT*Driv_PFC.out(i) + MSN.E ))/MSN.C); % + normrnd(0,NOISE);
            CN.u(i+1)= CN.u(i)+TAU*MSN.a*(MSN.b*(CN.v(i)-MSN.rv)-CN.u(i));
            if CN.v(i+1)>=MSN.vpeak
                CN.v(i)= MSN.vpeak;
                CN.v(i+1)= MSN.c;
                CN.u(i+1)= CN.u(i+1)+ MSN.d;
            end
            if CN.v(i) >= MSN.vpeak
                CN.out(i:n) = CN.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
            % GP Neuron
                    % Input from CN Neuron
                    % Output to MDN_A and MDN_B Neurons
            
                    %This part is less straightforward
    
            dGP = (-1)*CN.W_OUT*CN.out(i) + QIAF.beta + QIAF.gamma*(GP.v(i)- QIAF.rv)*(GP.v(i)-QIAF.vt);               
            GP.v(i+1) = GP.v(i) + dGP; % + normrnd(0,NOISE);
    
            if (GP.v(i+1) >= QIAF.vpeak);
                GP.v(i) = QIAF.vpeak;
                GP.v(i+1) = QIAF.vreset;
            end;
    
            if GP.v(i) >= QIAF.vpeak
                GP.out(i:n) = GP.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end

            % MDN_A Neuron
                    % Input from GP Neuron
                    % Input from pFC_A Neuron
                    % Output to pFC_A Neuron
            
            MDN_A.v(i+1)=((MDN_A.v(i) + TAU*(RSN.k*(MDN_A.v(i)-RSN.rv)*(MDN_A.v(i)-RSN.vt)-MDN_A.u(i)+ 10 + PFC_A.W_OUT2MDN*PFC_A.out(i) - GP.W_OUT*GP.out(i)))/RSN.C); % + normrnd(0,NOISE);
            MDN_A.u(i+1)=MDN_A.u(i)+TAU*RSN.a*(RSN.b*(MDN_A.v(i)-RSN.rv)-MDN_A.u(i));
            if MDN_A.v(i+1)>=RSN.vpeak
                MDN_A.v(i)= RSN.vpeak;
                MDN_A.v(i+1)= RSN.c;
                MDN_A.u(i+1)= MDN_A.u(i+1)+ RSN.d;
            end
            if MDN_A.v(i) >= RSN.vpeak
                MDN_A.out(i:n) = MDN_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
            
            % MDN_B Neuron
                    % Input from GP Neuron
                    % Input from pFC_A Neuron
                    % Output to pFC_A Neuron
                    
            MDN_B.v(i+1)=((MDN_B.v(i) + TAU*(RSN.k*(MDN_B.v(i)-RSN.rv)*(MDN_B.v(i)-RSN.vt)-MDN_B.u(i)+ 10 + PFC_B.W_OUT2MDN*PFC_B.out(i) - GP.W_OUT*GP.out(i)))/RSN.C); % + normrnd(0,NOISE);
            MDN_B.u(i+1)=MDN_B.u(i)+TAU*RSN.a*(RSN.b*(MDN_B.v(i)-RSN.rv)-MDN_B.u(i));
            if MDN_B.v(i+1)>=RSN.vpeak
                MDN_B.v(i)= RSN.vpeak;
                MDN_B.v(i+1)= RSN.c;
                MDN_B.u(i+1)= MDN_A.u(i+1)+ RSN.d;
            end
            
            if MDN_B.v(i) >= RSN.vpeak
                MDN_B.out(i:n) = MDN_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
            
            % AC_A Neuron
                    % Input from Rule Stimulus (arbitrary value - constant)
                    % Input from PFC_A Neuron
                    % Output to PFC_A
            
            AC_A.v(i+1)=((AC_A.v(i) + TAU*(RSN.k*(AC_A.v(i)-RSN.rv)*(AC_A.v(i)-RSN.vt)-AC_A.u(i)+ 10 + PFC_A.W_OUT2AC*PFC_A.out(i) + AC_A.rule_stim))/RSN.C); % + normrnd(0,NOISE);
            AC_A.u(i+1)=AC_A.u(i)+TAU*RSN.a*(RSN.b*(AC_A.v(i)-RSN.rv)-AC_A.u(i));
            if AC_A.v(i+1)>=RSN.vpeak
                AC_A.v(i)= RSN.vpeak;
                AC_A.v(i+1)= RSN.c;
                AC_A.u(i+1)= AC_A.u(i+1)+ RSN.d;
            end
            if AC_A.v(i) >= RSN.vpeak
                AC_A.out(i:n) = AC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
            
            % AC_B Neuron
                    % Input from Rule Stimulus (arbitrary value - constant)
                    % Input from PFC_B Neuron
                    % Output to PFC_B
            
            AC_B.v(i+1)=((AC_B.v(i) + TAU*(RSN.k*(AC_B.v(i)-RSN.rv)*(AC_B.v(i)-RSN.vt)-AC_B.u(i)+ 10 + PFC_B.W_OUT2AC*PFC_B.out(i) + AC_B.rule_stim))/RSN.C); % + normrnd(0,NOISE);
            AC_B.u(i+1)=AC_B.u(i)+TAU*RSN.a*(RSN.b*(AC_B.v(i)-RSN.rv)-AC_B.u(i));
            if AC_B.v(i+1)>=RSN.vpeak
                AC_B.v(i)= RSN.vpeak;
                AC_B.v(i+1)= RSN.c;
                AC_B.u(i+1)= AC_B.u(i+1)+ RSN.d;
            end
            if AC_B.v(i) >= RSN.vpeak
                AC_B.out(i:n) = AC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
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
        % Record "alpha" function, summing PMC A and PMC B output
        PMC.alpha(j,:) = PMC_A.out + PMC_B.out;
        trial_times(j) = toc(timeTrialStart);

        %% Determine decision neuron and reaction time, and record accuracy
        rt_start_time = tic;
        % Determine reacting neuron and latency
        [neuron_id_PFC, latency] = determine_reacting_neuron(PFC_A.out, PFC_B.out, PFC.DECISION_PT);
        PFC.rx_matrix(j,1:2) = [neuron_id_PFC, latency];
        PFC.rx_matrix(j,3) = r_group;
        [neuron_id_PMC, latency] = determine_reacting_neuron(PMC_A.out, PMC_B.out, PMC.DECISION_PT);
        PMC.rx_matrix(j,1:2) = [neuron_id_PMC, latency];
        PMC.rx_matrix(j,3) = r_group;
        % Determine accuracy
        
        if correct_rule == 1
            if (r_x <= 70 && neuron_id_PMC == 1) || (r_x >= 70 && neuron_id_PMC == 2)
            accuracy(j) = 1;
            else
            accuracy(j) = 0;
            end
        end
        
        if correct_rule == 2
            if (r_x >= 70 && neuron_id_PMC == 1) || (r_x <= 70 && neuron_id_PMC == 2)
            accuracy(j) = 1;
            else
            accuracy(j) = 0;
            end
        end
        
        if correct_rule == 3
            if (r_y <= 70 && neuron_id_PMC == 1) || (r_y >= 70 && neuron_id_PMC == 2)
            accuracy(j) = 1;
            else
            accuracy(j) = 0;
            end
        end
        
        if correct_rule == 4
            if (r_y >= 70 && neuron_id_PMC == 1) || (r_y <= 70 && neuron_id_PMC == 2)
            accuracy(j) = 1;
            else
            accuracy(j) = 0;
            end
        end
        
        rt_calc_times(j) = toc(rt_start_time);

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
        
        
        % Record Total Activations (per trials) for PFC, PMC, Head of Caudate, GPi, and Thalamus
        
        PFC_activations(j)  = trapz(PFC_A.out) + trapz(PFC_B.out);
        PMC_activations(j)  = trapz(PMC_A.out) + trapz(PMC_B.out);
        Caud_activations(j) = trapz(CN.out);
        GPI_activations(j)  = trapz(GP.out);
        Thal_activations(j) = trapz(MDN_A.out) + trapz(MDN_B.out);
        
        
        % Readjusting Saliences and Weights
            
            % If trial decision is correct
            
        if accuracy(j) == 1
            
            Salience_Matrix(rule) = Salience_Matrix(rule) + Delta_C;
            Rule_Weight_Matrix(rule) = Rule_Weight_Matrix(rule) + Persev_Const;
            
        end
        
            % If trial decision is incorrect
        
        if accuracy(j) == 0
            
            Salience_Matrix(rule) = Salience_Matrix(rule) + Delta_E;
            Rule_Weight_Matrix(rule) = Rule_Weight_Matrix(rule) + Persev_Const;
        end
        
        
            % Updating of Randomly Chosen Rule (Step Two in COVIS Weight Updating Algorithm)
        
        random_rule = randi(4);
        Rule_Weight_Matrix(random_rule) = Rule_Weight_Matrix(random_rule) + poissrnd(lambda_COVIS) ;
        
        
            % Calculate Rule Probabilities for Next Trial
      
       Total_Weights = 0;
            
       for q = 1:4
           Total_Weights = Total_Weights + Rule_Weight_Matrix(q) ;
       end
            
       for z = 1:4
           Rule_Probabilities_Matrix(z) = Rule_Weight_Matrix(z)/Total_Weights;
       end
       
       % Create Probability Space Matrix - This is used to select rule
        
        Probability_Space_Matrix(1) = Rule_Probabilities_Matrix(1);
        Probability_Space_Matrix(2) = Rule_Probabilities_Matrix(1)+Rule_Probabilities_Matrix(2);
        Probability_Space_Matrix(3) = Rule_Probabilities_Matrix(1)+Rule_Probabilities_Matrix(2)+Rule_Probabilities_Matrix(3);

        % Update Rule Log
    
        Rule_Log(j) = rule;
        
    end
    
        % Calculate Accuracy
        
        AveAccuracyS1  = mean(accuracy(FMRI_META.SES_1));
        AveAccuracyS4  = mean(accuracy(FMRI_META.SES_4));
        AveAccuracyS10 = mean(accuracy(FMRI_META.SES_10));
        AveAccuracyS20 = mean(accuracy(FMRI_META.SES_20));
    
        
        % Perform Convolutions
    
    % Set parameter values of hrf
    t1 = 1; n = 4; lamda = 2;
    
    % Define time axis
    seconds = LEARNING_TRIALS;
    t = 1:seconds;
    % Create hrf
    hrf = ((t-t1).^(n-1)).*exp(-(t-t1)/lamda)/((lamda^n)*factorial(n-1));
    
    % Compute convolution, dividing by 100 to set time unit at 0.01 seconds
    % Interpret stimvector in increments of 0.01 seconds
    
    Conv_PMC = conv(PMC_activations, hrf');
    Av_BOLD_Session1  = mean(Conv_PMC(FMRI_META.SES_1)) ;
    Av_BOLD_Session4  = mean(Conv_PMC(FMRI_META.SES_4)) ;
    Av_BOLD_Session10 = mean(Conv_PMC(FMRI_META.SES_10));
    Av_BOLD_Session20 = mean(Conv_PMC(FMRI_META.SES_20));
    
    x1 = Av_BOLD_Session1  ;
    x2 = Av_BOLD_Session4  ;
    x3 = Av_BOLD_Session10 ;
    x4 = Av_BOLD_Session20 ;
    x5 = AveAccuracyS1  ;
    x6 = AveAccuracyS4  ;
    x7 = AveAccuracyS10 ;
    x8 = AveAccuracyS20 ;
    
    % x2  = conv(PFC_activations, hrf');
    % x3  = conv(Caud_activations, hrf');
    % x4  = conv(GPI_activations, hrf');
    % x5  = conv(Thal_activations, hrf');
    
    
    
    
    Results = {x1,x2,x3,x4,x5,x6,x7,x8};
    
end
    
    