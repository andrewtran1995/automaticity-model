%AUTOMATICITYMODEL automaticity model for PMC neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Model of Automaticity in Rule-Based Learning        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
# Description
Run an automaticity model based upon a set of possible configurations:
 1. WALLIS
 2. MADDOX
 3. FMRI

# Programming Conventions
In general, variables that are written in all capital letters are meant
 to be constant values -- set once in the beginning of the program
 (variable initialization) and nowhere else.

Structures are used in the program by calling struct(...).
 These structures are meant to "bundle" or "group" together variables that
 have some kind of commonality, e.g., belonging to the same neuron, behavior,
 model, etc. The goal is to enforce readability by standardizing the names of grouped variables
 and making relationships between variables more apparent. Note that this
 has an effect on performance, but it should be negligible.

Note that the grid is set up column-major order, with points accessed as
 grid(y,x), with the origin situated at the top-left corner and axes
 increasing right and down for x and y, respectively.

# Tips/Tricks
If debugging, one can observe the workspace of the function by issuing the following
 command before execution: "dbstop if error".

# Function Signature
## Input Variables
arg_struct     - structure of n fields used to pass parameters that are
                 exposed in global optimization; if not specified in
                 non-codegen version, will be given default values based
                 on configuration; if called in codegen version, a full
                 structure must be provided
optional_parms - struct that may contain additional arguments for the
                 model, typically those that influence things outside of
                 how parameters are set
## Output Variables
opt_val_1      - return value signifying value of some cost function, used
                 for global optimization
opt_val_2      - return value signifying array used in FMRI group run

%}
function [opt_val_1, opt_val_2] = automaticityModel(arg_struct, optional_parms) %#codegen
    %% ======================================= %%
    %%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
    %  =======================================  %    
    % Code-generation declarations
    coder.extrinsic('getautoparams','displayautoresults');
    coder.extrinsic('tic','toc','struct2table');
    coder.varsize('chosen_rule');
    
    % Load configuration and config parameters
    MADDOX = 1; WALLIS = 2; FMRI = 3;
    CONFIGURATIONS = {'MADDOX', 'WALLIS', 'FMRI'};
    CONFIGURATION = FMRI;
    
    % Get parameters
    PARAMS = struct('PRE_LEARNING_TRIALS',0,'LEARNING_TRIALS',0,'POST_LEARNING_TRIALS',0,'PFC_DECISION_PT',0,'PMC_DECISION_PT',0,'MC_DECISION_PT',0,'HEB_CONSTS',0,'NMDA',0,'AMPA',0,'W_MAX',0,'NOISE_PFC',0,'NOISE_PMC',0,'NOISE_MC',0,'PMC_A_W_OUT',0,'PMC_B_W_OUT',0,'PFC_A_W_OUT_MDN',0,'PFC_B_W_OUT_MDN',0,'DRIV_PFC_W_OUT',0,'MDN_A_W_OUT',0,'MDN_B_W_OUT',0,'COVIS_DELTA_C',0,'COVIS_DELTA_E',0,'COVIS_PERSEV',0,'COVIS_LAMBDA',0);
    PARAMS = getautoparams(CONFIGURATIONS{CONFIGURATION});
    
    % Struct to contain meta-data of FMRI configuration
    FMRI_META = struct('NUM_TRIALS', 11520, 'GROUP_RUN', 0, ...
                       'SES_1',      1:480, 'SES_4',    1681:2160, ...
                       'SES_10', 5161:5640, 'SES_20', 11041:11520);

    % Model parameters (default values)
    VIS_INPUT_FROM_PARM   = 0;
    SUPPRESS_UI           = 0;
    OPTIMIZATION_CALC     = 0;
    FROST_ENABLED         = 1;
    COVIS_ENABLED         = 1;
    BUTTON_SWITCH_ENABLED = 1;
    PERF_TEST             = 0; % Enable/disable performance output
    
    % Override parameter values if they were specified as inputs
    % This should be turned into a function
    if nargin >= 1
        % Determine if parms are valid
        if not(areparamsvalid(PARAMS))
            disp(struct2table(PARAMS));
            error('Parameters not valid.');
        end
        if BUTTON_SWITCH_ENABLED && not(COVIS_ENABLED)
            error('Button switch cannot be enabled without COVIS');
        end
        PARAMS = absorbstruct(PARAMS, arg_struct);
    end

    % Override values with optional_parms if it was passed as an argument
    if nargin == 2
        if isfield(optional_parms, 'FMRI_META_GROUP_RUN')
            FMRI_META.GROUP_RUN = optional_parms.FMRI_META_GROUP_RUN;
        end
        if isfield(optional_parms, 'VIS_INPUT_FROM_PARM')
            VIS_INPUT_FROM_PARM = optional_parms.VIS_INPUT_FROM_PARM;
        end
    end
    
    %% Load visual stimulus matrix
    % %% Random Visual Input, 100 x 100 %%
    if 0
        loaded_input = load('datasets/randomVisualInput.mat');
        r_x_vals = loaded_input.r_x_mat;
        r_y_vals = loaded_input.r_y_mat;
    % %% Visual Input Matrix from optional_parms struct %%
    elseif VIS_INPUT_FROM_PARM
        r_x_vals = optional_parms.visualinput(:,1);
        r_y_vals = optional_parms.visualinput(:,2);
        r_groups = zeros(1, length(r_x_vals));
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
    % %% FMRI Visual Input, 100 X 100 %%
        loaded_input = load('datasets/fMRI_data.mat');
        r_x_vals = loaded_input.r_x_mat;
        r_y_vals = loaded_input.r_y_mat;
        r_groups = zeros(1, length(r_x_vals));
    end

    %% Initialize/configure constants (though some data structure specific constants are initialized below)
    % Set properties of grid
    STIM_GRID_SIZE = 100;                          % Length of side of square grid used for visual input; should be an even number
    BORDER_SIZE    = 20;                           % Width of border used to pad the grid such that visual stimulus on the edge still has an appropriate effect
    GRID_SIZE      = STIM_GRID_SIZE+2*BORDER_SIZE; % Total length of grid, i.e., the stimulus grid size and the border
    
    % Button switch must be initialized before TRIALS
    BUTTON_SWITCH = struct('TRIALS', 600, 'PMC_A_weights', ones(GRID_SIZE,GRID_SIZE,1,4), 'PMC_B_weights', ones(GRID_SIZE,GRID_SIZE,1,4));
    
    % Set behavior and number of trials
    PRE_LEARNING_TRIALS  = PARAMS.PRE_LEARNING_TRIALS;                                   % Number of control trials run before learning trials
    LEARNING_TRIALS      = PARAMS.LEARNING_TRIALS;                                       % Number of learning trials in automaticity experiment
    POST_LEARNING_TRIALS = PARAMS.POST_LEARNING_TRIALS;                                  % Number of trials where no learning is involved after learning trials
    if BUTTON_SWITCH_ENABLED
        TRIALS           = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS + BUTTON_SWITCH.TRIALS;
        LEARNING         = [zeros(1,PRE_LEARNING_TRIALS), ones(1,LEARNING_TRIALS), zeros(1,POST_LEARNING_TRIALS), ones(1,BUTTON_SWITCH.TRIALS)];
    else
        TRIALS           = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS;
        LEARNING         = [zeros(1,PRE_LEARNING_TRIALS), ones(1,LEARNING_TRIALS), zeros(1,POST_LEARNING_TRIALS)];
    end
    
    % Other parameters
    n = 1000;                     % Time period for one trial (in milliseconds)
    TAU = 1;                      % Tau
    LAMBDA = 20;                  % Lambda
    W_MAX = PARAMS.W_MAX;         % Maximum possible weight for Hebbian Synapses
    INIT_PMC_WEIGHT = 0.08;       % Initial weight for PMC neurons
    accuracy = zeros(TRIALS, 1);  % Boolean matrix indicating if correct PMC neuron reacted
    NOISE = struct('PFC', PARAMS.NOISE_PFC, 'PMC', PARAMS.NOISE_PMC, 'MC', PARAMS.NOISE_MC);
    
    % Performance parameters
    loop_times    = zeros(1, TRIALS); % Time needed for each loop (outer loop)
    trial_times   = zeros(1, TRIALS); % Time required to run each time loop (inner loop)
    rt_calc_times = zeros(1, TRIALS); % Time required to run reaction time calculation

    % Quantity of Visual Stimulus
    AREA = struct('LOWER_HALF', 1:GRID_SIZE/2, 'UPPER_HALF', GRID_SIZE/2+1:GRID_SIZE, ...
    			  'OUTER', [1:STIM_GRID_SIZE/4+BORDER_SIZE, STIM_GRID_SIZE*3/4+BORDER_SIZE+1:GRID_SIZE], ...
    			  'INNER', STIM_GRID_SIZE/4+BORDER_SIZE+1:STIM_GRID_SIZE*3/4+BORDER_SIZE, ...
    			  'ALL', 1:GRID_SIZE);
    VISUAL = struct('STIM', 50, ...
                    'r_x', 0, 'r_y', 0, 'r_group', 0, ...
                    'RULES', [ ...
                        struct('A_X', AREA.LOWER_HALF, 'A_Y', AREA.ALL,        'B_X', AREA.UPPER_HALF, 'B_Y', AREA.ALL,        'INVERSE', 2); ...
                        struct('A_X', AREA.UPPER_HALF, 'A_Y', AREA.ALL,        'B_X', AREA.LOWER_HALF, 'B_Y', AREA.ALL,        'INVERSE', 1); ...
                        struct('A_X', AREA.ALL,        'A_Y', AREA.LOWER_HALF, 'B_X', AREA.ALL,        'B_Y', AREA.UPPER_HALF, 'INVERSE', 4); ...
                        struct('A_X', AREA.ALL,        'A_Y', AREA.UPPER_HALF, 'B_X', AREA.ALL,        'B_Y', AREA.LOWER_HALF, 'INVERSE', 3); ...
                        struct('A_X', AREA.OUTER,      'A_Y', AREA.ALL,        'B_X', AREA.INNER,      'B_Y', AREA.ALL,        'INVERSE', 6); ...
                        struct('A_X', AREA.INNER,      'A_Y', AREA.ALL,        'B_X', AREA.OUTER,      'B_Y', AREA.ALL,        'INVERSE', 5); ...
                        struct('A_X', AREA.ALL,        'A_Y', AREA.OUTER,      'B_X', AREA.ALL,        'B_Y', AREA.INNER,      'INVERSE', 8); ...
                        struct('A_X', AREA.ALL,        'A_Y', AREA.INNER,      'B_X', AREA.ALL,        'B_Y', AREA.OUTER,      'INVERSE', 7) ...
                    ]);
    chosen_rule = 1;
    RULE = VISUAL.RULES(chosen_rule);

    % Radial Basis Function
    [X, Y] = meshgrid(1:GRID_SIZE, 1:GRID_SIZE);
    RBF = struct('RADIUS', 0.8, 'rbv', zeros(GRID_SIZE), 'X', X, 'Y', Y);

    %% COVIS Model
    COVIS_PARMS = struct('DELTA_C', PARAMS.COVIS_DELTA_C, 'DELTA_E', PARAMS.COVIS_DELTA_E, 'PERSEV', PARAMS.COVIS_PERSEV, ...
                         'LAMBDA', PARAMS.COVIS_LAMBDA, 'NUM_GUESS', 5, 'NUM_RULES', 4);
    COVIS_VARS = struct('correct_rule', VISUAL.RULES(2), 'rules', 1:COVIS_PARMS.NUM_RULES, 'saliences', ones(1,COVIS_PARMS.NUM_RULES), ...
                        'rule_weights', ones(1,COVIS_PARMS.NUM_RULES), 'rule_prob', ones(1,COVIS_PARMS.NUM_RULES), ...
                        'rule_log', ones(1,TRIALS));
    
    %% General settings for PFC, PMC neurons
    % Note that rx_matrix is big enough for both learning trials and no-learning trials to allow for comparisons
    % PFC general information
    PFC = struct( ...
        'V_SCALE', 1, ...                            % scaling factor for visual input into PFC neurons
        'W_LI', 2, ...                               % lateral inhibition between PFC A / PFC B
        'DECISION_PT', PARAMS.PFC_DECISION_PT, ...   % Integral value which determines which PFC neuron acts on a visual input
        'rx_matrix', zeros(TRIALS,3), ...            % Stores information about PFC neuron reacting during trial
        'activations', zeros(TRIALS,1) ...
    );

    % PMC general information
    PMC = struct( ...
        'V_SCALE', 1, ...                            % can use to scale PMC visual input value if it comes out way too high
        'W_LI', 2, ...                               % lateral inhibition between PMC A / PMC B
        'DECISION_PT', PARAMS.PMC_DECISION_PT, ...   % Integral value which determines which PMC neuron acts on a visual input
        'rx_matrix', zeros(TRIALS,3), ...            % Stores information about PMC neuron reacting during trial
        'alpha', zeros(TRIALS,n), ...
        'activations', zeros(TRIALS,1) ...
    );

    % MC general parameters
    MC = struct( ...
        'V_SCALE', 1, ...
        'W_LI', 2, ...
        'DECISION_PT', PARAMS.MC_DECISION_PT, ...
        'rx_matrix', zeros(TRIALS,3), ...
        'activations', zeros(TRIALS,1) ...
    );

    %% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses)
    % NMDA - upper threshold
    % AMPA - lower threshold
    % Strengthening occurs if integral_PMCAvoltage > Hebbian.NMDA
    % Weakening occurs if Hebbian.NMDA - integral_PMCAvoltage - Hebbian.AMPA > 0, i.e., only if integral_PMCAvoltage < Hebbian.NMDA + Hebbian.AMPA
    Hebbian = struct( ...
        'heb_coef', PARAMS.HEB_CONSTS, ...
        'anti_heb', PARAMS.HEB_CONSTS, ...
        'NMDA',     PARAMS.NMDA, ...
        'AMPA',     PARAMS.AMPA ...
    );

    %% Neuron constants
    % RSN (Regular Spiking Neuron): cortical regular spiking neuron
    RSN = struct('C', 100, 'rv', -60, 'vt', -40, 'k', 0.7, 'a', 0.03, 'b', -2,  'c', -50, 'd', 100, 'vpeak', 35, 'E', 60);

    % MSN (Regular Spiking Neuron): medium spiny neuron in caudate nucleus
    MSN = struct('C', 50,  'rv', -80, 'vt', -25, 'k', 1,   'a', 0.01, 'b', -20, 'c', -55, 'd', 150, 'vpeak', 40, 'E', 100);

    % QIAF (Quadratic Integrate and Fire Neuron): stimulate neurons in Globus Pallidus
    QIAF = struct( ...
        'beta', 11.83, ...
        'gamma', 0.117, ...
        'vt', -40, ...
        'rv', -60, ...
        'vpeak', 35, ...
        'vreset', -50 ...
    );

    %% Neuron Set-Up
    % # Initialize neuron structs to logically group variables
    % Many neurons are initialized as pairs, and are identical (except for
    % some exceptions, which are intialized manually)
    % # Commonly used struct fields
    % W_OUT  | weight of output from one neuron to another (e.g., PFC to PMC or PMC to PFC)
    % out    | output array
    % spikes | variables tracking spiking rate per trial
    % v      | voltage matrix (positive)
    % u      | voltage matrix (negative)
    
    % PFC - Primary Frontal Cortex
    PFC_A = struct( ...
        'W_OUT', 9, ...
        'W_OUT_MDN', PARAMS.PFC_A_W_OUT_MDN, ...
        'W_OUT_AC', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(RSN.rv,n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0 ...
    );
    PFC_B = PFC_A;
    PFC_B.W_OUT_MDN = PARAMS.PFC_B_W_OUT_MDN;

    % Create synpatic weight matrix
    % Use a simplified weights matrix for FMRI (performance reasons)
    if CONFIGURATION == FMRI
        if COVIS_ENABLED
            WEIGHTS_MATRIX = INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,1,4);
        else
            WEIGHTS_MATRIX = INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,1);
        end
    else
        if COVIS_ENABLED
            WEIGHTS_MATRIX = INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,TRIALS,4);
        else
            WEIGHTS_MATRIX = INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,TRIALS);
        end
    end
    
    % PMC - Premotor Cortex
    PMC_A = struct( ...
        'W_OUT', PARAMS.PMC_A_W_OUT, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(RSN.rv,n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0, ...
        'weights', WEIGHTS_MATRIX, ...
        'weights_avg', zeros(TRIALS,1) ...
    );
    PMC_B = PMC_A;
    PMC_B.W_OUT = PARAMS.PMC_B_W_OUT;

    % MC - Motor Cortex
    MC_A = struct( ...
        'W_OUT', 0, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(RSN.rv,n,1), ...
        'u', zeros(n,1), ...
        'pos_volt', zeros(n,1), ...
        'v_stim', 0 ...
    );
    MC_B = MC_A;

    %% FROST Model Neurons
    % Driving signal from PFC
    Driv_PFC = struct( ...
        'W_OUT', PARAMS.DRIV_PFC_W_OUT, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(RSN.rv,n,1), ...
        'u', zeros(n,1), ...
        'rule_stim', 0 ...
    );

    % Caudate nucleus
    CN = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(MSN.rv,n,1), ...
        'u', zeros(n,1), ...
        'activations', zeros(TRIALS,1) ...
    );

    % Globus pallidus
    GP = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(QIAF.rv,n,1),...
        'activations', zeros(TRIALS,1) ...
    );

    % MDN
    MDN = struct('activations', zeros(TRIALS,1));
    MDN_A = struct( ...
        'W_OUT', PARAMS.MDN_A_W_OUT, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(RSN.rv,n,1), ...
        'u', zeros(n,1) ...
    );
    MDN_B = MDN_A;
    MDN_B.W_OUT = PARAMS.MDN_B_W_OUT;

    % AC
    AC_A = struct( ...
        'W_OUT', 1, ...
        'out', zeros(n,1), ...
        'spikes', 0, ...
        'v', repmat(RSN.rv,n,1), ...
        'u', zeros(n,1), ...
        'rule_stim', 0.1 ...
    );
    AC_B = AC_A;

    %% ============================ %%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %  ============================  %
    start_time = tic;
    %% Pre-calculations (for performance reasons)
    % Calculate lambda values for individual trials
    t = (0:n)';
    LAMBDA_PRECALC = (t/LAMBDA).*exp((LAMBDA-t)/LAMBDA);
    
    %% Learning trials
    for j=1:TRIALS
        loopStart = tic;
        %% Initialize appropriate variables for each loop
        % Spiking rate of each neuron
        PFC_A.spikes = 0; PFC_B.spikes = 0; PMC_A.spikes = 0; PMC_B.spikes = 0; MC_A.spikes = 0; MC_B.spikes = 0;

        % Positive voltage values for integral (for Hebbian learning equation)
        PFC_A.pos_volt(:) = 0; PFC_B.pos_volt(:) = 0; PMC_A.pos_volt(:) = 0; PMC_B.pos_volt(:) = 0;

        % Neuron Voltage Matrices
        PFC_A.v(:) = RSN.rv; PFC_B.v(:) = RSN.rv; PMC_A.v(:) = RSN.rv; PMC_B.v(:) = RSN.rv; MC_A.v(:) = RSN.rv; MC_B.v(:) = RSN.rv;
        PFC_A.u(:) = 0;      PFC_B.u(:) = 0;      PMC_A.u(:) = 0;      PMC_B.u(:) = 0;      MC_A.u(:) = 0;      MC_B.u(:) = 0;

        % Neuron output vectors
        PFC_A.out(:) = 0; PMC_A.out(:) = 0; MC_A.out(:) = 0;
        PFC_B.out(:) = 0; PMC_B.out(:) = 0; MC_B.out(:) = 0;

        %% Initialize FROST components
        if FROST_ENABLED
        	Driv_PFC.spikes = 0; Driv_PFC.v(:) = RSN.rv;  Driv_PFC.u(:) = 0; Driv_PFC.out(:) = 0;
            CN.spikes       = 0; CN.v(:)       = RSN.rv;  CN.u(:)       = 0; CN.out(:)       = 0;
            GP.spikes       = 0; GP.v(:)       = RSN.rv;                     GP.out(:)       = 0;
            MDN_A.spikes    = 0; MDN_A.v(:)    = RSN.rv;  MDN_A.u(:)    = 0; MDN_A.out(:)    = 0;
            MDN_B.spikes    = 0; MDN_B.v(:)    = RSN.rv;  MDN_B.u(:)    = 0; MDN_B.out(:)    = 0;
            AC_A.spikes     = 0; AC_A.v(:)     = RSN.rv;  AC_A.u(:)     = 0; AC_A.out(:)     = 0;
            AC_B.spikes     = 0; AC_B.v(:)     = RSN.rv;  AC_B.u(:)     = 0; AC_B.out(:)     = 0;
        end
        %% Initialize COVIS components (choose a rule)
        if COVIS_ENABLED
            if j <= COVIS_PARMS.NUM_GUESS
                chosen_rule = randi(length(COVIS_VARS.rules));
            elseif accuracy(j-1) == 1
                chosen_rule = COVIS_VARS.rule_log(j-1);
            else
                chosen_rule = rand_discrete(COVIS_VARS.rule_prob);
            end
            RULE = VISUAL.RULES(chosen_rule);
        end
        %% Button Switch if enabled and correct trials
        if BUTTON_SWITCH_ENABLED && j == TRIALS - BUTTON_SWITCH.TRIALS + 1
            COVIS_VARS.correct_rule = VISUAL.RULES(RULE(1).INVERSE);
            if CONFIGURATION == FMRI
                BUTTON_SWITCH.PMC_A_weights(:,:,1,:) = PMC_A.weights(:,:,1,:);
                BUTTON_SWITCH.PMC_B_weights(:,:,1,:) = PMC_B.weights(:,:,1,:);
            else
                BUTTON_SWITCH.PMC_A_weights(:,:,1,:) = PMC_A.weights(:,:,j,:);
                BUTTON_SWITCH.PMC_B_weights(:,:,1,:) = PMC_B.weights(:,:,j,:);
            end
        end

        %% Set PMC weights, potentially dependent on COVIS
        % If first trial, set to initial weights
        if CONFIGURATION == FMRI || j==1
            if COVIS_ENABLED
                PMC_A_weights = PMC_A.weights(:,:,1,chosen_rule);
                PMC_B_weights = PMC_B.weights(:,:,1,chosen_rule);
            else
                PMC_A_weights = PMC_A.weights(:,:,1);
                PMC_B_weights = PMC_B.weights(:,:,1);
            end
        % Else, set weights to results of previous trial
        else
            if COVIS_ENABLED
                PMC_A_weights = PMC_A.weights(:,:,j-1,chosen_rule);
                PMC_B_weights = PMC_B.weights(:,:,j-1,chosen_rule);
            else
                PMC_A_weights = PMC_A.weights(:,:,j-1);
                PMC_B_weights = PMC_B.weights(:,:,j-1);
            end
        end

        %% Determine visual stimulus in range [1, GRID_SIZE] to pick random gabor for each trial, padded with
        % the BORDER_SIZE such that the visual stimulus is accounted for properly
        VISUAL.r_y = r_y_vals(j) + BORDER_SIZE;
        VISUAL.r_x = r_x_vals(j) + BORDER_SIZE;
        VISUAL.r_group = r_groups(j);

        %% Calculate visual stimulus effect using Radial Basis Function (RBF) implementation
        % Calculate RBF grid
        RBF.rbv(:,:) = exp( -(sqrt((VISUAL.r_y-RBF.Y).^2 + (VISUAL.r_x-RBF.X).^2))/RBF.RADIUS ) * VISUAL.STIM;
        % Sum RBF values depending on rule to find PFC_A and PFC_B v_stim values
        % Note that stim matrices are row-major order (e.g., indexed by y, then x)
        PFC_A.v_stim = sum(sum(RBF.rbv(RULE(1).A_Y, RULE(1).A_X)));
        PFC_B.v_stim = sum(sum(RBF.rbv(RULE(1).B_Y, RULE(1).B_X)));
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(sum(RBF.rbv(:,:).*PMC_A_weights));
        PMC_B.v_stim = sum(sum(RBF.rbv(:,:).*PMC_B_weights));
        % Scale v_stim values to prevent them from becoming too large
        PFC_A.v_stim = PFC_A.v_stim * PFC.V_SCALE;
        PFC_B.v_stim = PFC_B.v_stim * PFC.V_SCALE;
        PMC_A.v_stim = PMC_A.v_stim * PMC.V_SCALE;
        PMC_B.v_stim = PMC_B.v_stim * PMC.V_SCALE;

        %% Individual Time Trial Loop (iterating through n)
        timeTrialStart = tic;
        if FROST_ENABLED
            %% FROST Calculations
            for i=1:n-1
                % PFC A Neuron
                PFC_A.v(i+1)=(PFC_A.v(i) + TAU*(RSN.k*(PFC_A.v(i)-RSN.rv)*(PFC_A.v(i)-RSN.vt)-PFC_A.u(i)+ RSN.E + MDN_A.W_OUT*MDN_A.out(i) + AC_A.W_OUT*AC_A.out(i) + PFC_A.v_stim + (PMC_A.W_OUT*PMC_A.out(i)) - PFC.W_LI*PFC_B.out(i))/RSN.C) + normrnd(0,NOISE.PFC);
                PFC_A.u(i+1)=PFC_A.u(i)+TAU*RSN.a*(RSN.b*(PFC_A.v(i)-RSN.rv)-PFC_A.u(i));
                if PFC_A.v(i+1)>=RSN.vpeak
                    PFC_A.v(i)= RSN.vpeak;
                    PFC_A.v(i+1)= RSN.c;
                    PFC_A.u(i+1)= PFC_A.u(i+1)+ RSN.d;
                    PFC_A.out(i:n) = PFC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % PFC B Neuron
                PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + MDN_B.W_OUT*MDN_A.out(i) + AC_B.W_OUT*AC_A.out(i) + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out(i)) - PFC.W_LI*PFC_A.out(i))/RSN.C) + normrnd(0,NOISE.PFC);
                PFC_B.u(i+1)=PFC_B.u(i)+TAU*RSN.a*(RSN.b*(PFC_B.v(i)-RSN.rv)-PFC_B.u(i));
                if PFC_B.v(i+1)>=RSN.vpeak
                    PFC_B.v(i)= RSN.vpeak;
                    PFC_B.v(i+1)= RSN.c;
                    PFC_B.u(i+1)= PFC_B.u(i+1)+ RSN.d;
                    PFC_B.out(i:n) = PFC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % PMC_A Neuron
                PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + PMC_A.v_stim + (PFC_A.W_OUT*PFC_A.out(i)) - PMC.W_LI*PMC_B.out(i) )/RSN.C) + normrnd(0,NOISE.PMC);
                PMC_A.u(i+1)=PMC_A.u(i)+TAU*RSN.a*(RSN.b*(PMC_A.v(i)-RSN.rv)-PMC_A.u(i));
                if PMC_A.v(i+1)>=RSN.vpeak
                    PMC_A.v(i)= RSN.vpeak;
                    PMC_A.v(i+1)= RSN.c;
                    PMC_A.u(i+1)= PMC_A.u(i+1)+ RSN.d;
                    PMC_A.out(i:n) = PMC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % PMC_B Neuron
                PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + PMC_B.v_stim + (PFC_B.W_OUT*PFC_B.out(i)) - PMC.W_LI*PMC_A.out(i) )/RSN.C) + normrnd(0,NOISE.PMC);
                PMC_B.u(i+1)=PMC_B.u(i)+TAU*RSN.a*(RSN.b*(PMC_B.v(i)-RSN.rv)-PMC_B.u(i));
                if PMC_B.v(i+1)>=RSN.vpeak
                    PMC_B.v(i)= RSN.vpeak;
                    PMC_B.v(i+1)= RSN.c;
                    PMC_B.u(i+1)= PMC_B.u(i+1)+ RSN.d;
                    PMC_B.out(i:n) = PMC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end
                
                % MC_A Neuron
                MC_A.v(i+1)=(MC_A.v(i) + TAU*(RSN.k*(MC_A.v(i)-RSN.rv)*(MC_A.v(i)-RSN.vt)-MC_A.u(i)+ RSN.E + (PMC_A.W_OUT*PMC_A.out(i)) - MC.W_LI*MC_B.out(i) )/RSN.C) + normrnd(0,NOISE.MC);
                MC_A.u(i+1)=MC_A.u(i)+TAU*RSN.a*(RSN.b*(MC_A.v(i)-RSN.rv)-MC_A.u(i));
                if MC_A.v(i+1)>=RSN.vpeak
                    MC_A.v(i)= RSN.vpeak;
                    MC_A.v(i+1)= RSN.c;
                    MC_A.u(i+1)= MC_A.u(i+1)+ RSN.d;
                    MC_A.out(i:n) = MC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end
                
                % MC_B Neuron
                MC_B.v(i+1)=(MC_B.v(i) + TAU*(RSN.k*(MC_B.v(i)-RSN.rv)*(MC_B.v(i)-RSN.vt)-MC_B.u(i)+ RSN.E + (PMC_B.W_OUT*PMC_B.out(i)) - MC.W_LI*MC_A.out(i) )/RSN.C) + normrnd(0,NOISE.MC);
                MC_B.u(i+1)=MC_B.u(i)+TAU*RSN.a*(RSN.b*(MC_B.v(i)-RSN.rv)-MC_B.u(i));
                if MC_B.v(i+1)>=RSN.vpeak
                    MC_B.v(i)= RSN.vpeak;
                    MC_B.v(i+1)= RSN.c;
                    MC_B.u(i+1)= MC_B.u(i+1)+ RSN.d;
                    MC_B.out(i:n) = MC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end                

                % Driv_PFC Neuron
                    % Input from Rule Stimulus (arbitrary value - constant)
                    % Output to CN Neuron
                Driv_PFC.v(i+1)=((Driv_PFC.v(i) + TAU*(RSN.k*(Driv_PFC.v(i)-RSN.rv)*(Driv_PFC.v(i)-RSN.vt)-Driv_PFC.u(i) + RSN.E + Driv_PFC.rule_stim))/RSN.C);
                Driv_PFC.u(i+1)=Driv_PFC.u(i)+TAU*RSN.a*(RSN.b*(Driv_PFC.v(i)-RSN.rv)-Driv_PFC.u(i));
                if Driv_PFC.v(i+1)>=RSN.vpeak
                    Driv_PFC.v(i)= RSN.vpeak;
                    Driv_PFC.v(i+1)= RSN.c;
                    Driv_PFC.u(i+1)= Driv_PFC.u(i+1)+ RSN.d;
                    Driv_PFC.out(i:n) = Driv_PFC.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % CN Neuron
                    % Input from Driv_PFC Neuron
                    % Output to GP Neuron
                CN.v(i+1)=((CN.v(i) + TAU*(MSN.k*(CN.v(i)-MSN.rv)*(CN.v(i)-MSN.vt)- CN.u(i) + Driv_PFC.W_OUT*Driv_PFC.out(i) + MSN.E ))/MSN.C);
                CN.u(i+1)= CN.u(i)+TAU*MSN.a*(MSN.b*(CN.v(i)-MSN.rv)-CN.u(i));
                if CN.v(i+1)>=MSN.vpeak
                    CN.v(i)= MSN.vpeak;
                    CN.v(i+1)= MSN.c;
                    CN.u(i+1)= CN.u(i+1)+ MSN.d;
                    CN.out(i:n) = CN.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end
                
                % GP Neuron
                    % Input from CN Neuron
                    % Output to MDN_A and MDN_B Neurons
                    % This part is less straightforward
                dGP = (-1)*CN.W_OUT*CN.out(i) + QIAF.beta + QIAF.gamma*(GP.v(i)- QIAF.rv)*(GP.v(i)-QIAF.vt);               
                GP.v(i+1) = GP.v(i) + dGP;
                if (GP.v(i+1) >= QIAF.vpeak)
                    GP.v(i) = QIAF.vpeak;
                    GP.v(i+1) = QIAF.vreset;
                    GP.out(i:n) = GP.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end;

                % MDN_A Neuron
                    % Input from GP Neuron
                    % Input from pFC_A Neuron
                    % Output to pFC_A Neuron
                MDN_A.v(i+1)=((MDN_A.v(i) + TAU*(RSN.k*(MDN_A.v(i)-RSN.rv)*(MDN_A.v(i)-RSN.vt)-MDN_A.u(i)+ 10 + PFC_A.W_OUT_MDN*PFC_A.out(i) - GP.W_OUT*GP.out(i)))/RSN.C);
                MDN_A.u(i+1)=MDN_A.u(i)+TAU*RSN.a*(RSN.b*(MDN_A.v(i)-RSN.rv)-MDN_A.u(i));
                if MDN_A.v(i+1)>=RSN.vpeak
                    MDN_A.v(i)= RSN.vpeak;
                    MDN_A.v(i+1)= RSN.c;
                    MDN_A.u(i+1)= MDN_A.u(i+1)+ RSN.d;
                    MDN_A.out(i:n) = MDN_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % MDN_B Neuron
                    % Input from GP Neuron
                    % Input from pFC_A Neuron
                    % Output to pFC_A Neuron
                MDN_B.v(i+1)=((MDN_B.v(i) + TAU*(RSN.k*(MDN_B.v(i)-RSN.rv)*(MDN_B.v(i)-RSN.vt)-MDN_B.u(i)+ 10 + PFC_B.W_OUT_MDN*PFC_B.out(i) - GP.W_OUT*GP.out(i)))/RSN.C);
                MDN_B.u(i+1)=MDN_B.u(i)+TAU*RSN.a*(RSN.b*(MDN_B.v(i)-RSN.rv)-MDN_B.u(i));
                if MDN_B.v(i+1)>=RSN.vpeak
                    MDN_B.v(i)= RSN.vpeak;
                    MDN_B.v(i+1)= RSN.c;
                    MDN_B.u(i+1)= MDN_A.u(i+1)+ RSN.d;
                    MDN_B.out(i:n) = MDN_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % AC_A Neuron
                	% Input from Rule Stimulus (arbitrary value - constant)
                	% Input from PFC_A Neuron
                    % Output to PFC_A
                AC_A.v(i+1)=((AC_A.v(i) + TAU*(RSN.k*(AC_A.v(i)-RSN.rv)*(AC_A.v(i)-RSN.vt)-AC_A.u(i)+ 10 + PFC_A.W_OUT_AC*PFC_A.out(i) + AC_A.rule_stim))/RSN.C);
                AC_A.u(i+1)=AC_A.u(i)+TAU*RSN.a*(RSN.b*(AC_A.v(i)-RSN.rv)-AC_A.u(i));
                if AC_A.v(i+1)>=RSN.vpeak
                    AC_A.v(i)= RSN.vpeak;
                    AC_A.v(i+1)= RSN.c;
                    AC_A.u(i+1)= AC_A.u(i+1)+ RSN.d;
                    AC_A.out(i:n) = AC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % AC_B Neuron
                    % Input from Rule Stimulus (arbitrary value - constant)
                    % Input from PFC_B Neuron
                    % Output to PFC_B
                AC_B.v(i+1)=((AC_B.v(i) + TAU*(RSN.k*(AC_B.v(i)-RSN.rv)*(AC_B.v(i)-RSN.vt)-AC_B.u(i)+ 10 + PFC_B.W_OUT_AC*PFC_B.out(i) + AC_B.rule_stim))/RSN.C);
                AC_B.u(i+1)=AC_B.u(i)+TAU*RSN.a*(RSN.b*(AC_B.v(i)-RSN.rv)-AC_B.u(i));
                if AC_B.v(i+1)>=RSN.vpeak
                    AC_B.v(i)= RSN.vpeak;
                    AC_B.v(i+1)= RSN.c;
                    AC_B.u(i+1)= AC_B.u(i+1)+ RSN.d;
                    AC_B.out(i:n) = AC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end
            end
        else
            %% Non-FROST Calculation
            for i=1:n-1
                % PFC A Neuron
                PFC_A.v(i+1)=(PFC_A.v(i) + TAU*(RSN.k*(PFC_A.v(i)-RSN.rv)*(PFC_A.v(i)-RSN.vt)-PFC_A.u(i)+ RSN.E + PFC_A.v_stim + (PMC_A.W_OUT*PMC_A.out(i)) - PFC.W_LI*PFC_B.out(i))/RSN.C) + normrnd(0,NOISE.PFC);
                PFC_A.u(i+1)=PFC_A.u(i)+TAU*RSN.a*(RSN.b*(PFC_A.v(i)-RSN.rv)-PFC_A.u(i));
                if PFC_A.v(i+1)>=RSN.vpeak
                    PFC_A.v(i)= RSN.vpeak;
                    PFC_A.v(i+1)= RSN.c;
                    PFC_A.u(i+1)= PFC_A.u(i+1)+ RSN.d;
                    PFC_A.out(i:n) = PFC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % PFC B Neuron
                PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out(i)) - PFC.W_LI*PFC_A.out(i))/RSN.C) + normrnd(0,NOISE.PFC);
                PFC_B.u(i+1)=PFC_B.u(i)+TAU*RSN.a*(RSN.b*(PFC_B.v(i)-RSN.rv)-PFC_B.u(i));
                if PFC_B.v(i+1)>=RSN.vpeak
                    PFC_B.v(i)= RSN.vpeak;
                    PFC_B.v(i+1)= RSN.c;
                    PFC_B.u(i+1)= PFC_B.u(i+1)+ RSN.d;
                    PFC_B.out(i:n) = PFC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % PMC_A Neuron
                PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + PMC_A.v_stim + (PFC_A.W_OUT*PFC_A.out(i)) - PMC.W_LI*PMC_B.out(i) )/RSN.C) + normrnd(0,NOISE.PMC);
                PMC_A.u(i+1)=PMC_A.u(i)+TAU*RSN.a*(RSN.b*(PMC_A.v(i)-RSN.rv)-PMC_A.u(i));
                if PMC_A.v(i+1)>=RSN.vpeak
                    PMC_A.v(i)= RSN.vpeak;
                    PMC_A.v(i+1)= RSN.c;
                    PMC_A.u(i+1)= PMC_A.u(i+1)+ RSN.d;
                    PMC_A.out(i:n) = PMC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end

                % PMC_B Neuron
                PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + PMC_B.v_stim + (PFC_B.W_OUT*PFC_B.out(i)) - PMC.W_LI*PMC_A.out(i) )/RSN.C) + normrnd(0,NOISE.PMC);
                PMC_B.u(i+1)=PMC_B.u(i)+TAU*RSN.a*(RSN.b*(PMC_B.v(i)-RSN.rv)-PMC_B.u(i));
                if PMC_B.v(i+1)>=RSN.vpeak
                    PMC_B.v(i)= RSN.vpeak;
                    PMC_B.v(i+1)= RSN.c;
                    PMC_B.u(i+1)= PMC_B.u(i+1)+ RSN.d;
                    PMC_B.out(i:n) = PMC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end
                
                % MC_A Neuron
                MC_A.v(i+1)=(MC_A.v(i) + TAU*(RSN.k*(MC_A.v(i)-RSN.rv)*(MC_A.v(i)-RSN.vt)-MC_A.u(i)+ RSN.E + (PMC_A.W_OUT*PMC_A.out(i)) - MC.W_LI*MC_B.out(i) )/RSN.C) + normrnd(0,NOISE.MC);
                MC_A.u(i+1)=MC_A.u(i)+TAU*RSN.a*(RSN.b*(MC_A.v(i)-RSN.rv)-MC_A.u(i));
                if MC_A.v(i+1)>=RSN.vpeak
                    MC_A.v(i)= RSN.vpeak;
                    MC_A.v(i+1)= RSN.c;
                    MC_A.u(i+1)= MC_A.u(i+1)+ RSN.d;
                    MC_A.out(i:n) = MC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end
                
                % MC_B Neuron
                MC_B.v(i+1)=(MC_B.v(i) + TAU*(RSN.k*(MC_B.v(i)-RSN.rv)*(MC_B.v(i)-RSN.vt)-MC_B.u(i)+ RSN.E + (PMC_B.W_OUT*PMC_B.out(i)) - MC.W_LI*MC_A.out(i) )/RSN.C) + normrnd(0,NOISE.MC);
                MC_B.u(i+1)=MC_B.u(i)+TAU*RSN.a*(RSN.b*(MC_B.v(i)-RSN.rv)-MC_B.u(i));
                if MC_B.v(i+1)>=RSN.vpeak
                    MC_B.v(i)= RSN.vpeak;
                    MC_B.v(i+1)= RSN.c;
                    MC_B.u(i+1)= MC_B.u(i+1)+ RSN.d;
                    MC_B.out(i:n) = MC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
                end
            end
        end
        %% Record post-time-loop numbers
        % Count number of spikes
        PFC_A.spikes = nnz(PFC_A.v >= RSN.vpeak); PFC_B.spikes = nnz(PFC_B.v >= RSN.vpeak);
        PMC_A.spikes = nnz(PMC_A.v >= RSN.vpeak); PMC_B.spikes = nnz(PMC_B.v >= RSN.vpeak);
        % Record voltage value if positive. Else, do nothing.
        % For computation of integral
        PFC_A.pos_volt(PFC_A.v > 0) = PFC_A.v(PFC_A.v > 0);
        PFC_B.pos_volt(PFC_B.v > 0) = PFC_B.v(PFC_B.v > 0);
        PMC_A.pos_volt(PMC_A.v > 0) = PMC_A.v(PMC_A.v > 0);
        PMC_B.pos_volt(PMC_B.v > 0) = PMC_B.v(PMC_B.v > 0);
        % Record "alpha" function, summing PMC A and PMC B output
        PMC.alpha(j,:) = PMC_A.out + PMC_B.out;
        % Record total neuron activations
        PFC.activations(j) = trapz(PFC_A.out + PFC_B.out);
        CN.activations(j) = trapz(CN.out);
        GP.activations(j) = trapz(GP.out);
        MDN.activations(j) = trapz(MDN_A.out + MDN_B.out);
        PMC.activations(j) = trapz(PMC.alpha(j,:));
        MC.activations(j) = trapz(MC_A.out + MC_B.out);
        trial_times(j) = toc(timeTrialStart);

        %% Determine decision neuron and reaction time, and record accuracy
        rt_start_time = tic;
        % Determine reacting neuron and latency
        [neuron_id_PFC, latency] = determine_reacting_neuron(PFC_A.out, PFC_B.out, PFC.DECISION_PT);
        PFC.rx_matrix(j,:) = [neuron_id_PFC, latency, VISUAL.r_group];
        [neuron_id_PMC, latency] = determine_reacting_neuron(PMC_A.out, PMC_B.out, PMC.DECISION_PT);
        PMC.rx_matrix(j,:) = [neuron_id_PMC, latency, VISUAL.r_group];
        [neuron_id_MC, latency] = determine_reacting_neuron(MC_A.out, MC_B.out, MC.DECISION_PT);
        MC.rx_matrix(j,:) = [neuron_id_MC, latency, VISUAL.r_group];
        % Determine accuracy
        if COVIS_ENABLED
            accuracy(j) = double((any(VISUAL.r_x == COVIS_VARS.correct_rule.B_X) && any(VISUAL.r_y == COVIS_VARS.correct_rule.B_Y)) + 1) == neuron_id_MC;
        else
            accuracy(j) = double((any(VISUAL.r_x == RULE.B_X) && any(VISUAL.r_y == RULE.B_Y)) + 1) == neuron_id_MC;
        end
        rt_calc_times(j) = toc(rt_start_time);

        %% Weight change calculations
        if CONFIGURATION == FMRI
            k = 1;
        else
            k = j;
        end
        if COVIS_ENABLED
            el = chosen_rule;
        else
            el = 1;
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
            PMC_A_weights(:,:,1,1) = PMC_A_weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputA)*g_t_1_A.*(W_MAX - PMC_A_weights) - (Hebbian.anti_heb)*(integral_visinputA)*g_t_2_A.*PMC_A_weights);
            PMC_A.weights(:,:,k,el) = PMC_A_weights;

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_A.weights(:,:,k,el) = max(PMC_A.weights(:,:,k,el), 0);
            PMC_A.weights(:,:,k,el) = min(PMC_A.weights(:,:,k,el), W_MAX);

            %% Calculation of Hebbian Weight for PMC_B
            % Visual input to PMC_B neuron (presynaptic)
            integral_visinputB   = trapz(PFC_B.pos_volt);
            % Activation of PMC_B neuron   (post-synaptic)
            integral_PMCBvoltage = trapz(PMC_B.pos_volt);

            % Ensures g(t)-1 and g(2)-2 are never less than zero
            g_t_1_B = max(0, integral_PMCBvoltage - Hebbian.NMDA);
            g_t_2_B = max(0, Hebbian.NMDA - integral_PMCBvoltage - Hebbian.AMPA);

            % Determine new weights of visual PMC_B synapses
            PMC_B_weights(:,:,1,1) = PMC_B_weights + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputB)*g_t_1_B.*(W_MAX - PMC_B_weights) - (Hebbian.anti_heb)*(integral_visinputB)*g_t_2_B.*PMC_B_weights);
            PMC_B.weights(:,:,k,el) = PMC_B_weights;

            % Limit values of PMC_A.weights to be in range [0,W_MAX]
            PMC_B.weights(:,:,k,el) = max(PMC_B.weights(:,:,k,el), 0);
            PMC_B.weights(:,:,k,el) = min(PMC_B.weights(:,:,k,el), W_MAX); 
        % Else, if not learning, set new weights to previous weights
        else
            PMC_A.weights(:,:,k,el) = PMC_A_weights;
            PMC_B.weights(:,:,k,el) = PMC_B_weights;
        end
        
        % If COVIS is enabled and weight matrix has time dimension, update
        % all other weight matrices for this iteration
        if COVIS_ENABLED && CONFIGURATION ~= FMRI
            PMC_A.weights(:,:,j,1:COVIS_PARMS.NUM_RULES ~= chosen_rule) = PMC_A.weights(:,:,j-1,1:COVIS_PARMS.NUM_RULES ~= chosen_rule);
            PMC_B.weights(:,:,j,1:COVIS_PARMS.NUM_RULES ~= chosen_rule) = PMC_B.weights(:,:,j-1,1:COVIS_PARMS.NUM_RULES ~= chosen_rule);
        end
        
        % Record average weight for PMC_A and PMC_B
        PMC_A.weights_avg(j) = mean(mean(PMC_A.weights(:,:,k,el)));
        PMC_B.weights_avg(j) = mean(mean(PMC_B.weights(:,:,k,el)));
        
        %% COVIS Calculations - readjusting saliences, weights
        if COVIS_ENABLED
            % Step 1: Readjust saliences & weights
            if accuracy(j) == 1
                COVIS_VARS.saliences(chosen_rule) = COVIS_VARS.saliences(chosen_rule) + COVIS_PARMS.DELTA_C;
            else
                COVIS_VARS.saliences(chosen_rule) = COVIS_VARS.saliences(chosen_rule) + COVIS_PARMS.DELTA_E;
            end
            COVIS_VARS.rule_weights(chosen_rule) = COVIS_VARS.rule_weights(chosen_rule) + COVIS_PARMS.PERSEV;
            
            % Step 2: updating of randomly chosen rule
            random_rule = randi(4);
            COVIS_VARS.rule_weights(random_rule) = COVIS_VARS.rule_weights(random_rule) + poissrnd(COVIS_PARMS.LAMBDA);

            % Step 3: calculate rule probabilities for next trial
            COVIS_VARS.rule_prob = COVIS_VARS.rule_weights./sum(COVIS_VARS.rule_weights);

            COVIS_VARS.rule_log(j) = chosen_rule;
        end

        %% Print data to console
        if not(SUPPRESS_UI)
            consoleprogressbar('TRIALS COMPLETED', j, TRIALS);
        end
        loop_times(j) = toc(loopStart);
    end
    
    %% ========================================= %%
    %%%%%%%%%% OPTIMIZATION CALCULATIONS %%%%%%%%%%
    %  =========================================  %
    opt_val_1 = 0;
    opt_val_2 = zeros(4,4);
    % Calculate Sum of Squared Errors of Prediction (SSE)
    if OPTIMIZATION_CALC
        if CONFIGURATION == MADDOX
            opt_val_1 = 0;
        elseif CONFIGURATION == WALLIS
            opt_val_1 = 0;
        elseif CONFIGURATION == FMRI
            if ~FMRI_META.GROUP_RUN
                target = load('fmri/targetMeans1dCondition.mat');
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
                % Weight reaction time greater than accuracy
                target_diff = [target.means1dCondition(1,:) - output_acc;
                               (target.means1dCondition(2,:) - norm_output_rt)*20];
                opt_val_1 = sum(sum(target_diff.^2));
            else
                % Set parameter values of hrf
                t1 = 1; n = 4; lamda = 2;
                % Define time axis
                t = 1:LEARNING_TRIALS;
                % Create hrf
                hrf = ((t-t1).^(n-1)).*exp(-(t-t1)/lamda)/((lamda^n)*factorial(n-1));
                % Get value for opt_val_2
                opt_val_2 = [get_FMRI_corr_data(CN.activations, FMRI_META, hrf); ...
                             get_FMRI_corr_data(MDN.activations, FMRI_META, hrf); ...
                             get_FMRI_corr_data(PMC.activations, FMRI_META, hrf); ...
                             mean(accuracy(FMRI_META.SES_1)), mean(accuracy(FMRI_META.SES_4)), mean(accuracy(FMRI_META.SES_10)), mean(accuracy(FMRI_META.SES_20))];
            end
        end
    end
    %% =============================== %%
    %%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%
    %  ===============================  %
    if not(SUPPRESS_UI)
        displayautoresults(FROST_ENABLED, COVIS_ENABLED, BUTTON_SWITCH_ENABLED, BUTTON_SWITCH, COVIS_VARS, FMRI_META, CONFIGURATION, MADDOX, WALLIS, FMRI, TAU, n, RBF, BORDER_SIZE, VISUAL, TRIALS, PRE_LEARNING_TRIALS, LEARNING_TRIALS, POST_LEARNING_TRIALS, accuracy, PFC, PMC, PFC_A, PFC_B, PMC_A, PMC_B, Driv_PFC, CN, GP, MDN_A, MDN_B, AC_A, AC_B, PERF_TEST, start_time, loop_times, trial_times, rt_calc_times, chosen_rule);
    end
    return;
end

%% =============================== %%
%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%
%  ===============================  %
% Return what neuron reacts to the stimuli, and the latency
% Returns neuron_id = 0 for n1, neuron_id = 1 for n2
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
    elseif n1_latency < n2_latency
        neuron_id = 1;
        latency = n1_latency;
    % If latencies are equal (decision point never reached), take the
    % higher integral as the reacting neuron
    else
        neuron_id = double(trapz(n1) < trapz(n2));
        latency = length(n1);
    end
end

% Given discrete distribution, return index of chosen index
function [idx] = rand_discrete(distr)
    cum_distr = cumsum(distr);
    idx = find(rand<cum_distr, 1);
end

% Finds correlation between different neurons and accuracy
function [corr_vec] = get_FMRI_corr_data(activations, FMRI_META, hrf)
    % Compute convolution for each trial
    boldPMC = conv(activations, hrf');
    corr_vec = [mean(boldPMC(FMRI_META.SES_1)), ...
                mean(boldPMC(FMRI_META.SES_4)), ...
                mean(boldPMC(FMRI_META.SES_10)), ...
                mean(boldPMC(FMRI_META.SES_20))];
end

% Console progress bar
% Resources used:
% https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar
% https://www.mathworks.com/matlabcentral/fileexchange/30297-consoleprogressbar
function consoleprogressbar(str, iter, total)
    % Constants and variable initialization
    coder.extrinsic('num2str');
    coder.varsize('STR_CR');
    BAR_LENGTH = 40;
    persistent STR_CR_COUNT;
    if isempty(STR_CR_COUNT)
        STR_CR_COUNT = 0;
    end
    
    % Create progress bar
    percentage = iter/total;
    fullbar = repmat('.', [1, BAR_LENGTH]);
    completedbar = round(percentage * BAR_LENGTH);
    fullbar(1:completedbar) = '#';
    fullbar = ['[', fullbar, ']'];
    
    strout = [fullbar, ' ', num2str(iter), '/', num2str(total), ' ', str];
    
    % Print bar
    if iter == 1
        fprintf(strout);
    else
        fprintf([repmat('\b', 1, STR_CR_COUNT), strout]);
    end
    
    % Update carriage return
    STR_CR_COUNT = length(strout);
end