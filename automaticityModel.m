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

# General Notes
Variables written in all capital letters are generally meant
 to be constant values, initialized once in the beginning of the program
 and nowhere else.

Structures are initialized in the program by calling struct(...).
 These structures are meant to group variables that have some
 kind of commonality, e.g., belonging to the same neuron, behavior,
 model, etc. Note that this has a (negligible) effect on performance.

The grid is set up column-major order, with points accessed as
 grid(y,x), with the origin situated at the top-left corner and axes
 increasing right and down for x and y, respectively.

# MATLAB Code Generation
This file is compatible with MATLAB Code Generation, which greatly improves
 runtime. In order to keep it compatible, UI-related logic has been moved
 into a separate function. Note that this file's function allows optional
 arguments if being used in the non-codegen version. However, the codegen
 version of this function does not allow any optional arguments or
 structure fields: All arguments must be supplied as defined by the
 function definition.

# Tips/Tricks
If debugging, one can observe the workspace of the function by
 issuing the following command before execution: "dbstop if error".

# Function Signature
## Input Variables
arg_struct     - structure of n fields used to pass parameters that are
                 exposed in global optimization; if any fields are not
                 specified in the non-codegen version, they will be given
                 default values based on configuration
optional_parms - struct that contains additional arguments for the
                 model separate from specific model parameters
## Output Variables
opt_val_1      - return value signifying value of some cost function, used
                 for global optimization
opt_val_2      - return value (array) signifying value of cost function
                 used in FMRI group runs
%}
function [opt_val_1, opt_val_2] = automaticityModel(arg_struct, optional_parms) %#codegen
    %% ======================================= %%
    %%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
    %  =======================================  %
    % Code-generation declarations
    coder.extrinsic('getautoparams','displayautoresults');
    coder.extrinsic('tic','toc','struct2table','addpath');
    coder.varsize('chosen_rule');
    
    % Load dependencies
    addpath('classes');
    
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

    % Model/function behavior parameters (default values)
    VIS_INPUT_FROM_PARM   = 0;
    SUPPRESS_UI           = 0;
    OPTIMIZATION_CALC     = 0;
    FROST_ENABLED         = 1;
    COVIS_ENABLED         = 1;
    BUTTON_SWITCH_ENABLED = 1;
    PERF_OUTPUT           = 1;
    
    % Validate any supplied arguments
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

    % Override values with optional_parms if passed as an argument
    if nargin == 2
        if isfield(optional_parms, 'FMRI_META_GROUP_RUN')
            FMRI_META.GROUP_RUN = optional_parms.FMRI_META_GROUP_RUN;
        end
        if isfield(optional_parms, 'VIS_INPUT_FROM_PARM')
            VIS_INPUT_FROM_PARM = optional_parms.VIS_INPUT_FROM_PARM;
        end
    end
    
    %% Load visual stimulus matrix
    % Random Visual Input
    if 0
        loaded_input = load('datasets/randomVisualInput.mat');
        x_coordinates = loaded_input.x_coordinates;
        y_coordinates = loaded_input.y_coordinates;
    % Visual Input Matrix from optional_parms struct
    elseif VIS_INPUT_FROM_PARM
        x_coordinates = optional_parms.visualinput(:,1);
        y_coordinates = optional_parms.visualinput(:,2);
        coordinate_groups = zeros(length(x_coordinates), 1);
    elseif CONFIGURATION == MADDOX
    % MADDOX
        loaded_input = load('datasets/maddoxVisualInput.mat');
        x_coordinates = loaded_input.maddoxVisualInput(:, 1);
        y_coordinates = loaded_input.maddoxVisualInput(:, 2);
        coordinate_groups = loaded_input.maddoxVisualInput(:, 3);
    elseif CONFIGURATION == WALLIS
    % WALLIS
        loaded_input = load('datasets/wallisVisualInput.mat');
        x_coordinates = loaded_input.wallisVisualInput5(:,1);
        y_coordinates = loaded_input.wallisVisualInput5(:,2);
        coordinate_groups = zeros(length(x_coordinates), 1);
    elseif CONFIGURATION == FMRI
    % FMRI
        loaded_input = load('datasets/fMRI_data.mat');
        x_coordinates = loaded_input.x_coordinates;
        y_coordinates = loaded_input.y_coordinates;
        coordinate_groups = zeros(length(x_coordinates), 1);
    end

    %% Initialize/configure constants (though some data structure specific constants are initialized below)
    % Set properties of grid
    STIMULUS_GRID_SIZE = 100;                              % Length of side of square grid used for visual input; should be an even number
    BORDER_SIZE        = 20;                               % Width of border used to pad the grid such that visual stimulus on the edge still has an appropriate effect
    GRID_SIZE          = STIMULUS_GRID_SIZE+2*BORDER_SIZE; % Total length of grid, i.e., the stimulus grid size and the border

    % Button switch must be initialized before TRIALS
    BUTTON_SWITCH = struct('TRIALS', 600, 'PMC_A_weights', ones(GRID_SIZE,GRID_SIZE,1,4), 'PMC_B_weights', ones(GRID_SIZE,GRID_SIZE,1,4));

    % Set behavior and number of trials
    PRE_LEARNING_TRIALS  = PARAMS.PRE_LEARNING_TRIALS;  % Number of control trials run before learning trials
    LEARNING_TRIALS      = PARAMS.LEARNING_TRIALS;      % Number of learning trials in automaticity experiment
    POST_LEARNING_TRIALS = PARAMS.POST_LEARNING_TRIALS; % Number of trials where no learning is involved after learning trials
    if BUTTON_SWITCH_ENABLED
        TRIALS           = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS + BUTTON_SWITCH.TRIALS;
        IS_LEARNING      = [zeros(1,PRE_LEARNING_TRIALS), ones(1,LEARNING_TRIALS), zeros(1,POST_LEARNING_TRIALS), ones(1,BUTTON_SWITCH.TRIALS)];
    else
        TRIALS           = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS;
        IS_LEARNING      = [zeros(1,PRE_LEARNING_TRIALS), ones(1,LEARNING_TRIALS), zeros(1,POST_LEARNING_TRIALS)];
    end

    % Other parameters
    n = 1000;                     % Time period for one trial (in milliseconds)
    TAU = 1;
    LAMBDA = 20;
    W_MAX = PARAMS.W_MAX;         % Maximum weight for Hebbian Synapses
    accuracy = zeros(TRIALS, 1);  % Boolean matrix indicating if correct PMC neuron reacted
    NOISE = struct('PFC', PARAMS.NOISE_PFC, 'PMC', PARAMS.NOISE_PMC, 'MC', PARAMS.NOISE_MC);

    % Performance parameters
    loop_times    = zeros(TRIALS, 1); % Time needed for each loop (outer loop)
    trial_times   = zeros(TRIALS, 1); % Time required to run each time loop (inner loop)
    rt_calc_times = zeros(TRIALS, 1); % Time required to run reaction time calculation

    % Quantity of Visual Stimulus
    AREA = struct('LOWER_HALF', 1:GRID_SIZE/2, 'UPPER_HALF', GRID_SIZE/2+1:GRID_SIZE, ...
    			  'OUTER', [1:STIMULUS_GRID_SIZE/4+BORDER_SIZE, STIMULUS_GRID_SIZE*3/4+BORDER_SIZE+1:GRID_SIZE], ...
    			  'INNER', STIMULUS_GRID_SIZE/4+BORDER_SIZE+1:STIMULUS_GRID_SIZE*3/4+BORDER_SIZE, ...
    			  'ALL', 1:GRID_SIZE);
    VISUAL = struct('STIM', 50, ...
                    'x_coord', 0, 'y_coord', 0, 'coordinate_group', 0, ...
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
    CORRECT_RULE = 2;
    RULE = VISUAL.RULES(chosen_rule);

    % Radial Basis Function
    [X, Y] = meshgrid(1:GRID_SIZE, 1:GRID_SIZE);
    RBF = struct('RADIUS', 0.8, 'rbv', zeros(GRID_SIZE), 'X', X, 'Y', Y);

    %% COVIS Model
    COVIS_PARMS = struct('DELTA_C', PARAMS.COVIS_DELTA_C, 'DELTA_E', PARAMS.COVIS_DELTA_E, 'PERSEV', PARAMS.COVIS_PERSEV, ...
                         'LAMBDA', PARAMS.COVIS_LAMBDA, 'NUM_GUESS', 5, 'NUM_RULES', 4);
    COVIS_VARS = struct('correct_rule', VISUAL.RULES(CORRECT_RULE), 'rules', 1:COVIS_PARMS.NUM_RULES, 'saliences', ones(COVIS_PARMS.NUM_RULES,1), ...
                        'rule_weights', ones(COVIS_PARMS.NUM_RULES,1), 'rule_prob', ones(COVIS_PARMS.NUM_RULES,1), ...
                        'rule_log', ones(TRIALS,1));
    
    %% General settings for PFC, PMC, and MC neurons
    % Note that reactions is big enough for both learning trials and no-learning trials to allow for comparisons
    % PFC general information
    PFC = struct( ...                       
        'DECISION_PT', PARAMS.PFC_DECISION_PT, ...   % threshold which determines which PFC neuron acts on a visual input
        'reactions', zeros(TRIALS,3), ...            % stores information about PFC neuron reactions during trial
        'activations', zeros(TRIALS,1) ...
    );

    % PMC general information
    PMC = struct( ...                           
        'DECISION_PT', PARAMS.PMC_DECISION_PT, ...   % threshold which determines which PMC neuron acts on a visual input
        'reactions', zeros(TRIALS,3), ...            % stores information about PMC neuron reactions during trial
        'alpha', zeros(TRIALS,n), ...                % PMC_A.out + PMC_B.out
        'activations', zeros(TRIALS,1) ...
    );

    % MC general parameters
    MC = struct( ...
        'DECISION_PT', PARAMS.MC_DECISION_PT, ...
        'reactions', zeros(TRIALS,3), ...
        'activations', zeros(TRIALS,1), ...
        'PRIMARY_WEIGHT', 0.9, ...
        'SECONDARY_WEIGHT', 0.1 ...
    );

    %% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses)
    %{
        NMDA - upper threshold
        AMPA - lower threshold
        Strengthening occurs if the [voltage integral] > [NMDA]
        Weakening occurs if [NMDA] - [voltage integral] - [AMPA] > 0
    %}
    Hebbian = struct('heb_coef', PARAMS.HEB_CONSTS, 'anti_heb', PARAMS.HEB_CONSTS, 'NMDA', PARAMS.NMDA, 'AMPA', PARAMS.AMPA, ...
                     'heb_coef_mc', PARAMS.HEB_CONSTS, 'anti_heb_mc', PARAMS.HEB_CONSTS, 'NMDA_MC', PARAMS.NMDA, 'AMPA_MC', PARAMS.AMPA);

    %% Set up neurons
    PFC_A = PFCNeuron(n, TAU, LAMBDA, PARAMS.PFC_A_W_OUT_MDN);
    PFC_B = PFCNeuron(n, TAU, LAMBDA, PARAMS.PFC_B_W_OUT_MDN);

    PMC_A = PMCNeuron(n, TAU, LAMBDA, TRIALS, PARAMS.PMC_A_W_OUT, CONFIGURATION == FMRI, COVIS_ENABLED, GRID_SIZE);
    PMC_B = PMCNeuron(n, TAU, LAMBDA, TRIALS, PARAMS.PMC_B_W_OUT, CONFIGURATION == FMRI, COVIS_ENABLED, GRID_SIZE);

    MC_A = MCNeuron(n, TAU, LAMBDA, TRIALS);
    MC_B = MCNeuron(n, TAU, LAMBDA, TRIALS);

    Driv_PFC = Driv_PFCNeuron(n, TAU, LAMBDA, PARAMS.DRIV_PFC_W_OUT);

    CN = CNNeuron(n, TAU, LAMBDA, TRIALS);
    
    GP = GPNeuron(n, TAU, LAMBDA, TRIALS);

    MDN = struct('activations', zeros(TRIALS,1));
    MDN_A = MDNNeuron(n, TAU, LAMBDA, PARAMS.MDN_A_W_OUT);
    MDN_B = MDNNeuron(n, TAU, LAMBDA, PARAMS.MDN_B_W_OUT);

    AC_A = ACNeuron(n, TAU, LAMBDA);
    AC_B = ACNeuron(n, TAU, LAMBDA);

    %% ============================ %%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %  ============================  %
    start_time = tic;
    %% Learning trials
    for trial=1:TRIALS
        loopStart = tic;
        %% Initialize each neuron for the trial
        PFC_A = PFC_A.reset(); PFC_B = PFC_B.reset();
        PMC_A = PMC_A.reset(); PMC_B = PMC_B.reset();
        MC_A = MC_A.reset(); MC_B = MC_B.reset();
        
        if FROST_ENABLED
            Driv_PFC = Driv_PFC.reset();
            CN = CN.reset();
            GP = GP.reset();
            MDN_A = MDN_A.reset(); MDN_B = MDN_B.reset();
            AC_A = AC_A.reset(); AC_B = AC_B.reset();
        end
        %% Initialize COVIS components (choose a rule)
        if COVIS_ENABLED
            if trial <= COVIS_PARMS.NUM_GUESS
                chosen_rule = randi(length(COVIS_VARS.rules));
            elseif accuracy(trial-1) == 1
                chosen_rule = COVIS_VARS.rule_log(trial-1);
            else
                chosen_rule = rand_discrete(COVIS_VARS.rule_prob);
            end
            RULE = VISUAL.RULES(chosen_rule);
            COVIS_VARS.rule_log(trial) = chosen_rule;
        end
        %% Button Switch if enabled and correct trials
        if BUTTON_SWITCH_ENABLED && trial == TRIALS - BUTTON_SWITCH.TRIALS + 1
            [MC.SECONDARY_WEIGHT, MC.PRIMARY_WEIGHT] = deal(MC.PRIMARY_WEIGHT, MC.SECONDARY_WEIGHT);
            if CONFIGURATION == FMRI
                BUTTON_SWITCH.PMC_A_weights(:,:,1,:) = PMC_A.weights(:,:,1,:);
                BUTTON_SWITCH.PMC_B_weights(:,:,1,:) = PMC_B.weights(:,:,1,:);
            else
                BUTTON_SWITCH.PMC_A_weights(:,:,1,:) = PMC_A.weights(:,:,trial,:);
                BUTTON_SWITCH.PMC_B_weights(:,:,1,:) = PMC_B.weights(:,:,trial,:);
            end
        end

        %% Set PMC weights, potentially dependent on COVIS
        % If first trial, set to initial weights
        if CONFIGURATION == FMRI || trial==1
            if COVIS_ENABLED
                PMC_A_weights = PMC_A.weights(:,:,1,chosen_rule);
                PMC_B_weights = PMC_B.weights(:,:,1,chosen_rule);
            else
                PMC_A_weights = PMC_A.weights(:,:,1);
                PMC_B_weights = PMC_B.weights(:,:,1);
            end
            MC_A_weights = MC_A.weights(:,1);
            MC_B_weights = MC_B.weights(:,1);
        % Else, set weights to results of previous trial
        else
            if COVIS_ENABLED
                PMC_A_weights = PMC_A.weights(:,:,trial-1,chosen_rule);
                PMC_B_weights = PMC_B.weights(:,:,trial-1,chosen_rule);
            else
                PMC_A_weights = PMC_A.weights(:,:,trial-1);
                PMC_B_weights = PMC_B.weights(:,:,trial-1);
            end
            MC_A_weights = MC_A.weights(:,trial-1);
            MC_B_weights = MC_B.weights(:,trial-1);
        end

        %% Determine visual stimulus in range [1, GRID_SIZE] to pick random gabor for each trial, padded with
        % the BORDER_SIZE such that the visual stimulus is accounted for properly
        VISUAL.y_coord = y_coordinates(trial) + BORDER_SIZE;
        VISUAL.x_coord = x_coordinates(trial) + BORDER_SIZE;
        VISUAL.coordinate_group = coordinate_groups(trial);

        %% Calculate visual stimulus effect using Radial Basis Function (RBF) implementation
        % Calculate RBF grid
        RBF.rbv(:,:) = exp( -(sqrt((VISUAL.y_coord-RBF.Y).^2 + (VISUAL.x_coord-RBF.X).^2))/RBF.RADIUS ) * VISUAL.STIM;
        % Sum RBF values depending on rule to find PFC_A and PFC_B v_stim values
        % Note that stim matrices are row-major order (e.g., indexed by y, then x)
        PFC_A.v_stim = sum(sum(RBF.rbv(RULE(1).A_Y, RULE(1).A_X)));
        PFC_B.v_stim = sum(sum(RBF.rbv(RULE(1).B_Y, RULE(1).B_X)));
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(sum(RBF.rbv(:,:).*PMC_A_weights));
        PMC_B.v_stim = sum(sum(RBF.rbv(:,:).*PMC_B_weights));
        % Scale v_stim values to prevent them from becoming too large
        PFC_A.v_stim = PFC_A.v_stim * PFC_A.V_SCALE;
        PFC_B.v_stim = PFC_B.v_stim * PFC_B.V_SCALE;
        PMC_A.v_stim = PMC_A.v_stim * PMC_A.V_SCALE;
        PMC_B.v_stim = PMC_B.v_stim * PMC_B.V_SCALE;

        %% Individual Time Trial Loop (iterating through n)
        timeTrialStart = tic;
        if FROST_ENABLED
            %% FROST Calculations
            for i=1:n-1
                PFC_A = PFC_A.iterate_FROST(NOISE.PFC, PFC_B, PMC_A, MDN_A, AC_A);
                PFC_B = PFC_B.iterate_FROST(NOISE.PFC, PFC_A, PMC_B, MDN_B, AC_B);

                PMC_A = PMC_A.iterate(NOISE.PMC, PMC_B, PFC_A);
                PMC_B = PMC_B.iterate(NOISE.PMC, PMC_A, PFC_B);
                
                MC_A = MC_A.iterate(trial, NOISE.MC, MC_B, PMC_A, PMC_B, MC.PRIMARY_WEIGHT, MC.SECONDARY_WEIGHT);
                MC_B = MC_B.iterate(trial, NOISE.MC, MC_A, PMC_B, PMC_A, MC.PRIMARY_WEIGHT, MC.SECONDARY_WEIGHT);               

                Driv_PFC = Driv_PFC.iterate();

                CN = CN.iterate(Driv_PFC);
                
                GP = GP.iterate(CN);

                MDN_A = MDN_A.iterate(PFC_A, GP);
                MDN_B = MDN_B.iterate(PFC_B, GP);

                AC_A = AC_A.iterate(PFC_A);
                AC_B = AC_B.iterate(PFC_B);
            end
        else
            %% Non-FROST Calculation
            for i=1:n-1
                PFC_A = PFC_A.iterate(NOISE.PFC, PFC_B, PMC_A);
                PFC_B = PFC_B.iterate(NOISE.PFC, PFC_A, PMC_B);

                PMC_A = PMC_A.iterate(NOISE.PMC, PMC_B, PFC_A);
                PMC_B = PMC_B.iterate(NOISE.PMC, PMC_A, PFC_B);
                
                MC_A = MC_A.iterate(trial, NOISE.MC, MC_B, PMC_A, PMC_B, MC.PRIMARY_WEIGHT, MC.SECONDARY_WEIGHT);
                MC_B = MC_B.iterate(trial, NOISE.MC, MC_A, PMC_B, PMC_A, MC.PRIMARY_WEIGHT, MC.SECONDARY_WEIGHT);
            end
        end
        %% Record post-time-loop numbers
        % Record "alpha" function, summing PMC A and PMC B output
        PMC.alpha(trial,:) = PMC_A.out + PMC_B.out;
        % Record total neuron activations
        PFC.activations(trial) = trapz(PFC_A.out + PFC_B.out);
        CN.activations(trial) = trapz(CN.out);
        GP.activations(trial) = trapz(GP.out);
        MDN.activations(trial) = trapz(MDN_A.out + MDN_B.out);
        PMC.activations(trial) = trapz(PMC.alpha(trial,:));
        MC.activations(trial) = trapz(MC_A.out + MC_B.out);
        trial_times(trial) = toc(timeTrialStart);

        %% Determine decision neuron and reaction time, and record accuracy
        rt_start_time = tic;
        % Determine reacting neuron and latency
        [neuron_id_PFC, latency] = determine_reacting_neuron(PFC_A.out, PFC_B.out, PFC.DECISION_PT);
        PFC.reactions(trial,:) = [neuron_id_PFC, latency, VISUAL.coordinate_group];
        [neuron_id_PMC, latency] = determine_reacting_neuron(PMC_A.out, PMC_B.out, PMC.DECISION_PT);
        PMC.reactions(trial,:) = [neuron_id_PMC, latency, VISUAL.coordinate_group];
        [neuron_id_MC, latency] = determine_reacting_neuron(MC_A.out, MC_B.out, MC.DECISION_PT);
        MC.reactions(trial,:) = [neuron_id_MC, latency, VISUAL.coordinate_group];
        % Determine accuracy
        if COVIS_ENABLED
            accuracy(trial) = double((any(VISUAL.x_coord == COVIS_VARS.correct_rule.B_X) && any(VISUAL.y_coord == COVIS_VARS.correct_rule.B_Y)) + 1) == neuron_id_MC;
        else
            accuracy(trial) = double((any(VISUAL.x_coord == RULE.B_X) && any(VISUAL.y_coord == RULE.B_Y)) + 1) == neuron_id_MC;
        end
        rt_calc_times(trial) = toc(rt_start_time);

        %% Weight change calculations
        if CONFIGURATION == FMRI
            idx_weight = 1;
        else
            idx_weight = trial;
        end
        if COVIS_ENABLED
            idx_rule = chosen_rule;
        else
            idx_rule = 1;
        end
        if IS_LEARNING(trial)
            %% Calculate and store activation integrals (for performance reasons)
            integral_PFC_A = PFC_A.integralPosVolt(); % visual input to PMC_A neuron (presynaptic)
            integral_PMC_A = PMC_A.integralPosVolt(); % activation of PMC_A neuron (postsynaptic)
            integral_PFC_B = PFC_B.integralPosVolt(); % visual input to PMC_B neuron (presynaptic)
            integral_PMC_B = PMC_B.integralPosVolt(); % activation of PMC_B neuron (postsynaptic)
            integral_MC_A = MC_A.integralPosVolt(); % activation of MC_A neuron
            integral_MC_B = MC_B.integralPosVolt(); % activation of MC_B neuron

            %% Calculation of Hebbian Weight for PMC_A
            % Ensure g(t)-1 and g(2)-2 are never less than zero
            g_t_1_A = max(0, integral_PMC_A - Hebbian.NMDA);
            g_t_2_A = max(0, Hebbian.NMDA - integral_PMC_A - Hebbian.AMPA);

            % Determine new weights of visual PMC_A synapses (limit values to range of [0, W_MAX])
            PMC_A_weights(:,:,1,1) = PMC_A_weights + RBF.rbv(:,:).*(Hebbian.heb_coef*integral_PFC_A*g_t_1_A.*(W_MAX - PMC_A_weights) - Hebbian.anti_heb*integral_PFC_A*g_t_2_A.*PMC_A_weights);
            PMC_A.weights(:,:,idx_weight,idx_rule) = min(max(PMC_A_weights,0),W_MAX);

            %% Calculation of Hebbian Weight for PMC_B
            % Ensures g(t)-1 and g(2)-2 are never less than zero
            g_t_1_B = max(0, integral_PMC_B - Hebbian.NMDA);
            g_t_2_B = max(0, Hebbian.NMDA - integral_PMC_B - Hebbian.AMPA);

            % Determine new weights of visual PMC_B synapses (limit values to range of [0, W_MAX])
            PMC_B_weights(:,:,1,1) = PMC_B_weights + RBF.rbv(:,:).*(Hebbian.heb_coef*integral_PFC_B*g_t_1_B.*(W_MAX - PMC_B_weights) - Hebbian.anti_heb*integral_PFC_B*g_t_2_B.*PMC_B_weights);
            PMC_B.weights(:,:,idx_weight,idx_rule) = min(max(PMC_B_weights,0),W_MAX);
            
            %% Calculation of Hebbian Weights for MC_A
            g_t_1_MCA = max(0, integral_MC_A - Hebbian.NMDA_MC);
            g_t_2_MCA = max(0, Hebbian.NMDA_MC - integral_MC_A - Hebbian.AMPA_MC);
            
            MC_A.weights(:,idx_weight) = MC_A_weights + [MC.PRIMARY_WEIGHT; MC.SECONDARY_WEIGHT].*(integral_PMC_A*(Hebbian.heb_coef_mc*g_t_1_MCA.*(MC_A.W_MAX - MC_A_weights) - Hebbian.anti_heb_mc*g_t_2_MCA.*MC_A_weights));
            MC_A.weights(:,idx_weight) = min(max(MC_A.weights(:,idx_weight),0),MC_A.W_MAX);
            
            %% Calculation of Hebbian Weights for MC_B
            g_t_1_MCB = max(0, integral_MC_B - Hebbian.NMDA_MC);
            g_t_2_MCB = max(0, Hebbian.NMDA_MC - integral_MC_B - Hebbian.AMPA_MC);
            
            MC_B.weights(:,idx_weight) = MC_B_weights + [MC.PRIMARY_WEIGHT; MC.SECONDARY_WEIGHT].*(integral_PMC_B*(Hebbian.heb_coef_mc*g_t_1_MCB.*(MC_B.W_MAX - MC_B_weights) - Hebbian.anti_heb_mc*g_t_2_MCB.*MC_B_weights));
            MC_B.weights(:,idx_weight) = min(max(MC_B.weights(:,idx_weight),0),MC_B.W_MAX);
            
        % Else, if not learning, set new weights to previous weights
        else
            PMC_A.weights(:,:,idx_weight,idx_rule) = PMC_A_weights;
            PMC_B.weights(:,:,idx_weight,idx_rule) = PMC_B_weights;
            MC_A.weights(:,idx_weight) = MC_A_weights;
            MC_B.weights(:,idx_weight) = MC_B_weights;
        end
        
        % If COVIS is enabled and weight matrix has time dimension, update
        % all other weight matrices for this iteration
        if COVIS_ENABLED && CONFIGURATION ~= FMRI
            PMC_A.weights(:,:,trial,1:COVIS_PARMS.NUM_RULES ~= chosen_rule) = PMC_A.weights(:,:,trial-1,1:COVIS_PARMS.NUM_RULES ~= chosen_rule);
            PMC_B.weights(:,:,trial,1:COVIS_PARMS.NUM_RULES ~= chosen_rule) = PMC_B.weights(:,:,trial-1,1:COVIS_PARMS.NUM_RULES ~= chosen_rule);
        end
        
        % Record average weight for PMC_A and PMC_B
        PMC_A.weights_avg(trial) = mean(mean(PMC_A.weights(:,:,idx_weight,CORRECT_RULE)));
        PMC_B.weights_avg(trial) = mean(mean(PMC_B.weights(:,:,idx_weight,CORRECT_RULE)));
        
        %% COVIS Calculations - readjusting saliences, weights
        if COVIS_ENABLED && trial < PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS
            % Step 1: re-adjust saliences & weights
            if accuracy(trial) == 1
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
        end

        %% Print data to console
        if not(SUPPRESS_UI)
            consoleprogressbar('TRIALS COMPLETED', trial, TRIALS);
        end
        loop_times(trial) = toc(loopStart);
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
                norm_output_rt = [median(PMC.reactions(FMRI_META.SES_1,2)), ...
                                  median(PMC.reactions(FMRI_META.SES_4,2)), ...
                                  median(PMC.reactions(FMRI_META.SES_10,2)), ...
                                  median(PMC.reactions(FMRI_META.SES_20,2))]./1000;
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
        displayautoresults(FROST_ENABLED, COVIS_ENABLED, BUTTON_SWITCH_ENABLED, BUTTON_SWITCH, COVIS_VARS, FMRI_META, CONFIGURATION, MADDOX, WALLIS, FMRI, TAU, n, RBF, BORDER_SIZE, VISUAL, TRIALS, PRE_LEARNING_TRIALS, LEARNING_TRIALS, POST_LEARNING_TRIALS, accuracy, PFC, PMC, MC, PFC_A, PFC_B, PMC_A, PMC_B, MC_A, MC_B, Driv_PFC, CN, GP, MDN_A, MDN_B, AC_A, AC_B, PERF_OUTPUT, start_time, loop_times, trial_times, rt_calc_times, chosen_rule);
    end
    return;
end

%% =============================== %%
%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%
%  ===============================  %
% Return what neuron reacts to the stimuli, and the latency
% Returns neuron_id = 1 for n1, neuron_id = 2 for n2
function [neuron_id, latency] = determine_reacting_neuron(neuron_1, neuron_2, decision_pt)
    n1_latency = find(cumtrapz(neuron_1) >= decision_pt, 1);
    n2_latency = find(cumtrapz(neuron_2) >= decision_pt, 1);
    % n1_latency or n2_latency could be empty if the decision_pt was never reached
    % If so, set it to the maximum allowed value
    if isempty(n1_latency)
        n1_latency = length(neuron_1);
    end
    if isempty(n2_latency)
        n2_latency = length(neuron_2);
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
        neuron_id = double(trapz(neuron_1) < trapz(neuron_2));
        latency = length(neuron_1);
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