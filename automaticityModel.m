%AUTOMATICITYMODEL automaticity model for PMC neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Model of Automaticity in Rule-Based Learning        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
# Description
Run an automaticity model based on configurations defined
 by ModelConfig.

# General Notes
Variables written in all capital letters are generally meant
 to be constant values, initialized once in the beginning of the program.

Structures are meant to group variables that have some commonality.
 Note that this has a (negligible) effect on performance.

The grid is column-major order, with points accessed as
 grid(y,x) The origin is situated at the top-left corner and axes
 increase right and down for x and y, respectively.

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
                 default values based on config
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
    coder.extrinsic('struct2table','cell2struct','addpath','genpath');
    coder.varsize('chosen_rule');
    
    % Load dependencies
    addpath(genpath('.'));
    coder.extrinsic('getmodelparams','getconstants','displayautoresults');

    % Load config and config parameters
    config = ModelConfigButtonSwitch();

    % Get parameters
    PARAMS = struct('PRE_LEARNING_TRIALS',0,'LEARNING_TRIALS',0,'POST_LEARNING_TRIALS',0,'PFC_DECISION_PT',0,'PMC_DECISION_PT',0,'MC_DECISION_PT',0,'HEB_CONSTS',0,'NMDA',0,'AMPA',0,'W_MAX',0,'NOISE_PFC',0,'NOISE_PMC',0,'NOISE_MC',0,'PMC_A_W_OUT',0,'PMC_B_W_OUT',0,'PFC_A_W_OUT_MDN',0,'PFC_B_W_OUT_MDN',0,'DRIV_PFC_W_OUT',0,'MDN_A_W_OUT',0,'MDN_B_W_OUT',0,'COVIS_DELTA_C',0,'COVIS_DELTA_E',0,'COVIS_PERSEV',0,'COVIS_LAMBDA',0);
    PARAMS = getmodelparams(config);

    % Model/function behavior parameters (default values)
    VIS_INPUT_FROM_PARM   = 0;
    SUPPRESS_UI           = 0;
    OPTIMIZATION_CALC     = 0;
    
    % Validate any supplied arguments
    if nargin >= 1
        % Determine if parms are valid
        if not(areparamsvalid(PARAMS))
            disp(struct2table(PARAMS));
            error('Parameters not valid.');
        end
        PARAMS = absorbstruct(PARAMS, arg_struct);
    end

    % Override values with optional_parms if passed as an argument
    if nargin == 2
        if isfield(optional_parms, 'FMRI_META_GROUP_RUN') && isa(config, 'ModelConfigButtonSwitch')
            config.meta.optimization.GROUP_RUN = optional_parms.FMRI_META_GROUP_RUN;
        end
        if isfield(optional_parms, 'VIS_INPUT_FROM_PARM')
            VIS_INPUT_FROM_PARM = optional_parms.VIS_INPUT_FROM_PARM;
        end
    end
    
    %% Load visual stimulus matrix
%     if VIS_INPUT_FROM_PARM
%     % Visual Input Matrix from optional_parms struct
%         x_coordinates = optional_parms.visualinput(:,1);
%         y_coordinates = optional_parms.visualinput(:,2);
%         coordinate_groups = zeros(length(x_coordinates), 1);
% %     elseif config == ModelConfig.IMAGE_CORR
% %         loaded_input = load('datasets/imageVisualInput.mat');
% %         x_coordinates = loaded_input.visualInput.x;
% %         y_coordinates = loaded_input.visualInput.y;
% %         coordinate_groups = loaded_input.visualInput.groups;
%     elseif isa(config, 'ModelConfigButtonSwitch')
%         loaded_input = load('datasets/fMRI_data.mat');
%         x_coordinates = loaded_input.x_coordinates;
%         y_coordinates = loaded_input.y_coordinates;
%         coordinate_groups = zeros(length(x_coordinates), 1);
%     end
    if VIS_INPUT_FROM_PARM
        x_coordinates = optional_parms.visualinput(:,1);
        y_coordinates = optional_parms.visualinput(:,2);
        coordinate_groups = zeros(length(x_coordinates), 1);
    else
        [x_coordinates, y_coordinates, coordinate_groups] = config.loadCoords();
    end

    %% Initialize/configure constants (though some data structure specific constants are initialized below)
    % Set behavior and number of trials
    PRE_LEARNING_TRIALS  = PARAMS.PRE_LEARNING_TRIALS;  % Number of control trials run before learning trials
    LEARNING_TRIALS      = PARAMS.LEARNING_TRIALS;      % Number of learning trials in automaticity experiment
    POST_LEARNING_TRIALS = PARAMS.POST_LEARNING_TRIALS; % Number of trials where no learning is involved after learning trials
    if isa(config, 'ModelConfigButtonSwitch')
        TRIALS           = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS + config.meta.trialsAfterSwitch;
        IS_LEARNING      = [zeros(1,PRE_LEARNING_TRIALS), ones(1,LEARNING_TRIALS), zeros(1,POST_LEARNING_TRIALS), ones(1,config.meta.trialsAfterSwitch)];
    else
        TRIALS           = PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS;
        IS_LEARNING      = [zeros(1,PRE_LEARNING_TRIALS), ones(1,LEARNING_TRIALS), zeros(1,POST_LEARNING_TRIALS)];
    end
    config = config.setTrials(TRIALS);

    % Other parameters
    n = Neuron.n;                 % Time period for one trial (in milliseconds)
    W_MAX = PARAMS.W_MAX;         % Maximum weight for Hebbian Synapses
    accuracy = zeros(TRIALS, 1);  % Boolean matrix indicating if correct PMC neuron reacted

    % Get visual stimulus variable and establish correct and chosen rules.
    VISUAL = config.visual;
    chosen_rule = 1;
    correct_rule = 2;
    RULE = VISUAL.RULES(chosen_rule);

    % Radial Basis Function
    RBF = RadialBasisFunction(config.GRID_SIZE, VISUAL.STIM);

    %% COVIS Model
    config = config.initCOVISRules(PARAMS, chosen_rule, correct_rule, config.trials);
    
    %% General settings for neurons
    % Note that reactions is big enough for both learning trials and no-learning trials to allow for comparisons
    PFC = struct( ...                       
        'DECISION_PT', PARAMS.PFC_DECISION_PT, ...   % threshold which determines which PFC neuron acts on a visual input
        'reactions',   zeros(TRIALS,3), ...          % stores information about PFC neuron reactions during trial
        'activations', zeros(TRIALS,1) ...
    );

    PMC = struct( ...                           
        'DECISION_PT', PARAMS.PMC_DECISION_PT, ...   % threshold which determines which PMC neuron acts on a visual input
        'reactions',   zeros(TRIALS,3), ...          % stores information about PMC neuron reactions during trial
        'alpha',       zeros(TRIALS,n), ...          % PMC_A.out + PMC_B.out
        'activations', zeros(TRIALS,1) ...
    );

    MC = struct( ...
        'DECISION_PT',      PARAMS.MC_DECISION_PT, ...
        'reactions',        zeros(TRIALS,3), ...
        'activations',      zeros(TRIALS,1), ...
        'A_area',           zeros(TRIALS,1), ...
        'B_area',           zeros(TRIALS,1) ...
    );

    MDN = struct('activations', zeros(TRIALS,1));

    default_hebbian_consts = HebbianConst(PARAMS.HEB_CONSTS, PARAMS.HEB_CONSTS, PARAMS.NMDA, PARAMS.AMPA);
                 
    %% Set up neurons
    PFC_A = PFCNeuron(PARAMS.PFC_A_W_OUT_MDN);
    PFC_B = PFCNeuron(PARAMS.PFC_B_W_OUT_MDN);

    PMC_A = PMCNeuron(TRIALS, PARAMS.PMC_A_W_OUT, W_MAX, config.isCOVISEnabled, config.GRID_SIZE, default_hebbian_consts);
    PMC_B = PMCNeuron(TRIALS, PARAMS.PMC_B_W_OUT, W_MAX, config.isCOVISEnabled, config.GRID_SIZE, default_hebbian_consts);

    MC_A = MCNeuron(TRIALS, W_MAX, default_hebbian_consts);
    MC_B = MCNeuron(TRIALS, W_MAX, default_hebbian_consts);

    Driv_PFC = Driv_PFCNeuron(PARAMS.DRIV_PFC_W_OUT);

    CN = CNNeuron(TRIALS);
    
    GP = GPNeuron(TRIALS);

    MDN_A = MDNNeuron(PARAMS.MDN_A_W_OUT);
    MDN_B = MDNNeuron(PARAMS.MDN_B_W_OUT);

    AC_A = ACNeuron();
    AC_B = ACNeuron();

    %% ============================ %%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %  ============================  %
    %% Learning trials
    for trial=1:TRIALS
        config.trial = trial;
        %% Initialize each neuron for the trial
        PFC_A = PFC_A.reset(); PFC_B = PFC_B.reset();
        PMC_A = PMC_A.reset(); PMC_B = PMC_B.reset();
        MC_A = MC_A.reset(); MC_B = MC_B.reset();
        
        if config.isFROSTEnabled
            Driv_PFC = Driv_PFC.reset();
            CN = CN.reset();
            GP = GP.reset();
            MDN_A = MDN_A.reset(); MDN_B = MDN_B.reset();
            AC_A = AC_A.reset(); AC_B = AC_B.reset();
        end
        %% Initialize COVIS components (choose a rule)
        if config.isCOVISEnabled
            if trial <= config.COVISRules.GUESSES
                config.COVISRules.chosen = randi(config.COVISRules.NUM);
            elseif accuracy(trial-1) == 1 || trial > PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS
                config.COVISRules.chosen = config.COVISRules.log(trial-1);
            else
                config.COVISRules.chosen = rand_discrete(config.COVISRules.prob);
            end
            RULE = VISUAL.RULES(config.COVISRules.chosen);
            config.COVISRules.log(trial) = config.COVISRules.chosen;
        end
        %% Button Switch if enabled and correct trials
        if config.shouldButtonSwitch()
            config = config.doButtonSwitch(PMC_A, PMC_B);
        end

        %% Determine visual stimulus in range [1, GRID_SIZE] to pick random gabor for each trial, padded with
        % the BORDER_SIZE such that the visual stimulus is accounted for properly
        VISUAL.coord.x = x_coordinates(trial) + config.BORDER_SIZE + criterialnoise();
        VISUAL.coord.y = y_coordinates(trial) + config.BORDER_SIZE;
        VISUAL.coord.group = coordinate_groups(trial);
       

        %% Calculate visual stimulus effect using Radial Basis Function (RBF) implementation
        % Calculate RBF grid
        RBF = RBF.resolvestimulus(VISUAL.coord.x, VISUAL.coord.y);
        % Sum RBF values depending on rule to find PFC_A and PFC_B v_stim values
        % Note that stim matrices are row-major order (e.g., indexed by y, then x)
        PFC_A.v_stim = sum(sum(RBF.rbv(RULE(1).A_Y, RULE(1).A_X))) * PFC_A.V_SCALE;
        PFC_B.v_stim = sum(sum(RBF.rbv(RULE(1).B_Y, RULE(1).B_X))) * PFC_B.V_SCALE;
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(sum(RBF.rbv(:,:).*(PMC_A.weightsForTrial(config)))) * PMC_A.V_SCALE;
        PMC_B.v_stim = sum(sum(RBF.rbv(:,:).*(PMC_B.weightsForTrial(config)))) * PMC_B.V_SCALE;

        %% Individual Time Trial Loop (iterating through n)
        if config.isFROSTEnabled
            %% FROST Calculations
            for i=1:n-1
                PFC_A = PFC_A.iterate_FROST(PFC_B, MDN_A, AC_A);
                PFC_B = PFC_B.iterate_FROST(PFC_A, MDN_B, AC_B);

                PMC_A = PMC_A.iterate(PMC_B, PFC_A);
                PMC_B = PMC_B.iterate(PMC_A, PFC_B);
                
                MC_A = MC_A.iterate(trial, MC_B, PMC_A, PMC_B);
                MC_B = MC_B.iterate(trial, MC_A, PMC_B, PMC_A);               

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
                PFC_A = PFC_A.iterate(PFC_B);
                PFC_B = PFC_B.iterate(PFC_A);

                PMC_A = PMC_A.iterate(PMC_B, PFC_A);
                PMC_B = PMC_B.iterate(PMC_A, PFC_B);
                
                MC_A = MC_A.iterate(trial, MC_B, PMC_A, PMC_B);
                MC_B = MC_B.iterate(trial, MC_A, PMC_B, PMC_A);
            end
        end
        %% Record post-time-loop numbers
        % Record "alpha" function, summing PMC A and PMC B output
        PMC.alpha(trial,:) = PMC_A.out + PMC_B.out;
        MC.A_area(trial) = trapz(MC_A.out);
        MC.B_area(trial) = trapz(MC_B.out);
        % Record total neuron activations
        PFC.activations(trial) = trapz(PFC_A.out + PFC_B.out);
        CN.activations(trial) = trapz(CN.out);
        GP.activations(trial) = trapz(GP.out);
        MDN.activations(trial) = trapz(MDN_A.out + MDN_B.out);
        PMC.activations(trial) = trapz(PMC.alpha(trial,:));
        MC.activations(trial) = trapz(MC_A.out + MC_B.out);

        %% Determine decision neuron and reaction time, and record accuracy
        % Determine reacting neuron and latency
        [neuron_id_PFC, latency] = determine_reacting_neuron(PFC_A.out, PFC_B.out, PFC.DECISION_PT);
        PFC.reactions(trial,:) = [neuron_id_PFC, latency, VISUAL.coord.group];
        [neuron_id_PMC, latency] = determine_reacting_neuron(PMC_A.out, PMC_B.out, PMC.DECISION_PT);
        PMC.reactions(trial,:) = [neuron_id_PMC, latency, VISUAL.coord.group];
        [neuron_id_MC, latency] = determine_reacting_neuron(MC_A.out, MC_B.out, MC.DECISION_PT);
        MC.reactions(trial,:) = [neuron_id_MC, latency, VISUAL.coord.group];
        % Determine accuracy
        if config.isCOVISEnabled
            accuracy(trial) = determinecorrectneuron(VISUAL.coord.x, VISUAL.coord.y, VISUAL.RULES(config.COVISRules.correct)) == neuron_id_MC;
        else
            accuracy(trial) = determinecorrectneuron(VISUAL.coord.x, VISUAL.coord.y, RULE(1)) == neuron_id_MC;
        end

        %% Weight change calculations
        if IS_LEARNING(trial)
            PMC_A = PMC_A.doHebbianLearning(config, RBF.rbv, PFC_A);
            PMC_B = PMC_B.doHebbianLearning(config, RBF.rbv, PFC_B);
            
            MC_A = MC_A.doHebbianLearning(config, PMC_A);
            MC_B = MC_B.doHebbianLearning(config, PMC_B);
        % Else, if not learning, set new weights to previous weights
        else
            PMC_A.weights(:,:,config.weightIdx,config.chosenRuleIdx) = PMC_A.weightsForTrial(config);
            PMC_B.weights(:,:,config.weightIdx,config.chosenRuleIdx) = PMC_B.weightsForTrial(config);
            MC_A.weights(:,trial) = MC_A.previousweights(trial);
            MC_B.weights(:,trial) = MC_B.previousweights(trial);
        end
        
        % If COVIS is enabled and weight matrix has time dimension, update
        % all other weight matrices for this iteration
        if config.isCOVISEnabled && TRIALS <= PMCNeuron.LARGE_TRIAL_BOUNDARY
            PMC_A.weights(:,:,trial,1:config.COVISRules.NUM ~= chosen_rule) = PMC_A.weights(:,:,trial-1,1:config.COVISRules.NUM ~= chosen_rule);
            PMC_B.weights(:,:,trial,1:config.COVISRules.NUM ~= chosen_rule) = PMC_B.weights(:,:,trial-1,1:config.COVISRules.NUM ~= chosen_rule);
        end
        
        % Record average weight for PMC_A and PMC_B
        PMC_A.weights_avg(trial) = mean(mean(PMC_A.weights(:,:,config.weightIdx,config.correctRuleIdx)));
        PMC_B.weights_avg(trial) = mean(mean(PMC_B.weights(:,:,config.weightIdx,config.correctRuleIdx)));
        
        %% COVIS Calculations - readjusting saliences, weights
        % TODO: Do we not want to do COVIS after the button switch is executed?
        if config.isCOVISEnabled && trial <= PRE_LEARNING_TRIALS + LEARNING_TRIALS + POST_LEARNING_TRIALS
            config.COVISRules = config.COVISRules.processRuleAttempt(accuracy(trial) == 1);
        end

        %% Print data to console
        if not(SUPPRESS_UI)
            consoleprogressbar('TRIALS COMPLETED', trial, TRIALS);
        end
    end
    
    %% ========================================= %%
    %%%%%%%%%% OPTIMIZATION CALCULATIONS %%%%%%%%%%
    %  =========================================  %
    opt_val_1 = 0;
    opt_val_2 = zeros(4,4);
    % Calculate Sum of Squared Errors of Prediction (SSE)
    if OPTIMIZATION_CALC && isa(config, 'ModelConfigButtonSwitch')
        FMRI_META = config.meta.optimization;
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
    %% =============================== %%
    %%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%
    %  ===============================  %
    if not(SUPPRESS_UI)
        displayautoresults(config, Neuron.TAU, n, RBF, config.BORDER_SIZE, VISUAL, TRIALS, ...
            PRE_LEARNING_TRIALS, LEARNING_TRIALS, POST_LEARNING_TRIALS, accuracy, ...
            PFC, PMC, MC, PFC_A, PFC_B, PMC_A, PMC_B, MC_A, MC_B, Driv_PFC, CN, GP, MDN_A, MDN_B, AC_A, AC_B, ...
            chosen_rule, y_coordinates, x_coordinates...
        );
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
        neuron_id = double(trapz(neuron_1) < trapz(neuron_2)) + 1;
        latency = length(neuron_1);
    end
end

function [neuron_id] = determinecorrectneuron(x, y, rule)
    % If x and y are found in "B", the boolean will evaluate to true.
    equalsNeuron2 = any(x == rule.B_X) && any(y == rule.B_Y);
    % Add one to the result, since neuron IDs are 1-indexed, and cast the result
    % to a double for code-generation compatibility.
    neuron_id = double(equalsNeuron2 + 1);
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