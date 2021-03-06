function [config, neurons] = automaticityModel(parameter_overrides, optional_params) %#codegen
    %% Pre-processing
    % Code-generation declarations.
    coder.extrinsic('struct2table','cell2struct','addpath','genpath','getmodelparams','getconstants','displayautoresults','dispResults','dispStimulus','dispWeightsWithSlider','dispCOVISLog','saveOutput');
    coder.varsize('chosen_rule');

    % Get model parameters.
    config = ModelConfigButtonSwitch();
    PARAMS = struct('PRE_LEARNING_TRIALS',0,'LEARNING_TRIALS',0,'POST_LEARNING_TRIALS',0,'W_MAX',0);
    PARAMS = getmodelparams(config);

    % Model/function behavior parameters (default values).
    VIS_INPUT_FROM_PARM   = 0;
    SUPPRESS_UI           = 0;
    
    if nargin >= 1
        % Determine if parameters are valid.
        if not(areparamsvalid(PARAMS))
            disp(struct2table(PARAMS));
            error('Parameters not valid.');
        end
        % Accept any parameter overrides from "arg_struct".
        PARAMS = absorbstruct(PARAMS, parameter_overrides);
    end

    % Override values with optional_parms if passed as an argument
    if nargin == 2
        if isfield(optional_params, 'VIS_INPUT_FROM_PARM')
            VIS_INPUT_FROM_PARM = optional_params.VIS_INPUT_FROM_PARM;
        end
    end
    
    % Load visual stimulus.
    if VIS_INPUT_FROM_PARM
        x_coords = optional_params.visualinput(:,1);
        y_coords = optional_params.visualinput(:,2);
        coord_groups = repmat(Category.NONE, length(x_coords), 1);
    else
        [x_coords, y_coords, coord_groups] = config.loadCoords();
    end

    % Initialize constants/variables.
    config = config.setTrials(PARAMS.PRE_LEARNING_TRIALS, PARAMS.LEARNING_TRIALS, PARAMS.POST_LEARNING_TRIALS);
    TRIALS = config.trials;
    W_MAX = PARAMS.W_MAX;
    RBF = config.RBF; % Create local copy of RBF for performance reasons.

    % Initialize COVIS model.
    chosen_rule = 1;
    correct_rule = 2;
    RULE = config.visual.RULES(chosen_rule);
    config = config.initCOVISRules(PARAMS, chosen_rule, correct_rule, config.trials);
    
    % Initialize neurons.
    PFC = struct( ...
        'reactions',   zeros(TRIALS,2), ...          % stores information about PFC neuron reactions during trial
        'activations', zeros(TRIALS,1) ...
    );
    PFC_Stroop_A = PFCStroopNeuron();
    PFC_Stroop_B = PFCStroopNeuron();
    PFC_A = PFCNeuron();
    PFC_B = PFCNeuron();

    PMC = struct( ...
        'reactions',   zeros(TRIALS,2), ...          % stores information about PMC neuron reactions during trial
        'alpha',       zeros(TRIALS,Neuron.n), ...   % PMC_A.out + PMC_B.out
        'activations', zeros(TRIALS,1) ...
    );
    PMC_A = PMCNeuron(config, W_MAX);
    PMC_B = PMCNeuron(config, W_MAX);

    MC = struct( ...
        'reactions',   zeros(TRIALS,2), ...
        'activations', zeros(TRIALS,1), ...
        'A_area',      zeros(TRIALS,1), ...
        'B_area',      zeros(TRIALS,1) ...
    );
    MC_A = MCNeuron(TRIALS);
    MC_B = MCNeuron(TRIALS);

    MDN = struct( ...
        'activations', zeros(TRIALS,1) ...
    );
    MDN_A = MDNNeuron();
    MDN_B = MDNNeuron();

    Driv_PFC = Driv_PFCNeuron();

    CN = CNNeuron(TRIALS);
    
    GP = GPNeuron(TRIALS);

    AC_A = ACNeuron();
    AC_B = ACNeuron();

    %% ============================ %%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %  ============================  %
    %% Learning trials
    for trial=1:config.trials
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
            [config, RULE] = config.chooseCOVISRule();
        end
        %% Button Switch if enabled and correct trials
        if config.shouldButtonSwitch()
            config = config.doButtonSwitch(PMC_A, PMC_B);
        end

        %% Determine visual stimulus in range [1, GRID_SIZE] to pick random gabor for each trial, padded with
        % the BORDER_SIZE such that the visual stimulus is accounted for properly
        config.visual.coord = Coord( ...
            x_coords(trial) + config.BORDER_SIZE + config.hasCriterialNoise * criterialnoise(), ...
            y_coords(trial) + config.BORDER_SIZE, ...
            coord_groups(trial) ...
        );
       
        %% Calculate visual stimulus effect using Radial Basis Function (RBF) implementation
        % Calculate RBF grid
        RBF = RBF.resolvestimulus(config.visual.coord);
        % Sum RBF values depending on rule to find PFC_A and PFC_B v_stim values
        % Note that stim matrices are row-major order (e.g., indexed by y, then x)
        PFC_A.v_stim = sum(sum(RBF.rbv(RULE(1).A_Y, RULE(1).A_X))) * PFC_A.V_SCALE;
        PFC_B.v_stim = sum(sum(RBF.rbv(RULE(1).B_Y, RULE(1).B_X))) * PFC_B.V_SCALE;
        % Scale RBF values by PMC_A and PMC_B weights to find respective v_stim values
        PMC_A.v_stim = sum(sum(RBF.rbv(:,:).*(PMC_A.weightsForTrial(config)))) * PMC_A.V_SCALE;
        PMC_B.v_stim = sum(sum(RBF.rbv(:,:).*(PMC_B.weightsForTrial(config)))) * PMC_B.V_SCALE;

        %% Individual Time Trial Loop (iterating through n)
        if config.isFROSTEnabled && config.hasStroopInterference
            for i=1:Neuron.n-1
                PFC_Stroop_A = PFC_Stroop_A.iterate(MDN_A);
                PFC_Stroop_B = PFC_Stroop_B.iterate(MDN_B);
                
                PFC_A = PFC_A.iterateWithStroopInterference(PFC_B, PFC_Stroop_A, PFC_Stroop_B);
                PFC_B = PFC_B.iterateWithStroopInterference(PFC_A, PFC_Stroop_B, PFC_Stroop_A);

                PMC_A = PMC_A.iterate(PMC_B, PFC_A);
                PMC_B = PMC_B.iterate(PMC_A, PFC_B);
                
                MC_A = MC_A.iterate(trial, MC_B, PMC_A, PMC_B);
                MC_B = MC_B.iterate(trial, MC_A, PMC_B, PMC_A);               

                Driv_PFC = Driv_PFC.iterate();

                CN = CN.iterate(Driv_PFC);
                
                GP = GP.iterate(CN);

                MDN_A = MDN_A.iterate(PFC_Stroop_A, GP);
                MDN_B = MDN_B.iterate(PFC_Stroop_B, GP);

                AC_A = AC_A.iterate(PFC_Stroop_A);
                AC_B = AC_B.iterate(PFC_Stroop_B);
            end
        elseif config.isFROSTEnabled
            for i=1:Neuron.n-1
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
            for i=1:Neuron.n-1
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
        [neuron_id_PFC, latency] = determine_reacting_neuron(PFC_A.out, PFC_B.out, PFCNeuron.RESPONSE_THRESHOLD);
        PFC.reactions(trial,:) = [neuron_id_PFC, latency];
        [neuron_id_PMC, latency] = determine_reacting_neuron(PMC_A.out, PMC_B.out, PMCNeuron.RESPONSE_THRESHOLD);
        PMC.reactions(trial,:) = [neuron_id_PMC, latency];
        [neuron_id_MC, latency] = determine_reacting_neuron(MC_A.out, MC_B.out, MCNeuron.RESPONSE_THRESHOLD);
        MC.reactions(trial,:) = [neuron_id_MC, latency];
        % Determine accuracy
        if config.isCOVISEnabled
            config.accuracy(trial) = determinecorrectneuron(config.visual.coord, config.visual.RULES(config.COVISRules.correct)) == neuron_id_MC;
        else
            config.accuracy(trial) = determinecorrectneuron(config.visual.coord, RULE(1)) == neuron_id_MC;
        end

        %% Weight change calculations
        if config.isLearningTrial
            PMC_A = PMC_A.doHebbianLearning(config, RBF.rbv.*PMC_A.weightsForTrial(config).*1000);
            PMC_B = PMC_B.doHebbianLearning(config, RBF.rbv.*PMC_B.weightsForTrial(config).*1000);
            
            if config.isMCLearningEnabled
                MC_A = MC_A.doHebbianLearning(config, PMC_A.integralPosVolt);
                MC_B = MC_B.doHebbianLearning(config, PMC_B.integralPosVolt);
            end
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
        
        % Record average weight for PMC_A and PMC_B.
        PMC_A.weights_avg(trial) = mean(mean(PMC_A.weights(:,:,config.weightIdx,config.correctRuleIdx)));
        PMC_B.weights_avg(trial) = mean(mean(PMC_B.weights(:,:,config.weightIdx,config.correctRuleIdx)));
        
        % Calculate COVIS weights.
        if config.isCOVISEnabled && trial <= config.preLearningTrials + config.learningTrials + config.postLearningTrials
            config.COVISRules = config.COVISRules.processRuleAttempt(config.accuracy(trial) == 1);
        end

        if not(SUPPRESS_UI)
            consoleprogressbar('TRIALS COMPLETED', trial, TRIALS);
        end
    end
    
    %% Post-processing
    % Re-assign local variables back into config.
    config.RBF = RBF;
    
    % Create neuron struct for ease of output.
    neurons = struct( ...
        'PFC', PFC, ...
        'PFC_A', PFC_A, ...
        'PFC_B', PFC_B, ...
        'PFC_Stroop_A', PFC_Stroop_A, ...
        'PFC_Stroop_B', PFC_Stroop_B, ...
        'PMC', PMC, ...
        'PMC_A', PMC_A, ...
        'PMC_B', PMC_B, ...
        'MC', MC, ...
        'MC_A', MC_A, ...
        'MC_B', MC_B, ...
        'Driv_PFC', Driv_PFC, ...
        'CN', CN, ...
        'GP', GP, ...
        'MDN_A', MDN_A, ...
        'MDN_B', MDN_B, ...
        'AC_A', AC_A, ...
        'AC_B', AC_B ...
    );

    % Save results to an output file.
    saveOutput(config,neurons);
    
    % Display results unless specified otherwise.
    if not(SUPPRESS_UI)
        dispResults(config, neurons);
    end
    return;
end

%% Helper functions
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

function [neuron_id] = determinecorrectneuron(coord, rule)
    x = coord.x;
    y = coord.y;
    % If x and y are found in "B", the boolean will evaluate to true.
    equalsNeuronB = any(x == rule.B_X) && any(y == rule.B_Y);
    % Add one to the result, since neuron IDs are 1-indexed, and cast the result
    % to a double for code-generation compatibility.
    neuron_id = double(equalsNeuronB + 1);
end
