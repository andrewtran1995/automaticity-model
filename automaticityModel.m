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

function automaticityModel()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load visual stimulus matrix
    % load('randomVisualInput.mat');
    load('maddoxVisualInput.mat');
    r_x_vals = maddoxVisualInput(:, 1);
    r_y_vals = maddoxVisualInput(:, 2);
    r_groups = maddoxVisualInput(:, 3);

    %% Initialize/configure constants (though some data structure specific constants are initialized below)
    % Experiment parameters
    n = 1000;                 % Time period for one trial (in milliseconds)
    TAU = 1;
    TRIALS = 300;              % Number of trials in automaticity experiment
    GRID_SIZE = 140;          % Length of side of square grid for visual input; should always be an even number
    BORDER_SIZE = 20;         % Width of border used to pad the grid such that visual stimulus on the edge still has an appropriate effect
    LAMBDA = 20;              % Lambda Value
    W_MAX = 10;              % maximum possible weight for Hebbian Synapses
    DECISION_PT = 4;          % Integral value which determines which PMC neuron acts on a visual input
    INIT_PMC_WEIGHT = 0.08;   % Initial weight for PMC neurons

    % Quantity of Visual Stimulus
    Visual = struct( ...
        'stim', 50 ...
    );

    % Radial Basis Function
    RBF = struct( ...
        'RADIUS', 2, ...
        'rbv', zeros(GRID_SIZE) ...
    );

    % PFC scaling information
    PFC = struct( ...
        'V_SCALE', 1, ... % scaling factor for visual input into PFC neurons
        'W_LI', 2 ...       % lateral inhibition between PFC A / PFC B
    );

    % PMC scaling information
    PMC = struct( ...
        'V_SCALE', 1, ...     % can use to scale PMC visual input value if it comes out way too high
        'W_LI', 2 ...         % lateral inhibition between PMC A / PMC B
    );

    % Hebbian Constants (determine the subtle attributes of learning at the
    % Hebbian synapses)
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
        'W_OUT', 1, ...
        'out', zeros(1,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0, ...
        'weights', INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,TRIALS+1), ...
        'weights_avg', zeros(1,TRIALS) ...
    );

    PMC_B = struct( ...
        'W_OUT', 1, ...
        'out', zeros(1,n), ...
        'spikes', 0, ...
        'v', RSN.rv*ones(1,n), ...
        'u', zeros(1,n), ...
        'pos_volt', zeros(1,n), ...
        'v_stim', 0, ...
        'weights', INIT_PMC_WEIGHT*ones(GRID_SIZE,GRID_SIZE,TRIALS+1), ...
        'weights_avg', zeros(1,TRIALS) ...
    );

    % Matrix to store information about which matrix responds during a trial
    Reaction_Matrix = zeros(TRIALS, 3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% CALCULATIONS %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Precompute RBF.rbv matrices
%     parfor j=1:TRIALS
%         
%     end
    
    %% Learning trials
    trial_number = 0;

    for j=1:TRIALS
        %% Initialize appropriate variables for each loop
        trial_number = trial_number + 1;    % track number of current trial

        % variables tracking spiking rate in each neuron
        PMC_A.spikes = 0;
        PMC_B.spikes = 0;
        PFC_A.spikes = 0;
        PFC_B.spikes = 0;

        % variables keep track of positive voltage values for calculation of
        % integral (for Hebbian learning equation)
        PFC_A.pos_volt(:) = 0;
        PFC_B.pos_volt(:) = 0;
        PMC_A.pos_volt(:) = 0;
        PMC_B.pos_volt(:) = 0;

        % Re-initialize Neuron Voltage Matrices (set all v to RSN.rv; all u to 0)
        PFC_A.v(:) = RSN.rv;     PFC_A.u(:) = 0;
        PFC_B.v(:) = RSN.rv;     PFC_B.u(:) = 0;

        PMC_A.v(:) = RSN.rv;     PMC_A.u(:) = 0;
        PMC_B.v(:) = RSN.rv;     PMC_B.u(:) = 0;

        % Re-initialize Neuron Output Matrices
        PFC_A.out(:) = 0;
        PFC_B.out(:) = 0;

        PMC_A.out(:) = 0;
        PMC_B.out(:) = 0;

        % Determine visual stimulus in range [1, GRID_SIZE] to pick
        % random gabor for each trial, padded with the BORDER_SIZE such
        % that the visual stimulus is accounted for properly
        r_y = r_y_vals(j) + BORDER_SIZE;
        r_x = r_x_vals(j) + BORDER_SIZE;
        r_group = r_groups(j);

        %% Radial Basis Function (RBF) Implementation
        % Use temp variables x and y to iterate the entire grid, calculating visual input
        % Note: Recall that r_y corresponds to rows, and r_x corresponds to columns
        PFC_A.v_stim = 0;
        PFC_B.v_stim = 0;
        for y=1:GRID_SIZE
            for x=1:GRID_SIZE
                distance = sqrt(((r_y-y)^2) + ((r_x-x)^2));
                % Initialize RBF.rbv matrix while iterating to save computation time
                RBF.rbv(y,x) = exp(-(distance/RBF.RADIUS))*Visual.stim;
                % if r_x <= GRID_SIZE/2, must be on "left" side for PFC A
                if (x <= GRID_SIZE/2)
                    PFC_A.v_stim = PFC_A.v_stim + RBF.rbv(y,x);
                % else, r_x must be > GRID_SIZE/2 and is on the "right" side for PFC B
                else
                    PFC_B.v_stim = PFC_B.v_stim + RBF.rbv(y,x);
                end
            end
        end
        % Scale PFC_A.v_stim and PFC_B.v_stim to prevent them from becoming too large
        PFC_A.v_stim = PFC_A.v_stim * PFC.V_SCALE;
        PFC_B.v_stim = PFC_B.v_stim * PFC.V_SCALE;
        % Iterate through entire grid, calculating PMC_A.v_stim (matrix or value?)
        PMC_A.v_stim = 0;
        PMC_B.v_stim = 0;
        for y=1:GRID_SIZE
            for x=1:GRID_SIZE
                PMC_A.v_stim = PMC_A.v_stim + RBF.rbv(y,x)*PMC_A.weights(y,x,j);
                PMC_B.v_stim = PMC_B.v_stim + RBF.rbv(y,x)*PMC_B.weights(y,x,j);
            end
        end
        % Scale PMC_A.v_stim and PMC_B.v_stim to prevent them from becoming too large
        PMC_A.v_stim = PMC_A.v_stim * PMC.V_SCALE;
        PMC_B.v_stim = PMC_B.v_stim * PMC.V_SCALE;

        %% Individual Time Trial
        for i=1:n-1
            % Neuron Equations
            % PFC A Neuron
            PFC_A.v(i+1)=(PFC_A.v(i) + TAU*(RSN.k*(PFC_A.v(i)-RSN.rv)*(PFC_A.v(i)-RSN.vt)-PFC_A.u(i)+ RSN.E + PFC_A.v_stim + (PMC_A.W_OUT*PMC_A.out(i)) - PFC.W_LI*PFC_B.out(i))/RSN.C);   % no noise
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
            PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out(i)) - PFC.W_LI*PFC_A.out(i))/RSN.C);   % no noise
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
            PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + PMC_A.v_stim + (PFC_A.W_OUT*PFC_A.out(i)) - PMC.W_LI*PMC_B.out(i) )/RSN.C);   % no noise
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
            PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + PMC_B.v_stim + (PFC_B.W_OUT*PFC_B.out(i)) - PMC.W_LI*PMC_A.out(i) )/RSN.C);   % no noise
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

            % Record PFC v value if positive. Else, do nothing.
            % For computation of integral
            if PFC_A.v(i) > 0
                PFC_A.pos_volt(i) = PFC_A.v(i);
            end

            if PFC_B.v(i) > 0
                PFC_B.pos_volt(i) = PFC_B.v(i);
            end

            % this code creates a positive number matrix for the PMC voltages
            % only the positive values are retained and all other values
            % are converted to zero (for computation of integral)
            if PMC_A.v(i) > 0
                PMC_A.pos_volt(i) = PMC_A.v(i);
            end

            if PMC_B.v(i) > 0
                PMC_B.pos_volt(i) = PMC_B.v(i);
            end

        end

        %% Determine decision neuron and reaction time
        for i=1:n
            % If PMC_A meets the decision point sooner, indicate it in the
            % first column with a '0'
            if trapz(PMC_A.out(1:i)) >= DECISION_PT
                Reaction_Matrix(j,:) = [0, i, r_group];
                break;
            % Else, indicate PMC_B with '1'
            elseif trapz(PMC_B.out(1:i)) >= DECISION_PT
                Reaction_Matrix(j,:) = [1, i, r_group];
                break;
            else
                continue;
            end
        end

        %% Calculation of Hebbian Weight for PMC_A
        % Visual input to PMC_A neuron (presynaptic)
        integral_visinputA   = trapz(PFC_A.pos_volt);

        % Activation of PMC_A neuron   (post-synaptic)
        integral_PMCAvoltage = trapz(PMC_A.pos_volt);

        % Ensure g(t)-1 and g(2)-2 are never less than zero
        g_t_1_A = max(0, integral_PMCAvoltage - Hebbian.NMDA);
        g_t_2_A = max(0, Hebbian.NMDA - integral_PMCAvoltage - Hebbian.AMPA);

        % Determine new weights of visual PMC_A synapses
        PMC_A.weights(:,:,j+1) = PMC_A.weights(:,:,j) + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputA)*g_t_1_A.*(W_MAX - PMC_A.weights(:,:,j)) - (Hebbian.anti_heb)*(integral_visinputA)*g_t_2_A.*PMC_A.weights(:,:,j));

        % Limit values of PMC_A.weights to be in range [0,W_MAX]
        PMC_A.weights(:,:,j+1) = max(PMC_A.weights(:,:,j+1), 0);
        PMC_A.weights(:,:,j+1) = min(PMC_A.weights(:,:,j+1), W_MAX);

        %% Calculation of Hebbian Weight for PMC_B
        % Visual input to PMC_B neuron (presynaptic)
        integral_visinputB   = trapz(PFC_B.pos_volt);

        % Activation of PMC_B neuron   (post-synaptic)
        integral_PMCBvoltage = trapz(PMC_B.pos_volt);

        % Ensures g(t)-1 and g(2)-2 are never less than zero
        g_t_1_B = max(0, integral_PMCBvoltage - Hebbian.NMDA);
        g_t_2_B = max(0, Hebbian.NMDA - integral_PMCBvoltage - Hebbian.AMPA);

        % Determine new weights of visual PMC_B synapses
        PMC_B.weights(:,:,j+1) = PMC_B.weights(:,:,j) + RBF.rbv(:,:).*((Hebbian.heb_coef)*(integral_visinputB)*g_t_1_B.*(W_MAX - PMC_B.weights(:,:,j)) - (Hebbian.anti_heb)*(integral_visinputB)*g_t_2_B.*PMC_B.weights(:,:,j));

        % Limit values of PMC_A.weights to be in range [0,W_MAX]
        PMC_B.weights(:,:,j+1) = max(PMC_B.weights(:,:,j+1), 0);
        PMC_B.weights(:,:,j+1) = min(PMC_B.weights(:,:,j+1), W_MAX);

        % Record average weight for PMC_A and PMC_B
        PMC_A.weights_avg(j) = mean(mean(PMC_A.weights(:,:,j+1)));
        PMC_B.weights_avg(j) = mean(mean(PMC_B.weights(:,:,j+1)));

        %% Print data to console
        fprintf('~~~ TRIAL #: %d ~~~\n', trial_number);
        fprintf('r_y: %d\n', r_y);
        fprintf('r_x: %d\n', r_x);
        fprintf('PFC_A.v_stim: %d\n', PFC_A.v_stim);
        fprintf('PFC_B.v_stim: %d\n', PFC_B.v_stim);
        fprintf('PMC_A.v_stim: %d\n', PMC_A.v_stim);
        fprintf('PMC_B.v_stim: %d\n', PMC_B.v_stim);
%         fprintf('g_t_1_A: %d\n', g_t_1_A);
%         fprintf('g_t_2_A: %d\n', g_t_2_A);
%         fprintf('g_t_1_B: %d\n', g_t_1_B);
%         fprintf('g_t_2_B: %d\n', g_t_2_B);

    end

    %% Determine decision neuron and reaction time
%     parfor j=1:TRIALS
%         for i=1:n
%             % If PMC_A meets the decision point sooner, indicate it in the
%             % first column with a '0'
%             if trapz(PMC_A.out(1:i)) >= DECISION_PT
%                 Reaction_Matrix(j,:) = [0, i, r_group];
%                 break;
%             % Else, indicate PMC_B with '1'
%             elseif trapz(PMC_B.out(1:i)) >= DECISION_PT
%                 Reaction_Matrix(j,:) = [1, i, r_group];
%                 break;
%             else
%                 continue;
%             end
%         end
%     end
    
    % Delete first matrix (initialization matrix) of PMC_A.weights and PMC_B.weights
    % so that the trial number matches the index
    PMC_A.weights(:,:,1) = [];
    PMC_B.weights(:,:,1) = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Figure 1 - neuron information from last trial or throughout trials
    figure;

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
    surf(RBF.rbv(BORDER_SIZE:end-BORDER_SIZE-1,BORDER_SIZE:end-BORDER_SIZE-1,:));
    title(sprintf('Stimulus: (%d,%d); Weight: %d', r_y, r_x, Visual.stim));
%     axis([0 100 0 100]);

    subplot(rows,columns,10);
    x_axis = linspace(1, TRIALS, TRIALS);
    plot(x_axis, PMC_A.weights_avg, 'r', x_axis, PMC_B.weights_avg, 'b');
    legend('PMC_A', 'PMC_B', 'Location', 'southeast');
    title('PMC_A & PMC_B Weight Average');

    subplot(rows,columns,11);
    x_axis = linspace(1, TRIALS, TRIALS);
    PMC_A_Rx = Reaction_Matrix(:,1) == 0;
    PMC_B_Rx = ~PMC_A_Rx;
    scatter(x_axis(PMC_A_Rx), Reaction_Matrix(PMC_A_Rx,2), 10, 'r', 'filled');
    hold on;
    scatter(x_axis(PMC_B_Rx), Reaction_Matrix(PMC_B_Rx,2), 10, 'b', 'filled');
    legend('PMC_A', 'PMC_B');
    title('PMC_A & PMC_B Reaction Time');
    

    %% Figure 2
    figure;

    % Force slider to integer/discrete value:
    % https://www.mathworks.com/matlabcentral/answers/45769-forcing-slider-
    % values-to-round-to-a-valid-number
    rows = 1;
    columns = 2;

    PMC_A_trial_num = 1;
    PMC_A_no_border = PMC_A.weights(BORDER_SIZE:end-BORDER_SIZE, ...
                                    BORDER_SIZE:end-BORDER_SIZE, :);

    subplot(rows,columns,1);
    data3 = PMC_A_no_border(:,:,PMC_A_trial_num);
    colormap('hot');
    imagesc(data3);
    colorbar;
    title(sprintf('PMC_A Synaptic Heatmap, Trial %d\n', PMC_A_trial_num));
    slider_PMC_A = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', TRIALS, ...
        'Value', 1, ...
        'Position', [100 50 300 20]);
    set(slider_PMC_A, 'Callback', {@synaptic_slider_callback, 1, PMC_A_no_border, 'PMC_A'});

    PMC_B_trial_num = 1;
    PMC_B_no_border = PMC_B.weights(BORDER_SIZE:end-BORDER_SIZE, ...
                                    BORDER_SIZE:end-BORDER_SIZE, :);

    subplot(rows,columns,2);
    data4 = PMC_B_no_border(:,:,PMC_B_trial_num);
    colormap('hot');
    imagesc(data4);
    colorbar;
    title(sprintf('PMC_B Synaptic Heatmap, Trial %d\n', PMC_B_trial_num));
    slider_PMC_B = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', TRIALS, ...
        'Value', 1, ...
        'Position', [500 50 300 20]);
    set(slider_PMC_B, 'Callback', {@synaptic_slider_callback, 2, PMC_B_no_border, 'PMC_B'});
    
    %% Figure 3
    figure;
    
    PMC_S = Reaction_Matrix(:,3) == 'S';
    PMC_M = Reaction_Matrix(:,3) == 'M';
    PMC_L = Reaction_Matrix(:,3) == 'L';
    p1 = cdfplot(Reaction_Matrix(PMC_S, 2));
    set(p1, 'Color', 'r');
    hold on;
    p2 = cdfplot(Reaction_Matrix(PMC_M, 2));
    set(p2, 'Color', 'b');
    hold on;
    p3 = cdfplot(Reaction_Matrix(PMC_L, 2));
    set(p3, 'Color', 'g');
    legend('S', 'M', 'L', 'Location', 'southeast');
    title('CDF of RT by Grouping');
    
    %% Starts debug mode, allowing variables to be observed before the
    % function ends
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
