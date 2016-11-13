%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Model of Automaticity in Rule Based Learning %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% NOTES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This version of the RB_Automaticity model doesn't have the Radial Basis
% function (instead of a 100x100 visual grid it has a 2x2 grid)

% In general, variables that are written in all capital letters are meant
% to be constant values -- set once in the beginning of the program
% (variable initialization) and nowhere else

% Structures are used in the program by calling struct(...).
% These structures are meant to "bundle" or "group" together variables that
% have some kind of commonality, e.g., belonging to the same neuron, behavior,
% model, etc.
% The goal is to enforce readability by standardizing the names of grouped variables
% and making relationships between variables more apparent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call clearvars to remove all variables from currently active workspace
% clearvars;
% clf;

% Experiement parameters
n = 1000;                   % Time period for one trial (in milliseconds)
TAU = 1;
TRIALS = 100;                % Number of trials in automaticity experiment
GRID_SIZE = 100;            % Length of side of square grid for visual input; should always be an even number
LAMBDA = 20;                % Lambda Value
W_MAX = 10;                 % maximum possible weight for Hebbian Synapses

% Quantity of Visual Stimulus
Visual = struct( ...
    'stim', 50 ...
);

% Radial Basis Function
RBF = struct( ...
    'RADIUS', 1, ...
    'rbv', zeros(GRID_SIZE, GRID_SIZE) ...
);

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
    'v_stim', 0 ...
);

PFC_B = struct( ...
    'W_OUT', 9, ...
    'out', zeros(1,n), ...
    'spikes', 0, ...
    'v', RSN.rv*ones(1,n), ...
    'u', zeros(1,n), ...
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
    'v_stim', 0, ...
    'v_stim_matrix', zeros(GRID_SIZE,GRID_SIZE), ...
    'weights', 0.08*ones(GRID_SIZE,GRID_SIZE,TRIALS+1) ...
);

PMC_B = struct( ...
    'W_OUT', 1, ...
    'out', zeros(1,n), ...
    'spikes', 0, ...
    'v', RSN.rv*ones(1,n), ...
    'u', zeros(1,n), ...
    'v_stim', 0, ...
    'v_stim_matrix', zeros(GRID_SIZE,GRID_SIZE), ...
    'weights', 0.08*ones(GRID_SIZE,GRID_SIZE,TRIALS+1) ...
);

% Initializing PMC Synaptic weight matrix (all start with value of 0.1)
% PMCA_weights = 0.1*ones(GRID_SIZE);
% PMCB_weights = 0.1*ones(GRID_SIZE);

% Multi-dimensional array (contains important matrices)
% full_matrix = cat(3, PMCA_weights, PMCB_weights);

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

% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses);
Hebbian = struct( ...
    'heb_coef', 0.00000001, ...
    'anti_heb', 0.001, ...
    'NMDA', 380, ...
    'AMPA', 0 ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LOOP ON CALCULATIONS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Beginning of Learning Trials

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
    pos_PFCAvolt_values = zeros(1,n);
    pos_PFCBvolt_values = zeros(1,n);
    pos_PMCAvolt_values = zeros(1,n);
    pos_PMCBvolt_values = zeros(1,n);

    % Re-initialize Neuron Voltage Matrices
    PFC_A.v = RSN.rv*ones(1,n) ;     PFC_A.u = zeros(1,n);
    PFC_B.v = RSN.rv*ones(1,n) ;     PFC_B.u = zeros(1,n);

    PMC_A.v = RSN.rv*ones(1,n) ;     PMC_A.u = zeros(1,n);
    PMC_B.v = RSN.rv*ones(1,n) ;     PMC_B.u = zeros(1,n);

    % Re-initialize Neuron Output Matrices
    PFC_A.out = zeros(1,n);
    PFC_B.out = zeros(1,n);

    PMC_A.out = zeros(1,n);
    PMC_B.out = zeros(1,n);

    % Determine visual stimulus
%     r1 = randi([1 GRID_SIZE],1,1);      % this function generates new random integer between 1 and 2 everytime the function runs
%     r2 = randi([1 GRID_SIZE],1,1);      % the effect is to pick a new random gabor for each trial
    
    r1 = r_y_mat(j);
    r2 = r_x_mat(j);

    % % if r2 <= GRID_SIZE/2, must be on "left" side for PFC A
    % if (r2 <= GRID_SIZE/2)
    %     PFC_A.v_stim = 50;
    %     PFC_B.v_stim = 0;
    % % else, r2 must be > GRID_SIZE/2 and is on the "right" side for PFC B
    % else
    %     PFC_A.v_stim = 0;
    %     PFC_B.v_stim = 50;
    % end

    %% Radial Basis Function (RBF) Implementation
    % Use temp variables x and y to iterate the entire grid, calculating visual input
    % Note: Recall that r1 corresponds to rows, and r2 corresponds to columns
    PFC_A.v_stim = 0;
    PFC_B.v_stim = 0;
    for y=1:GRID_SIZE
        for x=1:GRID_SIZE
            distance = (((r1-y)^2) + ((r2-x)^2))^(1/2);
            % Initialize RBF.rbv matrix while iterating to save computation time
            RBF.rbv(y,x) = exp(-(distance/RBF.RADIUS))*Visual.stim;
            % if r2 <= GRID_SIZE/2, must be on "left" side for PFC A
            if (x <= GRID_SIZE/2)
                PFC_A.v_stim = PFC_A.v_stim + RBF.rbv(y,x);
            % else, r2 must be > GRID_SIZE/2 and is on the "right" side for PFC B
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
            PMC_A.v_stim_matrix(y,x) = RBF.rbv(y,x)*PMC_A.weights(y,x,j);
            PMC_A.v_stim = PMC_A.v_stim + RBF.rbv(y,x)*PMC_A.weights(y,x,j);
            PMC_B.v_stim_matrix(y,x) = RBF.rbv(y,x)*PMC_B.weights(y,x,j);
            PMC_B.v_stim = PMC_B.v_stim + RBF.rbv(y,x)*PMC_B.weights(y,x,j);
        end
    end
    % Scale PMC_A.v_stim and PMC_B.v_stim to prevent them from becoming too large
    PMC_A.v_stim = PMC_A.v_stim * PMC.V_SCALE;
    PMC_B.v_stim = PMC_B.v_stim * PMC.V_SCALE;

    %% Beginning of Individual Trial
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
            pos_PFCAvolt_values(i) = PFC_A.v(i);
        end

        if PFC_B.v(i) > 0
            pos_PFCBvolt_values(i) = PFC_B.v(i);
        end

        % this code creates a positive number matrix for the PMC voltages
        % only the positive values are retained and all other values
        % are converted to zero (for computation of integral)
        if (PMC_A.v(i) > 0)
            pos_PMCAvolt_values(i) = PMC_A.v(i);
        end

        if (PMC_B.v(i) > 0)
            pos_PMCBvolt_values(i) = PMC_B.v(i);
        end

    end


    %% Calculation of Hebbian Weight for PMC_A
    % Visual input to PMC_A neuron (presynaptic)
    integral_visinputA   = trapz(pos_PFCAvolt_values);

    % Activation of PMC_A neuron   (post-synaptic)
    integral_PMCAvoltage = trapz(pos_PMCAvolt_values);

    % Ensure g(t)-1 and g(2)-2 are never less than zero
    g_t_1_A = max(0, integral_PMCAvoltage - Hebbian.NMDA);
    g_t_2_A = max(0, Hebbian.NMDA - integral_PMCAvoltage - Hebbian.AMPA);

    % Iterate through grid to determine new weights of visual PMC_A synapses
    for y=1:GRID_SIZE
        for x=1:GRID_SIZE
            PMC_A.weights(y,x,j+1) = PMC_A.weights(y,x,j) + RBF.rbv(y,x)*(Hebbian.heb_coef)*(integral_visinputA)*g_t_1_A*(W_MAX - PMC_A.weights(y,x,j)) - RBF.rbv(y,x)*(Hebbian.anti_heb)*(integral_visinputA)*g_t_2_A*PMC_A.weights(y,x,j);
        end
    end
    % Limit values of PMC_A.weights to be in range [0,100]
    PMC_A.weights(:,:,j+1) = max(PMC_A.weights(:,:,j+1), 0);
    PMC_A.weights(:,:,j+1) = min(PMC_A.weights(:,:,j+1), 10);

    % PMC_A.weights(:,:,j) = full_matrix(:,:,1);

    %% Calculation of Hebbian Weight for PMC_B
    % Visual input to PMC_B neuron (presynaptic)
    integral_visinputB   = trapz(pos_PFCBvolt_values);

    % Activation of PMC_B neuron   (post-synaptic)
    integral_PMCBvoltage = trapz(pos_PMCBvolt_values);

    % Ensures g(t)-1 and g(2)-2 are never less than zero
    g_t_1_B = max(0, integral_PMCBvoltage - Hebbian.NMDA);
    g_t_2_B = max(0, Hebbian.NMDA - integral_PMCBvoltage - Hebbian.AMPA);

    % Iterate through grid to determine new weights of visual PMC_B synapses
    for y=1:GRID_SIZE
        for x=1:GRID_SIZE
            PMC_B.weights(y,x,j+1) = PMC_B.weights(y,x,j) + RBF.rbv(y,x)*(Hebbian.heb_coef)*(integral_visinputB)*g_t_1_B*(W_MAX - PMC_B.weights(y,x,j)) - RBF.rbv(y,x)*(Hebbian.anti_heb)*(integral_visinputB)*g_t_2_B*PMC_B.weights(y,x,j);
        end
    end
    % Limit values of PMC_A.weights to be in range [0,100]
    PMC_B.weights(:,:,j+1) = max(PMC_B.weights(:,:,j+1), 0);
    PMC_B.weights(:,:,j+1) = min(PMC_B.weights(:,:,j+1), 10);

    %% Print data to console
    fprintf('~~~ TRIAL #: %d ~~~\n', trial_number);
    
    
    fprintf('r1: %d\n', r1);
    fprintf('r2: %d\n', r2);
    fprintf('Visual.stim: %d\n', Visual.stim);
    fprintf('PFC_A.v_stim: %d\n', PFC_A.v_stim);
    fprintf('PFC_B.v_stim: %d\n', PFC_B.v_stim);
    fprintf('PMC_A.v_stim: %d\n', PMC_A.v_stim);
    fprintf('PMC_B.v_stim: %d\n', PMC_B.v_stim);
    fprintf('integral_PMCAvoltage: %d\n', integral_PMCAvoltage);
    fprintf('integral_PMCBvoltage: %d\n', integral_PMCBvoltage);
    
    
    % disp('PMC_A.weights(:,:,j+1)');
    % disp(PMC_A.weights(:,:,j+1));
    % disp('PMC_B.weights(:,:,j+1)');
    % disp(PMC_B.weights(:,:,j+1));

end

% Delete first matrix (initialization matrix) of PMC_A.weights and PMC_B.weights
% so that the trial number matches the index
PMC_A.weights(:,:,1) = [];
PMC_B.weights(:,:,1) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this code generates the spiking diagrams and the heatmaps
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
axis([0 n -1 20]);
title('PMC_A Neuron Output');

subplot(rows,columns,8);
plot(TAU*(1:n),PMC_B.out);
axis([0 n -1 20]);
title('PMC_B Neuron Output');

% subplot(rows,columns,9);
% data1 = full_matrix(:,:,1);
% colormap('hot');   % set colormap
% imagesc(data1);    % draw image and scale colormap to values range
% colorbar;          % show color scale
% title('PMC_A Synaptic Heatmap');
%
% subplot(rows,columns,10);
% data2 = full_matrix(:,:,2);
% colormap('hot');   % set colormap
% imagesc(data2);    % draw image and scale colormap to values range
% colorbar;          % show color scale
% title('PMC_B Synaptic Heatmap');

% Force slider to integer/discrete value: https://www.mathworks.com/matlabcentral/answers/45769-forcing-slider-values-to-round-to-a-valid-number

PMC_A_trial_num = 1;

subplot(rows,columns,9);
data3 = PMC_A.weights(:,:,PMC_A_trial_num);
colormap('hot');
imagesc(data3);
colorbar;
title(sprintf('PMC_A Synaptic Heatmap, Trial %d\n', PMC_A_trial_num));
slider_PMC_A = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', TRIALS, ...
    'Value', 1, ...
    'Position', [100 50 300 20]);
set(slider_PMC_A, 'Callback', 'subplot(rows,columns,9,''replace'');PMC_A_trial_num=round(get(slider_PMC_A,''value''));set(slider_PMC_A, ''Value'', PMC_A_trial_num);data3 = PMC_A.weights(:,:,PMC_A_trial_num);colormap(''hot'');imagesc(data3);colorbar;title(sprintf(''PMC_A Synaptic Heatmap, Trial %d\n'', PMC_A_trial_num));');

PMC_B_trial_num = 1;

subplot(rows,columns,10);
data4 = PMC_B.weights(:,:,PMC_B_trial_num);
colormap('hot');
imagesc(data4);
colorbar;
title(sprintf('PMC_B Synaptic Heatmap, Trial %d\n', PMC_B_trial_num));
slider_PMC_B = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', TRIALS, ...
    'Value', 1, ...
    'Position', [500 50 300 20]);
set(slider_PMC_B, 'Callback', 'subplot(rows,columns,10,''replace'');PMC_B_trial_num=round(get(slider_PMC_B,''value''));set(slider_PMC_B, ''Value'', PMC_B_trial_num);data4 = PMC_B.weights(:,:,PMC_B_trial_num);colormap(''hot'');imagesc(data4);colorbar;title(sprintf(''PMC_B Synaptic Heatmap, Trial %d\n'', PMC_B_trial_num));');

subplot(rows,columns,11);
surf(RBF.rbv(:,:,:));

% subplot(rows,columns,11);
% plot(TAU*(1:TRIALS), PMCA_weights_output);
% axis([0 TRIALS -1 12]);
%
% subplot(rows,columns,12);
% plot(TAU*(1:TRIALS), PMCB_weights_output);
% axis([0 TRIALS -1 12]);
