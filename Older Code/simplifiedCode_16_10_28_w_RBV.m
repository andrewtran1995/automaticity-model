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
clearvars;

% Experiement parameters
n = 1000;                   % Time period for one trial (in milliseconds)
TAU = 1;
TRIALS = 5;                 % Number of trials in automaticity experiment
GRID_SIZE = 2;              % Length of side of square grid for visual input; should always be an even number
LAMBDA = 20;                % Lambda Value

% Quantity of Visual Stimulus
Visual = struct( ...
    'stim', 50, ...
    'input_matrix', 50*ones(1,n) ...
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
    'W_OUT', 2, ...
    'out', zeros(1,n), ...
    'spikes', 0, ...
    'v', RSN.rv*ones(1,n), ...
    'u', zeros(1,n), ...
    'v_stim', 0 ...
);

PFC_B = struct( ...
    'W_OUT', 2, ...
    'out', zeros(1,n), ...
    'spikes', 0, ...
    'v', RSN.rv*ones(1,n), ...
    'u', zeros(1,n), ...
    'v_stim', 0 ...
);

PMC_A = struct( ...
    'W_OUT', 1, ...
    'out', zeros(1,n), ...
    'spikes', 0, ...
    'v', RSN.rv*ones(1,n), ...
    'u', zeros(1,n), ...
    'weights', 0.1*ones(GRID_SIZE,GRID_SIZE,TRIALS) ...
);

PMC_B = struct( ...
    'W_OUT', 1, ...
    'out', zeros(1,n), ...
    'spikes', 0, ...
    'v', RSN.rv*ones(1,n), ...
    'u', zeros(1,n), ...
    'weights', 0.1*ones(GRID_SIZE,GRID_SIZE,TRIALS) ...
);

% Initializing PMC Synaptic weight matrix (all start with value of 0.1)
PMCA_weights = 0.1*ones(GRID_SIZE);
PMCB_weights = 0.1*ones(GRID_SIZE);

% Multi-dimensional array (contains important matrices)
full_matrix = cat(3, PMCA_weights, PMCB_weights);

% PFC scaling information
PFC = struct( ...
    'V_SCALE', 1/5, ... % scaling factor for visual input into PFC neurons
    'W_LI', 1 ...       % lateral inhibition between PFC A / PFC B
);

% PMC scaling information
PMC = struct( ...
    'V_SCALE', 1, ...     % can use to scale PMC visual input value if it comes out way too high
    'W_LI', 1 ...         % lateral inhibition between PMC A / PMC B
);

% maximum possible weight for Hebbian Synapses
w_max = 100;

% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses);
Hebbian = struct( ...
    'heb_coef', 0.000000005, ...
    'anti_heb', 10, ...
    'NMDA', 200, ...
    'AMPA', 50 ...
);

% Radial Basis Function
RBF = struct( ...
    'RADIUS', 10, ...
    'matrix', zeros(GRID_SIZE) ...
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
    % pos_PFCAvolt_values = [0]; % seems unused?
    % pos_PFCBvolt_values = [0];
    pos_PMCAvolt_values = [0];
    pos_PMCBvolt_values = [0];

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
    r1 = randi([1 GRID_SIZE],1,1);      % this function generates new random integer between 1 and 2 everytime the function runs
    r2 = randi([1 GRID_SIZE],1,1);      % the effect is to pick a new random gabor for each trial

    % if r2 <= GRID_SIZE/2, must be on "left" side for PFC A
    if (r2 <= GRID_SIZE/2)
        PFC_A.v_stim = 50;
        PFC_B.v_stim = 0;
    % else, r2 must be > GRID_SIZE/2 and is on the "right" side for PFC B
    else
        PFC_A.v_stim = 0;
        PFC_B.v_stim = 50;
    end

    % Radial Basis Function (RBF) Implementation TODO
    % Use temp variables x and y to go through the entire grid, calculating visual input
    % Note: Recall that r1 corresponds to rows, and r2 corresponds to columns
    % PFC_A.v_stim = 0;
    % PFC_B.v_stim = 0;
    % for y=1:GRID_SIZE
    %     for x=1:GRID_SIZE
    %         distance = (((r1-y)^2) + ((r2-x)^2))^(1/2);
    %         rbv = exp(-(distance/RBF.RADIUS));
    %         % if r2 <= GRID_SIZE/2, must be on "left" side for PFC A
    %         if (r2 <= GRID_SIZE/2)
    %             PFC_A.v_stim = PFC_A.v_stim + rbv*Visual.stim;
    %         % else, r2 must be > GRID_SIZE/2 and is on the "right" side for PFC B
    %         else
    %             PFC_B.v_stim = PFC_B.v_stim + rbv*Visual.stim;
    %         end
    %     end
    % end
    % PFC_A.v_stim = PFC_A.v_stim * PFC.V_SCALE;
    % PFC_B.v_stim = PFC_B.v_stim * PFC.V_SCALE;

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
        PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + (Visual.stim*full_matrix(r1,r2,1)) + (PFC_A.W_OUT*PFC_A.out(i)) - PMC.W_LI*PMC_B.out(i) )/RSN.C);   % no noise
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
        PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + (Visual.stim*full_matrix(r1,r2,2)) + (PFC_B.W_OUT*PFC_B.out(i)) - PMC.W_LI*PMC_A.out(i) )/RSN.C);   % no noise
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

        % this code creates a positive number matrix for the PMC voltages
        % only the positive values are retained and all other values
        % are converted to zero (for computation of integral)
        if (PMC_A.v(i) > 0)
            pos_PMCAvolt_values = [pos_PMCAvolt_values;PMC_A.v(i)];
        else
            pos_PMCAvolt_values = [pos_PMCAvolt_values;0];
        end

        if (PMC_B.v(i) > 0)
            pos_PMCBvolt_values = [pos_PMCBvolt_values;PMC_B.v(i)];
        else
            pos_PMCBvolt_values = [pos_PMCBvolt_values;0];
        end

    end


    %% Calculation of Hebbian Weight for PMC_A
    % visual input to PMCA neuron (presynaptic)

    PMCA_vis_stim = full_matrix(r1,r2,1)*Visual.input_matrix;                           % visual input to PMCA neuron (presynaptic)

    integral_visinputA   = trapz(PMCA_vis_stim);

    integral_PMCAvoltage = trapz(pos_PMCAvolt_values);                               % activation of PMCA neuron   (post-synaptic)

    g_t_1_A = (integral_PMCAvoltage)-(Hebbian.NMDA);
    g_t_2_A = [Hebbian.NMDA - (integral_PMCAvoltage) - Hebbian.AMPA];


    % Ensures g(t)-1 and g(2)-2 are never less than zero

    if (g_t_1_A <= 0)
        g_t_1_A = 0;
    end

    if (g_t_2_A <= 0)
        g_t_2_A = 0;
    end

    % Equation determining weights at Visual-PMCA Synapses

    full_matrix(r1,r2,1) = full_matrix(r1,r2,1) + [[(Hebbian.heb_coef)*(integral_visinputA)*[g_t_1_A]*(w_max - full_matrix(r1,r2,1))] - [(Hebbian.anti_heb)*(integral_visinputA)*[g_t_2_A]*full_matrix(r1,r2,1)]];

    if (full_matrix(r1,r2,1) < 0)
        full_matrix(r1,r2,1) = 0;
    end

    if (full_matrix(r1,r2,1) > 100)
        full_matrix(r1,r2,1) = 100;
    end

    PMC_A.weights(:,:,j) = full_matrix(:,:,1);

    %% Calculation of Hebbian Weight for PMC_B

    PMCB_vis_stim = full_matrix(r1,r2,2)*Visual.input_matrix;                             % presynaptic

    integral_visinputB   = trapz(PMCB_vis_stim);

    integral_PMCBvoltage = trapz(pos_PMCBvolt_values);                                 % post-synaptic

    g_t_1_B = (integral_PMCBvoltage)-(Hebbian.NMDA);
    g_t_2_B = [Hebbian.NMDA - (integral_PMCBvoltage) - Hebbian.AMPA];


    % Ensures g(t)-1 and g(2)-2 are never less than zero

    if (g_t_1_B <= 0)
        g_t_1_B = 0;
    end

    if (g_t_2_B <= 0)
        g_t_2_B = 0;
    end

    % Equation determining weights at Visual-PMCB Synapses


    full_matrix(r1,r2,2) = full_matrix(r1,r2,2) + [[(Hebbian.heb_coef)*(integral_visinputB)*[g_t_1_B]*(w_max - full_matrix(r1,r2,2))] - [(Hebbian.anti_heb)*(integral_visinputB)*[g_t_2_B]*full_matrix(r1,r2,2)]];

    if (full_matrix(r1,r2,2) < 0)
        full_matrix(r1,r2,2) = 0;
    end

    if (full_matrix(r1,r2,2) > 100)
        full_matrix(r1,r2,2) = 100;
    end

    PMC_B.weights(:,:,j) = full_matrix(:,:,2);

    % Appends new Hebbian Weight value to heb_weight Matrix
    % this allows me to plot the Hebbian Values at the End

    % one weight value per trial


    %% Print data to console
    fprintf('~~~ TRIAL #: %d ~~~\n', trial_number);
    fprintf('r1: %d\n', r1);
    fprintf('r2: %d\n', r2);
    fprintf('Visual.stim: %d\n', Visual.stim);
    disp('full_matrix(r1,r2,1)');
    disp(full_matrix(:,:,1));
    disp('full_matrix(r1,r2,2)');
    disp(full_matrix(:,:,2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots of Final Voltages, Outputs, and Change of Hebbian Synapse over time

% creating gif of synaptic weight heatmap

%{

n = 1:0.5:5;
nImages = TRIALS;

figure;
for idx = 1:TRIALS
    y = x.^n(idx);
    plot(x,y,'LineWidth',3)
    title(['y = x^n,  n = ' num2str( n(idx)) ])
    drawnow
    frame = getframe(1);
    im{idx} = frame2im(frame);
end
close;

figure;
for idx = 1:nImages
    subplot(3,3,idx)
    imshow(im{idx});
end

%}

% this code generates the spiking diagrams and the heatmaps
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

% subplot(rows,columns,11);
% plot(TAU*(1:TRIALS), PMCA_weights_output);
% axis([0 TRIALS -1 12]);
%
% subplot(rows,columns,12);
% plot(TAU*(1:TRIALS), PMCB_weights_output);
% axis([0 TRIALS -1 12]);


% report important statistics

%{

disp('r1')
disp(r1)
disp('r2')
disp(r2)

disp('Synaptic Strength PMC_A')
disp(full_matrix(r1,r2,2))
disp('Synaptic Strength PMC_B')
disp(full_matrix(r1,r2,3))
disp('PMC_A Spikes')
disp(PMC_A.spikes)
disp('PMC_B Spikes')
disp(PMC_B.spikes)
disp('PFC_A Spikes')
disp(PFC_A.spikes)
disp('PFC_B Spikes')
disp(PFC_B.spikes)

%}
