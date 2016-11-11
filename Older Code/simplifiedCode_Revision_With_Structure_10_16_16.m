% Model of Automaticity in Rule Based Learning

% Could have class to hold RSN constants (must be included in separate .m class file)

% classdef RSN
%     properties
%         C = 100;
%     end
% end


% This version of the RB_Automaticity model doesn't have the Radial Basis
% function (instead of a 100x100 visual grid it has a 2x2 grid)

clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% VARIABLE INITIALIZATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Time period for one trial (in milliseconds)

n = 1000;
TAU = 1;

trials = 5;                         % number of trials in our automaticity experiment

tot_pfc_scale_factor = 1/5;          % scaling factor for visual input into PFC neurons


%{

% Total stimulus variable for PFCA and PFCB neurons

total_PFCA_stim = 0;
total_PFCB_stim = 0;

%}




% Quantity of Visual Stimulus

v_stim = 50;

vis_input_matrix = v_stim*ones(1,n);

% Initializing PMC Synaptic weight matrix (all start with value of 0.1)

PMCA_weights = 0.1*ones(2);
PMCB_weights = 0.1*ones(2);

% Multi-dimensional array (contains important matrices)

full_matrix = cat(3, PMCA_weights, PMCB_weights);    


%{

disp('row:')
disp(r1)
disp('column: ')
disp(r2)
disp('bar width value: ')                        
disp(full_matrix(r1,r2,1))
disp('frequency value: ')
disp(full_matrix(r1,r2,2))
disp('pFC weight')
disp(full_matrix(r1,r2,3))
disp('PMC weight')
disp(full_matrix(r1,r2,4))

%}


% PMC Scaling factor (I can use this to scale the PMC visual input value if it comes out way too high)

pmc_scaling_factor = 1;

% Neuron constants (set for a cortical regular spiking neuron) stored
% inside structure RSN
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

% RSN.C=100; RSN.rv=-60; RSN.vt=-40; RSN.k=0.7;       %RS -- regular spiking
% RSN.a=0.03; RSN.b=-2; RSN.c=-50; RSN.d=100;         %RS -- regular spiking
% RSN.vpeak=35; RSN.E=60;                             %RS -- regular spiking

% input weights (between PFC and PMC neurons)

w_IpmcA = 1;
w_IpfcA = 2;
w_IpmcB = 1;
w_IpfcB = 2;

% lateral inhibition connections between PFC/PFC and PMC/PMC

w_pFC_LI = 1;
w_PMC_LI = 1;


% maximum possible weight for Hebbian Synapses

w_max = 100;                           

% Lambda Value

lambda = 20; 

% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses);

heb_coef = 0.000000005;
anti_heb_const = 10;
nmda_const = 200;
ampa_const = 50;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LOOP ON CALCULATIONS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Beginning of Learning Trials

trial_number = 0;

for j=1:trials
    
    trial_number = trial_number + 1;    % tracks number of current trial
    
    
    % variables tracking spiking rate in each neuron
        
    PMCA_Spikes = 0;
    PMCB_Spikes = 0;
    PFCA_Spikes = 0;
    PFCB_Spikes = 0;

    % variables keep track of positive voltage values for calculation of
    % integral (for Hebbian learning equation)

    pos_PFCAvolt_values = [0];
    pos_PFCBvolt_values = [0];
    pos_PMCAvolt_values = [0]; 
    pos_PMCBvolt_values = [0]; 

    % Neuron Voltage Matrices

    v_PFC_A = RSN.rv*ones(1,n) ;     u_PFC_A = zeros(1,n);
    v_PFC_B = RSN.rv*ones(1,n) ;     u_PFC_B = zeros(1,n);
    v_PMC_A = RSN.rv*ones(1,n) ;     u_PMC_A = zeros(1,n);
    v_PMC_B = RSN.rv*ones(1,n) ;     u_PMC_B = zeros(1,n);

    % Neuron Output Matrices

    I_PFC_A = zeros(1,n);
    I_PFC_B = zeros(1,n);
    I_PMC_A = zeros(1,n);
    I_PMC_B = zeros(1,n);


    r1 = randi([1 2],1,1);                         % this function generates new random integer between 1 and 100 everytime the function runs
    r2 = randi([1 2],1,1);                         % the effect is to pick a new random gabor for each trial 


    if (r2 == 1)
        pfcA_stim = 50;
        pfcB_stim = 0;
    end

    if (r2 == 2)
        pfcA_stim = 0;
        pfcB_stim = 50;
    end
        
        
    % Beginning of Individual Trial

    for i=1:n-1
        
       
        % Neuron Equations   
        % PFC A Neuron

        v_PFC_A(i+1)=(v_PFC_A(i) + TAU*(RSN.k*(v_PFC_A(i)-RSN.rv)*(v_PFC_A(i)-RSN.vt)-u_PFC_A(i)+ RSN.E + pfcA_stim + (w_IpmcA*I_PMC_A(i)) - w_pFC_LI*I_PFC_B(i))/RSN.C);   % no noise
        u_PFC_A(i+1)=u_PFC_A(i)+TAU*RSN.a*(RSN.b*(v_PFC_A(i)-RSN.rv)-u_PFC_A(i));
        if v_PFC_A(i+1)>=RSN.vpeak;
            v_PFC_A(i)= RSN.vpeak;
            v_PFC_A(i+1)= RSN.c;
            u_PFC_A(i+1)= u_PFC_A(i+1)+ RSN.d;
        end 
        
        if (v_PFC_A(i) >= RSN.vpeak)
            PFCA_Spikes = PFCA_Spikes + 1;
            for j=i:n
               t= j-i;
               I_PFC_A(j)= I_PFC_A(j)+((t/lambda)*exp((lambda-t)/lambda));   
            end
        end
        
        % PFC B Neuron

        v_PFC_B(i+1)=(v_PFC_B(i) + TAU*(RSN.k*(v_PFC_B(i)-RSN.rv)*(v_PFC_B(i)-RSN.vt)-u_PFC_B(i)+ RSN.E + pfcB_stim + (w_IpmcB*I_PMC_B(i)) - w_pFC_LI*I_PFC_A(i))/RSN.C);   % no noise
        u_PFC_B(i+1)=u_PFC_B(i)+TAU*RSN.a*(RSN.b*(v_PFC_B(i)-RSN.rv)-u_PFC_B(i));
        if v_PFC_B(i+1)>=RSN.vpeak;
            v_PFC_B(i)= RSN.vpeak;
            v_PFC_B(i+1)= RSN.c;
            u_PFC_B(i+1)= u_PFC_B(i+1)+ RSN.d;
        end 
        
        if (v_PFC_B(i) >= RSN.vpeak)
            PFCB_Spikes = PFCB_Spikes + 1;
            for j=i:n
               t= j-i;
               I_PFC_B(j)= I_PFC_B(j)+((t/lambda)*exp((lambda-t)/lambda));   
            end
        end
        
        
        %PMC_A Neuron
        
        v_PMC_A(i+1)=(v_PMC_A(i) + TAU*(RSN.k*(v_PMC_A(i)-RSN.rv)*(v_PMC_A(i)-RSN.vt)-u_PMC_A(i)+ RSN.E + (v_stim*full_matrix(r1,r2,1)) + (w_IpfcA*I_PFC_A(i)) - w_PMC_LI*I_PMC_B(i) )/RSN.C);   % no noise
        u_PMC_A(i+1)=u_PMC_A(i)+TAU*RSN.a*(RSN.b*(v_PMC_A(i)-RSN.rv)-u_PMC_A(i));
        if v_PMC_A(i+1)>=RSN.vpeak;
            v_PMC_A(i)= RSN.vpeak;
            v_PMC_A(i+1)= RSN.c;
            u_PMC_A(i+1)= u_PMC_A(i+1)+ RSN.d;
        end 
        
        if (v_PMC_A(i) >= RSN.vpeak)
            PMCA_Spikes = PMCA_Spikes + 1;
            for j=i:n
               t= j-i;
               I_PMC_A(j)= I_PMC_A(j)+((t/lambda)*exp((lambda-t)/lambda));   
            end
        end
        
        
        
        %PMC_B Neuron
        
        v_PMC_B(i+1)=(v_PMC_B(i) + TAU*(RSN.k*(v_PMC_B(i)-RSN.rv)*(v_PMC_B(i)-RSN.vt)-u_PMC_B(i)+ RSN.E + (v_stim*full_matrix(r1,r2,2)) + (w_IpfcB*I_PFC_B(i)) - w_PMC_LI*I_PMC_A(i) )/RSN.C);   % no noise
        u_PMC_B(i+1)=u_PMC_B(i)+TAU*RSN.a*(RSN.b*(v_PMC_B(i)-RSN.rv)-u_PMC_B(i));
        if v_PMC_B(i+1)>=RSN.vpeak;
            v_PMC_B(i)= RSN.vpeak;
            v_PMC_B(i+1)= RSN.c;
            u_PMC_B(i+1)= u_PMC_B(i+1)+ RSN.d;
        end 
        
        if (v_PMC_B(i) >= RSN.vpeak)
            PMCB_Spikes = PMCB_Spikes + 1;
            for j=i:n
               t= j-i;
               I_PMC_B(j)= I_PMC_B(j)+((t/lambda)*exp((lambda-t)/lambda));   
            end
        end
        
        
        % this code creates a positive number matrix for the PMC voltages
        % only the positive values are retained and all other values 
        % are converted to zero (for computation of integral)
        
        if (v_PMC_A(i) > 0)
            pos_PMCAvolt_values = [pos_PMCAvolt_values;v_PMC_A(i)];
        else
            pos_PMCAvolt_values = [pos_PMCAvolt_values;0];
        end
        
        if (v_PMC_B(i) > 0)
            pos_PMCBvolt_values = [pos_PMCBvolt_values;v_PMC_B(i)];
        else
            pos_PMCBvolt_values = [pos_PMCBvolt_values;0];
        end
        
    end





    % Calculation of Hebbian Weight for PMC_A



    % visual input to PMCA neuron (presynaptic)

    PMCA_vis_stim = full_matrix(r1,r2,1)*vis_input_matrix;                           % visual input to PMCA neuron (presynaptic)

    integral_visinputA   = trapz(PMCA_vis_stim);                                     

    integral_PMCAvoltage = trapz(pos_PMCAvolt_values);                               % activation of PMCA neuron   (post-synaptic)

    g_t_1_A = (integral_PMCAvoltage)-(nmda_const);
    g_t_2_A = [nmda_const - (integral_PMCAvoltage) - ampa_const];


    % Ensures g(t)-1 and g(2)-2 are never less than zero

    if (g_t_1_A <= 0)
        g_t_1_A = 0;
    end

    if (g_t_2_A <= 0)
        g_t_2_A = 0;
    end

    % Equation determining weights at Visual-PMCA Synapses
          
            
    full_matrix(r1,r2,1) = full_matrix(r1,r2,1) + [[(heb_coef)*(integral_visinputA)*[g_t_1_A]*(w_max - full_matrix(r1,r2,1))] - [(anti_heb_const)*(integral_visinputA)*[g_t_2_A]*full_matrix(r1,r2,1)]];

    if (full_matrix(r1,r2,1) < 0)
        full_matrix(r1,r2,1) = 0;
    end

    if (full_matrix(r1,r2,1) > 100)
        full_matrix(r1,r2,1) = 100;
    end



    % Calculation of Hebbian Weight for PMC_B

    PMCB_vis_stim = full_matrix(r1,r2,2)*vis_input_matrix;                             % presynaptic

    integral_visinputB   = trapz(PMCB_vis_stim);                                

    integral_PMCBvoltage = trapz(pos_PMCBvolt_values);                                 % post-synaptic

    g_t_1_B = (integral_PMCBvoltage)-(nmda_const);
    g_t_2_B = [nmda_const - (integral_PMCBvoltage) - ampa_const];


    % Ensures g(t)-1 and g(2)-2 are never less than zero

    if (g_t_1_B <= 0)
        g_t_1_B = 0;
    end

    if (g_t_2_B <= 0)
        g_t_2_B = 0;
    end

    % Equation determining weights at Visual-PMCB Synapses


    full_matrix(r1,r2,2) = full_matrix(r1,r2,2) + [[(heb_coef)*(integral_visinputB)*[g_t_1_B]*(w_max - full_matrix(r1,r2,2))] - [(anti_heb_const)*(integral_visinputB)*[g_t_2_B]*full_matrix(r1,r2,2)]];

    if (full_matrix(r1,r2,2) < 0)
        full_matrix(r1,r2,2) = 0;
    end

    if (full_matrix(r1,r2,2) > 100)
        full_matrix(r1,r2,2) = 100;
    end


    % Appends new Hebbian Weight value to heb_weight Matrix
    % this allows me to plot the Hebbian Values at the End

    % one weight value per trial


    disp('Trial #')
    disp(trial_number)
    disp('r1')
    disp(r1)
    disp('r2')
    disp(r2)
    disp('v_stim')
    disp(v_stim)
    disp('full_matrix(r1,r2,1)')
    disp(full_matrix(:,:,1))
    disp('full_matrix(r1,r2,2)')
    disp(full_matrix(:,:,2))
    %disp('PMCA_vis_stim')
    %disp(PMCA_vis_stim)
    %disp('PMCA_vis_stim')
    %disp(PMCB_vis_stim)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DISPLAY RESULTS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots of Final Voltages, Outputs, and Change of Hebbian Synapse over time

% creating gif of synaptic weight heatmap

%{

n = 1:0.5:5;
nImages = trials;

figure;
for idx = 1:trials
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
plot(TAU*(1:n),v_PFC_A);
axis([0 n -100 100]);
title('PFC_A Neuron Voltage');

subplot(rows,columns,2);
plot(TAU*(1:n),v_PFC_B);
axis([0 n -100 100]);
title('PFC_B Neuron Voltage');

subplot(rows,columns,3);
plot(TAU*(1:n),v_PMC_A);
axis([0 n -100 100]);
title('PMC_A Neuron Voltage');

subplot(rows,columns,4);
plot(TAU*(1:n),v_PMC_B);
axis([0 n -100 100]);
title('PMC_B Neuron Voltage');

subplot(rows,columns,5);
plot(TAU*(1:n),I_PFC_A);
axis([0 n -1 10]);
title('PFC_A Neuron Output');    % I_PMC_B

subplot(rows,columns,6);
plot(TAU*(1:n),I_PFC_B);
axis([0 n -1 10]);
title('PFC_B Neuron Output');

subplot(rows,columns,7);
plot(TAU*(1:n),I_PMC_A);
axis([0 n -1 10]);
title('PMC_A Neuron Output');

subplot(rows,columns,8);
plot(TAU*(1:n),I_PMC_B);
axis([0 n -1 10]);
title('PMC_B Neuron Output');

subplot(rows,columns,9);
data1 = full_matrix(:,:,1);
colormap('hot');   % set colormap
imagesc(data1);    % draw image and scale colormap to values range
colorbar;          % show color scale
title('PMC_A Synaptic Heatmap');

subplot(rows,columns,10);
data2 = full_matrix(:,:,2);
colormap('hot');   % set colormap
imagesc(data2);    % draw image and scale colormap to values range
colorbar;          % show color scale
title('PMC_B Synaptic Heatmap');

% subplot(rows,columns,11);
% data3 = PMCA_weights_with_time(:,:,1);
% colormap('hot');
% imagesc(data3);
% colorbar;
% title('PMC_A Synaptic Heatmap with Slider');
% slider = uicontrol('Style', 'slider', ...
%     'Min', 1, 'Max', trials, ...
%     'Value', 1, ...
%     'Position', [1000 50 400 20]);
% set(slider, 'Callback', 'trial_number=nearest(get(slider,''value''));data3 = PMCA_weights_with_time(:,:,trial_number);colormap(''hot'');imagesc(data3);colorbar;')

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
disp(PMCA_Spikes)
disp('PMC_B Spikes')
disp(PMCB_Spikes)
disp('PFC_A Spikes')
disp(PFCA_Spikes)
disp('PFC_B Spikes')
disp(PFCB_Spikes)

%}