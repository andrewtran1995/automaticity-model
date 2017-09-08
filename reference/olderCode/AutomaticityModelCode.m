% Model of Automaticity in Rule Based Learning


trials = 200;                         % number of trials in our automaticity experiment

tot_pfc_scale_factor = 1/5;           % scaling factor for visual input into PFC neurons


% Initializing PMC stimulus matrices

PMCAstim = zeros(100);
PMCBstim = zeros(100);

% Initializing Radial Basis Function matrix 

pfc_rbv = zeros(100);

% Total stimulus variable for PFCA and PFCB neurons

total_PFCA_stim = 0;
total_PFCB_stim = 0;

% Quantity of Visual Stimulus

v_stim = 1;

% Initializing PMC Synaptic weight matrix (all start with value of 0.1)

PMCA_weights = 0.1*ones(100);
PMCB_weights = 0.1*ones(100);

% Multi-dimensional array (contains important matrices)

full_matrix = cat(3,pfc_rbv, PMCA_weights, PMCB_weights, PMCAstim, PMCBstim);      

%{

disp('row:')
disp(r1)
disp('column: ')
disp(r2)
disp('bar width value: ')                        % displays barwidth value of random trial
disp(full_matrix(r1,r2,1))
disp('frequency value: ')
disp(full_matrix(r1,r2,2))
disp('pFC weight')
disp(full_matrix(r1,r2,3))
disp('PMC weight')
disp(full_matrix(r1,r2,4))

%}


% Time period for one trial (in milliseconds)

n = 1000;
tau = 1;

% Radius for the radial basis function

radius = 10;


% PMC Scaling factor (I can use this to scale the PMC visual input value if it comes out way too high)

pmc_scaling_factor = 1;

% neuron constants (set for a cortical regular spiking neuron)

C_rsn=100; vr_rsn=-60; vt_rsn=-40; k_rsn=0.7;       %RS -- regular spiking
a_rsn=0.03; b_rsn=-2; c_rsn=-50; d_rsn=100;         %RS -- regular spiking
vpeak_rsn=35; E_rsn=60;                             %RS -- regular spiking

% input weights (between PFC and PMC neurons)

w_IpmcA = 1;
w_IpfcA = 100;
w_IpmcB = 1;
w_IpfcB = 100;

% lateral inhibition connections between PFC/PFC and PMC/PMC

w_pFC_LI = 1;
w_PMC_LI = 1;


% maximum possible weight for Hebbian Synapses

w_max = 100;                           

% Lambda Value

lambda = 20; 

% Hebbian Constants (determine the subtle attributes of learning at the Hebbian synapses);

heb_coef = 0.00000000005;
anti_heb_const = 10;
nmda_const = 200;
ampa_const = 50;





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

v_PFC_A = vr_rsn*ones(1,n) ;     u_PFC_A = zeros(1,n);
v_PFC_B = vr_rsn*ones(1,n) ;     u_PFC_B = zeros(1,n);
v_PMC_A = vr_rsn*ones(1,n) ;     u_PMC_A = zeros(1,n);
v_PMC_B = vr_rsn*ones(1,n) ;     u_PMC_B = zeros(1,n);

% Neuron Output Matrices

I_PFC_A = zeros(1,n);
I_PFC_B = zeros(1,n);
I_PMC_A = zeros(1,n);
I_PMC_B = zeros(1,n);


 r1 = randi([1 100],1,1);                         % this function generates new random integer between 1 and 100 everytime the function runs
 r2 = randi([1 100],1,1);                         % the effect is to pick a new random gabor for each trial 

 
 % Calculating visual stimulus to the PFCA and PFCB cells
 
 % left side A    (left 50 columns)
 % right side B   (right  50 columns)
 
 
 total_stim_valA = 0;
 total_stim_valB = 0;
 
 % PFCA stimulus
 
 for i=1:50
    for j=1:100

     
distance = (((r1-i)^2)+((r2-j)^2))^(1/2);
rbv = exp(-(distance/radius));

full_matrix(j,i,1) = rbv*v_stim;

total_stim_valA = total_stim_valA + full_matrix(j,i,1);

    end
 end

 
 
 
% PFCB stimulus

 for i=51:100
    for j=1:100
        
distance = (((r1-i)^2)+((r2-j)^2))^(1/2);
rbv = exp(-(distance/radius));

full_matrix(j,i,1) = rbv*v_stim;

total_stim_valB = total_stim_valB + full_matrix(j,i,1);

    end
 end
 
% final PFCA and PFCB stimulus values are scaled (because they would be too
% high otherwise)
 
 
scaled_total_stim_valA = total_stim_valA * tot_pfc_scale_factor;
scaled_total_stim_valB = total_stim_valB * tot_pfc_scale_factor;
 



total_PMCAstim_val = 0;
total_PMCBstim_val = 0;


% Calculate PMCA stimulus

for i=1:100
    for j=1:100

full_matrix(j,i,4) = full_matrix(j,i,1)*full_matrix(j,i,2)*v_stim;

total_PMCAstim_val = total_PMCAstim_val + full_matrix(j,i,4);

    end
end
 

% Calculate PMCB stimulus
 
for i=1:100
    for j=1:100
        
full_matrix(j,i,5) = full_matrix(j,i,1)*full_matrix(j,i,3)*v_stim;

total_PMCBstim_val = total_PMCBstim_val + full_matrix(j,i,5);

    end
end





% Beginning of Individual Trial

for i=1:n-1
    
   
% Neuron Equations  
       
%PFC A Neuron
    
    v_PFC_A(i+1)=(v_PFC_A(i) + tau*(k_rsn*(v_PFC_A(i)-vr_rsn)*(v_PFC_A(i)-vt_rsn)-u_PFC_A(i)+ E_rsn + scaled_total_stim_valA + (w_IpmcA*I_PMC_A(i)) - w_pFC_LI*I_PFC_B(i))/C_rsn);   % no noise
    u_PFC_A(i+1)=u_PFC_A(i)+tau*a_rsn*(b_rsn*(v_PFC_A(i)-vr_rsn)-u_PFC_A(i));
    if v_PFC_A(i+1)>=vpeak_rsn;
        v_PFC_A(i)= vpeak_rsn;
        v_PFC_A(i+1)= c_rsn;
        u_PFC_A(i+1)= u_PFC_A(i+1)+ d_rsn;
    end 
    
    if (v_PFC_A(i) >= vpeak_rsn)
        PFCA_Spikes = PFCA_Spikes + 1;
        for j=i:n
           t= j-i;
           I_PFC_A(j)= I_PFC_A(j)+((t/lambda)*exp((lambda-t)/lambda));   
        end
    end
    

%PFC B Neuron
    
    v_PFC_B(i+1)=(v_PFC_B(i) + tau*(k_rsn*(v_PFC_B(i)-vr_rsn)*(v_PFC_B(i)-vt_rsn)-u_PFC_B(i)+ E_rsn + scaled_total_stim_valB + (w_IpmcB*I_PMC_B(i)) - w_pFC_LI*I_PFC_A(i))/C_rsn);   % no noise
    u_PFC_B(i+1)=u_PFC_B(i)+tau*a_rsn*(b_rsn*(v_PFC_B(i)-vr_rsn)-u_PFC_B(i));
    if v_PFC_B(i+1)>=vpeak_rsn;
        v_PFC_B(i)= vpeak_rsn;
        v_PFC_B(i+1)= c_rsn;
        u_PFC_B(i+1)= u_PFC_B(i+1)+ d_rsn;
    end 
    
    if (v_PFC_B(i) >= vpeak_rsn)
        PFCB_Spikes = PFCB_Spikes + 1;
        for j=i:n
           t= j-i;
           I_PFC_B(j)= I_PFC_B(j)+((t/lambda)*exp((lambda-t)/lambda));   
        end
    end
    
    
%PMC_A Neuron
    
    v_PMC_A(i+1)=(v_PMC_A(i) + tau*(k_rsn*(v_PMC_A(i)-vr_rsn)*(v_PMC_A(i)-vt_rsn)-u_PMC_A(i)+ E_rsn + pmc_scaling_factor*total_PMCAstim_val + (w_IpfcA*I_PFC_A(i)) - w_PMC_LI*I_PMC_B(i) )/C_rsn);   % no noise
    u_PMC_A(i+1)=u_PMC_A(i)+tau*a_rsn*(b_rsn*(v_PMC_A(i)-vr_rsn)-u_PMC_A(i));
    if v_PMC_A(i+1)>=vpeak_rsn;
        v_PMC_A(i)= vpeak_rsn;
        v_PMC_A(i+1)= c_rsn;
        u_PMC_A(i+1)= u_PMC_A(i+1)+ d_rsn;
    end 
    
    if (v_PMC_A(i) >= vpeak_rsn)
        PMCA_Spikes = PMCA_Spikes + 1;
        for j=i:n
           t= j-i;
           I_PMC_A(j)= I_PMC_A(j)+((t/lambda)*exp((lambda-t)/lambda));   
        end
    end
    
    
    
%PMC_B Neuron
    
    v_PMC_B(i+1)=(v_PMC_B(i) + tau*(k_rsn*(v_PMC_B(i)-vr_rsn)*(v_PMC_B(i)-vt_rsn)-u_PMC_B(i)+ E_rsn + pmc_scaling_factor*total_PMCBstim_val + (w_IpfcB*I_PFC_B(i)) - w_PMC_LI*I_PMC_A(i) )/C_rsn);   % no noise
    u_PMC_B(i+1)=u_PMC_B(i)+tau*a_rsn*(b_rsn*(v_PMC_B(i)-vr_rsn)-u_PMC_B(i));
    if v_PMC_B(i+1)>=vpeak_rsn;
        v_PMC_B(i)= vpeak_rsn;
        v_PMC_B(i+1)= c_rsn;
        u_PMC_B(i+1)= u_PMC_B(i+1)+ d_rsn;
    end 
    
    if (v_PMC_B(i) >= vpeak_rsn)
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


integral_visinput   = trapz(total_PFCA_stim);                                     % need to rework this for the new visual stimulus arrangement

integral_PMCAvoltage = trapz(pos_PMCAvolt_values);

g_t_1_A = (integral_PMCAvoltage)-(nmda_const);
g_t_2_A = [nmda_const - (integral_PMCAvoltage) - ampa_const];


% Ensures g(t)-1 and g(2)-2 are never less than zero

if (g_t_1_A <= 0)
    g_t_1_A = 0;
end

if (g_t_2_A <= 0)
    g_t_2_A = 0;
end

% Equation determining weights of the Visual-PMC_A Synapses
      
      
for i=1:100
    for j=1:100
        
distance = (((r1-i)^2)+((r2-j)^2))^(1/2);
rbv = exp(-(distance/radius));
        
full_matrix(j,i,2) = full_matrix(j,i,2) + full_matrix(j,i,1)*[[(heb_coef)*(integral_visinput)*[g_t_1_A]*(w_max - full_matrix(j,i,2))] - [(anti_heb_const)*(integral_visinput)*[g_t_2_A]*full_matrix(j,i,2)]];

if (full_matrix(j,i,2) < 0)
    full_matrix(j,i,2) = 0;
end

if (full_matrix(j,i,2) > 100)
    full_matrix(j,i,2) = 100;
end

    end
end


% Calculation of Hebbian Weight for PMC_B


integral_visinput   = trapz(total_PFCB_stim);                                


integral_PMCBvoltage = trapz(pos_PMCBvolt_values);

g_t_1_B = (integral_PMCBvoltage)-(nmda_const);
g_t_2_B = [nmda_const - (integral_PMCBvoltage) - ampa_const];


% Ensures g(t)-1 and g(2)-2 are never less than zero

if (g_t_1_B <= 0)
    g_t_1_B = 0;
end

if (g_t_2_B <= 0)
    g_t_2_B = 0;
end

% Equation determining weights at Visual-PMC Synapses


for i=1:100
    for j=1:100
 
distance = (((r1-i)^2)+((r2-j)^2))^(1/2);
rbv = exp(-(distance/radius));

full_matrix(j,i,3) = full_matrix(j,i,3) + full_matrix(j,i,1)*[[(heb_coef)*(integral_visinput)*[g_t_1_B]*(w_max - full_matrix(j,i,3))] - [(anti_heb_const)*(integral_visinput)*[g_t_2_B]*full_matrix(j,i,3)]];

if (full_matrix(j,i,3) < 0)
    full_matrix(j,i,3) = 0;
end

if (full_matrix(j,i,3) > 100)
    full_matrix(j,i,3) = 100;
end

    end
end





% Ensures that weights are always between 0 and 10







% Appends new Hebbian Weight value to heb_weight Matrix
% this allows me to plot the Hebbian Values at the End

% one weight value per trial


if (trial_number > trials - 6)      % this reports statistics for the last five trials

disp('Trial #')
disp(trial_number)
disp('r1')
disp(r1)
disp('r2')
disp(r2)
disp('Strength of PFCA Stimulus')
disp(total_stim_valA)
disp('Strength of PFCB Stimulus')
disp(total_stim_valB)
disp('Strength of PMCA Stimulus')
disp(total_PMCAstim_val)
disp('Strength of PMCB Stimulus')
disp(total_PMCBstim_val)

end






end

data = full_matrix(:,:,2);

                   % matrix to draw
colormap('hot');   % set colormap
imagesc(data);        % draw image and scale colormap to values range
colorbar;          % show color scale



% Final Data Display


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
plot(tau*(1:n),v_PFC_A);
axis([0 n -100 100]);
title('PFC_A Neuron Voltage');
subplot(rows,columns,2);
plot(tau*(1:n),v_PFC_B);
axis([0 n -100 100]);
title('PFC_B Neuron Voltage');
subplot(rows,columns,3);
plot(tau*(1:n),v_PMC_A);
axis([0 n -100 100]);
title('PMC_A Neuron Voltage');
subplot(rows,columns,4);
plot(tau*(1:n),v_PMC_B);
axis([0 n -100 100]);
title('PMC_B Neuron Voltage');
subplot(rows,columns,5);
plot(tau*(1:n),I_PFC_A);
axis([0 n -1 10]);
title('PFC_A Neuron Output');    % I_PMC_B
subplot(rows,columns,6);
plot(tau*(1:n),I_PFC_B);
axis([0 n -1 10]);
title('PFC_B Neuron Output');
subplot(rows,columns,7);
plot(tau*(1:n),I_PMC_A);
axis([0 n -1 10]);
title('PMC_A Neuron Output');
subplot(rows,columns,8);
plot(tau*(1:n),I_PMC_B);
axis([0 n -1 10]);
title('PMC_B Neuron Output');
subplot(rows,columns,9);
data1 = full_matrix(:,:,2);
colormap('hot');   % set colormap
imagesc(data1);    % draw image and scale colormap to values range
colorbar;          % show color scale
title('PMC_A Synaptic Heatmap');
subplot(rows,columns,10);

data2 = full_matrix(:,:,3);
colormap('hot');   % set colormap
imagesc(data2);    % draw image and scale colormap to values range
colorbar;          % show color scale
title('PMC_B Synaptic Heatmap');



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