%function [ofb] = COVIS_Model()

%function [ofb, oresponder, ov_fb, oi_fb, ov_resp, oi_resp, orule_used, oZ, Aweights, Bweights, tempv, tempi, h_v, h_i] = COVIS_full()

% Combined 2D stimulus COVIS

clear all
clc
rand('twister',sum(100*clock))

%% General setup

n_sims = 1;

type = 1;

if type == 1
    
%      load input/unix_stim.dat;
%      stimuli = randrows(unix_stim);
    
% Load & transform stim data
    load input/Model_rb_stim.dat;
    data = Model_rb_stim;
    x = (data(:,2) - 0.25) * 30;
    y = (data(:,3) - (pi/9)) * (200/pi);
    %Data(:,1) = data(:,1);
    Data(:,1) = [2*ones(350,1); ones(350,1)];
    Data(:,2) = x;
    Data(:,3) = y;    
    stimuli = Data(randperm(size(Data,1)),:); % randomly reorder rows
    stimuli = stimuli(1:600,:);

    
else

    load input/ii_stim.dat
    stimuli = randrows(ii_stim);
    
end



%% Explicit system initial parameters 

% Rule switching and reward params
delta_c = .04; % positive, incremental salience parameter for correct response
delta_e = .01; % positive decremental salience parameter for error

gamma_exp = 5; % perseveration tendency ie. inverse of switching tendency
lambda = 10;   % selection parameter, Poisson distributed; assume lambda increases with DA levels in PFC

a = 1;         % rule selection gain parameter; amygdala
sigma_e = .5;  % standard deviation of criterial noise in the explicit system


% Various counters, etc
R = 1:(size(stimuli,2)-1);      % rule vector with values for each dimension of x, 1D rules
initial_salience = 1/length(R); % initial salience (Z) of rules
Y = zeros(length(R),1);         % rule weights vector
Z = ones(1, length(R)).* initial_salience; % initial salience vector

%rule_counter = 0; % counts number of times a rule is used
%P = zeros(size(stimuli,1),4); % keep track of selection probabilities
%rsw = zeros(640,1); % keep track of attempted rule switches
% b = ones(length(R),1); %Discriminant weights vector
% b(1) = 50; b(2) = 50; %Pre-define discriminant weights


% Initial rule guesses
n_guess = 5; % number of guessing trials
initial_rules = ceil(rand(1,n_guess) * length(R)); % Randomize initial guesses
initial_rules = initial_rules(randperm(n_guess));

% Random rule switches
switch_attempt = ceil(rand(1,size(stimuli,1)) * length(R));
switch_num = 0;




%% Procedural system initial params

% Stim params
visdim = 100;
visbase = .01;
visrange = .025;
rbf_std = 1.75;%2; %Value = .75 from Matt's code
epsilon_i = .33; % procedural system noise parameter

% Reward prediction/DA error params
predicted(1) = 0;
pr_alpha = .1;%.05;
DA_slope = .8;
DA_base = .2;
DA_max = 1;

% Learning rate params
alpha = .65; %.1 works
beta = 0.19;
gamma = .02;
ampa_thresh = .01;
nmda_thresh = .0022;
w_max = 1;

Aweights(:,:,1) = visbase + visrange*rand(visdim, visdim);
Bweights(:,:,1) = visbase + visrange*rand(visdim, visdim);


%% Competition parameters

odelta_e = 0.04; %.01 % Explicit system incorrect
odelta_c = 0.01; % .004 % Explicit system correct

initial_theta_v = .99;
initial_theta_i = 1-initial_theta_v;

theta_v(1) = initial_theta_v; 
theta_i(1) = initial_theta_i;

%%  Begin simulations

tic

for iteration = 1:n_sims
    
    if (mod(iteration,50)==0 || iteration==1)
        disp(['Simulation #' num2str(iteration)])
    end
    
    stimuli = stimuli(randperm(size(stimuli,1)),:);
    
    if iteration > 1 
        % Re-initialize parameters for each system
        clear Y Z P v_resp v_fb h_v rule_used rule_weights
        clear Aweights Bweights predicted actual DA i_fb i_resp h_i R1 R2
        clear resp fb responder theta_v theta_i

        % Explicit sys reset
        initial_salience = 1/length(R); %(Z_initial)
        Y = zeros(length(R),1); 
        Z = ones(1, length(R)).* initial_salience; 
        %b = ones(length(R),1); %Discriminant weights vector pre-defined
        %P = zeros(size(stimuli,1),4); 
        %b(1) = 50; b(2) = 50; %Pre-define discriminant weights
        initial_rules = ceil(rand(1,n_guess) * length(R)); % Randomize 'initial' rules during model guessing
        initial_rules = initial_rules(randperm(n_guess));
       
        %Initialize random rule switches
        switch_attempt = ceil(rand(1,size(stimuli,1))*length(R));
        switch_num = 0;
        
        % Implicit sys reset   
        predicted(1) = 0;
        Aweights(:,:,1) = visbase + visrange*rand(visdim, visdim);
        Bweights(:,:,1) = visbase + visrange*rand(visdim, visdim);

        %Overall sys reset
        theta_v(1) = initial_theta_v; 
        theta_i(1) = initial_theta_i;
    
    end

%% Begin experiment

for n = 1:size(stimuli,1)

    if (mod(n,25)==0 || n==1)
        disp(['Trial: ' num2str(n)]);
    end
    
    stim(n,:) = stimuli(n,2:3); %Select only x,y coords on each trial
    

%% Explicit system

    epsilon_e = sigma_e * randn; % Calculate the noisy criterion value
    
    % Calculate discriminant value and response (guess on early trials)
    
    if n <= n_guess
        current = initial_rules(n);
        rule_used(n) = R(current); % Record rule used 
        rule_weights(n,:) = Y;
        
        h_v(n) = 1;
        v_resp(n) = ceil(rand(1)*2); % Random response
        
    else
  
        % Calculate discriminant value
        h_v(n) = stim(n,R(current)) - mean(stim(:,R(current))) - epsilon_e;
              
        % Decide response: A = 1; B = 2
        if h_v(n) <= 0
            v_resp(n) = 1;
        else
            v_resp(n) = 2;
        end
    
    end
    
%     if n == 1
%         
%         %current = ceil(rand(1)*length(R));
%         
%         weights = initial_salience/sum(initial_salience);
% 
%         wt_int = weights(1);
% 
%         for i = 2:length(weights)
% 
%             wt_int(i) = wt_int(i-1)+weights(i);
% 
%         end
%         
%         current = min(find(wt_int > rand));
%         rule_used(n) = current; %Record rule used on trial 1
%         rule_weights(n,:) = Y;
%         
%     end
        
%     %Calculate discriminant function
%     h_v(n) = mean(abs(stim'.*R(:,current) - b(:,current).*R(:,current)));
    
    
%     h_v(n) = stim(current) - b(current) - epsilon_e;

%     % Decide response: A = 1; B = 2
%     if h_v(n) < 0 
%         resp(n) = 1;
%     elseif h_v(n) > 0
%         resp(n) = 2;
%     end

%     if all(stim'.*R(:,current) >= (b(:,current)+epsilon_e).*R(:,current))
% 
%         v_resp(n) = resp_key(current);
% 
%     else v_resp(n) = 3-resp_key(current);
% 
%     end


    % Give feedback and update salience parameters
    if v_resp(n) == stimuli(n,1)

        v_fb(n) = 1; %correct
        Z(current) = Z(current) + delta_c;

    else 

        v_fb(n) = 0; %error
        Z(current) = Z(current) - delta_e;

    end

    % Rule switching/selection
    if v_fb(n) ~= 1 %Use same rule if correct 
        
        % Step 1: select rule for next trial
        switch_num = switch_num + 1;
        
        new = switch_attempt(switch_num);

        %new = ceil(rand(1)*length(R)); %choose random rule otherwise
        
        %if new = current, force model to pick a different rule for
        %switching consideration:
        while new == current

            switch_num = switch_num + 1;

            new = switch_attempt(switch_num);

            %new = ceil(rand(1)*length(R));

        end

        % Step 2: Adjust weights (Y)
        % First, set salience param (Z) at each rule weight
        for i = 1:length(R)

            Y(i) = Z(i);

        end

        % Next, add gamma and X (poisson distribured RV) to current and
        % new rule
        Y(current) = Z(current) + gamma_exp;

        Y(new) = Z(new) + gen_pois(lambda);


        % Step 3: Choose rule with max normalized probability 
        weights = (Y.^a)/sum(Y.^a);

        wt_int = weights(1);

        for i = 2:length(weights)

            wt_int(i) = wt_int(i-1)+weights(i);

        end
        
        current = min(find(wt_int > rand));

    end

    % Record rule used and weights
    rule_used(n + 1) = current;
    rule_weights(n + 1,:) = Y;


%% Procedural system

    % Calculate visual activity to stimulus (i.e., rbf centered at stim)
    for i = 1:visdim

            for j = 1:visdim

                ydist = stim(n,2) - j;
                xdist = stim(n,1) - i;

                A(i,j,n) = exp(-(xdist*xdist + ydist*ydist)/(2*rbf_std^2));

                %Another way using matrix algebra...
                %dist = [[i,j]' - stim(n,2:3)'];
                %A(i,j,n) = exp(-(dist'*dist)/(2*rbf_std^2));

            end

    end
        
    % Calculate striatal response to vis activity

    % Reshape matrices into vectors
    wAvect = reshape(Aweights(:,:,n), 1, visdim^2);
    wBvect = reshape(Bweights(:,:,n), 1, visdim^2);
    actVect = reshape(A(:,:,n), 1, visdim^2);

    %Calculate sum of weights * visactivity + noise (R1, R2) 
    R1(n) = (wAvect*actVect') + epsilon_i*randn;
    R2(n) = (wBvect*actVect') + epsilon_i*randn;

    %Calculate striatal discriminant function
    h_i(n) = R1(n) - R2(n);

    %Decide response: A = 1; B = 2
    if h_i(n) < 0
        
        i_resp(n) = 2;
        
    elseif h_i(n) > 0
        
        i_resp(n) = 1;
        
    elseif h_i(n) == 0 %if difference is zero respond randomly
        
        if rand >= .5
            
            i_resp(n) = 1;
            
        else
            
            i_resp(n) = 2;
            
        end
        
    end
    
    %Give feedback 
    if i_resp(n) == stimuli(n,1)
        i_fb(n) = 1; %correct
    else 
        i_fb(n) = 0; %error
    end
    
    
    % Update weights 
    
    %Calculate actual reward
    if i_fb(n) == 1
        
        actual(n) = 1;
        
    elseif i_fb(n) == 0
        
        actual(n) = -1;
        
    elseif i_fb(n) == 2 %ie, no feedback
        
        actual(n) = 0;
        
    end
    

    %Calculate predicted reward, dopamine
%     predicted(n+1) = predicted(n) + pr_alpha*(actual(n) - predicted(n));
%     
%     DA(n) = 0.8*(actual(n) - predicted(n+1)) + 0.2;
    
    predicted(n+1) = predicted(n) + pr_alpha*(actual(n) - predicted(n));
    
    DA(n) = cap(DA_slope*(actual(n) - predicted(n)) + DA_base, DA_max);
    
%     if DA(n) < 0
%         DA(n) = 0;
%     elseif DA(n) >1
%         DA(n) = 1;
%     end


    % Update weights
%     for i = 1:visdim
%         
%         for j = 1:visdim
%             
%             Aweights(i,j,n+1) = pos(Aweights(i,j,n) ...
%                 + alpha * A(i,j,n) * pos(pos(R1(n)) - nmda_thresh) * pos(DA(n) - DA_base) * (w_max - Aweights(i,j,n)) ...
%                 - beta * A(i,j,n) * pos(pos(R1(n)) - nmda_thresh) * pos(DA_base - DA(n)) * Aweights(i,j,n) ...
%                 - gamma * A(i,j,n) * pos(pos(nmda_thresh - pos(R1(n))) - ampa_thresh) * Aweights(i,j,n));
%             
%             Bweights(i,j,n+1) = pos(Bweights(i,j,n) ...
%                 + alpha * A(i,j,n) * pos(pos(R2(n)) - nmda_thresh) * pos(DA(n) - DA_base) * (w_max - Bweights(i,j,n)) ...
%                 - beta * A(i,j,n) * pos(pos(R2(n)) - nmda_thresh) * pos(DA_base - DA(n)) * Bweights(i,j,n) ...
%                 - gamma * A(i,j,n) * pos(pos(nmda_thresh - pos(R2(n))) - ampa_thresh) * Bweights(i,j,n));
%             
%         end
%         
%     end
    
    tempAw = Aweights(:,:,n);
    tempAw = tempAw(:);
    tempBw = Bweights(:,:,n);
    tempBw = tempBw(:);

    for i = 1:length(tempAw)
        
        if R1(n) >= R2(n) % Implement simple lateral inhibition: only update weights for the responding unit
        
            tempAw(i) = pos(tempAw(i) ...
                + alpha * actVect(i) * pos(pos(R1(n)) - nmda_thresh) * pos(DA(n) - DA_base) * (w_max - tempAw(i)) ...
                - beta * actVect(i) * pos(pos(R1(n)) - nmda_thresh) * pos(DA_base - DA(n)) * tempAw(i) ...
                - gamma * actVect(i) * pos(pos(nmda_thresh - pos(R1(n))) - ampa_thresh) * tempAw(i));
            
        else

            tempBw(i) = pos(tempBw(i) ...
                + alpha * actVect(i) * pos(pos(R2(n)) - nmda_thresh) * pos(DA(n) - DA_base) * (w_max - tempBw(i)) ...
                - beta * actVect(i) * pos(pos(R2(n)) - nmda_thresh) * pos(DA_base - DA(n)) * tempBw(i) ...
                - gamma * actVect(i) * pos(pos(nmda_thresh - pos(R2(n))) - ampa_thresh) * tempBw(i));
            
        end
        
    end
    
    tempAw = reshape(tempAw,visdim,visdim);
    tempBw = reshape(tempBw,visdim,visdim);
    
    Aweights(:,:,n+1) = tempAw;
    Bweights(:,:,n+1) = tempBw;
    
    clear tempAw tempBw

    
    wAvect = [];
    wBvect = [];
    Actvect = [];
    
%% Calculate overall system response/responder and overall feedback

    %Normalize h_i and h_v on [0,1] scale
    h_i_old(n) = h_i(n);
    h_i_max(n) = max(abs(h_i_old)); %Keep track of max "distance" over experiment
    h_i(n) = h_i(n)/h_i_max(n);
    
    h_v_old(n) = h_v(n);
    h_v_max(n) = max(abs(h_v_old));
    h_v(n) = h_v(n)/h_v_max(n);
    
    if abs(theta_v(n)*h_v(n)) > abs(theta_i(n)*h_i(n))
        
        resp(n) = v_resp(n);
        
        responder(n) = 1; %Explicit system response
        
        fb(n) = v_fb(n);
        
    elseif abs(theta_v(n)*h_v(n)) < abs(theta_i(n)*h_i(n))
        
        resp(n) = i_resp(n);
        
        responder(n) = 2; %Procedural system response
        
        fb(n) = i_fb(n);
        
    end    
    
%     tempv(n) = theta_v(n)*h_v(n);
%     tempi(n) = theta_i(n)*h_i(n);

    % Update system weights
    theta_v(n+1) = theta_v(n) + (v_fb(n)-1)*odelta_e*theta_v(n) + v_fb(n)*odelta_c*(1-theta_v(n));
    theta_i(n+1) = 1 - theta_v(n+1);
    
end


%% Save data...

ofb(iteration,:) = fb;
oresp(iteration,:) = resp;
oresponder(iteration,:) = responder;
ostimuli(:,:,iteration) = stimuli;


ov_fb(iteration,:) = v_fb;
ov_resp(iteration,:) = v_resp;
oi_fb(iteration,:) = i_fb;
oi_resp(iteration,:) = i_resp;
oDA(iteration,:) = DA;

orule_used(iteration,:) = rule_used;
oZ(:,:,iteration) = Z;


end %end iteration

toc


figure(); surf(Aweights(:,:,600)); figure(); surf(Bweights(:,:,600));

if length(stimuli) == 600

% Calculate block accuracy for each system and overall system
for blk = 1:6
    
    for itr = 1:size(ofb,1)
        
        block_acc(itr,blk) = sum(ofb(itr,(blk-1)*100+1:(blk*100)))/100;
        
        block_acc_v(itr,blk) = sum(ov_fb(itr,(blk-1)*100+1:(blk*100)))/100;
        
        block_acc_i(itr,blk) = sum(oi_fb(itr,(blk-1)*100+1:(blk*100)))/100;
        
    end
    
end


% Plotting...

block_mean = mean(block_acc,1);
block_var = var(block_acc,1);

block_mean_v = mean(block_acc_v,1);
block_mean_i = mean(block_acc_i,1);

avgAw = squeeze(mean(Aweights,1));
avgBw = squeeze(mean(Bweights,1));


figure('position',[400 100 1000 800]);
subplot(5,1,1); plot(block_mean,'ko-'); hold on; plot(block_mean_v,'bo-'); plot(block_mean_i,'ro-'); hold off
    xlabel('Block'); ylabel('Acc'); legend('COVIS','RB Sys','II Sys','location','NorthWest'); title('Block Accuracy');
subplot(5,1,2); plot(mean(oDA),'k.');
    ylabel('Dopamine'); title('Average Dopamine Across Trials');
subplot(5,1,3); surf(Aweights(:,:,600)); %plot(avgAw(1,1:600),'bo'); hold on; plot(avgAw(2,1:600),'ro'); hold off;
    title('Category A weights from final simulation (trial 600)')
subplot(5,1,4); surf(Bweights(:,:,600)); %plot(avgBw(1,1:600),'bo'); hold on; plot(avgBw(2,1:600),'ro'); hold off;
    title('Category B weights from final simulation (trial 600)')
subplot(5,1,5); plot(mean(oresponder));
    xlabel('Trial'); ylabel('Responder'); title('Average Responder (1 = Explicit; 2 = Procedural)')
    
end