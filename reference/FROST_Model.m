%FROST

%Questions, Concerns, and Room for Improvement



%All Lambda's are set to 30 (more biologically accurate lambdas can be inserted)


%The weights of the synapses are all set fairly arbitrarily


%I increased the Thalamic resting firing rate above E_rsn

%why doesn't the signal degenerate?  Should it be degenerating?

%Greg seems to think that the memory signal should degenerate to zero


%signal from the driving region of pFC remains constant
%so why should the thalamic pFC loop degenerate?


%should the lateral inhibition only be active during the recall period?
%or do these cells display lateral inhibition all the time?





%instructions

%set # of trials
%set # of pFC units
%run program

%comment off approapraite regions in code below %{  }%


%output

%probability that stimulus retrieval is above threshold - for each data point



clear all;

%noise = 0.5;

trials = 10;

a1total = 0;
a2total = 0;
a3total = 0;
a4total = 0;
a5total = 0;
a6total = 0;
a7total = 0;
a8total = 0;
a9total = 0;
a10total = 0;

item =10;


for z=1:trials


% define the number of units and input strength and duration


%why do we get a bump at units = 6?

                                %number of pFC units
input_strength_PPC= 100;        %input strength to PPC
input_strength_drpfc= 100;


% constants

C_rsn=100; vr_rsn=-60; vt_rsn=-40; k_rsn=0.7;       %RS -- regular spiking
a_rsn=0.03; b_rsn=-2; c_rsn=-50; d_rsn=100;         %RS -- regular spiking
vpeak_rsn=35; E_rsn=60;                             %RS -- regular spiking

C_msn=50; vr_msn=-80; vt_msn=-25; k_msn=1;          %MSN -- striatum             
a_msn=0.01; b_msn=-20; c_msn=-55; d_msn=150;        %MSN -- striatum
vpeak_msn=40; E_msn=100;                            %MSN -- striatum

%Beta_iaf = 11.83; Gamma_iaf = 0.117; vr_iaf = -60;  %Quadratic I&F
%vt_iaf = -40; vpeak_iaf = 35; vreset_iaf = -50;     %Quadratic I&F

%time variables

T=12000; tau=1;     %if you use tau as milliseconds
n=round(T/tau);    %then total time n equals 3000 milliseconds
          
% Constants for GPi Quadratic Integrate and Fire Neurons

Beta_qiaf = 11.83 ; Gamma_qiaf = 0.117 ; Vr_qiaf = -60 ; Vt_qiaf = -40 ;
Vpeak_qiaf = 35 ; Vreset_qiaf = -50 ;

%data initialization

I=[input_strength_PPC*ones(1,0.2*n),zeros(1,0.8*n)];           %step-wise input to PPC (from visual area)
Ipfcdriv = [zeros(1,0.35*n),100*ones(1,0.25*n), zeros(1,0.4*n)];                     %step-wise input to CD  (from driving region of pFC)


% matrices for voltage data of each brain region

v_PPC = vr_rsn*ones(1,n) ;     u_PPC = zeros(1,n);
v_pFC = vr_rsn*ones(1,n) ;     u_pFC = zeros(1,n);
v_pFCdriv = vr_rsn*ones(1,n) ; u_pFCdriv = zeros(1,n);
v_CD = vr_msn*ones(1,n)      ; u_CD = zeros(1,n);
v_GP = Vr_qiaf*ones(1,n) ;
v_Thal = vr_rsn*ones(1,n);     u_Thal = zeros(1,n);


%alpha function matrices and constants

Ipfc = zeros(1,n);
Ilipfc = zeros(1,n);
Idpfc= zeros(1,n);
Ippc = zeros(1,n);
Icd  = zeros(1,n);
Igpi = zeros(1,n);
Ithal= zeros(1,n);

lambda = 54;
                        %lambda is set to be the same for all neurons
                        %this might not be a good idea
                        
                        
                            
% Weights of synaptic connections

w_Ippc = 16;            
w_Ipfc = 16;
w_IdpfctoThal = 0;
w_IdpfctoCD = 70;
w_Icd  = 5;
w_Igpi = 3.8;
w_Ithal= 35;
w_IthaltodpFC = 0;   %if this weight is too high the reverberation becomes constant (does not decay)
w_Ilipfc = 0.5;


%v_CD
%v_GP
%v_Thal
%other brain regions?



%explanation of constants

%C_msn - ?
%vr_msn - membrane voltage of MSN cell
%vt_msn - voltage threshold for AP
%k_msn - leakiness constant?
%a_msn - ?
%b_msn - ?
%c_msn - voltage reset value for MSN neuron after AP
%d_msn - value added to u(t) after AP (responsible for downswing after AP)
%vpeak_msn - peak value for AP
%E_msn - value sets the baseline firing rate (sometimes called Beta)

%all constants needed for PSP from PPC and pFC neurons
pFC_spikes = 0;
PPC_spikes = 0;
pFCdriv_spikes = 0;
CD_spikes = 0;
GP_spikes = 0;
Thal_spikes = 0;


for i=1:n-1

    
    %PPC Neurons
    
    v_PPC(i+1)=(v_PPC(i) + tau*(k_rsn*(v_PPC(i)-vr_rsn)*(v_PPC(i)-vt_rsn)-u_PPC(i) + E_rsn + I(i) + w_Ipfc*(Ipfc(i)))/C_rsn); 
    u_PPC(i+1)=u_PPC(i)+tau*a_rsn*(b_rsn*(v_PPC(i)-vr_rsn)-u_PPC(i));
    if v_PPC(i+1)>=vpeak_rsn
        v_PPC(i)=vpeak_rsn;
        v_PPC(i+1)=c_rsn;
        u_PPC(i+1)= u_PPC(i+1)+d_rsn;
    end
    
    %alpha functions for PPC PSPs (Post-synaptic input to pFC neurons)

% this input is excitatory

 
    if (v_PPC(i)>=vpeak_rsn)     
        PPC_spikes = PPC_spikes + 1;
        for j=i:n
           t= j-i;
           Ippc(j) = Ippc(j) + ((t/lambda)*exp((lambda-t)/lambda));
        end
    end
    
    %pFC WM unit Neurons
    
    
    
    
    v_pFC(i+1)=(v_pFC(i) + tau*(k_rsn*(v_pFC(i)-vr_rsn)*(v_pFC(i)-vt_rsn)-u_pFC(i)+ E_rsn + w_Ippc*(Ippc(i)) - (w_Ilipfc*(item-1)*(Ilipfc(i))) + w_Ithal*(Ithal(i)) )/C_rsn) + normrnd(0,.2);   %normrnd(0,0.027)
    u_pFC(i+1)=u_pFC(i)+tau*a_rsn*(b_rsn*(v_pFC(i)-vr_rsn)-u_pFC(i));
    if v_pFC(i+1)>=vpeak_rsn;
        v_pFC(i)=vpeak_rsn;
        v_pFC(i+1)=c_rsn;
        u_pFC(i+1)= u_pFC(i+1)+ d_rsn;
    end 
    
        
%alpha functions for pFC PSPs (Post-synaptic input to PPC neurons)

% this input is excitatory

    if (v_pFC(i) >= vpeak_rsn)     
        pFC_spikes = pFC_spikes + 1;
        for j=i:n
           t= j-i;
           Ipfc(j)= Ipfc(j)+((t/lambda)*exp((lambda-t)/lambda));   
        end
    end
    
    %alpha functions for pFC lateral inhibition

% this input is excitatory

    if (v_pFC(i) >= vpeak_rsn)     
        pFC_spikes = pFC_spikes + 1;
        for j=i:n
           t= j-i;
           Ilipfc(j)= Ilipfc(j)+((t/lambda)*exp((lambda-t)/lambda));   
        end
    end
    
    %pFC driving neurons
    
    v_pFCdriv(i+1)=v_pFCdriv(i) + tau*(k_rsn*(v_pFCdriv(i)-vr_rsn)*(v_pFCdriv(i)-vt_rsn) - u_pFCdriv(i) + E_rsn + Ipfcdriv(i))/C_rsn;
    u_pFCdriv(i+1)=u_pFCdriv(i)+tau*a_rsn*(b_rsn*(v_pFCdriv(i)-vr_rsn)-u_pFCdriv(i));
    if v_pFCdriv(i+1)>=vpeak_rsn;
        v_pFCdriv(i)=vpeak_rsn;
        v_pFCdriv(i+1)=c_rsn;
        u_pFCdriv(i+1)= u_pFCdriv(i+1)+ d_rsn;
    end
    
    %alpha functions of pFC driving neuron PSPs (Post-synaptic input to CD
%neurons

% this input is excitatory
    
    if (v_pFCdriv(i) >= vpeak_rsn)  
        pFCdriv_spikes = pFCdriv_spikes + 1;
        for j=i:n
           t= j-i;
           Idpfc(j) = Idpfc(j) + ((t/lambda)*exp((lambda-t)/lambda));
        end
    end
    
    %CD medium spiny neurons    
        
    v_CD(i+1)=v_CD(i) + tau*(k_msn*(v_CD(i)-vr_msn)*(v_CD(i)-vt_msn)-u_CD(i)+ E_msn + w_IdpfctoCD*(Idpfc(i)))/C_msn;
    u_CD(i+1)=u_CD(i)+tau*a_msn*(b_msn*(v_CD(i)-vr_msn)-u_CD(i));
    if v_CD(i+1)>=vpeak_msn;
        v_CD(i)=vpeak_msn;
        v_CD(i+1)=c_msn;
        u_CD(i+1)= u_CD(i+1)+ d_msn;  
    end
    
    %alpha functions of CD neuron PSPs (Post-synaptic input to GPi neurons)

% this input is inhibitory
    
    if (v_CD(i)>=vpeak_msn)    %alpha functions of CD neuron PSPs are calculated here  (25.5 value is arbitrary, but it returns accurate representation of spike number)
        CD_spikes = CD_spikes + 1;
        for j=i:n
           t= j-i;
           Icd(j) = Icd(j) + ((t/lambda)*exp((lambda-t)/lambda));
        end
    end
    
     % Model GP with Quadratic Integrate and Fire
     
    dGP = (-1)*w_Icd*Icd(i)+ Beta_qiaf + Gamma_qiaf*(v_GP(i)- Vr_qiaf)*(v_GP(i)-Vt_qiaf);
    v_GP(i+1) = v_GP(i) + dGP;
    if (v_GP(i+1) >= Vpeak_qiaf);
         v_GP(i) = Vpeak_qiaf;
         v_GP(i+1) = Vreset_qiaf;
    end;
    
    %alpha functions of GPi neuron PSPs (Post-synaptic input to Thalamic neurons)

    % this input is inhibitory
    
    if (v_GP(i)>= Vpeak_qiaf)    %alpha functions of CD neuron PSPs are calculated here  (25.5 value is arbitrary, but it returns accurate representation of spike number)
        GP_spikes = GP_spikes + 1;
        for j=i:n
           t= j-i;
           Igpi(j) = Igpi(j) + ((t/lambda)*exp((lambda-t)/lambda));
        end
    end
    
    %other pFC region - stimulated by the onset of the delay period (right
    %after the visual stimulus stops (attentional demands increase)
    
    %this means these neurons can be driven by a stepwise input stimulus
    %just like the PPC neurons
    
    
    %Thalamus
    
    v_Thal(i+1)=v_Thal(i) + tau*(k_rsn*(v_Thal(i)-vr_rsn)*(v_Thal(i)-vt_rsn)-u_Thal(i)+ E_rsn - w_Igpi*Igpi(i) + (w_Ipfc*Ipfc(i)))/C_rsn;    %Should I change E_rsn value here?
    u_Thal(i+1)=u_Thal(i)+tau*a_rsn*(b_rsn*(v_Thal(i)-vr_rsn)-u_Thal(i));                                                                                      %Does thalamus have high intrnsic firing rate?
    if v_Thal(i+1)>= vpeak_rsn
        v_Thal(i)= vpeak_rsn;
        v_Thal(i+1)= c_rsn;
        u_Thal(i+1)= u_Thal(i+1)+d_rsn;
    end
    
    if (v_Thal(i) >= vpeak_rsn)    %alpha functions of CD neuron PSPs are calculated here  (25.5 value is arbitrary, but it returns accurate representation of spike number)
    Thal_spikes = Thal_spikes + 1;
        for j=i:n
           t= j-i;
           Ithal(j) = Ithal(j) + ((t/lambda)*exp((lambda-t)/lambda));
        end
    end
    
end 


%disp(PPC_spikes);
%disp(pFC_spikes);
%disp(pFCdriv_spikes);
%disp(CD_spikes);
%disp(GP_spikes);

%v_pFC = max(v_pFC + 1 * randn(1, n), -100);      %adds noise to pFC signal

%do we need noise here?



thresh = 2.3;

a1 = 0;
for y=5750:6249
    a1 = a1 + Ipfc(y);
end

if ((a1/500)>=thresh)
    a1total = a1total + 1;
end

a2 = 0;
for y=6250:6749
    a2 = a2 + Ipfc(y);
end

if ((a2/500)>=thresh)
    a2total = a2total + 1;
end

a3 = 0;
for y=6750:7249
    a3 = a3 + Ipfc(y);
end

if ((a3/500)>=thresh)
    a3total = a3total + 1;
end

a4 = 0;
for y=7250:7749
    a4 = a4 + Ipfc(y);
end

if ((a4/500)>=thresh)
    a4total = a4total + 1;
end

a5 = 0;
for y=7750:8249
    a5 = a5 + Ipfc(y);
end

if ((a5/500)>=thresh)
    a5total = a5total + 1;
end

a6 = 0;
for y=8250:8749
    a6 = a6 + Ipfc(y);
end

if ((a6/500)>=thresh)
    a6total = a6total + 1;
end

a7 = 0;
for y=8750:9249
    a7 = a7 + Ipfc(y);
end

if ((a7/500)>=thresh)
    a7total = a7total + 1;
end

a8 = 0;
for y=9250:9749
    a8 = a8 + Ipfc(y);
end

if ((a8/500)>=thresh)
    a8total = a8total + 1;
end

a9 = 0;
for y=9750:10249
    a9 = a9 + Ipfc(y);
end

if ((a9/500)>=thresh)
    a9total = a9total + 1;
end

a10 = 0;
for y=10250:10749
    a10 = a10 + Ipfc(y);
end

if ((a10/500)>=thresh)
    a10total = a10total + 1;
end

end



plot1 = a1total/(trials*item);


fprintf('Probability of Recall of First Data Point');
disp(a1total/(trials));
fprintf('Probability of Recall of Second Data Point');
disp(a2total/(trials));
fprintf('Probability of Recall of Third Data Point');
disp(a3total/(trials));
fprintf('Probability of Recall of Fourth Data Point');
disp(a4total/(trials)); 
fprintf('Probability of Recall of Fifth Data Point');
disp(a5total/(trials));
fprintf('Probability of Recall of Sixth Data Point');
disp(a6total/(trials));
fprintf('Probability of Recall of Seventh Data Point');
disp(a7total/(trials));
fprintf('Probability of Recall of Eighth Data Point');
disp(a8total/(trials));
fprintf('Probability of Recall of Ninth Data Point');
disp(a9total/(trials));
fprintf('Probability of Recall of Tenth Data Point');
disp(a10total/(trials));
    

subplot(7,2,1);                  %how to graph multiple plots in one figure
plot(tau*(1:n),v_PPC);           %to read more about this: http://www.mathworks.com/help/matlab/examples/displaying-multiple-plots-in-a-single-figure.html?refresh=true
axis([0 13000 -100 100]);
title('PPC Neurons');
subplot(7,2,3);
plot(tau*(1:n),v_pFC);
axis([0 13000 -100 100]);
title('pFC Neurons');
subplot(7,2,5);
plot(tau*(1:n),v_pFCdriv);
axis([0 13000 -100 100]);
title('pFC CD Driving Neurons');
subplot(7,2,2);
plot(tau*(1:n),Ipfc);
axis([0 13000 -2 10]);
title('pFC Output');
subplot(7,2,4);
plot(tau*(1:n),Ipfc);
axis([0 13000 -2 10]);
title('PPC Output');
subplot(7,2,6);
plot(tau*(1:n),Idpfc);
axis([0 13000 -2 10]);
title('drpFC Output');
subplot(7,2,7);
plot(tau*(1:n),v_CD);
axis([0 13000 -100 0]);
title('CD Neurons');
subplot(7,2,8);
plot(tau*(1:n),Icd);
axis([0 13000 -2 10]);
title('CD Output');
subplot(7,2,9);
plot(tau*(1:n),v_GP);
axis([0 13000 -100 100]);
title('GPi Neurons');
subplot(7,2,10);
plot(tau*(1:n),Igpi);
axis([0 13000 -2 10]);
title('GPi Output');
subplot(7,2,11);
plot(tau*(1:n),v_Thal);
axis([0 13000 -100 100]);
title('Thalamic Neurons');
subplot(7,2,12);
plot(tau*(1:n),Ithal);
axis([0 13000 -2 10]);
title('Thalamic Output');
subplot(7,2,13);
plot(1,a1total/(trials),'b--o',2,a2total/(trials),'b--o',3,a3total/(trials),'b--o',4,a4total/(trials),'b--o',5,a5total/(trials),'b--o',6,a6total/(trials),'b--o',7,a7total/(trials),'b--o',8,a8total/(trials),'b--o',9,a9total/(trials),'b--o',10,a10total/(trials), 'b--o');
axis([0 10 0 1]);
title('Probability Plot');


%Adding Noise to the System?
%
%What is a good way to model noise?
%
%pFC = max(pFC + 0.027 * randn(num, tTotal), 0); what does this do exactly?






