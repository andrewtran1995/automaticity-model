function [ ] = run_trials( PFC, PMC, PFC_A, PFC_B, PMC_A, PMC_B, NOISE, LAMBDA_PRECALC, TAU, n ) %#codegen
%RUN_TRIALS Summary of this function goes here
    % Define RSN locally (to limit need to pass it as a parameter)
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

    for i=1:n-1
        % Neuron Equations
        % PFC A Neuron
        PFC_A.v(i+1)=(PFC_A.v(i) + TAU*(RSN.k*(PFC_A.v(i)-RSN.rv)*(PFC_A.v(i)-RSN.vt)-PFC_A.u(i)+ RSN.E + PFC_A.v_stim + (PMC_A.W_OUT*PMC_A.out(i)) - PFC.W_LI*PFC_B.out(i))/RSN.C) + normrnd(0,NOISE);
        PFC_A.u(i+1)=PFC_A.u(i)+TAU*RSN.a*(RSN.b*(PFC_A.v(i)-RSN.rv)-PFC_A.u(i));
        if PFC_A.v(i+1)>=RSN.vpeak
            PFC_A.v(i)= RSN.vpeak;
            PFC_A.v(i+1)= RSN.c;
            PFC_A.u(i+1)= PFC_A.u(i+1)+ RSN.d;
        end
        if PFC_A.v(i) >= RSN.vpeak
            PFC_A.out(i:n) = PFC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
        end

        % PFC B Neuron
        PFC_B.v(i+1)=(PFC_B.v(i) + TAU*(RSN.k*(PFC_B.v(i)-RSN.rv)*(PFC_B.v(i)-RSN.vt)-PFC_B.u(i)+ RSN.E + PFC_B.v_stim + (PMC_B.W_OUT*PMC_B.out(i)) - PFC.W_LI*PFC_A.out(i))/RSN.C) + normrnd(0,NOISE);
        PFC_B.u(i+1)=PFC_B.u(i)+TAU*RSN.a*(RSN.b*(PFC_B.v(i)-RSN.rv)-PFC_B.u(i));
        if PFC_B.v(i+1)>=RSN.vpeak
            PFC_B.v(i)= RSN.vpeak;
            PFC_B.v(i+1)= RSN.c;
            PFC_B.u(i+1)= PFC_B.u(i+1)+ RSN.d;
        end
        if PFC_B.v(i) >= RSN.vpeak
            PFC_B.out(i:n) = PFC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
        end

        % PMC_A Neuron
        PMC_A.v(i+1)=(PMC_A.v(i) + TAU*(RSN.k*(PMC_A.v(i)-RSN.rv)*(PMC_A.v(i)-RSN.vt)-PMC_A.u(i)+ RSN.E + PMC_A.v_stim + (PFC_A.W_OUT*PFC_A.out(i)) - PMC.W_LI*PMC_B.out(i) )/RSN.C) + normrnd(0,NOISE);
        PMC_A.u(i+1)=PMC_A.u(i)+TAU*RSN.a*(RSN.b*(PMC_A.v(i)-RSN.rv)-PMC_A.u(i));
        if PMC_A.v(i+1)>=RSN.vpeak
            PMC_A.v(i)= RSN.vpeak;
            PMC_A.v(i+1)= RSN.c;
            PMC_A.u(i+1)= PMC_A.u(i+1)+ RSN.d;
        end
        if PMC_A.v(i) >= RSN.vpeak
            PMC_A.out(i:n) = PMC_A.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
        end

        % PMC_B Neuron
        PMC_B.v(i+1)=(PMC_B.v(i) + TAU*(RSN.k*(PMC_B.v(i)-RSN.rv)*(PMC_B.v(i)-RSN.vt)-PMC_B.u(i)+ RSN.E + PMC_B.v_stim + (PFC_B.W_OUT*PFC_B.out(i)) - PMC.W_LI*PMC_A.out(i) )/RSN.C) + normrnd(0,NOISE);
        PMC_B.u(i+1)=PMC_B.u(i)+TAU*RSN.a*(RSN.b*(PMC_B.v(i)-RSN.rv)-PMC_B.u(i));
        if PMC_B.v(i+1)>=RSN.vpeak
            PMC_B.v(i)= RSN.vpeak;
            PMC_B.v(i+1)= RSN.c;
            PMC_B.u(i+1)= PMC_B.u(i+1)+ RSN.d;
        end
        if PMC_B.v(i) >= RSN.vpeak
            PMC_B.out(i:n) = PMC_B.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
        end
    end

end

