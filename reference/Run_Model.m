%Running Automaticity Model for Correlations

RUNS = 8;
PMC_BOLD_S1  = zeros(1,RUNS);
PMC_BOLD_S4  = zeros(1,RUNS);
PMC_BOLD_S10 = zeros(1,RUNS);
PMC_BOLD_S20 = zeros(1,RUNS);
Accuracy_S1  = zeros(1,RUNS);
Accuracy_S4  = zeros(1,RUNS);
Accuracy_S10 = zeros(1,RUNS);
Accuracy_S20 = zeros(1,RUNS);


for i=1:RUNS
Results = autoModel();

x1 =cell2mat(Results(1,1));
x2 =cell2mat(Results(1,2));
x3 =cell2mat(Results(1,3));
x4 =cell2mat(Results(1,4));
x5 =cell2mat(Results(1,5));
x6 =cell2mat(Results(1,6));
x7 =cell2mat(Results(1,7));
x8 =cell2mat(Results(1,8));

PMC_BOLD_S1(i)  = x1;
PMC_BOLD_S4(i)  = x2;
PMC_BOLD_S10(i) = x3;
PMC_BOLD_S20(i) = x4;
Accuracy_S1(i)  = x5;
Accuracy_S4(i)  = x6;
Accuracy_S10(i) = x7;
Accuracy_S20(i) = x8;

end

Correlation_PMC_S1  = corr(PMC_BOLD_S1 , Accuracy_S1) ;
Correlation_PMC_S4  = corr(PMC_BOLD_S4 , Accuracy_S4) ;
Correlation_PMC_S10 = corr(PMC_BOLD_S10, Accuracy_S10);
Correlation_PMC_S20 = corr(PMC_BOLD_S20, Accuracy_S20);

PMC_Correlation_Data = [Correlation_PMC_S1 Correlation_PMC_S4 Correlation_PMC_S10 Correlation_PMC_S20];

figure;
plot(PMC_Correlation_Data);
title('PMC Correlations (Average Bold vs. Accuracy)');
axis([0 4 -1.1 1.1]);



%{
figure
subplot(2,3,1)
plot(x1)
title('Convolution of PMC')
subplot(2,3,2)
%}
