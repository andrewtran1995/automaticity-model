NUM_SAMPLES = 1e6;

figure;
rows = 3;
columns = 2;

%% Normal distribution
m = 0.5;
v = 0.25;
x = normrnd(m,v,1,NUM_SAMPLES);

subplot(rows,columns,1);
semilogy(1:NUM_SAMPLES, sort(x));
title('Normal');

%% Log-normal distribution (using lognrnd)
m = 0.5;
v = 0.25;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

subplot(rows,columns,2);
semilogy(1:NUM_SAMPLES, sort(lognrnd(mu,sigma,1,NUM_SAMPLES)));
title('Log-Normal');

%% Uniform distribution
a = 1e-10;
b = 1;

subplot(rows,columns,3);
semilogy(1:NUM_SAMPLES, sort(a + (b-a)*rand(1,NUM_SAMPLES)));
title('Uniform');

%% Uniform distribution (log-space)
range = 10;
x = zeros(1,NUM_SAMPLES);
for i=1:NUM_SAMPLES
    x(i) = randi(range) * 10^(-randi(range));
end

subplot(rows,columns,4);
semilogy(1:NUM_SAMPLES, sort(x));
title('Uniform (Log-Space)');

%% Coefficient and Exponent
N = randi(9);
M = randi(15);
a = M-1;
b = M+1;
x = zeros(1,NUM_SAMPLES);
for i=1:NUM_SAMPLES
    x(i) = N*10^(-(a + (b-a)*rand));
end

subplot(rows,columns,5);
semilogy(1:NUM_SAMPLES, sort(x));
title('Coefficient & Exponent');