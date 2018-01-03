% algorithm for determining sampling range for hebbian coefficient

% steps
%    select coefficient between 9 and 1    coeff
%    select exponent    between 1 and 15   exp
%    determine range    always one exponent above and one  below   (if coef=5 and exp=9 then range is between 5x10^-8 and 5x10^-10)
%    take uniform sample of size n from within the specified range

n=100;
samplematrix = zeros(n,1);

N = randi(9);
M = randi(15)+3;

% Determines Sample Range (one exponent above and one below)

A = M-1;
B = M+1;

mu_1 =  N*(10^(-M));

for i=1:n
samplematrix(i) = N*(10^(-(A + (B-A)*rand)));
end

X = (1:n);

sortedsamplematrix = sort(samplematrix);

N
M
A
B
mu_1
max(sortedsamplematrix)
min(sortedsamplematrix)

figure;
semilogy(X, sortedsamplematrix);

