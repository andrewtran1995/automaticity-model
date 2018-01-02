% alternate way to represent numbers to get a continuous distribution
% through a logarithmic space

n = 1000;

functionmatrix = zeros(n,1);
Nmatrix = zeros(n,1);
Mmatrix = zeros(n,1);

% number generation based on coefficient/exponent approach
for i=1:n
    
N = randi(10);
M = randi(10);

function1 =  N*(10^(-M));

functionmatrix(i) = function1;
Nmatrix(i) = N;
Mmatrix(i) = M;

end

sortedfunctionmatrix = sort(functionmatrix);    % sorts matrix by size

X = (1:n);

figure;
semilogy(X,sortedfunctionmatrix);     % graphs with logarithmic scale on y-axis
