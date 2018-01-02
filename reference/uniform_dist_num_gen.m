% uniform continuous distribution

n = 1000;

amatrix = zeros(n,1);
bmatrix = zeros(n,1);
rmatrix = zeros(n,1);


for i=1:n
    
a = 10e-1;
b = 1e-10;

r = a + (b-a)*rand; %  online this was a recommended way to generate random uniformly distributed numbers between two endpoints     

amatrix(i) = a;
bmatrix(i) = b;
rmatrix(i) = r;

end

X = (1:n);

sortedrmatrix = sort(rmatrix);

figure;
semilogy(X, sortedrmatrix);