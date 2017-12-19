function [Endo, Exo, SigTrue]  = Generate_Data(n, T, lag)
% Generate the data for the simulation experiment in Section 3. 
 
a = 0.3;  b = 0.9;
B0 = sign(randn(n)).*(a + (b-a)*rand(n));
BG = double(rand(n)<0.2);
B = tril(BG,-1).*B0 + 0.5*eye(n);
B = B + B'; w = eig(B); delta = (n*min(w)-max(w))/(1-n);
SigTrue = B + delta*eye(n);

Exo = [];
Endo = randn(T,n)*chol(SigTrue);
if lag ~= 0; 
    A0 = sign(randn(n,1)).*(a + (b-a)*rand(n,1));
    A = diag(A0);
    ch_coeff = sum(abs(eig(A))>1); 
    while ch_coeff > 0
        A0 = sign(randn(n,1)).*(a + (b-a)*rand(n,1));
        A = diag(A0);
        ch_coeff = sum(abs(eig(A))>1); 
    end    
    T0 = 0.5*T; 
    Endo = randn(T0+T, n);
    U = randn(T0+T,n)*chol(SigTrue);
    % Now generate data from VAR
    for t = lag+1:T0+T   
        Endo(t,:) = A*Endo(t-1,:)' + U(t,:)';
    end
    Endo = Endo(T0+1:end,:);
end