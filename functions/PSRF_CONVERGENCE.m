function R = PSRF_CONVERGENCE(X)
% Compute the potential scale reduction factor (PSRF) convergence 
% diagnostics of the MCMC (Gelman and Rubin, 1992)

[N,M]=size(X);

if M == 1;
    K = fix(N/3);
    Y = zeros(K,2);
    for i=1:2
        if i==1;
            Y(:,i) = X(1:K,:);    
        else
            Y(:,i) = X(end-K+1:end,:);
        end
    end
    Z = Y;
else
    Z = X;
end
[N,M]=size(Z);

W = 0;
for n = 1:M
    x = Z(1:N,n) - repmat(mean(Z(1:N,n)),N,1);
    W = W + sum(x.*x);
end
W = W / ((N-1) * M);
    
% Calculate variances B (in fact B/n) of the means.
Bpn = 0;
MX = mean(Z);
m = mean(MX);
for n=1:M
    x = mean(Z(1:N,n)) - m;
    Bpn = Bpn + x.*x;
end
Bpn = Bpn / (M-1);
    
% Calculate reduction factors
S = (N-1)/N * W + Bpn;
R = (M+1)/M * S ./ W - (N-1)/M/N;
% V = R .* W;
R = real(sqrt(R));  