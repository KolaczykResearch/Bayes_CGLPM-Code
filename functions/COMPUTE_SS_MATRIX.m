function [SS, T, lambda_psi] = COMPUTE_SS_MATRIX(Endo, Exo, lag)
% Compute the unnormalized covariance matrix of Y conditional on X.

% PRELIMINARY DATA PROCESSING
[nobs_1, n] = size(Endo);
[nobs_2, nz] = size(Exo);
if nobs_2 > nobs_1;
    Exo = Exo(1:nobs_1,:);
end
Data = [Endo, Exo];
Dt = Data;
nx = n + nz;

D = [];
for j = 0:lag    
    D = [D   Dt(lag+1-j:end-j,:)]; %#ok
end

[T, ~] = size(D);
D = (D - ones(T,1)*nanmean(D))./(ones(T,1)*nanstd(D));
Y = D(:, 1:n);
X = D(:, nx+1:end);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lag == 0;
    SS = Y'*Y;
    lambda_psi = 0;
else
    [~, lambda_psi] = GRID_SEARCH(Endo, Exo);
    YX = Y'*X;          XX = X'*X;          YY = Y'*Y; 
    k = size(XX,2);     A0 = zeros(n,k);    
    Psi_inv = lambda_psi*eye(k); 
    S_xx = XX + Psi_inv;
    S_yx = YX + A0*Psi_inv;
    S_yy = YY + A0*Psi_inv*A0'; 
    SS = S_yy - S_yx*(S_xx\S_yx');
end 