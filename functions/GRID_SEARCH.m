function [Score, lambda_psi] = GRID_SEARCH(Endo, Exo)
% Implements a grid search of the ridge parameter (lambda_psi)

lag = 1;
grid_rate = 0.1:0.1:10;
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

% Data Arrangement
[~,nvar] = size(Y);
k = size(X,2);

train_sz = round(0.8*nobs_1);
test_sz = T - train_sz;
Yr = Y(1:train_sz,:);       % dependent variable used for training the model
Yf = Y(train_sz+1:end,:);    % Out-of-sample observations: responses 

Xr = X(1:train_sz,:);       % Predictors used for training the model
Xf = X(train_sz+1:end,:);    % % Out-of-sample observations : Predictors 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A0 = zeros(nvar,k);
YX = Yr'*Xr;
XX = Xr'*Xr;
% Tk = sqrt(T*k);
Tk = k;
% Cross-validation setup (alpha hyperparameter selection with grid search)
gvalues = grid_rate*Tk;
num_grids = length(gvalues);
Score = zeros(num_grids,1);

for i=1:num_grids
    Psi_inv = gvalues(i)*eye(k);
    S_xx = XX + Psi_inv;
    S_yx = YX + A0*Psi_inv;
    A_post = S_yx/(S_xx); 
    
    %Store Predictions
    Y_pred = zeros(test_sz,nvar);
    for j=1:test_sz
        Y_pred(j,:)  = real(Xf(j,:)*A_post');
    end
    U = Yf - Y_pred;
    SFE = U.^2;
    MSFE = mean(SFE);
    Score(i) = mean(MSFE);
end

% Take first differences
Ds = diff(Score); 
tol = 1e-1;  
[id1, ~] = find(abs(Ds) <tol);
id = id1(1)+1;
lambda_psi = real(grid_rate(id(1)+1)*Tk);