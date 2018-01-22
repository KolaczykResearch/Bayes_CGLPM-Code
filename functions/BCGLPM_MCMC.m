function M = BCGLPM_MCMC(Endo, Exo, lag, nsimu, h, v_0)
% MCMC Sampling the Bayesian covariance graphical structure (G) and the 
% latent positions (U) of the nodes.

% SET NECESSARY PARAMETERS
dimU = 2;      % Dimension of co-ordinate system
burnin = 0.3*nsimu; % number of discarded samples

% INITIALIZE
[SS, T, lambda_psi] = COMPUTE_SS_MATRIX(Endo, Exo, lag);
[~,Sig] = cov2corr(SS);
Sig0 = Sig;
n = size(Sig,1);
S = Sig*T;
C = inv(Sig);
vii = 1;   

% Predefined parameters for GIG simulation 
G_a = 1; 
G_q = 1-T/2;     

V = ones(n);
V((1:n+1:n*n))=vii.^2;
ind_noi_all = zeros(n-1,n);
for i = 1:n
    cid = 1:n;
    cid(cid==i)=[];
    ind_noi_all(:,i) = cid;
end

G =  zeros(n);
G_mat = zeros(n);
Pr_prob = zeros(n);
Pst_prob = zeros(n);
Pst_mat = zeros(n);

nsave = nsimu - burnin;

NLogL = zeros(nsimu,1);
Edge = zeros(nsimu,1);
Sigma_mat = zeros(n);
Omega_mat = zeros(n);
U_mat = zeros(n,dimU);
L_save = zeros(nsimu, dimU);
Th_save = zeros(nsimu,1);

fprintf('\n#######')
fprintf(' SAMPLING BAYESIAN COVARIANCE GRAPH MODEL ')
fprintf('#########\n')

% starting values for latent position model
%%% THETA
num_off_diag = n*(n-1)/2;
t2_theta = 100;
v_theta = 1/(1/t2_theta + num_off_diag);
theta_0 = -0.5;
theta = theta_0;
th_div_t2 = (theta_0/t2_theta);

%%% LAMBDA
t2_lambda = n ;
lam_var = 2*t2_lambda./(2+t2_lambda); 
L = zeros(dimU);

%%% INITIALIZE U
X0 = randn(n,dimU);
[tmp.vec, tmp.val] = eig(X0'*X0);
U = X0*(tmp.vec*diag(sqrt(1./diag(tmp.val)))*(tmp.vec)');

Z = theta + U*L*U';

v0 = v_0^2;
v1 = h^2*v0;
v0_vc = v0*ones(n-1, 1);
v1_vc = v1*ones(n-1, 1);

tic;
for iter = 1:nsimu
    
    % UPDATE EDGE INCLUSION PRIOR PROBABILITIES
    Pr_mat = normcdf(Z);
    
    for i = 1:n
        ind_noi = ind_noi_all(:,i);  %ind_noi = negative i
        v_ind_noi = V(ind_noi,i);
        
        C11 = C(ind_noi,ind_noi);
        C12 = C(ind_noi,i);
        S11 = S(ind_noi,ind_noi); 
        S12 = S(ind_noi,i);         
        
        invSig11 = C11 - C12*C12'./C(i,i);
        invSig11S12 = invSig11*S12;
        S11invSig11 = S11*invSig11;
        
        W1 = invSig11*S11invSig11;
        
        %====== v from generalized inverse gaussian (GIG)
        u_s = Sig(ind_noi,i);   
        G_b = u_s'*W1*u_s - 2*u_s'*invSig11S12 + S(i,i);  
        gam = gigrnd(G_q, G_a, G_b, 1);     % v ~ GIG(q, a, b)
        
        %===== Sample mu
        inv_Dv = diag(1./v_ind_noi);
        B = W1/gam + V(i,i)*invSig11 + inv_Dv;
        B_chol = chol(B);
        mu_mn = B_chol\(B_chol'\invSig11S12)./gam;  
        mu = mu_mn + B_chol\randn(n-1,1);
        
        pii = Pr_mat(ind_noi,i); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w_1 = (1./sqrt(v1_vc)).*exp(-0.5*(mu.^2)./(v1_vc)).*pii;
        w_0 = (1./sqrt(v0_vc)).*exp(-0.5*(mu.^2)./(v0_vc)).*(1-pii);
        ost = w_1./(w_1 + w_0);
        
        dwr = rand(n-1,1);
        gij = (dwr < ost);
        
        % UPDATE V
        v_t = v0_vc;
        v_t(gij) = v1_vc(gij);
        V(ind_noi,i) = v_t;
        V(i,ind_noi) = v_t;
        
        % UPDATE G
        G(i,ind_noi) = gij;
        Pst_mat(i,ind_noi) = ost;
        
        % Updating Covariance and Precision matrix 
        Sig(ind_noi,i) = mu;
        Sig(i,ind_noi) = mu; 
        Sig(i,i) = gam + mu'*invSig11*mu;
        
        invSig11mu = invSig11*mu;
        C(ind_noi,ind_noi) = invSig11 + invSig11mu*invSig11mu'./gam;
        C12 = -invSig11mu/gam;
        C(ind_noi,i) = C12;
        C(i,ind_noi) = C12';
        C(i,i) = 1/gam;
    end
    G = double((G + G')>0);
    
    % SAMPLE Z FROM NORMAL DENSITIES WITH MEAN ULU'
    EZ = theta + U*L*U';
    Z = Sample_Z_FC(G, EZ); 
    
    % SAMPLE THETA FROM NORMAL CONDITIONAL DISTRIBUTION
    Et = Z - U*L*U';
    Up_Et = triu(Et,1);
    theta_pst = v_theta*(th_div_t2 + sum(Up_Et(:)));
    theta = theta_pst + randn*v_theta;
    
    % SAMPLE L FROM NORMAL CONDITIONAL DISTRIBUTION
    ZTh = Z - theta;
    lam_mn = lam_var*diag((U'*ZTh*U)./2);
    L = diag(lam_mn + randn(dimU,1)*sqrt(lam_var)); 
    
    % SAMPLE U FROM FULL CONDITIONAL UNDER BvMF(U|ZTh/2, L, C=0)
    U = Matrix_Bingham_VMF_Gibbs(U, ZTh, L);  %*** 

    if(mod(iter,100)==0)
      fprintf('  %4.0f out of %4.0f ', iter, nsimu);
      toc;
    end
    % COMPUTE LOG PREDICTIVE DENSITY
    Zul = (ZTh - U*L*U')'*(ZTh - U*L*U');
    score = T*( -log(det(C)) + sum(diag(C*Sig0))) + 0.5*sum(diag(Zul));
    NLogL(iter) = score;
    Edge(iter) = nnz(G)/2;
    L_save(iter,:) = diag(L);
    Th_save(iter) = theta;
    
    if iter > burnin
        %rl = iter - burnin;
        G_mat = G_mat + G;
        U_mat = U_mat + U;
        Sigma_mat = Sigma_mat + Sig;
        Omega_mat = Omega_mat + C;
        Pr_prob = Pr_prob + Pr_mat;
        Pst_prob = Pst_prob + Pst_mat;
    end
end
Time = toc;
M.U = U_mat./nsave;
M.L = mean(L_save(burnin+1:end,:));
M.L_pst = L_save;
M.Lambda_psi = lambda_psi;
M.Theta = mean(Th_save(burnin+1:end));
M.Theta_pst = Th_save;
M.PostG = G_mat./nsave;
M.Pr_prob = Pr_prob./nsave;
M.Pst_prob = Pst_prob./nsave;
M.Sigma = Sigma_mat./nsave;
M.Omega = Omega_mat./nsave;
M.Time = Time;
M.NLogL = NLogL;
M.Edge = Edge;
M.PSRF_G = PSRF_CONVERGENCE(M.NLogL(burnin+1:end));