function M = GRAPH_ACCURACY(AdjTrue, M)

tau = [0.25 0.5 0.75];
PostG = M.PostG;
n = size(AdjTrue,2);
nt = length(tau);
TP = zeros(nt,1);
FP = zeros(nt,1);
TN = zeros(nt,1);
FN = zeros(nt,1);
indmx = reshape([1:n^2],n,n); %#ok
upperind = indmx(triu(indmx,1)>0); 
for j = 1:length(tau)
    t = tau(j); 
    Pred_Adj = double(PostG > t); 
    TP(j) = sum(Pred_Adj(upperind)~=0 & AdjTrue(upperind)~=0);  % True positive 
    FP(j) = sum(Pred_Adj(upperind)~=0 & AdjTrue(upperind)==0);  % False positve
    TN(j) = sum(Pred_Adj(upperind)==0 & AdjTrue(upperind)==0);  % True negative
    FN(j) = sum(Pred_Adj(upperind)==0 & AdjTrue(upperind)~=0);  % False negative  
end
NN = TP + FP + TN + FN;
ACC  = 1e2*(TP + TN)./(NN); 
SP = TN./(TN + FP); % Specificity (SP%)
SE = TP./(TP + FN); % Sensitivity (SE%)
TPR = SE;
FPR = 1 - SP;

% add trivial end-points
tp = [0 ; TPR ; 1];
fp = [0 ; FPR ; 1];

[~, id] = sort(fp);
fpr = fp(id);
tpr = tp(id);

M.TP  = TP(2);
M.FP  = FP(2);
M.FN  = FN(2);
M.TN  = TN(2);
M.ACC = ACC(2);
M.AUC = 1e2.*trapz(fpr,tpr);