function Z = Sample_Z_FC(G,EZ)
% Sampling the latent Z following Hoff's rstiefel R package

n = size(G,1);
indmx = reshape(1:n^2,n,n); 
upperind = indmx(triu(indmx,1)>0); 
y = G(upperind);        
ez = EZ(upperind);       
ly = length(y);
lb = -Inf(ly,1) ;       
ub = Inf(ly,1);
lb(y==1) = 0;           
ub(y==0) = 0;

temp_lb = normcdf(lb, ez, 1);   
temp_ub = normcdf(ub, ez, 1);
prob = unifrnd(temp_lb, temp_ub, ly, 1);
prob(prob==1)=0.9999;  prob(prob==0)=0.0001;
z = norminv(prob,ez,1);
Z = zeros(n);
Z(upperind) = z ; 
Z = Z + Z';
Z(1:size(Z,2)+1:end) = normrnd(diag(EZ),sqrt(2));