function Proc = Procrustean_Analysis(U_0, U_1)
% Procrustes transformation of U_1 to U_0 (with U_0 as the target)
X = U_0;    
Y = U_1;
[nx, mx] = size(X);
[ny, my] = size(Y);

nz = min(nx, ny);
mz = max(mx, my);
if nx ~= ny
    X = X(1:nz,:); 
    Y = Y(1:nz,:);
end

% Center at the origin.
muX = mean(X,1);
muY = mean(Y,1);
Xc = X - repmat(muX, nz, 1);
Yc = Y - repmat(muY, nz, 1);
normX = sqrt(trace(Xc*Xc')); 
normY = sqrt(trace(Yc*Yc')); 

% Scale to equal (unit) norm.
X0 = Xc / normX;        Y0 = Yc / normY;
if my < mz
    Y0 = [Y0 zeros(nz, mz-my)];
end

% The optimum rotation matrix of Y.
[Uz, Lz, Vz] = svd(X0'*Y0);
A = Vz*Uz';

% The optimum scaling of Y.
b = trace(Lz) * normX / normY;

% symmetric Procrustes statistic (or distance between X and b*Y*A+c), ranges from 0-1
% (smaller values indicate better fit)
d = 1 - trace(Lz).^2;       
%d = 1 - (trace(Lz)^2)/(trace(Y0'*Y0)*trace(X0*X0'));

c = muX - b*muY*A;
Z = normX*trace(Lz)*Y0*A + repmat(muX, nz, 1);

numerator = sum(sum((X-Z).^2));
denominator = sum(sum(bsxfun(@minus,X,mean(X)).^2));
ratio = numerator/denominator;

Proc.Translation = c;
Proc.Scale = b;
Proc.Rotation = A;
Proc.Z = Z;
Proc.Distance = d;
Proc.ratio = ratio;