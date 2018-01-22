function U = Matrix_Bingham_VMF_Gibbs(U, Z, L) 
% Perform Gibbs sampling of the latent positions (U) from a matrix 
% Bingham-von Mises-Fisher distribution.

[~,R] = size(U);
A = Z./2;
pr = randperm(R);
for i = 1:R
    r = pr(i);
    cid = 1:R;  cid(r)=[];
    U_r = U(:,cid);
    % Obtain null space of U_r: Nr = null(U_r');
    Nr = null(U_r');
    Al = L(r,r)*Nr'*A*Nr;
    x0 = Nr'*U(:,r);
    x = Vector_Bingham_VMF_Gibbs(Al,x0);
    U(:,r) = Nr*x;
end

end

function xs = Vector_Bingham_VMF_Gibbs(Al,x0)

% Do Gibbs sampling for a vector Bingham-von Mises-Fisher distribution.
% Generates nSamples starting from x0.  If x0 is not given, a default is
% used.  If nSamples isn't provided it defaults to 1.
%
% The vector Bingham-von Mises-Fisher distribution has pdf:
%    p(x) \propto exp(c.'*x + 0.5*x.'(A + A.')x) for ||x|| = 1
%
% This function implements the Gibbs sampler described by Peter Hoff,
% "Simulation of the matrix Bingham-von Mises-Fisher distribution..." (2009)
c = 0;
nDims = size(Al,1);
[E,L] = eig(Al);
Et = E.';
L = diag(L);
d = Et*c;

y = Et*x0;
s = sign(y);
y2 = y.^2;

for i = randperm(nDims)
    q = y2./(1-y2(i));
    subi = [1:(i-1) (i+1):nDims];
    aa = L(i) - (q(subi).'*L(subi));
    bb = sum(s(subi).*sqrt(q(subi)).*d(subi)');
    theta = sampleThetaGibbs(y2(i), (nDims-3)/2, aa, bb, d(i));

    y2(i) = theta;
    y2(subi) = (1-theta).*q(subi);    
    y2 = y2./sum(y2); % prevent drift
end
xs = Et.'*(s.*sqrt(y2));

end

function theta = sampleThetaGibbs(theta,k,a,b,c)
    d1=1/2;
    if k > min([a,b,c])
        d2 = k;
    else
        d2=1+min([k, max( [k-a,-1/2] )] );
    end
    
    betaBlockSize = 1;  % original=20
    th = betarnd(d1,d2,[betaBlockSize,1]);
    logpth = computeLogPTheta(th,k,a,b,c);
    logptheta = computeLogPTheta(theta,k,a,b,c);
    logbth = -betaNegLikelihood(th,d1,d2,false);
    logbtheta = -betaNegLikelihood(theta,d1,d2,false);
    for i = 1:betaBlockSize
        acc = log(rand) < (logpth(i) - logptheta) + (logbtheta - logbth(i));
        if acc
            theta = th(i);
            logptheta = logpth(i);
            logbtheta = logbth(i);
        end
    end
end

function logptheta = computeLogPTheta(th,k,a,b,c)
    logptheta = -0.5*log(th) + k*log(1-th) + th.*a + (1-th).^(1/2) .* b + ...
        log(exp(-c*th.^0.5) + exp(c*th.^0.5));
end

function [f, g] = betaNegLikelihood(x,a,b,normalized)
%[f g] = betaNegLikelihood(x,a,b)
%
% Evaluate the negative log likelihood of x under a beta distribution 
% with parameters a and b.
%
% E[X] = a/(a+b)
% E[(X-E[X])^2] = a*b/((a+b)^2*(a+b+1))
%
% f = -log( x^(a-1) * (1-x)^(b-1) / beta(a,b) )
% g = df/dx
    if nargin == 2 && length(a) == 2
        b = a(2);
        a = a(1);
    end

    f = zeros(size(x));
    posX = x > 0 & x < 1;
    f(posX) = -(a-1)*log( x(posX) ) - (b-1)*log( (1-x(posX)) );
    if normalized
        f(posX) = f(posX) + betaln(a,b);
    end
    f(~posX) = inf;
    if nargout > 1
        g = zeros(size(x));
        g(posX) = -(a-1)./x(posX) + (b-1)./(1-x(posX));
    end
end