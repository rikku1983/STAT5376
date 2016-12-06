function [Hval, Hgrad, fmat] = ...
                  PCA_PLS(avec, xmat, ymat, gamval, lambda, Pmat1, Pmat2)
%  PCA_PLS evaluates the loss function and its gradient for a parameter
%  cascading approach to principal components analysis of 
%  multivariate N by n data YMAT using K << n factors.
%  Outer loss function is H(A|Y) = ||Y - F(A)A||^2
%  Inner loss function is J(F|A,Y) = ||Y - FA||^2 + lambda ||F'RF||
%  In this version, d
%  ||Y - F(A)A||^2 = (1-\gamval)*||xmat - F(A)A1||^2 +     
%                       \gamval *||ymat'[I-F(A)(F'(A)F(A))^{-1}F'(A)]||^2
%  for 0 <= \gamval <= 1
%  fitting ymat on the basis of scores that also depend on the fit to xmat
%
%  Arguments:
%  AVEC    ... Vectorized K by n principal component coefficient matrix
%  XMAT    ... N by n matrix of independent variable values
%  YMAT    ... N by p matrix of dependent variable values
%  P       ... Number of columns for ymat
%  GAMVAL  ... weight in [0,1] placed on fitting ymat
%  LAMBDA  ... smoothing parameter for score estimation.
%  PMAT1   ... Order K roughness penalty matrix
%
%  Last modified 12 June 2012

%  Set default values for penalized least squares estimation of PC scores

if nargin < 8, Pmat2  = []; end
if nargin < 7, Pmat1  = []; end
if nargin < 6, lambda = 0;  end
if nargin < 5, gamval = 1;  end

if gamval < 0 || gamval > 1
    error('GAMMAVAL is not within [0,1[.');
end

%  get dimensions of data and convert AVEC to K by n matrix AMAT

[N1,p] = size(ymat);
[N2,n] = size(xmat);
if N1 == N2
    N = N1;
else
    error('XMAT and YMAT do not have same number of rows.');
end
    
nK   = length(avec);
K    = nK/n;
amat = reshape(avec, K, n);

%  Compute N by K matrix FMAT containing smoothed PC scores

Mmat = amat*amat';
if lambda > 0
    Mmat = Mmat + lambda*Pmat1;
end
Mmatinv = inv(Mmat);
fmat = xmat*amat'*Mmatinv;

%  Evaluate the outer loss function H(A)

xhat  = fmat * amat;        %  fit to the data
rmat  = xmat - xhat;        %  residuals
Hval1 = sum(sum(rmat.^2));  % error sum of squares
fmatinvCP = inv(fmat'*fmat);
Q2    = eye(N) - fmat*fmatinvCP*fmat';
Hval2 = trace(ymat'*Q2*ymat);
Hval  = (1-gamval)*Hval1 + gamval*Hval2;

%  Evaluate the gradient if required

if nargout > 1
    Hgrad = zeros(nK,1);
    %  loop through the gradient elements
    m = 0;
    for q=1:n
        for p=1:K
            m = m + 1;
            %  Evaluate the derivative of Mmat wrt Amat(p,q)
            DMmatpq = zeros(K);
            DMmatpq(:,p) = amat(:,q);
            DMmatpq(p,:) = DMmatpq(p,:) + amat(:,q)';
            %  Evaluate the derivative of Fmat wrt Amat(p,q)
            DMmatinvpq = -Mmatinv*DMmatpq*Mmatinv;
            Dffacpq      = amat'*DMmatinvpq;
            Dffacpq(q,:) = Dffacpq(q,:) + Mmatinv(p,:);
            Dfmatpq      = xmat*Dffacpq;
            %  Evaluate the derativative of H wrt Amat(p,q)
            DHfacpq      = Dfmatpq*amat;
            DHfacpq(:,q) = DHfacpq(:,q) + fmat(:,p);
            Hgrad(m)     = -2*(1-gamval)*sum(sum(rmat.*DHfacpq)) + ...
              gamval*trace( ...
              ymat'*(fmat*fmatinvCP*(Dfmatpq'*fmat+fmat'*Dfmatpq)*fmatinvCP*fmat' ...
                   - Dfmatpq*fmatinvCP*fmat' - fmat*fmatinvCP*Dfmatpq')*ymat);
        end
    end   
end