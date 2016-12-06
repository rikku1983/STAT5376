function [Fval, Fgrad, Fhess, Fcross, fvec, mu] = ...
              PCA_MDA_GLM(Y, DCell, mu0, Amat, fvec0, Wtvec, PenStruct)
%PCA_MDA_GLM does a parameter cascaded principal component analysis of a 
%   data matrix with regularization where the loss function for each
%   variable varies among the GLM family.
%  Arguments:
%
%  Y       ... an n by 2 matrix of data to be fitted.  The first column
%              contains the values of the observations and only in the
%              binomial case, the second colum contains the values M, the
%              binomial sample sizes.
%  DCELL   ... a cell array of length n containing information about the
%              distribution of each variable.  Each entry is a struct
%              object with these fields:
%                family:  a string indicating which of the three GLM family 
%                  is assumed:
%                    'normal' or 'gaussian' or 'Gaussian'
%                    'binomial' or 'binary' or 'Bernoulli'
%                    'poisson'
%                devFn:   deviance function
%                stdFn:   variance function
%                linkFn:  link function data expectation mu to canonical 
%                         parameter eta
%                DlinkFn: link function derivative
%                IlinkFn: inverse of link function mapping eta to mu
%  MU0     ... an n-vector of initial values for data expectations
%  AMAT    ... an k by n matrix of values of factor coefficients 
%              
%  FVEC0   ... starting values for regression coefficients
%  WTVEC   ... a vector of prior weights, such as the inverses of the
%              relative variance of each observation.
%  PENSTRUCT ... A struct object containing specifications of penalty terms
%                and their associated smoothing parameter vectors
%
%  Returns:
%  BVEC     ... Final estimate of coefficients
%  DEVIANCE ... Deviance values
%
%  Last modified 16 June 2014 by Jim Ramsay

%--------------------------------------------------------------------------
%                    Check arguments
%--------------------------------------------------------------------------

if nargin < 7,  PenStruct = [];  end

if isempty(PenStruct)
    FRowMat = [];
    FRowPar = 0;
else
    FRowMat = PenStruct.FRowMat;
    FRowPar = PenStruct.FRowPar;
    %  check dimensions of roughness penalty matrices
    if ~isempty(FRowMat) && any(size(FRowMat) ~= K)  
        error('FRowMat is not of order K');  
    end
end

%--------------------------------------------------------------------------
%                   Initialize eta 
%--------------------------------------------------------------------------

[K,n] = size(Amat);

% compute eta = E(y) from mu

mu   = mu0;
eta  = zeros(n,1);
Deta = zeros(n,1);
stdm = zeros(n,1);
for j=1:n
    DStructj = DCell{j};
    linkFnj  = DStructj.linkFn;
    eta(j)   = linkFnj(mu0(j));
end

%--------------------------------------------------------------------------
%                        Set up for iterations
%--------------------------------------------------------------------------

iter     = 0;
iterLim  = 100;
seps     = sqrt(eps);
convcrit = 1e-6;
sqrtwt   = sqrt(Wtvec);

%  set up starting value Bvec0 if required

if isempty(fvec0)
    fvec0 = zeros(nbasis,ncurve);
end
fvec = fvec0;

%--------------------------------------------------------------------------
%                       Start of GLM iteration loop
%--------------------------------------------------------------------------

while iter <= iterLim
    iter = iter+1;
    
    % Compute adjusted dependent variable for least squares fit
    
    for j=1:n
        DStructj = DCell{j};
        DlinkFnj = DStructj.DlinkFn;
        stdFnj   = DStructj.stdFn;
        muj      = mu(j);
        Deta(j)  = DlinkFnj(muj);
        if strcmp(DStructj.family,'binomial')
            stdm(j) = stdFnj(muj,Y(j,2));
        else
            stdm(j) = stdFnj(muj);
        end
    end
    Zvec = eta + (Y(:,1) - mu).*Deta;
  
    % Compute IRLS weights the inverse of the variance function
    
    sqrtw = sqrtwt./(abs(Deta).*stdm);
    
    % Compute coefficient estimates for this iteration - the IRLS step
    
    fvec_old = fvec;
    Yw       = Zvec.*sqrtw;
    Xmatw    = Amat'.*(sqrtwt*ones(1,K));
    Mmat     = Xmatw'*Xmatw;
    if FRowPar > 0 && ~isempty(FRowMat)
        Mmat = Mmat + FRowPar*FRowMat;
    end
    fvec     = Mmat\Xmatw'*Yw;
    eta      = Amat'*fvec;
    for j=1:n
        DStructj = DCell{j};
        IlinkFnj = DStructj.IlinkFn;
        mu(j) = IlinkFnj(eta(j));
    end
    
    % Force mean in bounds, in case the linkFn function is faulty    
    
    for j=1:n
        DStructj = DCell{j};
        switch DStructj.family
            case 'binomial'
                muLims = [eps 1-eps];
                if mu(j) < muLims(1) || muLims(2) < mu(j)
                    mu(j) = max(min(mu(j),muLims(2)),muLims(1));
                end
            case 'poisson'
                muLims = realmin.^.25;
                if mu(j) < muLims(1)
                    mu(j) = max(mu(j),muLims(1));
                end
        end
    end

    % Check stopping conditions
    
    crit = max(abs(fvec-fvec_old));
    disp([iter,crit])
    if crit < convcrit*max(abs(fvec_old)) 
        break; 
    end
    
end

%--------------------------------------------------------------------------
%                    end of GLM iteration loop
%--------------------------------------------------------------------------

if iter > iterLim
    warning(['smooth_GLM:','Iteration'],'Iteration limit reached.');
end

Fval   = 0;
Fgrad  = zeros(K,1);
Fhess  = zeros(K,K);
Fcross = zeros(n*K,K);
Gradt  = zeros(n,1);
Hesst  = zeros(n,n);
for j=1:n;
    DStructj = DCell{j};
    devFnj   = DStructj.devFn;
    if strcmp(DStructj.family,'binomial')
        devj = devFnj(mu(j),Y(j,1),Y(j,2));
    else
        devj = devFnj(mu(j),Y(j,1));
    end
    Fval     = Fval + Wtvec(j).*devj;
    DStructj = DCell{j};
    switch DStructj.family
        case 'normal'
            Gradt(j)   = Y(j,1) - mu(j);
            Hesst(j,j) = 1;
        case 'binomial'
            Gradt(j)   = Y(j,2).*(Y(j,1) - mu(j));
            Hesst(j,j) = mu(j).*(1 - mu(j));
        case 'poisson'
            Gradt(j)   = Y(j,1).*log(mu(j)) - mu(j);
            Hesst(j,j) = mu(j);
    end
end
Fgrad  = -2.*Wtvec(j).*Amat*Gradt;
Fhess  = +2.*Wtvec(j).*Amat*Hesst*Amat';
Fcross = -2.*Wtvec(j).*kron(Gradt,eye(K));

