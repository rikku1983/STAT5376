% -------------------------------------------------------------------
%                     Test the GLM version of PCA
% -------------------------------------------------------------------


addpath('../fdaM')

%  define dimensions

N = 50;
n = 5;
K = 2;

%  ---------------  define true Amat  --------------------

%  random Amat_pop

Amat_pop = randn(K,n);

%  ---------------  define true Fmat  --------------------

%  random Fmat_pop

Fmat_pop = randn(N,K);

%  --------------   generate data  Ymat   ----------------

%  Ymat is model value plus error

Ymat0 = Fmat_pop*Amat_pop;
sigerr = 0.1;
Ymat = Ymat0 + sigerr.*randn(N,n);

%  center the data

Ymat = Ymat - ones(N,1)*mean(Ymat);

Y     = Ymat(1,:)';
Wtvec = ones(n,1);
fvec0 = zeros(K,1);
mu0   = Y;
Amat  = Amat_pop;
PenStruct = [];

DCell = cell(n,1);
%  set up for normal
DStruct.family = 'normal';
DStruct.devFn   = @(mu,Ymat) (Ymat - mu).^2;
DStruct.stdFn   = @(mu)  ones(size(mu));
DStruct.linkFn  = @(mu)  mu;
DStruct.DlinkFn = @(mu)  ones(size(mu));
DStruct.IlinkFn = @(eta) eta;
for j=2:n
    DCell{j} = DStruct;
end

%  set up for binomial

DStruct.family = 'binomial';
DStruct.devFn   = @(mu,Ymat,M) 2*M.*(Ymat.*log((Ymat+(Ymat==0))./mu) + ...
    (1-Ymat).*log((1-Ymat+(Ymat==1))./(1-mu)));
DStruct.stdFn   = @(mu,M)  sqrt(mu.*(1-mu)./M);
DStruct.linkFn  = @(mu)  log(mu./(1-mu));
DStruct.DlinkFn = @(mu)  1./(mu.*(1-mu));
DStruct.IlinkFn = @(eta) 1./(1 + exp(-constrain(eta,log(eps),-log(eps))));
DCell{1} = DStruct;

Y = Ymat(1,:)';
Y = [Y ones(n,1)];
Y(1,1) = 0.5;
mu0(1) = 0.5;

%  set up for poisson

DStruct.family = 'poisson';
DStruct.devFn   = @(mu,Ymat) 2*(Ymat.*(log((Ymat+(Ymat==0))./mu)) - ...
    (Ymat - mu));
DStruct.stdFn   = @(mu)  sqrt(mu);
DStruct.linkFn  = @(mu)  log(mu);
DStruct.DlinkFn = @(mu)  1./mu;
DStruct.IlinkFn = @(eta) exp(constrain(eta,log(eps),-log(eps)));
DCell{1} = DStruct;

Y = Ymat(1,:)';
Y(1,1) = 2;
mu0(1) = 1;

[Fval, Fgrad, Fhess, Fcross, fvec, mu] = ...
              PCA_MDA_GLM(Y, DCell, mu0, Amat, fvec0, Wtvec, PenStruct);
