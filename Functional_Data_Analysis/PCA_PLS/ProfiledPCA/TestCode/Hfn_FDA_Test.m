addpath('../fdaM')

%% -------------------------------------------------------------------
%                 Test function Hfn_FDA, Multivariate data
% -------------------------------------------------------------------

%  define dimensions

N = 200;
n = 11;
K = 2;

%  Compute Grassman map

GrssMnPts = linspace(0,1,n*K)';
GrssBasis = create_fourier_basis([0,1],(n-K)*K);
GrssMnMap = eval_basis(GrssMnPts,GrssBasis);
GrssMnMap = GrssMnMap(:,1:(n-K)*K);
[GrssMnMap,S,V] = svd(GrssMnMap,0);

%  ---------------  define true Amat_pop  --------------------

%  nonrandom Amat_pop ... two straight lines

Amat_pop = ones(K,n);
Amat_pop(2,:) = -5:5;

%  ---------------  define true Fmat_pop  --------------------

%  fixed Fmat_pop  ... points on a circle

tvec = (1:N)'*2*pi/N;
Fmat_pop = [sin(tvec),cos(tvec)];

%  -----------  normalize Fmat_pop and Amat_pop  -------------

Fmatssq  = sum(Fmat_pop.^2);
Fmat_pop = sqrt(N).*Fmat_pop./(ones(N,1)*sqrt(Fmatssq));
Amat_pop = (sqrt(Fmatssq')*ones(1,n)).*Amat_pop./sqrt(N);

%  ----------------  generate data  Xmat  --------------------

%  Xmat is model value plus error

Xmat = Fmat_pop*Amat_pop;
sigerr = 0.5;
Xmat = Xmat + sigerr.*randn(N,n);

%  -----------  generate initial values for Amat  ------------

%  starting Amat from eigenanalysis

[V,D] = eig(cov(Xmat));
[Dsort,Isort] = sort(diag(D),'descend');
Vsort = V(:,Isort);
avec0 = reshape(Vsort(:,1:K)',n*K,1);

avec0 = avec0 + randn(n*K,1);

%  assess starting values

Amat_est0  = reshape(avec0, K, n);

%  ---------  define roughness penalty struct PenStruct  ---------

PenStruct.FRowMat = eye(K);  
% PenStruct.FColMat = eye(N);
PenStruct.FColMat = eye(N) - Fmat_pop*inv(Fmat_pop'*Fmat_pop)*Fmat_pop';
PenStruct.ARowMat = eye(K);
PenStruct.AColMat = eye(n);

radius    = 1;
anglewrd  = 0;
GrssMnWrd = 1;

if anglewrd
    if GrssMnWrd
        parvec0 = cartes2sphere(GrssMnMap'*avec0);
    else
        parvec0 = cartes2sphere(parvec0);
    end
else
    if GrssMnWrd
        parvec0 = GrssMnMap'*avec0;
    else
        parvec0 = avec0;
    end
end
        
[Hval0, Hgrad0, Fmat_est0] = ...
    Hfn_FDA(parvec0, Xmat, [], PenStruct, ...
            radius, anglewrd, GrssMnWrd, GrssMnMap);

sqrt(mean(Hgrad0.^2))

% canonical correlations between true and estimated Amat

[U1, S1, V1] = svd(Amat_pop');
[U2, S2, V2] = svd(Amat_est0');
CCR = svd(U1(:,1:K)'*U2(:,1:K));
disp(['Canonical correlations between true and estimated A:   ', ...
      num2str(CCR')])

%  plot initial factor scores

figure(1)
plot(Fmat_est0(:,1), Fmat_est0(:,2), 'bo', ...
     Fmat_pop(:,1),  Fmat_pop(:,2),  'ro')

%  ---  set options for optimization by function fminunc ----

options = optimset('LargeScale',  'off',    ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on',     ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  ---------  Minimize the fitting criterion  ---------------

PenStruct.FRowPar = 0;    %  smoothing parameter lambda for Rmat
PenStruct.FColPar = 2e0;  %  smoothing parameter lambda for Smat
PenStruct.ARowPar = 0;

radius    = n*K;
anglewrd  = 0;
GrssMnWrd = 0;

if anglewrd
    if GrssMnWrd
        parvec0 = cartes2sphere(GrssMnMap'*avec0);
    else
        parvec0 = cartes2sphere(parvec0);
    end
else
    if GrssMnWrd
        parvec0 = GrssMnMap'*avec0;
    else
        parvec0 = avec0;
    end
end
        
tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, Xmat, [], PenStruct, ...
                 radius, anglewrd, GrssMnWrd, GrssMnMap);
toc

[Hval, Hgrad, Fmat0_est] = ...
    Hfn_FDA(parvec, Xmat, [], PenStruct, ...
            radius, anglewrd, GrssMnWrd, GrssMnMap);
disp(['RMSQ gradient = ', num2str(sqrt(mean(Hgrad.^2)))])

if anglewrd
    if GrssMnWrd
        avec = GrssMnMap*sphere2cartes(parvec);
    else
        avec = sphere2cartes(parvec);
    end
else
    if GrssMnWrd
        avec = GrssMnMap*parvec;
    else
        avec = parvec;
    end
end

Amat_est  = reshape(avec, K, n);

% canonical correlations between true and estimated Amat

[U1, S1, V1] = svd(Amat_pop');
[U2, S2, V2] = svd(Amat_est');
CCR = svd(U1(:,1:K)'*U2(:,1:K));
disp(['Canonical correlations between true and estimated A:   ', ...
      num2str(CCR')])

%  plot factor scores

figure(2)
sclmat0 = inv(Fmat_pop'*Fmat_pop)*Fmat_pop'*Fmat0_est;
subplot(1,2,1)
% Fmatscl_est = Fmat_est*sclmat;
Fmat0scl_pop = Fmat_pop*sclmat0;
plot(Fmat0_est(:,1),    Fmat0_est(:,2),    'bo', ...
     Fmat0scl_pop(:,1), Fmat0scl_pop(:,2), 'b.')
xlabel('\fontsize{13} PC I')
ylabel('\fontsize{13} PC II')
% axis([-2.5,2.5,-2.5,2.5])
axis('square')

sclmat0 = inv(Fmat_pop'*Fmat_pop)*Fmat_pop'*Fmat0_est;
subplot(1,2,2)
Fmat0scl_pop = Fmat_pop*sclmat0;
plot(Fmat0_est(:,1),    Fmat0_est(:,2), 'bo', ...
     Fmat0scl_pop(:,1), Fmat0scl_pop(:,2), 'b.')
xlabel('\fontsize{13} PC I')
% ylabel('\fontsize{13} PC II')
axis('square')


%% -------------------------------------------------------------------
%                 Test function Hfn_FDA, Discrete functional data
% -------------------------------------------------------------------

%  define dimensions

N = 200;
n = 11;
K = 2;

%  Compute Grassman map

GrssMnPts = linspace(0,1,n*K)';
GrssBasis = create_fourier_basis([0,1],(n-K)*K);
GrssMnMap = eval_basis(GrssMnPts,GrssBasis);
GrssMnMap = GrssMnMap(:,1:(n-K)*K);
[GrssMnMap,S,V] = svd(GrssMnMap,0);

%  define Abasis

nAbasis = 9;
Abasis = create_bspline_basis([0,1], nAbasis);
Phimat = inprod(Abasis, Abasis);

%  -----------  generate data  -----------------

%  cmat0 random matrix

cmat0 = randn(nAbasis,K);

%  cmat a fixed matrix

cmat0 = ones(nAbasis,K);
cmat0(:,2) = (-4:4)';

%  functional data object for factor loadings

Afd0  = fd(cmat0,Abasis);

%  Fmat random

Fmat0 = randn(N,K);

%  Fmat patterned for K = 2

fvec  = (1:N)'*2*pi/N;
Fmat0 = [sin(fvec),cos(fvec)];

%  generate Xfd

nfine = 101;
tfine = linspace(0,1,nfine)';

Xcoef = cmat0*Fmat0';
Xfd0  = fd(Xcoef,Abasis);
Xmat0 = eval_fd(tfine, Xfd0);

%  generate random fourier functions for error

T = 1;
nfourierbasis = 21;
fourierbasis    = create_fourier_basis([0,T],nfourierbasis);
fourierbasismat = eval_basis(tfine,fourierbasis);
basiswt = exp(-(1:nfourierbasis)/1)';
ecoef  = randn(nfourierbasis,N).*(basiswt*ones(1,N));
efd    = fd(ecoef, fourierbasis);
emat   = eval_fd(tfine, efd);

sigerr = 5;
Xmat = Xmat0 + sigerr.*emat;

nXbasis = 101;
Xbasis  = create_bspline_basis([0,1],nXbasis);
XLfd    = int2Lfd(2);
Xlambda = 1e-4;
XfdPar  = fdPar(Xbasis, XLfd, Xlambda);
Xfd     = smooth_basis(tfine, Xmat, XfdPar);

plotfit_fd(Xmat, tfine, Xfd)

%  Generate random Ymat  

Ymat = randn(N,1);

%  generate initial values for cvec

%  true value

cvec0 = reshape(cmat0,nAbasis*K,1);

%  random normal

cvec0 = randn(nAbasis*K,1);

%  from eigenanalysis

% [V,D] = eig(cov(Xmat));
% [Dsort,Isort] = sort(diag(D),'descend');
% Vsort = V(:,Isort);
% avec0 = reshape(Vsort(:,1:K)',n*K,1);

%  smoothing parameters

Hlambda = 0;
Rmat    = eval_penalty(Abasis, 2);
Flambda = 0;
Smat    = [];

%  Load PenStruct with these variables

PenStruct.FRowMat = [];
PenStruct.FRowPar = 0;
PenStruct.FColMat = [];
% PenStruct.FColMat = eye(N) - Fmat0*inv(Fmat0'*Fmat0)*Fmat0';
PenStruct.FColPar = 0;
% PenStruct.FColPar = 1e-8;

PenStruct = [];

AfdPar = fdPar(Afd0, 2, Hlambda, 1, Rmat);

%  evaluate function and gradient

parvec0 = cartes2sphere(cvec0);

gamval =0;
tic;
[Hval0, Hgrad0] = Hfn_FDA(parvec0, Xfd, AfdPar, Ymat, ...
                          gamval, PenStruct, GrssMnMap);
toc

%  set options for optimization

options = optimset('LargeScale',  'off',  ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on', ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  -------------   fit functional PCA model  ----------------------

cvec0 = randn(nAbasis*K,1);

parvec0 = cartes2sphere(cvec0);

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, Xfd, AfdPar, Ymat, gamval, ...
                 PenStruct);
toc

[Hval, Hgrad, Fmat] = Hfn_FDA(parvec, Xfd, AfdPar, Ymat, gamval, ...
                              PenStruct);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

cvec    = sphere2cartes(parvec);
cmatest = reshape(cvec, K, nAbasis);
Afdest  = fd(cmatest',Abasis);

figure(1)
subplot(1,1,1)
plot(Afdest)

%  plot factor scores

figure(2)
subplot(1,1,1)
plot(Fmat(:,1),Fmat(:,2),'o')


% %  Histograms of columns of F
% 
% figure(2)
% for p=1:K
%     subplot(K,1,p)
%     hist(Fmat(:,p))
% end

%% -------------------------------------------------------------------
%                 Test function Hfn_FDA, Functional data
% -------------------------------------------------------------------

%  define dimensions

N = 20;  %  number of functions
K = 2;    %  number of dimensions of model

%  define Abasis

nabasis = 7;
Abasis = create_bspline_basis([0,1], nabasis);

%  Compute Grassman map, which in this case is 
%    NABASIS*K by (NABASIS-K)*K

GrssMnPts = linspace(0,1,nabasis*K)';
GrssBasis = create_fourier_basis([0,1],(nabasis-K)*K);
GrssMnMap = eval_basis(GrssMnPts,GrssBasis);
GrssMnMap = GrssMnMap(:,1:(nabasis-K)*K);
[GrssMnMap,S,V] = svd(GrssMnMap,0);


%  -----------  generate data  -----------------

%  random matrix of rank n

% Xcoef = randn(nybasis,N);
% Xfd = fd(Xcoef, Xbasis);

%  random coef0 matrix

coef0 = randn(nabasis,K);

%  designed coef0 matrix

coef0 = ones(nabasis,K);
coef0(:,2) = -3:3;

Afd0 = fd(coef0, Abasis);

%  Fmat random

Fmat0 = randn(N,K);

%  Fmat patterned for K = 2

tvec  = (1:N)'*2*pi/N;
Fmat0 = [sin(tvec),cos(tvec)];
Fmat_pop = Fmat0;

Ymat0 = Fmat0(:,1);
Ymat  = Ymat0 + randn(N,1);

%  generate data for PCA in N by NFINE matrix XMAT 

nfine  = 101;
tfine  = linspace(0,1,nfine)';
Amat0  = eval_fd(tfine, Afd0);
Phimat = eval_basis(tfine, Abasis);

Xmat0  = Fmat0*Amat0';
sigerr = 0.1;
Xmat   = Xmat0 + sigerr.*randn(N,nfine);
Xbasis = create_bspline_basis([0,1],11);
Xfd    = smooth_basis(tfine, Xmat', Xbasis);

%  generate inner product matrices  (not needed at this time)

% Wmat = inprod(Abasis,Abasis);
% Zmat = inprod(Xfd,   Abasis);

%  generate initial values for cvec

%  true value

cvec0 = reshape(coef0,nabasis*K,1);

%  random normal

cvec0 = randn(nabasis*K,1);

%  Load PenStruct with these variables

PenStruct.FRowMat = [];
PenStruct.FRowPar = 0;
PenStruct.FColMat = [];
PenStruct.FColMat = eye(N) - Fmat0*inv(Fmat0'*Fmat0)*Fmat0';
PenStruct.FColPar = 0;
PenStruct.FColPar = 1e-2;

%  fdPar object for factor coefficient functions

Hlambda = 0;
Rmat    = eval_penalty(Abasis, 2);
Cmat0   = reshape(cvec0, nabasis, K);
Afd0    = fd(Cmat0, Abasis);
AfdPar  = fdPar(Afd0, 2, Hlambda, 1, Rmat);

%  evaluate function and gradient (disable angular transformation first)

parvec0 = cvec0;
parvec0 = cartes2sphere(cvec0);
parvec0 = GrssMnMap'*cvec0;
parvec0 = cartes2sphere(GrssMnMap'*cvec0);

Ymat      = [];
gamval    = 0;
radius    = 1;
anglewrd  = 0;
GrssMnWrd = 1;
[Hval0, Hgrad0] = Hfn_FDA(parvec0, Xfd, AfdPar, ...
                          Ymat, gamval, PenStruct, ...
                          radius, anglewrd, GrssMnWrd, GrssMnMap);

h0 = 1e-1;
toler = 1e-6;

grdhat = grdHessm_FDA(@Hfn_FDA, parvec0, h0, toler, Xfd, AfdPar, ...
                      [], gamval, PenStruct, radius, anglewrd, ...
                      GrssMnWrd, GrssMnMap);

[Hgrad0, grdhat]
max(abs(Hgrad0-grdhat))

%  set options for optimization

options = optimset('LargeScale',  'off',  ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on', ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6,     ...
                   'RelLineSrchBnd', 1e1);

%  -------------   fit functional PCA model  ----------------------

cvec0 = randn(nabasis*K,1);

gamval = 1.0;

Ymat      = [];
gamval    = 0;
radius    = 1;
anglewrd  = 1;
GrssMnWrd = 1;

parvec0 = cvec0;
parvec0 = cartes2sphere(cvec0);
parvec0 = GrssMnMap'*cvec0;
parvec0 = cartes2sphere(GrssMnMap'*cvec0);

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, Xfd, AfdPar, ...
                 Ymat, gamval, PenStruct, ...
                 radius, anglewrd, GrssMnWrd, GrssMnMap);
toc

[Hval, Hgrad, Fmat_est] = Hfn_FDA(parvec, Xfd, AfdPar, ...
                                  Ymat, gamval, PenStruct, ...
                                  radius, anglewrd, GrssMnWrd, GrssMnMap);
                              
% [Hval, Hgrad, Fmat] = Hfn_FDA(cvec, Xfd, AfdPar, ...
%                                    Ymat, gamval, PenStruct);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

cvec = sphere2cartes(parvec);
cvec = GrssMnMap*sphere2cartes(parvec);

coefest = reshape(cvec, K, nabasis);
Afdest  = fd(coefest',Abasis);
Amatest = eval_fd(tfine, Afdest);

figure(1)
subplot(1,1,1)
plot(Afdest)

%  plot factor scores

figure(2)
subplot(1,1,1)
plot(Fmat_est(:,1),Fmat_est(:,2),'o')

% canonical correlations between true and estimated Fmat

[U1, S1, V1] = svd(Fmat_pop,0);
[U2, S2, V2] = svd(Fmat_est,0);
CCR = svd(U1'*U2);
disp(['Canonical correlations between true and estimated A:   ', ...
      num2str(CCR')])

