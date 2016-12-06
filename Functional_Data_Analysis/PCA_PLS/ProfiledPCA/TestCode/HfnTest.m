% -------------------------------------------------------------------
%                     test function Hfn_MDA
% -------------------------------------------------------------------

addpath('../fdaM')

%  define dimensions

N = 50;
n = 5;
K = 2;

%  ---------------  define true Amat  --------------------

%  random Amat_pop

Amat_pop = randn(K,n);

%  nonrandom Amat_pop ... two straight lines

Amat_pop = ones(K,n);
Amat_pop(2,:) = -2:2;

%  ---------------  define true Fmat  --------------------

%  random Fmat_pop

Fmat_pop = randn(N,K);

%  fixed Fmat_pop  ... points on a circle

tvec = (1:N)'*2*pi/N;
Fmat_pop = [sin(tvec),cos(tvec)];

%  fixed Fmat_pop  ... two horizontal lines

Fmat_pop = [ones(N,1), -ones(N,1)];

%  --------------   generate data  Ymat   ----------------

%  random matrix of rank n

% Ymat = randn(N,n);

%  Ymat is model value plus error

Ymat = Fmat_pop*Amat_pop;
sigerr = 0.1;
Ymat = Ymat + sigerr.*randn(N,n);

%  center the data

Ymat = Ymat - ones(N,1)*mean(Ymat);

%  ------  generate initial values for A  -------------

%  true value plus random noise

avec0 = reshape(Amat_pop,n*K,1) + sigerr*randn(n*K, 1);

%  true value

avec0 = reshape(Amat_pop,n*K,1);

%  starting Amat from eigenanalysis

[V,D] = eig(cov(Ymat));
[Dsort,Isort] = sort(diag(D),'descend');
Vsort = V(:,Isort);
avec0 = reshape(Vsort(:,1:K)',n*K,1);

%  smoothing parameters

lambda = 1;

%  Rmat is projection on to complement of one(K,1)

Rmat = eye(K) - ones(K,1)*ones(K,1)'./K;

%  evaluate function and gradient

[Hval0, Hgrad0] = Hfn_MDA(avec0, Ymat, lambda, Rmat);

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;
Hgradhat = grdHessm_MDA(@Hfn_MDA, avec0, h0, toler, Ymat, lambda, Rmat);

%  compare

[Hgrad0, Hgradhat]

%  set options for optimization

options = optimset('LargeScale',  'on',  ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on', ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  fit covariance surface

lambda = 1e-1;

tic;
avec = fminunc(@Hfn_MDA, avec0, options, Ymat, lambda, Rmat);
toc

Amatest  = reshape(avec, K, n);

[Hval, Hgrad, Fmat] = Hfn_MDA(avec, Ymat, lambda, Rmat);

sqrt(mean(Hgrad.^2))

disp('A:')
disp([Amat',Amatest'])

%  Histograms of columns of F

figure(1)
for p=1:K
    subplot(K,1,p)
    hist(Fmat(:,p))
end

%  plot factor scores

figure(2)
plot(Fmat(:,1),Fmat(:,2),'o')

% -------------------------------------------------------------------
%                 Test function Hfn_FDA, Discrete data
% -------------------------------------------------------------------

addpath('../fdaM')

%  define dimensions

N = 200;
n = 11;
K = 2;

%  define Abasis

tvec = linspace(0,1,n);
nbasis = 9;
Abasis = create_bspline_basis([0,1], nbasis);
Phimat = eval_basis(tvec, Abasis);

%  -----------  generate data  -----------------

%  random matrix of rank n

% Ymat = randn(N,n);

%  cmat random matrix

cmat0 = randn(K,nbasis);
Amat0 = cmat*Phimat';

%  cmat a fixed matrix

cmat0 = ones(K,nbasis);
cmat0(2,:) = -4:4;
Amat0 = cmat0*Phimat';

%  Fmat random

Fmat0 = randn(N,K);

%  Fmat patterned for K = 2

fvec  = (1:N)'*2*pi/N;
Fmat0 = [sin(fvec),cos(fvec)];

%  generate Ymat 

Ymat0 = Fmat0*Amat0;
sigerr = 0.5;
Ymat = Ymat0 + sigerr.*randn(N,n);

%  generate initial values for cvec

%  true value

cvec0 = reshape(cmat,nbasis*K,1);

%  random normal

cvec0 = randn(nbasis*K,1);

%  from eigenanalysis

% [V,D] = eig(cov(Ymat));
% [Dsort,Isort] = sort(diag(D),'descend');
% Vsort = V(:,Isort);
% avec0 = reshape(Vsort(:,1:K)',n*K,1);

%  smoothing parameters

Hlambda = 0e-0;
Rmat    = eval_penalty(Abasis, 2);
Flambda = 0;
Smat    = [];

%  Load fitstruct with these variables

fitstruct.tfine   = tvec; 
fitstruct.Ymat    = Ymat;
fitstruct.Phimat  = Phimat;
fitstruct.Hlambda = Hlambda; 
fitstruct.Flambda = Flambda; 
fitstruct.Rmat    = Rmat;
fitstruct.Smat    = Smat;

%  evaluate function and gradient

[Hval0, Hgrad0] = Hfn_FDA(cvec0, fitstruct);

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;
Hgradhat = grdHessm_FDA(@Hfn_FDA, cvec0, h0, toler, fitstruct);

%  compare

[Hgrad0, Hgradhat]

%  set options for optimization

options = optimset('LargeScale',  'on',  ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on', ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  -------------   fit functional PCA model  ----------------------

Flambda = 0;
Smat   = eye(K);

cvec0 = randn(nbasis*K,1);

tic;
cvec = fminunc(@Hfn_FDA, cvec0, options, fitstruct);
toc

[Hval, Hgrad, Fmat] = Hfn_FDA(cvec, fitstruct);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

cmatest = reshape(cvec, K, nbasis);
Afdest  = fd(cmatest',Abasis);
Amatest = eval_fd(tvec, Afdest);

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

% -------------------------------------------------------------------
%                 Test function Hfn_FDA, Functional data
% -------------------------------------------------------------------

addpath('../fdaM')

%  define dimensions

N = 200;
nybasis = 13;
Ybasis = create_bspline_basis([0,1], nybasis);
K = 2;

%  define Abasis

nabasis = 7;
Abasis = create_bspline_basis([0,1], nabasis);

%  -----------  generate data  -----------------

%  random matrix of rank n

% Ycoef = randn(nybasis,N);
% Yfd = fd(Ycoef, Ybasis);

%  random coef0 matrix

coef0 = randn(nabasis,K);

%  designed coef0 matrix

coef0 = ones(nabasis,K);
coef0(:,2) = -3:3;

Afd = fd(coef0, Abasis);

%  Fmat random

Fmat0 = randn(N,K);

%  Fmat patterned for K = 2

tvec  = (1:N)'*2*pi/N;
Fmat0 = [sin(tvec),cos(tvec)];

%  generate Ymat 

nfine = 101;
tfine = linspace(0,1,nfine)';
Amat  = eval_fd(tfine, Afd);
Ymat0 = Fmat0*Amat';
sigerr = 0.0;
Ymat = Ymat0 + sigerr.*randn(N,nfine);
Ybasis = create_bspline_basis([0,1],11);
Yfd = smooth_basis(tfine, Ymat', Ybasis);

%  generate inner product matrices

Wmat = inprod(Abasis,Abasis);
Xmat = inprod(Yfd, Abasis);

%  generate initial values for cvec

%  true value

cvec0 = reshape(coef0,nabasis*K,1);

%  random normal

cvec0 = randn(nabasis*K,1);

%  from eigenanalysis

[U,S,V] = svd(Ymat');
[Sort,Isort] = sort(diag(S),'descend');
Usort = U(:,Isort);
Ufd  = smooth_basis(tfine, Usort, Ybasis);
cvec0 = reshape(getcoef(Vfd)',n*K,1);

%  smoothing parameters

Hlambda = 0e-0;
Rmat    = eval_penalty(Abasis, 2);
Flambda = 0;
Smat    = [];

%  Load fitstruct with these variables

fitstruct.tfine   = tvec; 
fitstruct.Ymat    = Ymat;
fitstruct.Xmat    = Xmat;
fitstruct.Phimat  = Wmat;
fitstruct.Hlambda = Hlambda; 
fitstruct.Flambda = Flambda; 
fitstruct.Rmat    = Rmat;
fitstruct.Smat    = Smat;

%  evaluate function and gradient

[Hval0, Hgrad0] = Hfn_FDA(cvec0, fitstruct);

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;
Hgradhat = grdHessm_FDA(@Hfn_FDA, cvec0, h0, toler, fitstruct);

%  compare

[Hgrad0, Hgradhat]

%  set options for optimization

options = optimset('LargeScale',  'on',  ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on', ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  -------------   fit functional PCA model  ----------------------

Flambda = 0;
Smat   = eye(K);

cvec0 = randn(nabasis*K,1);

tic;
cvec = fminunc(@Hfn_FDA, cvec0, options, fitstruct);
toc

[Hval, Hgrad, Fmat] = Hfn_FDA(cvec, fitstruct);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

coefest = reshape(cvec, K, nabasis);
Afdest  = fd(coefest',Abasis);
Amatest = eval_fd(tvec, Afdest);

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

% -------------------------------------------------------------------
%                 test function PCA_Reg
% -------------------------------------------------------------------

addpath('../fdaM')

%  test function 

nfine  = 101;
tfine  = linspace(0,1,nfine)';

%  unwarped version:  sine and cosine

sinfine = sin(2*pi*tfine);
cosfine = cos(2*pi*tfine);

%  set up warping function

cwarp = 1.2;
hfine = (1 - exp(cwarp*tfine))./(1 - exp(cwarp));

plot(tfine, hfine, 'b-', tfine, tfine, 'b--')

%  factor scores for two principal components

fvec = ones(2,1);

%  warped version to act as a target

y0fine = fvec(1).*sin(2*pi*hfine) + fvec(2).*cos(2*pi*hfine);

%  set up a basis for principal component functions

PCbasis = create_bspline_basis([0,1],11);

%  set up single principal component function as smooth of
%  unwarped spline

PCfd = smooth_basis(tfine, [sinfine, cosfine], PCbasis);

%  set up a basis for W functions defining warping functions h

Wbasis = create_bspline_basis([0,1],5);

%  compute W coefficients that exactly correct the warp

%  smooth t as a function of h to get inverse function
hinvfd = smooth_basis(hfine, tfine, PCbasis);
plot(hinvfd)
%  derivative of inverse warping function
Dhinvfine = eval_fd(tfine, hinvfd, 1);
%  log-derivative = W
wfine = log(Dhinvfine);
plot(tfine, wfine, 'b-', [0,1], [0,0], 'b--')
%  set up true W function
Wfd_pop = smooth_basis(tfine, wfine, Wbasis);
plotfit_fd(wfine, tfine, Wfd_pop)
%  coefficient vector for true W function
dvec_pop = getcoef(Wfd_pop);

%  set Kmat, periodic and crit

Kmat = [];
periodic = 0;
crit = 2;

%  Load fitstruct with these variables

fitstruct.tfine    = tfine;
fitstruct.yfine   = y0fine;
fitstruct.PCfd     = PCfd;
fitstruct.Wbasis   = Wbasis;
fitstruct.crit     = crit;
fitstruct.periodic = periodic;
fitstruct.Kmat     = Kmat;

%  criterion value for random coefficient

dvec0 = randn(5,1);

Fval = PCA_Reg(dvec0, fitstruct);
disp(Fval)

h0 = 1;
toler = 1e-6;

[Fval, Fgrad] = PCA_Reg(dvec0, fitstruct);

Fgradhat = grdHessm_Reg(@PCA_Reg, dvec0, h0, toler, fitstruct); 

[Fgrad, Fgradhat]

[Fval, Fgrad, Fhess] = PCA_Reg(dvec0, fitstruct);

[Fgradhat, Fhesshat] = grdHessm_Reg(@PCA_Reg, dvec0, h0, toler, fitstruct); 

Fhess 
Fhesshat

%  check full gradient and hessian version

dcvec0 = [dvec0; randn(22,1)];

Fval = PCA_Reg_Full(dcvec0, fitstruct);

[Fval, Fgrad] = PCA_Reg_Full(dcvec0, fitstruct);

Fgradhat = grdHessm_Reg(@PCA_Reg_Full, dcvec0, h0, toler, fitstruct); 

[Fgrad, Fgradhat]

[Fval, Fgrad, Fhess] = PCA_Reg_Full(dcvec0, fitstruct);

[Fgradhat, Fhesshat] = grdHessm_Reg(@PCA_Reg_Full, dcvec0, ...
                                    h0, toler, fitstruct); 

% Fhess   (1:5,1:5) 
% Fhesshat(1:5,1:5)

Fhess   (1:5,6:16) 
Fhesshat(1:5,6:16)

Fhess   (1:5,17:27) 
Fhesshat(1:5,17:27)

%  test cross-derivatives for PCA_Reg

cvec0 = dcvec(6:27);
Cmat = reshape(cvec0,11,2);
PCfd = putcoef(PCfd, Cmat);
fitstruct.PCfd = PCfd;

[Fval, Fgrad, Fhess, Fcross] = PCA_Reg(dvec, fitstruct);

Fcross  (1:5,1:11) 
Fhesshat(1:5,6:16)

fitstruct.dvec = dvec0;
cvec0 = reshape(getcoef(PCfd)',nAbasis*K,1);

%  set options for optimization

options = optimset('LargeScale',  'off',   ...
                   'Display',     'iter', ...
                   'Diagnostics', 'on',   ...
                   'GradObj',     'on',   ...
                   'Hessian',     'on',   ...               
                   'TolFun',      1e-6);

%  -------------   fit y0fine  ----------------------

tic;
dvec = fminunc(@PCA_Reg, dvec0, options, fitstruct);
toc

[Fval, Fgrad, Fhess, fitstruct, PCreg, Wfd] = PCA_Reg(dvec, fitstruct);
disp(Fval)
disp(Fgrad)

PCfine = eval_fd(tfine, PCreg);

figure(1)
plot(tfine, PCfine, 'b-', tfine, [sinfine, cosfine], 'r-')

figure(2)
plot(Wfd)

% -------------------------------------------------------------------
%                 test function Hfn_Reg
% -------------------------------------------------------------------

addpath('../fdaM')

%  define dimensions

N = 50;
K = 2;
nfine  = 101;
tfine  = linspace(0,1,nfine)';

%  define Abasis

nAbasis = 7;
Abasis = create_bspline_basis([0,1], nAbasis);

%   ----  generate coefficient matrix for coefficient fns.  ------

%  Acoefmat random matrix

Acoefmat0 = randn(nAbasis,K);

%  Acoefmat a fixed matrix

Acoefmat0 = ones(nAbasis,K);
Acoefmat0(:,2) = -3:3;
          
PCfd0 = fd(Acoefmat0, Abasis);

plot(PCfd0)

Amat0 = eval_fd(tfine, PCfd0);

%  define Wbasis

% nWbasis = 3;
% Wbasis = create_monomial_basis([0,1],nWbasis);

nWbasis = 5;
Wbasis = create_bspline_basis([0,1],nWbasis);

%   ---------------  generate coefficient scores  ------------------

%  Fmat patterned for K = 2

tvec  = (1:N)'*2*pi/N;
Fmat0 = [sin(tvec)';cos(tvec)'];

plot(Fmat0(1,:),Fmat0(2,:),'o')

%   ---------------  generate warping functions  ------------------

%  set up warping function

sigwarp = 0.3;
Wwarp = randn(1,N)*sigwarp;

hfinemat = (1 - exp(tfine*Wwarp))./(1 - exp(ones(nfine,1)*Wwarp));

plot(tfine, hfinemat, '-')

%  ----------------------  generate data  -------------------------

Ymat0 = zeros(nfine,N);
for i=1:N
    Ymat0(:,i) = eval_fd(hfinemat(:,i), PCfd0)*Fmat0(:,i);
end

figure(1)
subplot(2,1,1)
plot(tfine, Ymat0', '-')

%  generate noisy Ymat with smooth error process

T = 1;
nfourierbasis = 21;
fourierbasis  = create_fourier_basis([0,T],nfourierbasis);
fourierbasismat = eval_basis(tfine,fourierbasis);
basiswt = exp(-(1:nfourierbasis)/1)';
ecoef  = randn(nfourierbasis,N).*(basiswt*ones(1,N));
efd    = fd(ecoef, fourierbasis);
emat   = eval_fd(tfine, efd);

sigerr = 0.1;
Ymat = Ymat0 + sigerr.*emat;

subplot(2,1,2)
plot(tfine, Ymat', '-')

nybasis = 53;
Ybasis = create_bspline_basis([0,1], nybasis);
Yfd = smooth_basis(tfine, Ymat, Ybasis);

%  ----------  generate initial values for parameters  ---------------

Dmat0 = zeros(nWbasis,N);

%  true value

Cmat0 = Acoefmat0';
cvec0 = reshape(Cmat0',nAbasis*K,1);

%  random value

cvec0 = randn(K*nAbasis,1);

%  from eigenanalysis

pcastruct = pca_fd(Yfd, K);
harmmat   = eval_fd(tfine,pcastruct.harmfd);
PCfd      = smooth_basis(tfine, harmmat, Abasis);
coefmat   = getcoef(PCfd);
cvec0     = reshape(coefmat,nAbasis*K,1);

cvec = cvec0;  %  used for stepping through function

%  ------------------  smoothing parameters  -----------------------

Hlambda = 0;
Rmat    = eval_penalty(Abasis, 2);
Flambda = 0;
Smat    = [];
Wlambda = 0;
Kmat    = [];

%  set options for optimization

Joptions = optimset('LargeScale',  'off',  ...
                    'Display',     'off',   ...
                    'Diagnostics', 'off', ...
                    'GradObj',     'on',     ...
                    'Hessian',     'off',    ...               
                    'TolFun',      1e-6);
               
%  Load fitstruct with these variables

global fitstruct

fitstruct.tfine    = tfine; 
fitstruct.Ymat     = Ymat;
fitstruct.PCfd     = PCfd;
fitstruct.Hlambda  = Hlambda; 
fitstruct.Flambda  = Flambda; 
fitstruct.Wlambda  = Wlambda;
fitstruct.Rmat     = Rmat;
fitstruct.Smat     = Smat;
fitstruct.Kmat     = Kmat;
fitstruct.Dmat     = Dmat0;
fitstruct.Abasis   = Abasis;
fitstruct.Wbasis   = Wbasis;
fitstruct.crit     = 2;
fitstruct.periodic = 0;
fitstruct.Kmat     = [];
fitstruct.options  = Joptions;

%  ------------  check gradient and hessian of PCA_Reg  --------------

%  evaluate function, gradient, and hessian

dvec = Dmat0(:,1);

tic;
[Fval, Fgrad, Fhess] = PCA_Reg(dvec);
toc

h0 = 1e0;
toler = 1e-6;

tic;
[Fgradhat, Fhesshat] = grdHessm_Reg(@PCA_Reg, dvec, h0, toler); 
toc

[Fgrad, Fgradhat]

Fhess

Fhesshat

Joptions = optimset('LargeScale',  'off',  ...
                    'Display',     'off',   ...
                    'Diagnostics', 'off', ...
                    'GradObj',     'on',     ...
                    'Hessian',     'off',    ...               
                    'TolFun',      1e-6);
fitstruct.options = Joptions;

i=1;
fitstruct.yfine = Ymat(:,i);
dvec0 = Dmat0(:,i);

tic;
dvec = fminunc(@PCA_Reg, dvec0);
toc

tic;
[Fval, Fgrad, Fhess] = PCA_Reg(dvec);
toc

Fval
Fgrad'
Fhess
eig(Fhess)'

%  ------------  check gradient of Hfn_reg  --------------------------

%  evaluate function and gradient

tic;
[Hval, Hgrad] = Hfn_Reg(cvec0);
toc

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;

tic;
Hgradhat = grdHessm_Reg(@Hfn_Reg, cvec0, h0, toler); 
toc

[Hgrad, Hgradhat]

Hgradmat    = reshape(Hgrad,K,nAbasis)
Hgradhatmat = reshape(Hgradhat,nAbasis,K)'

%  -------------   fit functional PCA model  ----------------------

%  set smoothing for factor scores

Flambda = 0;
Smat    = [];
fitstruct.Flambda = Flambda;
fitstruct.Smat    = Smat;

%  set smoothing for warping functions

Wlambda = 1e-2;
Kmat    = eval_penalty(Wbasis,2);
fitstruct.Wlambda = Wlambda;
fitstruct.Kmat    = Kmat;

Hoptions = optimset('LargeScale',   'off',     ...
                    'Display',      'iter',   ...
                    'Diagnostics',  'on',     ...
                    'GradObj',      'on',     ...
                    'Hessian',      'off',    ...               
                    'TolFun',       1e-5);

tic;
cvec = fminunc(@Hfn_Reg, cvec0, Hoptions, fitstruct);
toc

[Hval, Hgrad, Fmat, fitstruct] = Hfn_Reg(cvec, fitstruct);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

Acoefmatest = reshape(cvec, K, nAbasis);
Afdest  = fd(Acoefmatest',Abasis);
Amatest = eval_fd(tfine, Afdest);

figure(1)
subplot(1,1,1)
plot(tfine, Amatest, 'b-', tfine, Amat0, 'b--')

%  plot factor scores

figure(2)
subplot(1,1,1)
plot(Fmat(1,:),Fmat(2,:),'o')

%  plot warping functions

Dmat = fitstruct.Dmat;
Wfd  = fd(Dmat, Wbasis);

plotdeform_fd(Wfd)

Yhat = fitstruct.Yhat;

figure(3)
for i=1:N
    subplot(2,1,1)
    plotdeform_fd(Wfd(i))
    hold on
    d0i = hfinemat(:,i) - tfine;
    plot(tfine, d0i, 'b--')
    hold off
    axis([0,1,-0.1,0.1])
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} d(t) = h(t) - t')
    title(['\fontsize{13} Row ',num2str(i)])
    subplot(2,1,2)
    plot(tfine, Yhat(:,i), 'b-', tfine, Ymat0(:,i), 'b--', ...
         [0,1], [0,0], 'r:')
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} y(t)')    
    pause
end


