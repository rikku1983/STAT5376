% -------------------------------------------------------------------
%                     Demo for function Hfn_MDA
%              This fits model Fmat*Amat to data Ymat
%                       by parameter cascading
% -------------------------------------------------------------------

addpath('../fdaM')

%  define dimensions

N = 200;  %  number of cases
n =  5;  %  number of variables
K =  2;  %  number of principal components

%  ---------------  define true Amat_pop  --------------------

%  nonrandom Amat_pop ... two straight lines

Amat_pop = ones(K,n);
Amat_pop(2,:) = -2:2;

%  ---------------  define true Fmat_pop  --------------------

%  fixed Fmat_pop  ... points on a circle

tvec = (1:N)'*2*pi/N;
Fmat_pop = [sin(tvec),cos(tvec)];

%  -----------  normalize Fmat_pop and Amat_pop  -------------

Fmatssq  = sum(Fmat_pop.^2);
Fmat_pop = sqrt(N).*Fmat_pop./(ones(N,1)*sqrt(Fmatssq));
Amat_pop = (sqrt(Fmatssq')*ones(1,n)).*Amat_pop./sqrt(N);

%  ----------------  generate data  Ymat  --------------------

%  Ymat is model value plus error

Ymat = Fmat_pop*Amat_pop;
sigerr = 0.5;
Ymat = Ymat + sigerr.*randn(N,n);

%  -----------  generate initial values for Amat  ------------

%  starting Amat from eigenanalysis

[V,D] = eig(cov(Ymat));
[Dsort,Isort] = sort(diag(D),'descend');
Vsort = V(:,Isort);
avec0 = reshape(Vsort(:,1:K)',n*K,1);

avec0 = avec0 + randn(n*K,1);

%  assess starting values

Amat_est0  = reshape(avec0, K, n);

[Hval0, Hgrad0, Fmat_est0] = Hfn_MDA(avec0, Ymat);

sqrt(mean(Hgrad0.^2))

%  normalize Fmat_est0 and Amat_est0

Fmatssq0  = sum(Fmat_est0.^2);
Fmat_est0 = sqrt(N).*Fmat_est0./(ones(N,1)*sqrt(Fmatssq0));
Amat_est0 = (sqrt(Fmatssq0')*ones(1,n)).*Amat_est0./sqrt(N);

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

%  ---------  define roughness penalty struct PenStruct  ---------

PenStruct.FRowMat = eye(K);  
PenStruct.FColMat = eye(N);
% PenStruct.FColMat = eye(N) - Fmat_pop*inv(Fmat_pop'*Fmat_pop)*Fmat_pop'./N;
PenStruct.ARowMat = eye(n);
PenStruct.AColMat = eye(K);
PenStruct.ANrmPar = 1000;

%  ---  set options for optimization by function fminunc ----

options = optimset('LargeScale',  'off',    ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on',     ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  ---------  Minimize the fitting criterion  ---------------

PenStruct.FRowPar = 0;  %  smoothing parameter lambda for Rmat
PenStruct.FColPar = 0;  %  smoothing parameter lambda for Smat
% PenStruct.FColPar = 1000;  %  smoothing parameter lambda for Smat
PenStruct.ARowPar = 0;

parvec0 = cartes2sphere(avec0);

tic;
parvec = fminunc(@Hfn_MDA, parvec0, options, Ymat, PenStruct);
toc

avec = sphere2cartes(parvec);

Amat_est  = reshape(avec, K, n);

[Hval, Hgrad, Fmat_est] = Hfn_MDA(parvec, Ymat, PenStruct);

disp(['RMSQ gradient = ', num2str(sqrt(mean(Hgrad.^2)))])

%  normalize Fmat_est and Amat_est

Fmatssq  = sum(Fmat_est.^2);
Fmat_est = sqrt(N).*Fmat_est./(ones(N,1)*sqrt(Fmatssq));
Amat_est = (sqrt(Fmatssq')*ones(1,n)).*Amat_est./sqrt(N);

% canonical correlations between true and estimated Amat

[U1, S1, V1] = svd(Amat_pop');
[U2, S2, V2] = svd(Amat_est');
CCR = svd(U1(:,1:K)'*U2(:,1:K));
disp(['Canonical correlations between true and estimated A:   ', ...
      num2str(CCR')])

%  plot factor scores

figure(2)
subplot(1,2,1)
plot(Fmat_pop(:,1), Fmat_pop(:,2), 'b*')
xlabel('\fontsize{16} PC I')
ylabel('\fontsize{16} PC II')
axis([-2.5,2.5,-2.5,2.5])
axis('square')
subplot(1,2,2)
plot(Fmat_est(:,1), Fmat_est(:,2), 'bo')
xlabel('\fontsize{16} PC I')
ylabel('\fontsize{16} PC II')
axis('square')



 
