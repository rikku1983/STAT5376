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

%  normalize Fmat_pop and Amat_pop

Fmatssq  = sum(Fmat_pop.^2);
Fmat_pop = sqrt(N).*Fmat_pop./(ones(N,1)*sqrt(Fmatssq));
Amat_pop = (sqrt(Fmatssq')*ones(1,n)).*Amat_pop./sqrt(N);

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

%  evaluate function and gradient

[Hval0, Hgrad0] = Hfn_MDA(avec0, Ymat, lambda_1, Rmat, lambda_2, Smat);

disp(['Hval0 = ',num2str(Hval0)])

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;
Hgradhat = grdHessm_MDA(@Hfn_MDA, avec0, h0, toler, Ymat, ...
                        lambda_1, Rmat, lambda_2, Smat);

%  compare

[Hgrad0, Hgradhat]

%  set options for optimization

options = optimset('LargeScale',  'off',  ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on', ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  fit covariance surface

lambda = 1e1;

tic;
avec = fminunc(@Hfn_MDA, avec0, options, Ymat, lambda, Rmat);
toc

Amat_est  = reshape(avec, K, n);

[Hval, Hgrad, Fmat_est] = Hfn_MDA(avec, Ymat, lambda, Rmat);

sqrt(mean(Hgrad.^2))

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

figure(1)
plot(Fmat_est(:,1), Fmat_est(:,2), 'bo', ...
     Fmat_pop(:,1), Fmat_pop(:,2), 'ro')

