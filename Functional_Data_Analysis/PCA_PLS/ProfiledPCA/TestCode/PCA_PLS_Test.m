% -------------------------------------------------------------------
%                     test function PCA_PLS
% -------------------------------------------------------------------

addpath('../fdaM')

%  define dimensions

N = 50;
n = 8;
p = 1;
K = 2;

%  ---------------  define true Amat  --------------------

%  random Amat1_pop and Amat2_pop 

Amat_pop = randn(K,n);

%  nonrandom Amat_pop ... two straight lines + constant

Amat_pop = ones(K,n);
Amat_pop(2,:) = linspace(-1,1,n)';

%  ---------------  define true Fmat  --------------------

%  random Fmat_pop

Fmat_pop = randn(N,K);

%  fixed Fmat_pop  ... points on a circle

tvec = (1:N)'*2*pi/N;
Fmat_pop = [sin(tvec),cos(tvec)];

%  normalize Fmat_pop and Amat_pop

Fmatssq  = sum(Fmat_pop.^2);
Fmat_pop = sqrt(N).*Fmat_pop./(ones(N,1)*sqrt(Fmatssq));
Amat_pop = (sqrt(Fmatssq')*ones(1,n)).*Amat_pop./sqrt(N);

%  --------------   generate data  Ymat   ----------------

%  random matrix of rank n

%  xmat and ymat are model value plus error

xmat_pop = Fmat_pop*Amat_pop;
sigerr1 = 0.2;
xmat = xmat_pop + sigerr1.*randn(N,n);

sigerr2 = 1;
Bmat_pop1 = 0.5*ones(K,1)./K;
Bmat_pop2 = 2.0*ones(n,1)./n;
Gmat_pop = (eye(N) - Fmat_pop*inv(Fmat_pop'*Fmat_pop)*Fmat_pop')*xmat_pop;
ymat = Fmat_pop*Bmat_pop1 + Gmat_pop*Bmat_pop2 + sigerr2.*randn(N,p);

%  errorless example n = 8

tvec   = (1:N)'*2*pi/N;
circle = [sin(tvec),cos(tvec)];
poly   = zeros(n,4);
svec   = linspace(-1,1,N)';
xmat   = zeros(N,8);
wtvec  = [4,4,1,1];
m2 = 0;
for j=1:4
    m1 = m2 + 1;
    m2 = m2 + 2;
    polyk = svec.^(j-1) * ones(1,2);
    xmat(:,m1:m2) = circle.*polyk.*wtvec(j);
end
beta = [1;1;4;4];
sigerr2 = 4;
ymat = xmat(:,1:4)*beta + sigerr2*randn(N,1);

%  center the data

xmat = xmat - ones(N,1)*mean(xmat);
ymat = ymat - ones(N,1)*mean(ymat);

%  ------  generate initial values for A  -------------

%  true value plus random noise

avec0 = reshape(Amat_pop,n*K,1) + 0.1*randn(n*K, 1);

%  true value

avec0 = reshape(Amat_pop,n*K,1);

%  starting Amat from eigenanalysis

[V,D] = eig(cov(Ymat));
[Dsort,Isort] = sort(diag(D),'descend');
Vsort = V(:,Isort);
avec0 = reshape(Vsort(:,1:K)',n*K,1);

%  smoothing parameters

lambda = 0;
Pmat1  = [];
gamval = 0.5;

%  evaluate function and gradient

[Hval0, Hgrad0] = PCA_PLS(avec0, xmat, ymat, gamval, lambda, Pmat1);

disp(['Hval0 = ',num2str(Hval0)])

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;
Hgradhat = grdHessm_MDA2(@PCA_PLS, avec0, h0, toler, xmat, ymat, ...
                         gamval, lambda, Pmat1);

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

gamval = 0.5;

%  compute starting values from de Jong and Kiers paper

Gmat = gamval*xmat*xmat' + (1 - gamval)*ymat*ymat';

[V,D] = eig(Gmat);       
D = diag(D);
[Dsort, Isort] = sort(D,'descend');
Vsort = V(:,Isort);
Fmat = Vsort(:,1:K);
amat = Fmat'*xmat;
avec0 = reshape(amat,n*K,1);


tic;
avec = fminunc(@PCA_PLS, avec0, options, xmat, ymat, gamval, lambda, Pmat1);
toc

Amat_est  = reshape(avec, K, n);

disp(Amat_est)

[Hval, Hgrad, Fmat_est] = PCA_PLS(avec, xmat, ymat, gamval, lambda, Pmat1);

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

Bmat = Fmat_est\ymat;
disp(Bmat')
Yhat2 = Fmat_est*Bmat;

figure(2)
plot(Yhat2, ymat, 'o')

disp(['Correlation for ymat = ',num2str(corr(Yhat2,ymat))])

%  PLS analysis

[T,P,Q,B,B0] = pls1(xmat, ymat, K);
[XL,YL,XS,YS,betapls] = plsregress(xmat, ymat, K);

[U2, S2, V2] = svd(P);
CCR = svd(U1(:,1:K)'*U2(:,1:K));
disp(['Canonical correlations between true and estimated A:   ', ...
      num2str(CCR')])

yhat = xmat*B + B0;
yhat = xmat*betapls(2:9);

disp(['Correlation for PLS = ',num2str(corr(yhat,ymat))])


yhat = xmat*inv(xmat'*xmat)*xmat'*ymat;
disp(['Correlation for regression = ',num2str(corr(yhat,ymat))])



