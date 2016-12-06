%%  -----------------------------------------------------------------------
%                    Berkeley Growth Data
%  -----------------------------------------------------------------------

addpath('../fdaM')

% load results from previous analysis

load BerkeleyGrowth

%%  ----------------------------------------------------------------------------
%                   Set up the data
%  ----------------------------------------------------------------------------

% -----------------  Female data  --------------------

index = 1:ncasef;
N = length(index);

nfine = 1001;
tfine = linspace(3,18,nfine)';

Yfine = zeros(nfine,N);
m = 0;
for i = index
    m = m + 1;
    beta = betaf(:,i);
    Yfine(:,m) = beta(2).*eval_mon(tfine, Wfdf(i), 2);
end

%  sum of squares for Y

SSY = sum(sum(Yfine.^2));

%  plot data

figure(1)
subplot(1,1,1)
Ymeanfine = mean(Yfine(:,index),2);
phdl = plot(tfine, Yfine(:,index), 'b-', [3,18], [0,0], 'b:');
set(phdl, 'LineWidth', 1)
% lhdl = line(tfine, Ymeanfine);
% set(lhdl, 'LineWidth', 2, 'LineStyle', '--', 'color', 'b')
xlabel('\fontsize{13} Age (years)')
ylabel('\fontsize{13} Acceleration (cm/yr^2)')
% title('\fontsize{13} Female acceleration curves from monotone smooth')
axis([3,18,-4,3])

%%  -----------------------------------------------------------------
%                PCA using function Hfn_FDA
%  -----------------------------------------------------------------

%  define dimensions

agerng = [3,18];

%  define Ybasis

nybasis = 53;
Ybasis = create_bspline_basis(agerng, nybasis);
lambda = 1e-8;
YfdPar = fdPar(Ybasis, 2, lambda);

%  set up the data

Yfd = smooth_basis(tfine, Yfine, YfdPar);

plotfit_fd(Yfine, tfine, Yfd)

% residuals of the order of 1e-5

%  define Abasis

nAbasis = 23;
Abasis = create_bspline_basis(agerng, nAbasis);

%  Compute Grassman map, which in this case is 
%    NABASIS*K by (NABASIS-K)*K

GrssMnPts = linspace(0,1,nAbasis*K)';
GrssBasis = create_fourier_basis([0,1],(nAbasis-K)*K);
GrssMnMap = eval_basis(GrssMnPts,GrssBasis);
GrssMnMap = GrssMnMap(:,1:(nAbasis-K)*K);
[GrssMnMap,S,V] = svd(GrssMnMap,0);


%  ----------  generate initial values for parameters  ---------------

K = 1;
K = 2;
K = 3;
K = 4;

%  cvec0 from eigenanalysis

eigval  = eig(Yfine'*Yfine./N);
eigval  = sort(eigval, 'descend');
varprop = eigval./sum(eigval);
disp(cumsum(varprop(1:5)))
%     0.6716
%     0.8703
%     0.9263
%     0.9653
%     0.9789

pcastruct = pca_fd(Yfd, K, Abasis, 0);
harmmat   = eval_fd(tfine,pcastruct.harmfd);
PCfd0     = smooth_basis(tfine, harmmat, Abasis);
coefmat   = getcoef(PCfd0);
cvecPCA     = reshape(coefmat,nAbasis*K,1);
% cvec0     = randn(nAbasis*K,1);

pcaFmat = pcastruct.harmscr;

pcastruct_1 = pcastruct;
save pcastruct_1 pcastruct_1
pcastruct_2 = pcastruct;
save pcastruct_2 pcastruct_2
pcastruct_3 = pcastruct;
save pcastruct_3 pcastruct_3
pcastruct_4 = pcastruct;
save pcastruct_4 pcastruct_4

Hlambda = 0;
Rmat    = [];
AfdPar = fdPar(PCfd0, 2, Hlambda, 1, Rmat);

PenStruct.FRowMat = [];
PenStruct.FRowPar = 0;
PenStruct.FColMat = [];
PenStruct.FColMat = Pmat;
PenStruct.FColPar = 0;
PenStruct.FColPar = 1e2;

%  set options for optimization

options = optimset('LargeScale',  'off',    ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on',     ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  carry out the pca

% cvec0 = randn(length(cvecPCA),1);
cvec0 = cvecPCA;

Ymat      = [];
gamval    = 0;
radius    = 1;
anglewrd  = 1;
GrssMnWrd = 1;

parvec0 = cvec0;
% parvec0 = cartes2sphere(cvec0);
% parvec0 = GrssMnMap'*cvec0;
parvec0 = cartes2sphere(GrssMnMap'*cvec0);  %  Default

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, Yfd, AfdPar, ...
                 [], 0, PenStruct, 1, 0, 0);
toc

[Hval, Hgrad, Fmat] = Hfn_FDA(parvec, Yfd, AfdPar, [], 0, PenStruct, 1, 0, 0);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

cvec = parvec;
% cvec = sphere2cartes(parvec,1);
% cvec = GrssMnMap*parvec;
cvec = GrssMnMap*sphere2cartes(parvec);  %  Default

Hfn_FDA_cvec_1 = cvec;
save Hfn_FDA_cvec_1 Hfn_FDA_cvec_1 
Hfn_FDA_cvec_2 = cvec;
save Hfn_FDA_cvec_2 Hfn_FDA_cvec_2 
Hfn_FDA_cvec_3 = cvec;
save Hfn_FDA_cvec_3 Hfn_FDA_cvec_3 
Hfn_FDA_cvec_4 = cvec;
save Hfn_FDA_cvec_4 Hfn_FDA_cvec_4 

Afdest  = fd(reshape(cvec, nAbasis, K),Abasis);
Amat    = eval_fd(tfine, Afdest);

%  compute approximation to Y and residuals if required

Yhatmat = Fmat*Amat';
Rmat    = Ymat - Yhatmat';

%  plot principal component functions

nfine = 51;
tfine = linspace(3,18,nfine)';
Amatest = eval_fd(tfine,Afdest);

figure(1)
subplot(1,1,1)
phdl=plot(tfine, -Amatest(:,1), 'b-', tfine, -Amatest(:,2), 'b--', ...
          [3,18],[0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} D^2 Height (cm/sec^2)')
% axis([3,18,-0.2,0.4])
legend('\fontsize{16} Component I', '\fontsize{16} Component II', ...
       'Location', 'NorthWest')

%  plot factor scores

figure(2)
subplot(1,2,1)
phdl=plot(Fmat(:,2), -Fmat(:,1), 'o', ...
          [-4,4], [0,0], 'b:', ...
          [0,0], [-1,6], 'b:');
set(phdl, 'LineWidth', 2)
axis([-4, 4, -1, 6])
xlabel('\fontsize{16} Component I Scores')
ylabel('\fontsize{16} Component II Scores)')
title('\fontsize{16} \lambda_2 = 0')

subplot(1,2,2)
phdl=plot(Fmat2(:,1), -Fmat2(:,2), 'o', ...
          [-10,6], [0,0], 'b:', ...
          [0,0], [-1,7], 'b:');
set(phdl, 'LineWidth', 2)
axis([-10, 6, -1, 7])
xlabel('\fontsize{16} Component I Scores')
ylabel('\fontsize{16} Component II Scores)')
title('\fontsize{16} \lambda_2 = 100')

 %  compute angles to columns of Fmat
 
Fangl = zeros(N,1);
for i=1:N
    Fangl(i) = atan(Fmat(i,2)/Fmat(i,1));
end
Fangl = Fangl + 3*pi/2;

 %  compute radial distances to columns of Fmat
 
Fdist = zeros(N,1);
for i=1:N
    Fdist(i) = sqrt(sum(Fmat(i,:).^2));
end

%  smooth radial distances over angles

rngangl   = [min(Fangl),max(Fangl)];
anglbasis = create_bspline_basis(rngangl,56);
lambda    = 1e0;
anglfdPar = fdPar(anglbasis, 2, lambda);
anglfd    = smooth_basis(Fangl, Fdist, anglfdPar);
Fanglfine = linspace(rngangl(1),rngangl(2),101)';
figure(3)
plotfit_fd(Fdist, Fangl, anglfd)
xlabel('\fontsize{16} Angle')
ylabel('\fontsize{16} Distance')

%  plot radial distances against angles along with fit

Fdisthat = eval_fd(Fangl, anglfd);

figure(3)
phdl = plot(Fangl, Fdist, 'b.', Fangl, Fdisthat, 'bo');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Angle')
ylabel('\fontsize{16} Distance')
axis([2, 5, 1.5, 5.5])

%  set up matrix of fitted values for Fmat

Fhat = Pmat*Fmat;

Fhat = zeros(N,2);
for i=1:N
    Fhat(i,:) = Fdisthat(i).*Fmat(i,:)./Fdist(i);
end

%  compute projection operator on to the complement of cols of Fhat

Pmat = eye(N) - Fhat*inv(Fhat'*Fhat)*Fhat';

%  plot factors scores plus fitted values

figure(2)
subplot(1,1,1)
phdl=plot(Fmat(:,1), Fmat(:,2), 'b.', Fhat(:,1), Fhat(:,2), 'bo');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Component I Scores')
ylabel('\fontsize{16} Component II Scores)')

%  compute RSQR for the fit

SSEF0 = sum(sum(Fmat.^2));
SSEF1 = sum(sum((Pmat*Fmat).^2));
RSQR  = (SSEF0-SSEF1)/SSEF0;
disp(['RSQR = ',num2str(RSQR)])

% RSQR = 0.97869

%  3D case
figure(2)
subplot(1,1,1)
plot3(Fmat(:,1),Fmat(:,2),Fmat(:,3),'o')
grid on

% canonical correlations between pca_fd and Hfn_FDA Fmats

[U1, S1, V1] = svd(Fmat,0);
[U2, S2, V2] = svd(pcaFmat,0);
CCR = svd(U1'*U2);
disp(['Canonical correlations:   ', num2str(CCR')])

% K = 3:  Canonical correlations:   1     0.98938    0.076982
% K = 2:  Canonical correlations:   0.99978     0.25219
% K = 1:  Canonical correlations:   0.19302
% canonical correlations between pca_fd coefmat and Hfn_FDA cmatest

[U1, S1, V1] = svd(coefmat,0);
[U2, S2, V2] = svd(cmatest,0);
CCR = svd(U1'*U2);
disp(['Canonical correlations:   ', num2str(CCR')])

% K = 3:  Canonical correlations:   1.00000     0.99962     0.36231
% K = 2:  Canonical correlations:   0.99999     0.64604
% K = 1:  Canonical correlations:   0.32922

%  compute approximation to data

for i=1:N
    subplot(2,1,1)
    plot(tfine, Yfine(:,i), 'b--', tfine, Yhatmat(i,:)', 'b-', ...
         [3,18], [0,0], 'r:', [PGSctrmean, PGSctrmean], [-4,3], 'r:')
    axis([3,18,-4,3])
    title(['Curve ',num2str(i)])
    subplot(2,1,2)
    plot(tfine, Rfine(:,i), 'r-', [3,18], [0,0], 'b--')
    pause
end

%  display proportion of variance accounted for by both methods

SSE = sum(sum(Rfine.^2));
SSE_FDA_2 = SSE;

Rsq_FDA = (SSY - SSE)/SSY;

disp(['FDA    R^2 = ', num2str(Rsq_FDA)])
disp(['pca_fd R^2 = ', num2str(sum(pcastruct.varprop))])

%  reanalysis with hyperspherical factor score model

PenStruct.FColMat = Pmat;
PenStruct.FColPar = 1e0;  %  smoothing parameter lambda for Smat

parvec0 = parvec;

tic;
parvec_pen = fminunc(@Hfn_FDA, parvec0, options, Yfd, AfdPar, ...
                 [], gamval, PenStruct, 1, 1);
toc
[Hval_pen, Hgrad_pen, Fmat_pen] = ...
    Hfn_FDA(parvec_pen, Yfd, AfdPar, [], gamval, PenStruct, 1, 1);

figure(3)
subplot(1,1,1)
plot(Fmat_pen(:,1), Fmat_pen(:,2), 'o')


% cvec = sphere2cartes(parvec,1);
cvec = parvec;

cvec0 = cvec;

Hfn_FDA_cvec_1 = cvec;
save Hfn_FDA_cvec_1 Hfn_FDA_cvec_1 
Hfn_FDA_cvec_2 = cvec;
save Hfn_FDA_cvec_2 Hfn_FDA_cvec_2 
Hfn_FDA_cvec_3 = cvec;
save Hfn_FDA_cvec_3 Hfn_FDA_cvec_3 
Hfn_FDA_cvec_4 = cvec;
save Hfn_FDA_cvec_4 Hfn_FDA_cvec_4 

%  K = 1:   71 iterations, Elapsed time is 11.221869 seconds.
%  K = 2:   54 iterations, Elapsed time is 9.181508 seconds.
%  K = 3:  112 iterations, Elapsed time is 18.945006 seconds.
%  K = 4:  135 iterations, Elapsed time is 21.248448 seconds.

[Hval, Hgrad, Fmat] = Hfn_FDA(parvec, Yfd, AfdPar, [], gamval, PenStruct);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

% K = 4:

% Hfn_FDA R^2 = 0.96462
% pca_fd  R^2 = 0.9646

% K = 3:  
% Hfn_FDA R^2 = 0.9263
% pca_fd  R^2 = 0.9263

% K = 2:  
% Hfn_FDA R^2 = 0.87027
% pca_fd  R^2 = 0.8703

%  K = 1:  
% Hfn_FDA R^2 = 0.67159
% pca_fd  R^2 = 0.6716

%  comparative plots of three solutions

load Hfn_FDA_cvec_1  
load Hfn_FDA_cvec_2  
load Hfn_FDA_cvec_3  
load Hfn_FDA_cvec_4  

K = 1;
AfdPar = fdPar(fd(zeros(nAbasis,K),Abasis),2, Hlambda, 1, Rmat);
parvec1 = cartes2sphere(Hfn_FDA_cvec_1);
[Hval_1, Hgrad_1, Fmat_1] = Hfn_FDA(parvec1, Yfd, AfdPar, ...
                                    0, [], PenStruct);
K = 2;
AfdPar = fdPar(fd(zeros(nAbasis,K),Abasis),2, Hlambda, 1, Rmat);
parvec2 = cartes2sphere(Hfn_FDA_cvec_2);
[Hval_2, Hgrad_2, Fmat_2] = Hfn_FDA(parvec2, Yfd, AfdPar, ...
                                    0, [], PenStruct);
K = 3;
AfdPar = fdPar(fd(zeros(nAbasis,K),Abasis),2, Hlambda, 1, Rmat);
parvec3 = cartes2sphere(Hfn_FDA_cvec_3);
[Hval_3, Hgrad_3, Fmat_3] = Hfn_FDA(parvec3, Yfd, AfdPar, ...
                                    0, [], PenStruct);
K = 4;
AfdPar = fdPar(fd(zeros(nAbasis,K),Abasis),2, Hlambda, 1, Rmat);
parvec4 = cartes2sphere(Hfn_FDA_cvec_4);
[Hval_4, Hgrad_4, Fmat_4] = Hfn_FDA(parvec4, Yfd, AfdPar, ...
                                    0, [], PenStruct);

K = 1;
Afdest_1  = fd(reshape(Hfn_FDA_cvec_1, nAbasis, K),Abasis);
Amat_1    = eval_fd(tfine, Afdest_1);
Yhatmat_1 = Fmat_1*Amat_1';
Rfine_1   = Yfine - Yhatmat_1';
RMSE_1fine = sqrt(mean(Rfine_1.^2,2));
K = 2;
Afdest_2  = fd(reshape(Hfn_FDA_cvec_2, nAbasis, K),Abasis);
Amat_2    = eval_fd(tfine, Afdest_2);
Yhatmat_2 = Fmat_2*Amat_2';
Rfine_2   = Yfine - Yhatmat_2';
RMSE_2fine = sqrt(mean(Rfine_2.^2,2));
K = 3;
Afdest_3  = fd(reshape(Hfn_FDA_cvec_3, nAbasis, K),Abasis);
Amat_3    = eval_fd(tfine, Afdest_3);
Yhatmat_3 = Fmat_3*Amat_3';
Rfine_3   = Yfine - Yhatmat_3';
RMSE_3fine = sqrt(mean(Rfine_3.^2,2));
K = 4;
Afdest_4  = fd(reshape(Hfn_FDA_cvec_4, nAbasis, K),Abasis);
Amat_4    = eval_fd(tfine, Afdest_4);
Yhatmat_4 = Fmat_4*Amat_4';
Rfine_4   = Yfine - Yhatmat_4';
RMSE_4fine = sqrt(mean(Rfine_4.^2,2));

load PGSctr
PGSctrmean = mean(PGSctr);

%  plot RMSE for k = 1,...,4

figure(5)
subplot(1,1,1)
plot(tfine, RMSE_1fine, 'b--', tfine, RMSE_2fine, 'b-',  ...
     tfine, RMSE_3fine, 'b:',  tfine, RMSE_4fine, 'b-.', ... 
     [PGSctrmean,PGSctrmean],[0,1.8], 'b:')
axis([3,18,0,1.8])
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} RMSE(t) (cm/year^2)') 
legend('\fontsize{13} K = 1', '\fontsize{13} K = 2', ...
       '\fontsize{13} K = 3', '\fontsize{13} K = 4')

%  plot RMSE for k = 2,...,4

figure(5)
subplot(1,1,1)
phdl = plot(tfine, RMSE_2fine, 'b-.',  ...
            tfine, RMSE_3fine, 'b--',  tfine, RMSE_4fine, 'b-', ... 
            [PGSctrmean,PGSctrmean],[0,0.7], 'b:');
set(phdl, 'LineWidth', 2)
axis([3,18,0,0.7])
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} RMSE(t) (cm/year^2)') 
legend('\fontsize{13} K = 2', ...
       '\fontsize{13} K = 3', '\fontsize{13} K = 4')

%%  -----------------------------------------------------------------
%                PLS fit to age(18) using function Hfn_FDA
%  -----------------------------------------------------------------

%  data to be fit  ... height at age 18

Ymat = hgtfmat(31,:)';

%  define dimensions

agerng = [3,17];

%  define Abasis

nAbasis = 13;
Abasis  = create_bspline_basis(agerng, nAbasis);

%  dimension of subspace

K = 1;

%  define Ybasis

nxbasis = 53;
Xbasis = create_bspline_basis(agerng, nybasis);
lambda = 1e-8;
XfdPar = fdPar(Xbasis, 2, lambda);

%  set up the data

index = 1:ncasef;
N = length(index);

nfine = 1001;
tfine = linspace(3,17,nfine)';

Xfine = zeros(nfine,N);
m = 0;
for i = index
    m = m + 1;
    beta = betaf(:,i);
    Xfine(:,m) = beta(2).*eval_mon(tfine, Wfdf(i), 2);
end

Xfd = smooth_basis(tfine, Xfine, XfdPar);

%  set up initial coefficients

pcastruct = pca_fd(Xfd, K, Abasis, 0);
harmmat   = eval_fd(tfine,pcastruct.harmfd);
PCfd0     = smooth_basis(tfine, harmmat, Abasis);
coefmat   = getcoef(PCfd0);
cvecPCA   = reshape(coefmat,nAbasis*K,1);

%  set up regulatization and AfdPar

Hlambda = 0;
Rmat    = [];
AfdPar = fdPar(PCfd0, 2, Hlambda, 1, Rmat);

%  define PenStruct

PenStruct.FRowMat = [];
PenStruct.FRowPar = 0;
PenStruct.FColMat = [];
PenStruct.FColPar = 0;

%  set options for optimization

options = optimset('LargeScale',  'off',  ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on', ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  carry out the pca

cvec0 = cvecPCA;

parvec0 = cartes2sphere(cvec0);

Hlambda = 1e-4;
Rmat    = [];
AfdPar  = fdPar(PCfd0, 2, Hlambda, 1, Rmat);

gamval  = 0;
gamval  = 1;

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, Xfd, AfdPar, ...
                 gamval, Ymat, PenStruct);
toc

cvec = sphere2cartes(parvec,1);

[Hval, Hgrad, Fmat] = Hfn_FDA(parvec, Xfd, AfdPar, ...
                              gamval, Ymat, PenStruct);
pcaFmat = Fmat;

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

Afdest  = fd(reshape(cvec, nAbasis, K),Abasis);
Amat    = eval_fd(tfine, Afdest);
Yhatmat = Fmat*Amat';
Rfine   = Yfine - Yhatmat';

figure(1)
subplot(1,1,1)
plot(Afdest)
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} a_k')

figure(2)
plot(Fmat(:,1),Fmat(:,2), 'o')
xlabel('\fontsize{13} PC I score')
ylabel('\fontsize{13} PC II score')

% canonical correlations between pca_fd and Hfn_FDA Fmats

[U1, S1, V1] = svd(Fmat,0);
[U2, S2, V2] = svd(pcaFmat,0);
CCR = svd(U1'*U2);
disp(['Canonical correlations:   ', num2str(CCR')])

%  compute approximation to Y

Bmat = Fmat\Ymat;
disp(Bmat')
Yhat2 = Fmat*Bmat;

figure(3)
plot(Yhat2, Ymat, 'o', Yhat2, Yhat2, 'b-')
hold on
for i=1:54
  text(Yhat2(i)+1, Ymat(i), num2str(i))
end
hold off
xlabel('\fontsize{13} Yhat')
ylabel('\fontsize{13} Y')

figure(4)
plot(Xfd(16))
hold on
plot(mean(Xfd))
hold off

disp(['Rsqrd = ', num2str(corr(Yhat2, Ymat)^2)])

%  -----------------------------------------------------------------
%           Registration to PCA using function Hfn_Reg
%  -----------------------------------------------------------------

agerng = [3,18];

%  define Ybasis

nybasis = 53;
Ybasis = create_bspline_basis(agerng, nybasis);
lambda = 1e-8;
YfdPar = fdPar(Ybasis, 2, lambda);

%  set up the data

Yfd = smooth_basis(tfine, Yfine, YfdPar);

%  define Wbasis

%  four order 3 basis functions, no roughness penalty

nWbasis = 4;
nWorder = 3;
breaks = [3,12,18];
Wbasis  = create_bspline_basis([3,18], nWbasis, nWorder, breaks);
Wlambda = 0e0;
WfdPar  = fdPar(Wbasis, 1, Wlambda);

%  seven order 3 basis functions, roughness penalty

nWbasis = 7;
Wbasis  = create_bspline_basis([3,18], nWbasis, nWorder);
Wlambda = 1e0;
Kmat    = eval_penalty(Wbasis, 1);
WfdPar  = fdPar(Wbasis, 1, Wlambda, 1, Kmat);

%  define Ybasis

%  see pca-fd analysies for set up of yfd

% %  landmark registration to provide starting values for Dmat%  
%   This was found to be a bad idea ... trapped the solution in local
%   optima for K = 3 and even K = 2.
%

load PGSctr
PGSctrmean = mean(PGSctr);

wbasisLM = create_bspline_basis([3,18], 4, 3, [3,PGSctrmean,18]);
WfdParLM = fdPar(wbasisLM,1,1e-1);

[YregfdLM, warpfdLM, WfdLM] = ...
       landmarkreg(Yfd, PGSctr, PGSctrmean, WfdParLM, 1);
YregmatLM = eval_fd(tfine, YregfdLM);

WmatLM = eval_fd(tfine, WfdLM);
WfdLM  = smooth_basis(tfine, WmatLM, Wbasis);

Dmat0LM = getcoef(WfdLM);

%  carry out a pca_fd of WfdLM to get a starting value for cvec

K = 2;

%  define Abasis

nAbasis = 23;
Abasis  = create_bspline_basis(agerng, nAbasis);
Hlambda = 0;
Rmat    = [];
AfdPar  = fdPar(Abasis, 2, Hlambda, 1, Rmat);

pcastructLM = pca_fd(YregfdLM, K, AfdPar, 0);

disp(pcastructLM.varprop)

disp(pcastructLM.values(1:6)./sum(pcastructLM.values))

harmfd = pcastructLM.harmfd;

cmatLM = getcoef(harmfd);


% 
% disp('Std dev of cols of Dmat0LM: ')
% disp(num2str(sqrt(var(Dmat0LM'))))
% 
% % Std dev of cols of Dmat0LM: 
% % 0.00051893   0.0040383    0.061792     0.18687     0.34598     0.41423     0.41564

figure(4)
index = 1:N;
index = 15;
for i=index
    subplot(2,1,1)
    plotdeform_fd(WfdLM(i))
    axis([3,18,-3,3])
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} W(t)')
    title(['\fontsize{13} Case ',num2str(i)])
    subplot(2,1,2)
    plot(tfine, YregmatLM(:,i), 'b-', tfine, Yfine(:,i), 'b--', ...
         [agerng(1),agerng(2)], [0,0], 'r:', ...
         [PGSctrmean, PGSctrmean],[-6,3],'r:')
    axis([3,18,-6,3])
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} y(t)')    
    if length(index) > 1
        pause
    end
end

%  ------------------  smoothing parameters  -----------------------

Flambda = 0;
Smat    = [];

Hlambda = 0;
Rmat    = [];

AfdPar = fdPar(Abasis, 2, Hlambda, 1, Rmat);

Wlambda = 1e0;
% WfdPar  = fdPar(Wbasis, 1, Wlambda, 1, Kmat);
WfdPar  = fdPar(Wbasis, 1, Wlambda);

%  set options for optimization

Joptions = optimset('LargeScale',  'off',  ...
                    'Display',     'off',  ...
                    'Diagnostics', 'off',  ...
                    'GradObj',     'on',   ...
                    'Hessian',     'off',  ...  
                    'MaxIter',     50,     ...
                    'TolFun',      1e-6);
               
Hoptions = optimset('LargeScale',   'off',    ...
                    'Display',      'iter',   ...
                    'Diagnostics',  'on',     ...
                    'GradObj',      'on',     ...
                    'Hessian',      'off',    ...
                    'MaxIter',      100,      ...
                    'TolFun',       1e-5);

%  Initialize global matrix Dmat

global Dmat

%  -------------   fit functional PCA model  ----------------------

matlabpool open

%  Initialize Dmat, registration functions play role opposite to 
%  what they do in landmark registration

WfdPar = putlambda(WfdPar,1e-1);

Dmat = -Dmat0LM;  %  Used only for K = 1

cvec0 = cvec;
cvec0 = Hfn_FDA_cvec_3;
cvec0 = Hfn_FDA_cvec_2;
cvec0 = Hfn_FDA_cvec_1;
cvec0 = reshape(cmatLM,nAbasis*K,1);
AfdPar = fdPar(fd(reshape(cvec0,nAbasis,K),Abasis), 2, Hlambda, 1, Rmat);

i = 1;
dvec = Dmat(:,i);
yfine = Yfine(:,i);

[Fval, Fgrad, Fhess, Fcross, fvec, yhat, Hi, Dc_Hi, Wmat] = ...
                   PCA_Reg_trapz(dvec, tfine, yfine, AfdPar, WfdPar, ...
                           PenStruct, periodic, crit);
                       
Dmat = zeros(nWbasis,N);
WfdPar = putlambda(WfdPar,1e-1);
tic;
[Hval, Hgrad, Fmat, Yhat] = Hfn_Reg(cvec0, tfine, Yfine, Joptions, ...
                                    AfdPar, WfdPar);
toc
disp(['Hval for Hfn_Reg = ',num2str(Hval)])
Wfd  = fd(Dmat, Wbasis);
figure(1)
plotdeform_fd(Wfd)
axis([3,18,-3,3])
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} Deformation d(t) = h(t) - t')

Dmat = zeros(nWbasis,N);
WfdPar = putlambda(WfdPar,1e1);
tic;
[Hval, Hgrad, Fmat, Yhat] = Hfn_Reg_trapz(cvec0, tfine, Yfine, ...
        AfdPar, WfdPar, PenStruct, Joptions);
toc
disp(['Hval Hfn_Reg_trapz = ',num2str(Hval)])
Wfd  = fd(Dmat, Wbasis);
figure(2)
plotdeform_fd(Wfd)
axis([3,18,-3,3])
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} Deformation d(t) = h(t) - t')

Dmat = zeros(nWbasis,N);
WfdPar = putlambda(WfdPar,1e-1);
tic;
cvec = fminunc(@Hfn_Reg, cvec0, Hoptions, ...
               tfine, Yfine, Joptions, AfdPar, WfdPar);
toc

Dmat = zeros(nWbasis,N);
WfdPar = putlambda(WfdPar,1e1);
tic;
cvec = fminunc(@Hfn_Reg_trapz, cvec0, Hoptions, ...
               tfine, Yfine, AfdPar, WfdPar, PenStruct, Joptions);
toc

% K = 3:  Elapsed time is 391.138224 seconds,  77 iterations
% K = 2:  Elapsed time is 360.625407 seconds,  78 iterations
% K = 1:  Elapsed time is 173.648638 seconds,  32 iterations

Hfn_Reg_cvec_1 = cvec;
save Hfn_Reg_cvec_1 Hfn_Reg_cvec_1 
Hfn_Reg_cvec_2 = cvec;
save Hfn_Reg_cvec_2 Hfn_Reg_cvec_2 
Hfn_Reg_cvec_3 = cvec;
save Hfn_Reg_cvec_3 Hfn_Reg_cvec_3 

[Hval, Hgrad, Fmat, Yhat] = Hfn_Reg(cvec, tfine, Yfine, Joptions, ...
                                    AfdPar, WfdPar);

Hfn_Reg_Dmat_1 = Dmat;
save Hfn_Reg_Dmat_1 Hfn_Reg_Dmat_1 
Hfn_Reg_Dmat_2 = Dmat;
save Hfn_Reg_Dmat_2 Hfn_Reg_Dmat_2 
Hfn_Reg_Dmat_3 = Dmat;
save Hfn_Reg_Dmat_3 Hfn_Reg_Dmat_3 
      
%  for K = 2, Hval = 19.87 on the final 49th iteration represents a 
%  failure of the Quasi-Newton optimization.  Restarting with the current
%  value of cvec achieved the optimal value

cvec0 = cvec;

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

% K = 3:  Hval = 12.6742
% Hgrad RMSE = 0.00080082
% K = 2:  Hval = 18.9191
% Hgrad RMSE = 0.00099839
%  K = 1: Hval = Hval = 46.0143
% Hgrad RMSE = 0.00045058

Acoefmatest = reshape(Hfn_Reg_cvec_2, nAbasis, K);
Afdest  = fd(Acoefmatest,Abasis);
Amatest = eval_fd(tfine, Afdest);

figure(1)
subplot(1,1,1)
phdl = plot(tfine, Amatest, 'b-', [3,18], [0,0], 'r:', ...
            [PGSctrmean, PGSctrmean],[-2,2],'r:');
set(phdl, 'LineWidth', 2)
axis([3,18,-2,2])
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} A(t)')    

%  plot factor scores

figure(2)
subplot(1,1,1)
plot(Fmat(1,:),Fmat(2,:),'o')

%  plot warping functions

Wfd  = fd(Dmat, Wbasis);

figure(3)
plotdeform_fd(Wfd)
axis([3,18,-3,3])
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} Deformation d(t) = h(t) - t')

%  plot fit to data ... Yhat extracted from fitstruct_final

MSE = zeros(N,1);
SSE = 0;
figure(4)
index = 1:N;
for i=index
    subplot(3,1,1)
    plot(Wfd(i))
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} W(t)')
    title(['\fontsize{13} Row ', num2str(i), ...
           ',  fvec: ',num2str(Fmat(:,i)')])
    subplot(3,1,2)
    plotdeform_fd(Wfd(i))
%     axis([3,18,-1.5,1.5])
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} d(t) = h(t) - t')
    MSE(i) = mean((Yfine(:,i) - Yhat(:,i)).^2);
    SSE = SSE + sum((Yfine(:,i) - Yhat(:,i)).^2);
    subplot(3,1,3)
    plot(tfine, Yfine(:,i), 'b-', tfine, Yhat(:,i), 'b--', ...
         [agerng(1),agerng(2)], [0,0], 'r:', ...
         [PGSctrmean, PGSctrmean],[-6,3],'r:')
    axis([3,18,-6,3])
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} y(t)') 
    title(['\fontsize{13} RMSE = ',num2str(sqrt(MSE(i)))])
    if length(index) > 1
        pause
    end
end

%  plot MSE(t)

Rfine = Yfine - Yhat;
Rfine_Reg_2 = Rfine;

RMSEfine = sqrt(mean(Rfine.^2,2));
figure(4)
subplot(1,1,1)
plot(tfine, RMSEfine, 'b-')
axis([3,18,0,max(RMSEfine)])
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} RMSE(t)') 
title(['\fontsize{16} K = ',num2str(K)])

SSE = sum(sum(Rfine.^2));

RMSE = sqrt(mean(MSE));
disp(['RMSE = ',num2str(RMSE)])

figure(5)
hist(MSE)

HfnRegRsq = (SSY - SSE)/SSY;

disp(['Hfn_FDA R^2 = ', num2str(HfnFDARsq)])
disp(['pca_fd  R^2 = ', num2str(sum(pcastruct.varprop))])
disp(['Hfn_Reg R^2 = ', num2str(HfnRegRsq)])

% K = 3, Wlambda = 1:
% Hfn_FDA R^2 = 0.8959
% pca_fd  R^2 = 0.8959
% Hfn_Reg R^2 = 0.97912  (5 unregistered factors)
% K = 2, Wlambda = 1:
% Hfn_FDA R^2 = 0.87027
% pca_fd  R^2 = 0.8703
% Hfn_Reg R^2 = 0.96883  (4 unregistered factors)
% K = 1, Wlambda = 1:
% Hfn_FDA R^2 = 0.9263
% pca_fd  R^2 = 0.9263
% Hfn_Reg R^2 = 0.92418  (3 unregistered factors)

% 13 July 12  What has been learned:
% Starting Dmat from the LM registration values is a bad idea ...
%     possibly because the factor scores do a lot of the registration, 
%     especially for K = 3.  2-4 complete failures found, but not when
%     starting from 0.
% Setting Wlambda = 1 gives better results than Wlambda = -1 with 7
% basis functions.  Need this level of smoothing.  Try with fewer 
% basis functions.
% 16 July 12  What has been learned:
% Much better results by starting Hfn_Reg from Hfn_FDA cvec.
% K = 2 produced an optimization failure, but restarting with the
% terminal cvec then converged properly
% The discrepancy between the pca_fd results and Hfn_FDA remains a
% puzzle.
% Plot Rfine = Yfine - Yhat

%  comparative plots of three solutions

load Hfn_Reg_cvec_1  
load Hfn_Reg_cvec_2  
load Hfn_Reg_cvec_3  
load Hfn_Reg_Dmat_1  
load Hfn_Reg_Dmat_2 
load Hfn_Reg_Dmat_3  

K = 1;
AfdPar = fdPar(fd(zeros(nAbasis,K),Abasis),2, Hlambda, 1, Rmat);
Dmat = Hfn_Reg_Dmat_1;
[Hval_1, Hgrad_1, Fmat_1, Yhat_1] = ...
    Hfn_Reg(Hfn_Reg_cvec_1, tfine, Yfine, Joptions, AfdPar, WfdPar);
K = 2;
AfdPar = fdPar(fd(zeros(nAbasis,K),Abasis),2, Hlambda, 1, Rmat);
Dmat = Hfn_Reg_Dmat_2;
[Hval_2, Hgrad_2, Fmat_2, Yhat_2] = ...
    Hfn_Reg(Hfn_Reg_cvec_2, tfine, Yfine, Joptions, AfdPar, WfdPar);
K = 3;
AfdPar = fdPar(fd(zeros(nAbasis,K),Abasis),2, Hlambda, 1, Rmat);
Dmat = Hfn_Reg_Dmat_3;
[Hval_3, Hgrad_3, Fmat_3, Yhat_3] = ...
    Hfn_Reg(Hfn_Reg_cvec_3, tfine, Yfine, Joptions, AfdPar, WfdPar);

Rfine_1 = Yfine - Yhat_1;
RMSE_1fine = sqrt(mean(Rfine_1.^2,2));
Rfine_2 = Yfine - Yhat_2;
RMSE_2fine = sqrt(mean(Rfine_2.^2,2));
Rfine_3 = Yfine - Yhat_3;
RMSE_3fine = sqrt(mean(Rfine_3.^2,2));

figure(5)
subplot(1,1,1)
phdl = plot(tfine, RMSE_1fine, 'b--', tfine, RMSE_2fine, 'b-', ...
            tfine, RMSE_3fine, 'b:', ...
            [PGSctrmean,PGSctrmean],[0,0.7], 'b:');
set(phdl, 'LineWidth', 2)
axis([3,18,0,0.7])
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} RMSE(t) (cm/year^2)') 
legend('\fontsize{13} K = 1', '\fontsize{13} K = 2', '\fontsize{13} K = 3')

%  -----------------------------------------------------------------
%           Registration to PCA using function pcaregfn
%  -----------------------------------------------------------------

%  basis for warping function:

nwbasis   =  7;
wbasis    = create_bspline_basis(agerng, nwbasis);
warpfdPar = fdPar(wbasis);

%  set up the PCA

harmfdPar = fdPar(ybasis, 2, 0);
centerfns = 0;

%  set up cell object for regression coefficient function

betaLfdobj   = 2;
betalambda   = 1e3;
betaestimate = 1;
betapenmat   = eval_penalty(dbasis,betaLfdobj);
betafdPar    = fdPar(dbasis, betaLfdobj, betalambda, betaestimate, ...
                     betapenmat);

%  ----------------------  do the analysis  --------------------

yfd0 = accffd;

plotwrd =  0;
iterlim =  7;
conv    = 1e-6;
nharm   =  2;

tic
[yfd, dcoefs, ycoefregs, tcoefregs, RMSE, RSQR, D2int, ...
        MSamp, MSpha, MSres, C] = ...
    pcaregfn(t, yfd0, nharm, betafdPar, warpfdPar, harmfdPar, ...
               iterlim, centerfns, plotwrd, conv);
toc

%  set up warping functions

tcoef = squeeze(tcoefregs(:,:,nharm,iterlim));
wfd = fd(tcoef, wbasis);

%  plot registered and unregistered curves

figure(1)

subplot(2,1,1)
ymat0 = eval_fd(t, yfd0);
lhdl = plot(t, ymat0, 'k-');
set(lhdl, 'LineWidth', 1)
ylabel('\fontsize{16} D^2 Height')
title('\fontsize{16} Unregistered')
axis([3,18,-4,3])

subplot(2,1,2)
ymat = eval_fd(t, yfd);
lhdl = plot(t, ymat, 'k-');
set(lhdl, 'LineWidth', 1)
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} D^2 Height')
title('\fontsize{16} Registered')
axis([3,18,-4,3])

print -dps2 'c:\Matlab7\Reg2PCA\figs\BerkAcc.ps'


%  compute principal components

[harmfd,  eigvals,  harmscr,  varprop] = ...
                   pca_fd(yfd, nharm, harmfdPar, 0);
               
disp(['Proportions of variance: ',num2str(varprop)])

%  plot harmonics

figure(2)
subplot(1,1,1)
harmmat = -eval_fd(t, harmfd);
lhdl = plot(t, harmmat(:,1), 'k-', t, harmmat(:,2), 'k--', ...
     [3,18], [0,0], 'k:');
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Components of D^2 Height')
axis([3,18,-0.8,0.6])
legend('\fontsize{16} 86.8%', '\fontsize{16}   6.0%')

print -dps2 'c:\Matlab7\Reg2PCA\figs\BerkPC.ps'

%  plot progress of iterations

figure(3)
subplot(3,1,1)
plot(0:iterlim, D2int, 'o-')
ylabel('\fontsize{12} Tot. Curv.')
subplot(3,1,2)
plot(0:iterlim, RSQR, 'o-')
ylabel('\fontsize{12} R^2')
subplot(3,1,3)
plot(0:iterlim, MSres, 'o-')
ylabel('\fontsize{12} MS_{res}')

% subplot(5,1,1)
% plot(0:iterlim, MSamp, 'o-')
% ylabel('\fontsize{12} MS_{amp}')
% subplot(5,1,3)
% plot(0:iterlim, MSpha, 'o-')
% ylabel('\fontsize{12} MS_{phase}')

wmat = eval_fd(t, wfd);
dmat = wmat - t*ones(1,ncasef);
Dwfd = deriv(wfd);
Wmat = log(eval_fd(t,Dwfd));

figure(3)
subplot(2,1,1)
lhdl = plot(t, dmat, 'k-');
set(lhdl, 'LineWidth', 1)
ylabel('\fontsize{19} Deformation h(t) - t')
axis([3,18,-2,2])

subplot(2,1,2)
lhdl = plot(t, Wmat, 'k-');
set(lhdl, 'LineWidth', 1)
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Log Dh(t)')
axis([3,18,-0.4,0.4])

print -dps2 'c:\Matlab7\Reg2PCA\figs\BerkWarp.ps'

Wfd = data2fd(Wmat, t, wbasis);

figure(5)
subplot(1,1,1)
plot(Wfd)

[Wharmfd,  Weigvals,  Wharmscr,  Wvarprop] = ...
                   pca_fd(Wfd, 2, harmfdPar, 0);
               
disp(['Proportions of variance: ', num2str(Wvarprop)])


save Berkeley

