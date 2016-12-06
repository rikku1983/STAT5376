% -------------------------------------------------------------------
%                 test function Hfn_Reg
% -------------------------------------------------------------------

%  Last modified 28 January 2013

addpath('../fdaM')

%  define dimensions

N = 50;
K = 2;
nfine  = 101;
tfine  = linspace(0,1,nfine)';

%  define Abasis

nAbasis = 7;
Abasis = create_bspline_basis([0,1], nAbasis);

%   generate coefficient matrix for coefficient fns.

Acoef0 = randn(nAbasis,K);
Afd0   = fd(Acoef0, Abasis);
AfdPar = fdPar(Afd0);

%  plot coefficient functions

figure(1)
plot(Afd0)
xlabel('\fontsize{13} Time t')
ylabel('\fontsize{13} Coefficient \alpha(t)')

%  define Wbasis

% nWbasis = 3;
% Wbasis = create_monomial_basis([0,1],nWbasis);

nWbasis = 4;
nWorder = 3;
Wbasis = create_bspline_basis([0,1],nWbasis,nWorder);

%  define WfdPar

Wfd = fd(zeros(nWbasis,1),Wbasis);
Wlambda = 1e0;
Kmat    = eval_penalty(Wbasis,1);
WfdPar  = fdPar(Wfd, 1, Wlambda, 1, Kmat);

%   ---------------  generate coefficient scores  ------------------

%  Fmat patterned for K = 2

tvec  = (1:N)'*2*pi/N;
Fmat0 = [sin(tvec)'; cos(tvec)'];

figure(2)
plot(Fmat0(1,:),Fmat0(2,:),'o')
xlabel('\fontsize{13} True Component score I')
ylabel('\fontsize{13} True Component score II')

%   ---------------  generate warping functions  ------------------

%  set up warping function

sigwarp = 0.3;
Wwarp = randn(1,N)*sigwarp;

hfinemat = (1 - exp(tfine*Wwarp))./(1 - exp(ones(nfine,1)*Wwarp));

figure(3)
plot(tfine, hfinemat, '-')
xlabel('\fontsize{13} Time t')
ylabel('\fontsize{13} True warped time h(t)')

%  ----------------------  generate data  -------------------------

Ymat0 = zeros(nfine,N);
for i=1:N
    Ymat0(:,i) = eval_fd(hfinemat(:,i), Afd0)*Fmat0(:,i);
end

figure(4)
subplot(1,1,1)
plot(tfine, Ymat0', '-')
xlabel('\fontsize{13} Time t')
ylabel('\fontsize{13} Data x_i(t)')

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

figure(5)
subplot(1,1,1)
plot(tfine, Ymat', '-')
xlabel('\fontsize{13} Time t')
ylabel('\fontsize{13} Data x_i(t)')

nybasis = 53;
Ybasis = create_bspline_basis([0,1], nybasis);
Yfd = smooth_basis(tfine, Ymat, Ybasis);

%  ----------  generate initial values for parameters  ---------------

%  coefficients for Wfd's

Dmat0 = zeros(nWbasis,N);

%  true value

Cmat0 = Acoef0;
Afd0  = fd(Cmat0, Abasis);
Amat0 = eval_fd(tfine, Afd0);
cvec0 = reshape(Cmat0',nAbasis*K,1);

%  random value

cvec0 = randn(K*nAbasis,1);

%  from eigenanalysis

pcastruct = pca_fd(Yfd, K, AfdPar, 0);  %  something wrong ... fix this!!
PCfd0     = pcastruct.harmfd;
harmmat   = eval_fd(tfine,PCfd0);
Afd0      = smooth_basis(tfine, harmmat, AfdPar);
coefmat   = getcoef(Afd0);
cvec0     = reshape(coefmat,nAbasis*K,1);

Fmat0 = harmmat\Ymat;

Hlambda = 0;
Rmat    = [];
AfdPar = fdPar(Afd0, 2, Hlambda, 1, Rmat);

cvec = cvec0;  %  used for stepping through function

%  ------------------  smoothing parameters  -----------------------

Flambda = 0;
Smat    = [];

%  set options for optimization

Joptions = optimset('LargeScale',  'off',  ...
                    'Display',     'off',   ...
                    'Diagnostics', 'off', ...
                    'GradObj',     'on',     ...
                    'Hessian',     'off',    ...               
                    'TolFun',      1e-6);
               
%  ------------  check gradient and hessian of PCA_Reg  --------------

%  evaluate function, gradient, and hessian

dvec  = Dmat0(:,1);
yfine = Ymat(:,1);

tic;
[Fval, Fgrad, Fhess] = PCA_Reg(dvec, tfine, yfine, AfdPar, WfdPar);
toc

h0 = 1e0;
toler = 1e-6;

tic;
[Fgradhat, Fhesshat] = grdHessm_PCA_Reg(@PCA_Reg, dvec, h0, toler, ...
                                    tfine, yfine, AfdPar, WfdPar); 
toc

[Fgrad, Fgradhat]

Fhess

Fhesshat

%  ------------  check gradient and hessian of PCA_Reg_trapz  --------------

%  evaluate function, gradient, and hessian

dvec  = Dmat0(:,7);
yfine = Ymat(:,7);

tic;
[Fval, Fgrad, Fhess] = PCA_Reg_trapz(dvec, tfine, yfine, AfdPar, WfdPar);
toc

h0 = 1e0;
toler = 1e-6;

tic;
[Fgradhat, Fhesshat] = grdHessm_PCA_Reg(@PCA_Reg_trapz, dvec, h0, toler, ...
                                    tfine, yfine, AfdPar, WfdPar); 
toc

[Fgrad, Fgradhat]

Fhess

Fhesshat

%  compare Fcross for the two functions

tic;
[Fval, Fgrad, Fhess, Fcross] = ...
    PCA_Reg(dvec, tfine, yfine, AfdPar, WfdPar);
toc

Fcross

tic;
[Fval, Fgrad, Fhess, Fcross] = ...
    PCA_Reg_trapz(dvec, tfine, yfine, AfdPar, WfdPar);
toc

Fcross

%  ------------  check gradient of Hfn_reg  --------------------------

Hoptions = optimset('LargeScale',   'off',    ...
                    'Display',      'iter',   ...
                    'Diagnostics',  'on',     ...
                    'GradObj',      'on',     ...
                    'Hessian',      'off',    ...
                    'MaxIter',      100,      ...
                    'TolFun',       1e-5);

%  evaluate function and gradient

global Dmat

Dmat0 = zeros(nWbasis,N);

Dmat = Dmat0;

matlabpool open

tic;
[Hval, Hgrad] = Hfn_Reg(cvec0, tfine, Ymat, Joptions, AfdPar, WfdPar);
toc
% Elapsed time is 2.746344 seconds.
Hval
%     0.3258

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;

tic;
Hgradhat = grdHessm_Hfn_Reg(@Hfn_Reg, cvec0, h0, toler, ...
                            tfine, Ymat, Joptions, AfdPar, WfdPar); 
toc/60
% Elapsed time is 5.9 minutes.

[Hgrad, Hgradhat]

%    -0.0106   -0.0106
%    -0.0185   -0.0189
%    -0.0512   -0.0512
%    -0.0626   -0.0600
%    -0.0195   -0.0171
%    -0.0232   -0.0214
%    -0.0043   -0.0038
%    -0.0025   -0.0030
%    -0.0055   -0.0059
%    -0.0144   -0.0160
%    -0.0174   -0.0192
%    -0.0052   -0.0045
%    -0.0079   -0.0079
%    -0.0038   -0.0038

%  ------------  check gradient of Hfn_reg_trapz  --------------------------

PenStruct.FRowMat = [];
PenStruct.FRowPar = 0;
PenStruct.FColMat = [];
% PenStruct.FColMat = eye(N) - Fmat0*inv(Fmat0'*Fmat0)*Fmat0';
PenStruct.FColPar = 0;
% PenStruct.FColPar = 1e-8;
PenStruct.ANrmPar = 0e0;

%  evaluate function and gradient

Dmat = Dmat0;

tic;
[Hval, Hgrad] = Hfn_Reg_trapz(cvec0, tfine, Ymat, AfdPar, WfdPar, ...
                              PenStruct, Joptions);
toc
% Elapsed time is 2.217231 seconds.
Hval
%     0.3254

%  approximate gradient by differencing

h0 = 1e-1;
toler = 1e-6;

tic;
Hgradhat = grdHessm_Hfn_Reg_trapz(@Hfn_Reg_trapz, cvec0, h0, toler, ...
                         tfine, Ymat, AfdPar, WfdPar, PenStruct, Joptions); 
toc/60

% Elapsed time is 8.4 minutes.

[Hgrad, Hgradhat]

%    -0.0106   -0.0106
%    -0.0183   -0.0189
%    -0.0511   -0.0507
%    -0.0627   -0.0600
%    -0.0196   -0.0171
%    -0.0232   -0.0216
%    -0.0044   -0.0038
%    -0.0025   -0.0030
%    -0.0055   -0.0059
%    -0.0144   -0.0164
%    -0.0175   -0.0189
%    -0.0053   -0.0046
%    -0.0079   -0.0074
%    -0.0038   -0.0035

matlabpool close

Hgradmat    = reshape(Hgrad,nAbasis,K)'
Hgradhatmat = reshape(Hgradhat,nAbasis,K)'

% Hgradmat =
% 
%    -0.0106   -0.0183   -0.0511   -0.0627   -0.0196   -0.0232   -0.0044
%    -0.0025   -0.0055   -0.0144   -0.0175   -0.0053   -0.0079   -0.0038
% 
% 
% Hgradhatmat =
% 
%    -0.0106   -0.0189   -0.0507   -0.0600   -0.0171   -0.0216   -0.0038
%    -0.0030   -0.0059   -0.0164   -0.0189   -0.0046   -0.0074   -0.0035

%  -------------   fit functional PCA model  ----------------------

%  set smoothing for factor scores

Flambda = 0;
Smat    = [];

%  -------------  fit data using Hfn_Reg  ------------------------

global Dmat

Dmat0 = zeros(nWbasis,N);

Dmat = Dmat0;

parvec0 = cartes2sphere(cvec0);

tic;
parvec = fminunc(@Hfn_Reg, parvec0, Hoptions, ...
               tfine, Ymat, Joptions, AfdPar, WfdPar);
toc

% Elapsed time is 93.154438 seconds in 55 iterations
% Elapsed time is 1085 seconds in 87 iterations (on ThinkPad)

[Hval, Hgrad, Fmat, Yhat] = Hfn_Reg(parvec, tfine, Ymat, Joptions, ...
                                    AfdPar, WfdPar);

disp(['Hval = ',num2str(Hval)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])

%  Dell results:
% Hval = 0.010716
% Hgrad RMSE = 3.7928e-006
%  ThinkPad results:
% Hval = 0.010716
% Hgrad RMSE = 1.2658e-005

% %  -------------  fit data using Hfn_Reg_trapz  ------------------------
% 
% Dmat = Dmat0;
% 
% tic;
% cvec = fminunc(@Hfn_Reg_trapz, cvec0, Hoptions, ...
%                tfine, Ymat, AfdPar, WfdPar, PenStruct, Joptions);
% toc
% 
% % Elapsed time is 111.178188 seconds in 55 iterations.
% 
% [Hval, Hgrad, Fmat, Yhat] = Hfn_Reg_trapz(cvec, tfine, Ymat, ...
%                                     AfdPar, WfdPar, PenStruct, Joptions);
% 
% disp(['Hval = ',num2str(Hval)])
% disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad.^2)))])
% 
% % Hval = 0.010706
% % Hgrad RMSE = 3.6648e-006

%  ------------------  Set up Afdest and Wfdest  -----------------------

cvec = sphere2cartes(parvec, 1);
Acoefmatest = reshape(cvec, K, nAbasis);
Afdest  = fd(Acoefmatest',Abasis);
Amatest = eval_fd(tfine, Afdest);

Wfdest = fd(Dmat, Wbasis);

%  ------------------  display results  ---------------------------------

%  load previous results if desired

load Hfn_Reg_Test

%  plot coefficient functions

figure(6)
subplot(1,1,1)
plot(tfine, Amatest, 'b-', tfine, Amat0, 'r--')
xlabel('\fontsize{13} Time t')
ylabel('\fontsize{13} Coefficients \alpha_k(t)')
legend('\fontsize{12} Estimated I', '\fontsize{12} Estimated II', ...
       '\fontsize{12} True I',      '\fontsize{12} True II')

%  plot factor scores

figure(7)
subplot(1,1,1)
plot(Fmat(1,:),Fmat(2,:),'bo',Fmat0(1,:),Fmat0(2,:),'ro')
xlabel('\fontsize{13} Component score I')
ylabel('\fontsize{13} Component score II')
legend('\fontsize{12} Estimated', '\fontsize{12} True')

%  plot deformation functions

figure(8)
plotdeform_fd(Wfdest)
xlabel('\fontsize{13} Time t')
ylabel('\fontsize{13} Deformation d_i(t) = h_i(t) - t')

figure(9)
for i=1:N
    subplot(2,1,1)
    plotdeform_fd(Wfdest(i))
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} d(t) = h(t) - t')
    title(['\fontsize{13} Row ',num2str(i)])
    subplot(2,1,2)
    plot(tfine, Yhat(:,i), 'b-', tfine, Ymat(:,i), 'b--', ...
         [0,1], [0,0], 'r:')
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} y(t)')    
    pause
end


save Hfn_Reg_Test 


