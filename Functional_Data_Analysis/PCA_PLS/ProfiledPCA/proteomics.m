addpath('../fdaM')

%  Last modified 16 June 2014

%  The proteomics data are mass spectrometry data with the independent
%  variable being time in minutes and the output total ion count (TIC)
%  Minutes varies from 20 to 220 in the raw data.  
%  There are 15 spectra.

%  These analyses are designed to show HfnReg in action as a continuous 
%  registration tool.  It is compared with landmark registration.

%%  ----------------------------------------------------------------------------
%              Set up and display the data
%  ----------------------------------------------------------------------------

%  load the mat file for these analyses

load proteomics  % last saved 25 March 2014

load RawData.txt

n = 15;
traw = RawData(:,1);

xmat  = RawData(:,2:16);

%  plot the raw data (ranges from 6.8 to 9)

plot(traw, xmat)

%  log base 10 of the raw data (ranges from 0.83 to 0.95)

lxmat = log10(xmat);

%  plot the log raw data

plot(traw, lxmat)

%  plot a small segment of the data

indplt  = find(traw >= 100 & traw <= 119);
nindplt = length(indplt);

%  line plots

plot(traw(indplt), lxmat(indplt,:))

%  point plot

plot(traw(indplt), lxmat(indplt,1), '.')

%  Run minutes from 0 to 150 by subtracting 20 minutes and then using only
%  the first 150 minutes. The last part of the spectrum contains no useful
%  information and is eliminated in this way.

t = traw - 20;
t = t(t <= 150);
n = length(t);  %  15small
lxmat = lxmat(1:n,:);

%  set up basis for data
lxbasis = create_bspline_basis([0,150],n+2);
%  define smoothing
lxfdPar = fdPar(lxbasis,2,1e-4);
%  smooth the data
lxfd    = smooth_basis(t, lxmat, lxfdPar);
%  evaluate the smooth
plotfit_fd(lxmat, t, lxfd)
%  evaluate the first derivative
Dlxmat = eval_fd(t, lxfd, 1);

plot(t(indplt), Dlxmat(indplt,1), '-', ...
    [t(indplt(1)), t(indplt(nindplt))], [0,0], 'r:')

%  display all spectra in a single plot

figure(1)
for i=1:15
    subplot(3,5,i)
    plot(t, lxmat(:,i), '-')
    title(['Curve ', num2str(i)])
end

%%  ------------------------------------------------------------
%             Landmark registration using 3 landmarks
%  ------------------------------------------------------------

%  select three landmarks for each spectrum

peaktime = zeros(15,3);
figure(2)
subplot(1,1,1)
recind = 1:15;
recind = 4;
for j=1:3
    for i=recind
        plot(t, lxmat(:,i), '-')
        title(['Curve ', num2str(i)])
        j = 1;
        disp(['Peak ',num2str(j)])
        [xj,yj] = ginput(1);
        peaktime(i,j) = xj;
    end
end

peaktime = ...
  [32.5806   94.5161  135.3763 ...
   22.0430   83.9785  125.9140 ...
   16.0215   78.3871  121.1828 ...
   33.0108   94.9462  136.0215 ...
   30.4301   91.7204  133.2258 ...
   18.6022   80.5376  122.9032 ...
   33.2258   94.7312  135.5914 ...
   26.7742   88.4946  129.5699 ...
   14.0860   78.3871  121.6129 ...
   33.6559   94.7005  136.0215 ...
   28.9247   89.7849  131.5054 ...
   14.9462   77.7419  121.6129 ...
   33.2258   94.7312  135.5914 ...
   20.1075   82.2581  124.4086 ...
   14.0860   79.2473  123.5484];

peaktime = reshape(peaktime, 3, 15)';

peaktimemean = mean(peaktime);

%    24.7814   86.9442  128.9391

%  register using interpolation

[lxreg1fd, warp1fd] = landmarkreg(lxfd, peaktime, peaktimemean);

%  register with B-spline basis

wbasis = create_bspline_basis([0,150],9);
Wlambda = 1e-8;
WfdParobj = fdPar(wbasis, 2, Wlambda);

%  register with B-spline basis and strictly monotone warping

wbasis = create_bspline_basis([0,150],9);
Wlambda = 1e-8;
WfdParobj = fdPar(wbasis, 2, Wlambda);

[lxreg1fd, warp1fd] = landmarkreg(lxfd, peaktime, peaktimemean, ...
                                WfdParobj, 0);

%  display warping and deformation functions

figure(3)
subplot(2,1,1)
warp1mat = eval_fd(t, warp1fd);
plot(t, warp1mat, '-', [0,150], [0,150], 'r:')
axis([0,150,0,150])
title('Warping functions')
subplot(2,1,2)
plot(t, warp1mat-t*ones(1,15), '-', [0,150], [0,0], 'r:')
title('Deformation functions')

%  display unregistered and registered functions

lxreg1mat = eval_fd(t,lxreg1fd);

figure(4)
subplot(2,1,1)
plot(t,lxmat)
axis([0,150,0.83,0.95])
title('Un-registered functions')
subplot(2,1,2)
plot(t,lxreg1mat)
axis([0,150,0.83,0.95])
title('Registered functions')

%  regress lxreg1fd over [0,60] on non-responder/responder dummy variable

Yvec = zeros(15,1);
Yvec(10:15) = 1;

index = find(t >= 50 & t <= 100);
tshort = t(index)-50;
nshort = length(tshort);
lxbasisshort = create_bspline_basis([0,60],nshort+2);
%  define smoothing
lxfdParshort = fdPar(lxbasisshort,2,1e-4);
%  smooth the data
lxfdshort = smooth_basis(tshort, lxreg1mat(index,:), lxfdParshort);

beta1basis = create_constant_basis([0,60]);
beta1fd = fd(1,beta1basis);
nbetabasis = 103;
beta2basis = create_bspline_basis([0,60],nbetabasis);
beta2fd    = fd(zeros(nbetabasis,1), beta2basis);
betacell{1} = fdPar(beta1fd);
betacell{2} = fdPar(beta2fd, 2, 1e-5);

xfdcell{1} = ones(15,1);
xfdcell{2} = lxfdshort;

fRegressStruct = fRegress(Yvec, xfdcell, betacell);

betahat = fRegressStruct.betahat;
beta1fd = getfd(betahat{1});
beta2fd = getfd(betahat{2});

figure(1)
plot(beta1fd)
figure(2)
plot(beta2fd)

xreg1mat = 10.^lxreg1mat;
beta2mat = eval_fd(tshort, beta2fd);

figure(3)
subplot(2,1,1)
plot(tshort, beta2mat, 'b-', [0,50], [0,0], 'r:')
subplot(2,1,2)
plot(tshort, xreg1mat(index,:))

yhatfd = fRegressStruct.yhat;
figure(4)
plot(yhatfd)

%%  ------------------------------------------------------------
%    Continuous registration of landmark-registered curves
%  ------------------------------------------------------------

lxreg0fd = mean(lxreg1fd);

wbasis = create_bspline_basis([0,150],9);
Wlambda = 1e-8;
Wfd0Par = fdPar(wbasis, 2, Wlambda);

[lxreg2fd, warp2fd] = register_fd(lxreg0fd, lxreg1fd, Wfd0Par, 0, 2, 1e-6);

% Note:  register_fd does much better here than register_fd_QD. which
% does not decrease the function in the first iteration for most curves.

%  display warping and deformation functions

figure(3)
subplot(2,1,1)
warp2mat = eval_fd(t, warp2fd);
plot(t, warp2mat, '-', [0,150], [0,150], 'r:')
axis([0,150,0,150])
title('Warping functions')
subplot(2,1,2)
plot(t, warp2mat-t*ones(1,15), '-', [0,150], [0,0], 'r:')
title('Deformation functions')

%  display unregistered and registered functions

lxreg2mat = eval_fd(t,lxreg2fd);

figure(4)
subplot(2,1,1)
plot(t,lxreg1mat)
axis([0,150,0.83,0.95])
title('Un-registered functions')
subplot(2,1,2)
plot(t,lxreg2mat)
axis([0,150,0.83,0.95])
title('Registered functions')

%  ------------------------------------------------------------
%    Registration of landmark-registered curves to PCA over
%    the interval of 40 to 90
%  ------------------------------------------------------------

%  select a subset of the data

index = find(t >= 40 & t <= 90);

%  plot the landmark registered data

figure(1)
subplot(1,1,1)
phdl = plot(t(index), lxreg1mat(index,:));
set(phdl, 'LineWidth', 1)
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} log_{10} intensity y(t)') 

%  define the functions over [0,1]

tsmall   = t(index);
rangeval = [40,90];
nsmall   = length(tsmall);
Ymat0     = lxreg1mat(index,:);  %  the log registered data

%  center the data

Ymat0 = Ymat0 - ones(nshort,1)*mean(Ymat0);

%  plot the centered landmark registered data

figure(1)
subplot(1,1,1)
phdl = plot(tsmall, Ymat0);
set(phdl, 'LineWidth', 1)
xlabel('\fontsize{16} Time')
ylabel('\fontsize{16} Centered log_{10} intensity y(t)') 

%  smooth these data

%  set up a basis for smoothing
lxsmallbasis = create_bspline_basis(rangeval,nsmall+2);
%  define smoothing
lxsmallfdPar = fdPar(lxsmallbasis,2,1e-10);
%  smooth the data
lxreg1smallfd = smooth_basis(tsmall, Ymat0, lxsmallfdPar);

%  define dimensions

N = 15;
K = 3;

%  define Abasis

nAbasis = 103;
Abasis = create_bspline_basis(rangeval, nAbasis);
Acoef0 = randn(nAbasis,K);
Afd0   = fd(Acoef0, Abasis);
AfdPar = fdPar(Afd0);

nfine  = max(10*nAbasis+1,501);
tfine  = linspace(rangeval(1),rangeval(2),nfine)';
Ymat0  = eval_fd(tfine, lxreg1smallfd);

%  define Wbasis

nWbasis = 13;
nWorder = 3;
Wbasis = create_bspline_basis(rangeval,nWbasis,nWorder);
Wfd = fd(zeros(nWbasis,1),Wbasis);
Wlambda = 1e0;
Kmat    = eval_penalty(Wbasis,1);
WfdPar  = fdPar(Wfd, 1, Wlambda, 1, Kmat);

%  coefficients for Wfd's

Dmat0 = zeros(nWbasis,N);

%  initial values for coefficients

cvec0 = randn(K*nAbasis,1);

%  from eigenanalysis

pcastruct = pca_fd(lxreg1smallfd, K, AfdPar, 0);
PCfd0     = pcastruct.harmfd;
plot(pcastruct.values)
harmmat   = eval_fd(tfine,PCfd0);
Afd0      = smooth_basis(tfine, harmmat, AfdPar);
coefmat   = getcoef(Afd0);
cvec0     = reshape(coefmat,nAbasis*K,1);

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

Hoptions = optimset('LargeScale',   'off',    ...
                    'Display',      'iter',   ...
                    'Diagnostics',  'on',     ...
                    'GradObj',      'on',     ...
                    'Hessian',      'off',    ...
                    'MaxIter',      100,      ...
                    'TolFun',       1e-5);

%  compute initial function and gradient values

                       
cvec0   = randn(K*nAbasis,1);
% parvec0 = cartes2sphere(cvec0);
parvec0 = cvec0;

matlabpool open local

global Dmat
Dmat0 = zeros(nWbasis,N);
Dmat = Dmat0;

Wlambda = 1e-4;
WfdPar = putlambda(WfdPar, Wlambda);

tic;
parvec = fminunc(@Hfn_Reg, parvec0, Hoptions, ...
                 lxreg1smallfd, Joptions, AfdPar, WfdPar);
toc

parvec0 = parvec;

[Hval, Hgrad, Fmat, Yhat] = Hfn_Reg(parvec, ...
                           lxreg1smallfd, Joptions, AfdPar, WfdPar);

disp(['RMS gradient = ',num2str(sqrt(mean(Hgrad.^2)))])

% cvec = sphere2cartes(parvec);
cvec = parvec;

Acoefmatest = reshape(cvec, nAbasis, K);
Afdest  = fd(Acoefmatest,Abasis);
Amatest = eval_fd(tfine, Afdest);
AfdPar = fdPar(Afdest);

%  plot coefficient functions

figure(3)
subplot(3,1,1)
phdl = plot(tfine, Amatest(:,1), 'b-', rangeval, [0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} a_1(t)')    
subplot(3,1,2)
phdl = plot(tfine, Amatest(:,2), 'b-', rangeval, [0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} a_2(t)')    
subplot(3,1,3)
phdl = plot(tfine, Amatest(:,3), 'b-', rangeval, [0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} a_3(t)')    

%  3D plot of coefficient functions

figure(4)
subplot(1,1,1)
phdl=plot3(Amatest(:,1),Amatest(:,2),Amatest(:,3),'-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} a_1(t)')    
ylabel('\fontsize{13} a_2(t)')    
zlabel('\fontsize{13} a_3(t)')  
grid on

%  plot factor scores

figure(5)
subplot(1,1,1)
plot(Fmat(1,:),Fmat(2,:),'o');
set(phdl, 'LineWidth', 2)
for i=1:15
    thdl=text(Fmat(1,i)+0.002,Fmat(2,i),num2str(i));
    set(thdl, 'LineWidth', 2)
end
xlabel('\fontsize{13} f_1')    
ylabel('\fontsize{13} f_2')   

%  3D plot of factor scores

figure(6)
phdl=plot3(Fmat(1,:),Fmat(2,:),Fmat(3,:),'o');
set(phdl, 'LineWidth', 2)
grid on
for i=1:15
    thdl=text(Fmat(1,i)+0.005,Fmat(2,i)+0.005,Fmat(3,i)+0.005,num2str(i));
    set(thdl, 'LineWidth', 2)
end
xlabel('\fontsize{13} f_1')    
ylabel('\fontsize{13} f_2')    
zlabel('\fontsize{13} f_3')  

%  plot warping functions

Wfd  = fd(Dmat, Wbasis);

%  deformation function plot

figure(7)
subplot(1,1,1)
plotdeform_fd(Wfd)
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} Deformation d(t) = h(t) - t')
 
%  plot the fit to the data

%  step through cases

MSE = zeros(N,1);
SSE = 0;
figure(8)
subplot(1,1,1)
index = 1:N;
ptindex = 1:10:1031;
for i=index
    MSE(i) = mean((Ymat0(:,i) - Yhat(:,i)).^2);
    SSE = SSE + sum((Ymat0(:,i) - Yhat(:,i)).^2);
    phdl = plot(tfine(ptindex), Ymat0(ptindex,i), 'b-', ...
                tfine(ptindex), Yhat(ptindex,i), 'b:');
    set(phdl, 'LineWidth', 2)
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} y(t)') 
%     title(['\fontsize{16} Record ', num2str(i), ...
%            ',  RMSE = ',num2str(sqrt(MSE(i)))])
    if length(index) > 1
        pause
    end
end
MSE1 = mean(MSE);

%  plot the first and last fits

figure(9)
subplot(2,1,1)
i=1;
phdl = plot(tfine(ptindex), Ymat0(ptindex,i), 'b-', ...
    tfine(ptindex), Yhat(ptindex,i), 'b:');
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{16} y_1(t)')
subplot(2,1,2)
i=15;
phdl = plot(tfine(ptindex), Ymat0(ptindex,i), 'b-', ...
    tfine(ptindex), Yhat(ptindex,i), 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} t')
ylabel('\fontsize{16} y_{15}(t)')

%  plot the first and last fits and deformations

figure(10)
subplot(2,2,1)
i=1;
phdl = plot(tfine(ptindex), Ymat0(ptindex,i), 'b-', ...
    tfine(ptindex), Yhat(ptindex,i), 'b:');
set(phdl, 'LineWidth', 1)
ylabel('\fontsize{16} y_1(t)')
axis([40,90,-0.02,0.06])
subplot(2,2,2)
plotdeform_fd(Wfd(i))
xlabel('')
ylabel('\fontsize{16} d_1(t) = h_1(t) - t')
axis([40,90,-1,1])
i=15;
subplot(2,2,3)
phdl = plot(tfine(ptindex), Ymat0(ptindex,i), 'b-', ...
    tfine(ptindex), Yhat(ptindex,i), 'b:');
set(phdl, 'LineWidth', 1)
xlabel('\fontsize{16} Time t (minutes)')
ylabel('\fontsize{16} y_{15}(t)')
axis([40,90,-0.02,0.06])
subplot(2,2,4)
plotdeform_fd(Wfd(i))
xlabel('\fontsize{16} Time t (minutes)')
ylabel('\fontsize{16} d_{15} = h_{15}(t) - t')
axis([40,90,-1,1])

%  display root mean squared errors

disp(['RMSE = ',num2str(sqrt(MSE1))])
% RMSE = 0.0038227

%  fit the unregistered data with three components

[U,S,V] = svd(Ymat0,0);
Yhat0 = U(:,1:3)*S(1:3,1:3)*V(:,1:3)';
Rmat0 = Ymat0 - Yhat0;
MSE0 = mean(mean(Rmat0.^2));

disp(['RMSE0 = ',num2str(sqrt(mean(MSE0)))])
% RMSE0 = 0.0052297

RSQ = (MSE0 - MSE1)/MSE0;

disp(['RSQ = ',num2str(RSQ)])
% RSQ = 0.4657

save proteomics





% Wlambda = 1e0

% Wlamda = 1e-4;
% RMSE = 0.0053136

