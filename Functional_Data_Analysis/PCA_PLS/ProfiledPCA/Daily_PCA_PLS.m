addpath('../fdaM')
addpath('../fdaM/examples/weather')

%  This code registers the daily weather data to its principal components
%  using function HfnReg.

%  Last modified 15 November 2014

%%  -----------------------------------------------------------------------
%                 set up the data for analysis and display them
%  ------------------------------------------------------------------------

%  load the daily weather data

load daily

%  load the results of the last analysis of these data

load Daily_PCA_PLS

%  convert data to display winter in the center of the plot
%  and start of plot on July 1st, end on June 30

winterindex = [182:365,1:181];
tempav = tempav(winterindex,:);
precav = precav(winterindex,:);

%  two-character place names for plotting

place2 = [ ...
    'Aa'; ...
    'Be'; ...
    'Cy'; ...
    'Cn'; ...
    'Cl'; ...
    'Dn'; ...
    'En'; ...
    'Fn'; ...
    'Hx'; ...
    'Ik'; ...
    'It'; ...
    'Ks'; ...
    'Ln'; ...
    'Ml'; ...
    'Oa'; ...
    'PA'; ...
    'PG'; ...
    'PR'; ...
    'Qc'; ...
    'Ra'; ...
    'Re'; ...
    'Sl'; ...
    'Se'; ...
    'SJ'; ...
    'Sy'; ...
    'TP'; ...
    'TB'; ...
    'To'; ...
    'UC'; ...
    'Vr'; ...
    'Va'; ...
    'We'; ...
    'Wg'; ...
    'Yh'; ...
    'Ye'];

%  set up indices that order the stations from east to west to north

geogindex = [24,  9, 25, 34,  4,  8, 22,  1,  2, 19, 23, 14, 15, 28, ...
             13, 27, 33, 26,  5, 20, 16, 29,  7,  3, 12, 30, 31, 17, ...
             18, 32,  6, 35, 11, 10, 21];

%  sort place names by east-west

place2 = place2(geogindex,:);

%  -------------  set up fourier basis  ---------------------------
%  Here it was decided that 65 basis functions captured enough of
%  the detail in the temperature data: about one basis function
%  per week.  However, see below for smoothing with a saturated
%  basis (365 basis functions) where smoothing is defined by the
%  GCV criterion.

daytime = daytime/365;
dayrange = [0,1];

%  low resolution basis

nbasis     = 13;
daybasis13 = create_fourier_basis(dayrange, nbasis);

%  medium resolution basis

nbasis     = 65;
daybasis65 = create_fourier_basis(dayrange, nbasis);

%  harmonic acceleration linear differential operator

Lbasis  = create_constant_basis(dayrange);  %  create a constant basis
Lcoef   = [0,(2*pi/1)^2,0];      %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

%  ---------  create fd objects for temp. and prec. ------------

daytempfd = smooth_basis(daytime, tempav, daybasis13);
daytempfd_fdnames{1} = 'Proportion of year';
daytempfd_fdnames{2} = 'Station';
daytempfd_fdnames{3} = 'Deg C';
daytempfd = putnames(daytempfd, daytempfd_fdnames);

% figure(1) 
% subplot(1,1,1)
% plotfit_fd(tempav, daytime, daytempfd)

dayprecfd = smooth_basis(daytime, precav, daybasis13);
dayprecfd_fdnames{1} = 'Day';
dayprecfd_fdnames{2} = 'Station';
dayprecfd_fdnames{3} = 'mm';
dayprecfd = putnames(dayprecfd, dayprecfd_fdnames);

for i=1:35
    ind = find(precav(:,i) == 0);
    if length(ind) > 0
        precav(ind,i) = 0.1;
    end
end
lnprecav = log10(precav);
daylnprecfd = smooth_basis(daytime, lnprecav, daybasis13);
daylnprecfd_fdnames{1} = 'Day';
daylnprecfd_fdnames{2} = 'Station';
daylnprecfd_fdnames{3} = 'log_{10} mm';
daylnprecfd = putnames(daylnprecfd, daylnprecfd_fdnames);

figure(1)
dayprecmat = eval_fd(daytime, daytempfd);
phdl = plot(daytime, dayprecmat, 'b-', [0,1], [0,0], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Year (July 1 to June 30)')
ylabel('\fontsize{16} Temperature (deg C)')

%%  -----------------------------------------------------------------------
%    PLS/PCA predicting mean precipitation from daily temperature data
%  ------------------------------------------------------------------------

%  set up mean daily precipitation in tenths of a millimetre

meanprec = mean(precav)'*10;

%  set options for optimization using Matlab function fminunc

options = optimset('LargeScale',  'off',    ...
                   'Display',     'iter',   ...
                   'Diagnostics', 'on',     ...
                   'GradObj',     'on',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-6);

%  ----------  generate initial values for parameters  ---------------

K = 2;

%  define Abasis

nAbasis = 23;
Abasis = create_bspline_basis(dayrange, nAbasis);

%  Compute Grassman map, which in this case is 
%    NABASIS*K by (NABASIS-K)*K

GrssMnPts = linspace(0,1,nAbasis*K)';
GrssBasis = create_fourier_basis([0,1],(nAbasis-K)*K);
GrssMnMap = eval_basis(GrssMnPts,GrssBasis);
GrssMnMap = GrssMnMap(:,1:(nAbasis-K)*K);
[GrssMnMap,S,V] = svd(GrssMnMap,0);

%  cvec0 from eigenanalysis

%  carry out PCA
pcastruct = pca_fd(daytempfd, K, Abasis, 0);
%  set up matrix of harmonic values
harmmat   = eval_fd(daytime,pcastruct.harmfd);
%  convnert harmonic values to functional data objects by smoothing
PCfd0     = smooth_basis(daytime, harmmat, Abasis);
%  extract the coefficients
coefmat   = getcoef(PCfd0);
%  reshape these into vector format for use with fminunc
cvecPCA   = reshape(coefmat,nAbasis*K,1);
%  extract the factor scores
pcaFmat = pcastruct.harmscr;
%  dipslay proportions of variance
disp(pcastruct.varprop)

%     0.7556
%     0.2333

Hlambda = 0;  %  smoothing parameter for outer criterion
%  set up functional parameter object for harmonics
Rmat    = []; 
AfdPar  = fdPar(PCfd0, 2, Hlambda, 1, Rmat);

%  define the penalty structure (all null)

PenStruct.FRowMat = [];
PenStruct.FRowPar = 0;
PenStruct.FColMat = [];
PenStruct.FColPar = 0;

%%  carry out the pca with gamma = 0 to give straight PCA solution

gamval = 0;

cvec0 = cvecPCA;  %  initial coefficient values

parvec0 = cvec0;  %  unaltered coefficient vector
% parvec0 = cartes2sphere(cvec0);  %  convert to spherical coordinates
% parvec0 = GrssMnMap'*cvec0;  %  convert to Grassmann coordinates
% parvec0 = cartes2sphere(GrssMnMap'*cvec0);  % convert to Grassman and  
%                                             then to spherical coordinates

%  run with unaltered coordiantes

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, daytempfd, AfdPar, ...
                 meanprec, gamval, PenStruct, 1, 0, 0);
toc

%  converged in 1 iteration, 0.02 seconds

%  run with spherical coordinates

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, daytempfd, AfdPar, ...
                 meanprec, gamval, PenStruct, 1, 1, 0);
toc

%  run with Grassmann coordinates

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, daytempfd, AfdPar, ...
                 meanprec, gamval, PenStruct, 1, 0, 1);
toc

%  Grassmann coordinates converged in 99 iterations, 0.41 seconds

%  run with Grassmann/spherical coordinates

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, daytempfd, AfdPar, ...
                 meanprec, gamval, PenStruct, 1, 1, 1);
toc

%  Grassmann/spherical coordinates failed to converge in 401 iterations, 
%  took 8.74 seconds

%  evaluate solution at final value

[Hval0, Hgrad0, Fmat0, meanprechat0] = ...
    Hfn_FDA(parvec, daytempfd, AfdPar, meanprec, gamval, PenStruct, ...
            1, 0, 0);

disp(['Hval = ',num2str(Hval0)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad0.^2)))])

% Hval = 65.2386
% Hgrad RMSE = 1.2234e-05

cvec = parvec;  % final coefficient vector
% cvec = sphere2cartes(parvec,1);
% cvec = GrssMnMap*parvec;
% cvec = GrssMnMap*sphere2cartes(parvec);  %  Default

%  set up the estimated functional data object for harmonics

Afdest0 = fd(reshape(cvec, nAbasis, K),Abasis);
Amat0   = eval_fd(daytime, Afdest0);  %  matrix of values at days

%  plot principal component functions

figure(2)
subplot(1,1,1)
phdl=plot(daytime, Amat0(:,1), 'b-', ...
          daytime, Amat0(:,2), 'r-', ...
          dayrange,[0,0], 'r:');
% phdl=plot(daytime, Amat0(:,1), 'b-', ...
%           daytime, Amat0(:,2),  'b--', ...
%           daytime, Amat0(:,3),  'b:', ...
%           dayrange,[0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Temperature (deg C)')
% axis([3,18,-0.2,0.4])
legend('\fontsize{16} Component I', '\fontsize{16} Component II', ...
       'Location', 'SouthWest')
title('\fontsize{16} Principal Components Analysis')

%  plot factor scores

figure(3)
h = 0.02;
subplot(1,1,1)
phdl=plot(Fmat0(:,1), Fmat0(:,2), '.');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(Fmat0(i,1)+h, Fmat0(i,2), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
xlabel('\fontsize{16} Component I  Scores')
ylabel('\fontsize{16} Component II Scores')
title('\fontsize{16} Principal Components Analysis')

%  plot fit to meanprec

disp(['PLS correlation = ',num2str(corr(meanprec, meanprechat0))])

% PLS correlation = 0.39657

figure(4)
h = 1.0;
subplot(1,1,1)
phdl=plot(meanprechat0, meanprec, '.', ...
         [-20,40], [-20,40], 'b--');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(meanprechat0(i)+h, meanprec(i), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
xlabel('\fontsize{16} Fit to mean precipitation (mm*10)')
ylabel('\fontsize{16} Mean precipitation (mm*10)')
title('\fontsize{16} PCA correlation = 0.40')

%%  carry out the pca with gamma = 1 to give straight PLS solution

gamval = 1.0;

cvec0 = cvecPCA;

parvec0 = cvec0;
% parvec0 = cartes2sphere(cvec0);
% parvec0 = GrssMnMap'*cvec0;
% parvec0 = cartes2sphere(GrssMnMap'*cvec0);  %  Default

tic;
parvec = fminunc(@Hfn_FDA, parvec0, options, daytempfd, AfdPar, ...
                 meanprec, gamval, PenStruct, 1, 0, 0);
toc

%  converged in 85 iterations in 0.42 seconds

[Hval1, Hgrad1, Fmat1, meanprechat1] = ...
    Hfn_FDA(parvec, daytempfd, AfdPar, meanprec, gamval, PenStruct, ...
            1, 0, 0);

disp(['Hval = ',num2str(Hval1)])
disp(['Hgrad RMSE = ',num2str(sqrt(mean(Hgrad1.^2)))])

% Hval = 0.0025479
% Hgrad RMSE = 1.2979e-07

cvec = parvec;
% cvec = sphere2cartes(parvec,1);
% cvec = GrssMnMap*parvec;
% cvec = GrssMnMap*sphere2cartes(parvec);  %  Default

%  PLS harmonic object

Afdest1 = fd(reshape(cvec, nAbasis, K),Abasis);
Amat1   = eval_fd(daytime, Afdest1);

%  plot PLS component functions

figure(4)
subplot(1,1,1)
phdl=plot(daytime*365, Amat1(:,1), 'b-', ...
          daytime*365, Amat1(:,2), 'r-', ...
          dayrange*365,[0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Temperature (deg C)')
axis([0,365,-150,150])
legend('\fontsize{16} Component I', '\fontsize{16} Component II', ...
       'Location', 'SouthEast')
title('\fontsize{16} Partial Least Squares')

%  Plot both PCA and PLS component functions.  
%  In these two-panel plots the sign of the second harmonic has been
%    reversed in order to display the results in a way that leads
%    to a simpler descriptionm.

ind = 3:5:363;
figure(2)

subplot(1,2,1)
phdl=plot(daytime(ind),  Amat0(ind,1), 'b-', ...
          daytime(ind),  Amat0(ind,2), 'b--', ...
          dayrange,[0,0], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} PCA Components')
axis('square')

subplot(1,2,2)
phdl=plot(daytime(ind),  Amat1(ind,1), 'b-', ...
          daytime(ind),  Amat1(ind,2), 'b--', ...
          dayrange,[0,0], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} PLS Components')
axis('square')

%  plot factor scores

figure(5)
h = 0.02;
subplot(1,1,1)
phdl=plot(Fmat1(:,1), Fmat1(:,2), 'b.');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(Fmat1(i,1)+h, Fmat1(i,2), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
xlabel('\fontsize{16} Component I  Scores')
ylabel('\fontsize{16} Component II Scores')
title('\fontsize{16} Partial Least Squares')

%  plot both PCA and PLS factor scores

figure(3)

subplot(1,2,1)
h = 0.05;
phdl=plot(Fmat0(:,1), Fmat0(:,2), 'bo', ...
          [0.5,3.0], [0,0], 'b:');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(Fmat0(i,1)+h, Fmat0(i,2), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
xlabel('\fontsize{16} PCA I  Scores')
ylabel('\fontsize{16} PCA II Scores')
axis('square')

subplot(1,2,2)
h = 0.05;
phdl=plot(Fmat1(:,1), Fmat1(:,2), 'bo', ...
          [-0.2,0.8], [0,0], 'b:', [0,0], [-0.35,0.04], 'b:');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(Fmat1(i,1)+h, Fmat1(i,2), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
xlabel('\fontsize{16} PLS I  Scores')
ylabel('\fontsize{16} PLS II Scores')
axis('square')

%  plot fit to meanprec

disp(['PLS correlation = ',num2str(corr(meanprec, meanprechat1))])

% PLS correlation = 0.78136
 
figure(6)
h = 1.0;
phdl = plot(meanprechat1, meanprec, 'b.', ...
            [-5,45], [-5,45], 'b--');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(meanprechat1(i)+h, meanprec(i), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
xlabel('\fontsize{16} Fit to mean precipitation (mm*10)')
ylabel('\fontsize{16} Mean precipitation (mm*10)')
title('\fontsize{16} PLS correlation = 0.78')
   
%  plot both fits

figure(4)

subplot(1,2,1)
h = 2.0;
phdl=plot(meanprechat0, meanprec, 'o', ...
         [-20,40], [-20,40], 'b--');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(meanprechat0(i)+h, meanprec(i), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
axis([-15,45,0,50])
axis('square')
xlabel('\fontsize{16} PCA fit (mm*10)')
ylabel('\fontsize{16} Mean precipitation (mm*10)')

subplot(1,2,2)
h = 2.0;
phdl = plot(meanprechat1, meanprec, 'bo', ...
            [-5,45], [-5,45], 'b--');
set(phdl, 'LineWidth', 2)
hold on
for i=1:35
    thdl=text(meanprechat1(i)+h, meanprec(i), place2(i,:));
    set(thdl, 'FontSize',13)
end
hold off
axis([-15,45,0,50])
axis('square')
xlabel('\fontsize{16} PLS fit (mm*10)')
ylabel('\fontsize{16} Mean precipitation  (mm*10)')

% canonical correlations between PCA and PLS Fmats

Fmat0ctr = Fmat0 - ones(35,1)*mean(Fmat0);
Fmat1ctr = Fmat0 - ones(35,1)*mean(Fmat1);

[U1, S1, V1] = svd(Fmat0ctr,0);
[U2, S2, V2] = svd(Fmat1ctr,0);
CCR = svd(U1'*U2);
disp(['Canonical correlations:   ', num2str(CCR')])

% Canonical correlations:   1.000     0.175

%%  compute prediction over a grid of gamma values

gamvec = linspace(0.4,0.5,11);
ngam   = length(gamvec);
corrstore = zeros(ngam,1);
parvstore = zeros(ngam,length(cvec0));
for igam=1:ngam
    gami = gamvec(igam);
    if igam == 1
        parvec0i = parvec0;
    else
        parvec0i = parvstore(igam-1,:)';
    end
    parvec = fminunc(@Hfn_FDA, parvec0i, options, daytempfd, AfdPar, ...
        meanprec, gami, PenStruct, 1, 0, 0);
    [Hvali, Hgradi, Fmati, meanprechati] = ...
        Hfn_FDA(parvec, daytempfd, AfdPar, meanprec, gami, PenStruct, ...
                1, 0, 0);
    parvstore(igam,:) = parvec';
    corrstore(igam) = corr(meanprec, meanprechati);
end
