function [Hval, Hgrad, Fmat, Yhat] = Hfn_Reg_trapz(cvec, tfine, Yfine, ...
        AfdPar, WfdPar, PenStruct, Joptions, periodic, crit, dbgwrd)
%  HFN_FDA evaluates the loss function and its gradient for a parameter
%  cascading approach to principal components analysis of 
%  functional N by n data YMAT using K << n factors.
%  Outer loss function is H(A|Y) = ||Y - F(A)A||^2
%  Inner loss function is J(F|A,Y) = ||Y = FA|| + lambda ||F'RF||
%
%  Arguments:
%
%  CVEC    ... Vectorized K by nbasis PC basis coefficient matrix
%  A single struct variable FITSTRUCT that contains: 
%
%  Last modified 26 July 2012

global Dmat

if nargin < 10, dbgwrd = 0;   end
if nargin <  9, crit = 2;     end
if nargin <  8, periodic = 0; end
if nargin <  7
    Joptions = optimset('LargeScale',  'off',  ...
                        'Display',     'off',  ...
                        'Diagnostics', 'off',  ...
                        'GradObj',     'on',   ...
                        'Hessian',     'off',  ...  
                        'MaxIter',     50,     ...
                        'TolFun',      1e-6);
end
if nargin < 6 || isempty(PenStruct)  
    PenStruct.FRowMat = [];
    PenStruct.FRowPar = 0;
    ANrmPar = 0;
else
    ANrmPar = PenStruct.ANrmPar;
end

%  update AfdPar

Afd          = getfd(AfdPar);
Acoef        = getcoef(Afd);
[nAbasis, K] = size(Acoef);

Cmat = reshape(cvec, nAbasis, K);

Afd    = putcoef(Afd, Cmat);
AfdPar = putfd(AfdPar, Afd);

%  optimize wrt warping functions

ncvec  = length(Acoef(:));

%  get dimensions of data and convert AVEC to K by n matrix AMAT

[nfine,N] = size(Yfine);
nWbasis   = size(Dmat,1);

%  check dimensions of Dmat

if size(Dmat,1) ~= nWbasis
    error('First dimension of Dmat not equal to no. basis functions.');
end
if size(Dmat,2) ~= N
    error('Second dimension of Dmat not equal to no. cases.');
end

%  screen Dmat for large values

lrgdvec = 0;
for i=1:N
    if any(abs(Dmat(:,i)) > 10)
        disp(['Case ', num2str(i), ',  dvec0: ', num2str(Dmat(:,i)')])
        lrgdvec = 1;
    end
end
if lrgdvec
    error('One of more coefficients outside of [-10,10].');
end

%  optimize wrt warping functions

Hval  = 0;
Hgrad = zeros(ncvec,1);
Fmat  = zeros(K,N);
Yhat  = zeros(nfine,N);
DmatTmp = Dmat;

% dbgwrd = 1;
parfor i=1:N
    yfine = Yfine(:,i);
    dvec0 = DmatTmp(:,i);
%     tic;
%     dvec1  = fminunc(@PCA_Reg, dvec0, Joptions, ...
%                     tfine, yfine, AfdPar, WfdPar, ...
%                     Flambda, Smat, periodic, crit);
%     toc
%     tic;
    dvec  = fminunc(@PCA_Reg_trapz, dvec0, Joptions, ...
                    tfine, yfine, AfdPar, WfdPar, ...
                    PenStruct, periodic, crit);
%     toc
    if dbgwrd
         disp(['Case ',num2str(i),'  dvec: ',num2str(dvec')])
    end
    DmatTmp(:,i) = dvec;
    [Fval, Fgrad, Fhess, Fcross, fvec, yhat, Hi, Dc_Hi, Wmat] = ...
            PCA_Reg_trapz(dvec, tfine, yfine, AfdPar, WfdPar, ...
                    PenStruct, periodic, crit);
    Fmat(:,i) = fvec;
    Yhat(:,i) = yhat;
    Hval      = Hval  + Hi/N;
    Hgrad     = Hgrad + Dc_Hi/N;
    if ANrmPar > 0
        ANrmDif = trace(Cmat'*Wmat*Cmat)/nAbasis - 1;
        Hval  = Hval  + ANrmPar*ANrmDif^2;
        Hgrad = Hgrad + ...
            4*ANrmPar*ANrmDif*reshape(Cmat'*Wmat,K*nAbasis,1)/nAbasis;
    end
    if dbgwrd
        disp(['Fgrad: ',num2str(Fgrad')])
%         disp('Fhess:')
%         disp(Fhess)
        eigval = eig(Fhess);
        disp(['Eigenvalues: ',num2str(eigval')])
    end
end

Dmat = DmatTmp;

%  extract H function value and gradient

% Hval  = fitstruct.Hval;
% Hgrad = fitstruct.Hgrad;

if dbgwrd
    disp(cvec')
    disp(Hval)
    disp(Hgrad')
end



