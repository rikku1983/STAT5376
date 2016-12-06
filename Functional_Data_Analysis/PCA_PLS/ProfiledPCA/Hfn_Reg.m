function [Hval, Hgrad, Fmat, Yhat] = Hfn_Reg(parvec, ...
                           Xfd, Joptions, AfdPar, WfdPar, ...
                           Flambda, Smat, periodic, crit, anglewrd, dbgwrd)
%  HFN_FDA evaluates the loss function and its gradient for a parameter
%  cascading approach to principal components analysis of 
%  functional N by n data YMAT using K << n factors.
%  Outer loss function is H(A|Y) = ||Y - F(A)A||^2
%  Inner loss function is J(F|A,Y) = ||Y - FA|| + lambda ||F'RF||
%
%  Arguments:
%
%  PARVEC  ... Vectorized K by N BASIS PC basis coefficient matrix,
%                converted to hyperspherical coordinates if 
%                ANGLEWRD ~= 0
%  XFD     ... Functional data object for the observed functions
%  AFDPAR  ... A functional parameter object for the K principal component 
%              functions.  The coefficient matrix is the transpose of
%              coefficient matrix CMAT defining these functions
%  WFDPAR  ... A functional parameter object for the function W_i(t) 
%              defining the initial warping functions h_i(t).
%  FLAMBDA ... roughness penalty parameter 
%  SMAT    ... Roughness penalty matrix for scores
%  PERIODIC... 1 implies the data are periodic and a shift is estimated.
%              In this case the initial coefficient in DVEC is the
%              time-shift
%              0 implies data are not periodic.  In this case the initial
%              coefficient in DVEC is not used.
%  CRIT    ... 1, 2, or 3 defining the fitting criterion. 
%              2 implies the minimum eigenvalue of the cross-product
%              matrix
%  ANGLEWRD .. If nonzero, PARVEC is transformed into hyperspherical 
%              coordinates.  Defaults to 0
%  DBGWRD  ... Level of output on each iteration
%
%  Returns:
%
%  HVAL   ...  Value of objective function H
%  HGRAD  ...  Gradient of H with respect to PARVEC
%  FMAT   ...  Matrix of component scores
%  YHAT   ...  Approximation to the data by the model
%
%  Last modified 9 June 2014

global Dmat

if nargin < 12, dbgwrd = 0;   end
if nargin < 11, anglewrd = 0; end
if nargin < 10, crit = 2;     end
if nargin <  9, periodic = 0; end
if nargin <  8, Smat = [];    end
if nargin <  7, Flambda = 0;  end

%  update AfdPar

Afd            = getfd(AfdPar);
PCcoef         = getcoef(Afd);
[nAbasis, K]   = size(PCcoef);

%  convert PARVEC to CVEC,  and CVEC to K by n matrix CMAT

if anglewrd
    [cvec,Dcvec] = sphere2cartes(parvec,1);
else
    cvec = parvec;
end

Cmat   = reshape(cvec, nAbasis, K);
Afd    = putcoef(Afd,Cmat);
AfdPar = putfd(AfdPar, Afd);

%  check WfdPar for consistency with periodic argument

Wbasis = getbasis(getfd(WfdPar));
if periodic
   if ~strcmp(getbasistype(Wbasis), 'fourier')
      error('PERIODIC is true, Wbasis not of fourier type.');
   end
end

%  optimize wrt warping functions

ncvec  = length(PCcoef(:));

%  get dimensions of data and convert AVEC to K by n matrix AMAT

[nXbasis,N] = size(getcoef(Xfd));
nWbasis   = size(Dmat,1);

%  screen Dmat for large values

climit = 50;
lrgdvec = 0;
for i=1:N
    if any(abs(Dmat(:,i)) > climit)
        disp(['Case ', num2str(i), ',  dvec0: ', num2str(Dmat(:,i)')])
        lrgdvec = 1;
    end
end
if lrgdvec
    error('One of more coefficients outside of limit.');
end

%  optimize wrt warping functions

nfine   = max(10*nAbasis+1,501);

Hval  = 0;
Hgrad = zeros(ncvec,1);
Fmat  = zeros(K,N);
Yhat  = zeros(nfine,N);
DmatTmp = Dmat;

% dbgwrd = 1;
parfor i=1:N
    Xfdi  = Xfd(i);
    dvec0 = DmatTmp(:,i);
    dvec  = fminunc(@PCA_Reg, dvec0, Joptions, ...
                    Xfd(i), AfdPar, WfdPar, ...
                    Flambda, Smat, periodic, crit);
    if dbgwrd
         disp(['Case ',num2str(i),'  dvec: ',num2str(dvec')])
    end
    DmatTmp(:,i) = dvec;
    [Fval, Fgrad, Fhess, Fcross, fvec, yhat, Hi, Dc_Hi] = ...
            PCA_Reg(dvec, Xfdi, AfdPar, WfdPar);
    Fmat(:,i) = fvec;
    Yhat(:,i) = yhat;
    Hval      = Hval  + Hi/N;
    Hgrad     = Hgrad + Dc_Hi/N;
    if dbgwrd
        disp(['Fgrad: ',num2str(Fgrad')])
        eigval = eig(Fhess);
        disp(['Eigenvalues: ',num2str(eigval')])
    end
end

Dmat = DmatTmp;

%  convert gradient to angular coordinates

if anglewrd
    Hgrad = Dcvec'*Hgrad;
end

if dbgwrd
    disp(cvec')
    disp(Hval)
    disp(Hgrad')
end



