function [Fval, Fgrad, Fhess, Fcross, fvec, yhat, Hi, Dc_Hi] = ...
                   PCA_Reg(dvec, Xfd, AfdPar, WfdPar, ...
                           Flambda, Smat, periodic, crit)
%  PCA_REG  evaluates the function and gradient for the inner objective
%  function in PCA with registration.  It is optimized in Hfn_Reg 
%  for each record through a call to Matlab optimizer fminunc 
%  Inner loss function is J(F|A,Y) = ||Y - FA|| + lambda ||F'RF||
%  
%  Arguments:
%  DVEC    ... A vector of coefficients defining the function W_i(t) that
%              in turn defines the time-warping function h_i(t).
%  TFINE   ... Column vector of NFINE times of observation for the data 
%              values
%  YFINE   ... Column vector of NFINE data values corresponding to
%              observation tines in TFINE
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
%
%  Returned objects:
%  FVAL   ... The function value
%  FGRAD  ... The gradient wrt to DVEC
%  FHESS  ... The second derivative wrt DVEC
%  FCROSS ... The second cross-derivative wrt to CVEC and DVEC
%  FVEC   ... The princpal component scores for this record
%  YHAT   ... The approximation to the record's data
%  HI     ... The integrated squared residual
%  DCHI   ... The derivative wrt C of HI
%
%  Last modified 29 January 2013 by Jim Ramsay


if nargin <  9, crit = 2;     end
if nargin <  8, periodic = 0; end
if nargin <  7, Smat = [];    end
if nargin <  6, Flambda = 0;  end

Afd     = getfd(AfdPar);
Abasis  = getbasis(Afd);
nAbasis = getnbasis(Abasis);

Wfd     = getfd(WfdPar);
Wbasis  = getbasis(Wfd);
Kmat    = getpenmat(WfdPar);
Wlambda = getlambda(WfdPar);

trng   = getbasisrange(Wbasis);
tlo    = trng(1);
thi    = trng(2);
width  = thi - tlo;
nfine  = max(10*nAbasis+1,501);
tfine  = linspace(tlo,thi,nfine)';
yfine  = eval_fd(tfine, Xfd);

ndvec  = length(dvec);
K      = size(getcoef(Afd),2);

if periodic
    shift  = dvec(1);
    dvec(1) = 0;
end

%  W function for input coefficients

nWbasis = getnbasis(Wbasis);
Wfd     = fd(dvec, Wbasis);

%  Warping function values defined by these coefficients

ffine        = monfn(tfine, Wfd);
fmax         = ffine(nfine);
hfine        = tlo + width.*ffine./fmax;
hfine(1)     = tlo;
hfine(nfine) = thi;

%  carry out shift if period and shift ~= 0

if periodic && shift ~= 0
   yfine = shifty(tfine, yfine, shift);
end

%  evaluate Amat ... the principal component values 

% Proposed new code, but needs checking, 29 Jan 13
% Amat = eval_fd(hfine, Afd);
% Afd  = smooth_basis(hfine, Amat, Abasis);
% Cmat = getcoef(Afd);

Phimat  = eval_basis(hfine,Abasis);
Cmat    = getcoef(Afd);
Amat    = Phimat*Cmat;

%  compute fvec ... the K factor scores for this case

Wmat    = Phimat'*Phimat;
CWCtmat = Cmat'*Wmat*Cmat;
Mmat    = CWCtmat;
if ~isempty(Smat)
    Mmat = Mmat + Flambda.*Smat;
end
Mmatinv = inv(Mmat);
fvec    = Mmatinv*Amat'*yfine;

%  compute predicted values 

Fmat = ones(nfine,1)*fvec';
yhat = sum(Fmat.*Amat,2);

%  compute criterion value according to the choice of criterion

aa = mean(yfine.^2);
bb = mean(yfine.*yhat);
cc = mean(yhat.^2);

if crit == 1
    %  least squares criterion
    res  = yfine - yhat;
    Fval = aa - 2*bb + cc;
elseif crit == 2
    % minimum eigenvalue criterion
    ff   = aa - cc;
    dd   = sqrt(ff^2 + 4*bb^2);
%     Fval = sum(ffine);
%     Fval = sum(hfine);
%     Fval = sum(Amat(:,1));
%     Fval = sum(fvec);
%     Fval = sum(yhat);
    Fval = aa + cc;
    Fval = Fval - dd;
else
    error('Inadmissible value for CRIT.');
end

if ~isempty(Kmat) && Wlambda > 0
    Fval  = Fval + Wlambda.*dvec'*Kmat*dvec;
end

Dd_yhat = [];

%  ------------------------------------------------------------------------
%                         Compute gradient if needed
%  ------------------------------------------------------------------------

if nargout > 1

    %  derivative of warping function with respect dvec
    
    Dd_ffine = mongrad(tfine, Wfd);
    Dd_fmax  = Dd_ffine(nfine,:);
    Dd_h     = width.*(fmax.*Dd_ffine - ffine*Dd_fmax)./fmax^2;
    if periodic
        Dd_h(:,1) = 1;
    else
        Dd_h(:,1) = 0;
    end
    
    %  compute derivative of Amat wrt dvec
    
    ndvec = length(dvec); 
    onecoef = ones(1,ndvec);
    Dd_Amat = zeros(nfine,K,ndvec);
    D_A     = eval_fd(hfine, Afd, 1);
    for k=1:K
        Dd_Amat(:,k,:) = (D_A(:,k)*onecoef).*Dd_h;
    end
    
    %  compute derivative of Mmat wrt dvec
    %  compute derivative of fvec and yhat wrt dvec
    
    Dd_Mmat    = zeros(K,K,ndvec);
    Dd_Mmatinv = Dd_Mmat;
    Dd_fvec    = zeros(K,ndvec);
    Dd_yhat    = zeros(nfine,ndvec);
    for r=1:ndvec
        Dd_Amatr    = squeeze(Dd_Amat(:,:,r));
        Dd_Mmatr    = Amat'*Dd_Amatr + Dd_Amatr'*Amat;
        Dd_Mmatinvr = -Mmatinv*Dd_Mmatr*Mmatinv;
        Dd_Mmat(:,:,r)    = Dd_Mmatr;
        Dd_Mmatinv(:,:,r) = Dd_Mmatinvr;
        Dd_fvec(:,r) = (Dd_Mmatinvr*Amat' + Mmatinv*Dd_Amatr')*yfine;
        Dd_yhat(:,r) = Amat*Dd_fvec(:,r) + Dd_Amatr*fvec;
    end
    
    %  compute gradient according to the choice of criterion
    
    if crit == 1
        %  least squares criterion
        Fgrad = -2.*Dd_yhat'*res./nfine;
    elseif crit == 2
        % minimum eigenvalue criterion
        Dd_bb   =    Dd_yhat'*yfine./nfine;
        Dd_cc   = 2.*Dd_yhat'*yhat./nfine;
        gg      = 4.*bb.*Dd_bb - (aa-cc).*Dd_cc;
        Dd_dd   = gg./dd;
%         Fgrad   = sum(Dd_ffine)';
%         Fgrad   = sum(Dd_h)';
%         Fgrad   = sum(squeeze(Dd_Amat(:,1,:)))';
%         Fgrad   = sum(Dd_fvec)';
%         Fgrad   = sum(Dd_yhat)';
        Fgrad   = zeros(ndvec,1);
        Fgrad   = Dd_cc;
        Fgrad   = Fgrad - Dd_dd;
    else
        error('Inadmissible value for CRIT.');
    end
    
    if ~isempty(Kmat) && Wlambda > 0
        Fgrad = Fgrad + 2.*Wlambda.*Kmat*dvec;
    end
    
    if ~periodic
        Fgrad(1) = 0;
    else
        Fgrad(1) = 1;
    end
    
end

%  ------------------------------------------------------------------------
%                   Compute hessian wrt dvec if needed
%  ------------------------------------------------------------------------

if nargout > 2
    
    %  second derivative of warping function with respect dvec
    
    D2d_h    = zeros(nfine,ndvec,ndvec);
    %  Problem:  monhess only works with bsplines ... can't use monomials
    D2h      = monhess_Reg(tfine, Wfd);
    D2d_fmax = D2h(nfine,:);
     m = 0;
     for r=1:ndvec
         for s=1:r
             m = m + 1;
             D2h(:,m) = (width./fmax.^3).* ...
                (2.*ffine.*Dd_fmax(r).*Dd_fmax(s) -    ...
                 fmax.*(Dd_ffine(:,r).*Dd_fmax(s)  + ...
                        Dd_ffine(:,s).*Dd_fmax(r)) + ...
                 fmax.^2.*D2h(:,m) - ffine.*fmax.*D2d_fmax(m));
            D2d_h(:,r,s) = D2h(:,m);
            D2d_h(:,s,r) = D2d_h(:,r,s);
        end
    end

    if ~periodic
        D2d_h(1,:) = 0;
        D2d_h(:,1) = 0;
        D2d_h(1,1) = 1;
    end

    %  compute second derivative of Amat wrt dvec
    
    D2d_Amat = zeros(nfine,K,ndvec,ndvec);
    D2_A     = eval_fd(hfine, Afd, 2);
    for k=1:K
        for r=1:ndvec
            for s=1:r
                D2d_Amat(:,k,r,s) = ...
                    D_A(:,k).*squeeze(D2d_h(:,r,s)) + ...
                    D2_A(:,k).*Dd_h(:,r).*Dd_h(:,s);
                D2d_Amat(:,k,s,r) = D2d_Amat(:,k,r,s);
            end
        end
    end
    
    %  compute second derivative of yhat wrt dvec
    
    D2d_fvec = zeros(K,ndvec,ndvec);
    D2d_yhat = zeros(nfine,ndvec,ndvec);
    for r=1:ndvec
        Dd_Amatr    = squeeze(Dd_Amat(:,:,r));
        Dd_Mmatr    = squeeze(Dd_Mmat(:,:,r));
        Dd_Mmatinvr = squeeze(Dd_Mmatinv(:,:,r));
        Facmatr     = Mmatinv*Dd_Mmatr*Mmatinv;
        for s=1:r
            Dd_Amats      = squeeze(Dd_Amat(:,:,s));
            Dd_Mmats      = squeeze(Dd_Mmat(:,:,s));
            Dd_Mmatinvs   = squeeze(Dd_Mmatinv(:,:,s));
            D2d_Amatrs    = squeeze(D2d_Amat(:,:,r,s));
            D2d_Mmatrs    = Amat'*D2d_Amatrs + Dd_Amatr'*Dd_Amats;
            D2d_Mmatrs    = D2d_Mmatrs + D2d_Mmatrs';
            D2d_Mmatinvrs = Mmatinv*Dd_Mmats*Facmatr;
            D2d_Mmatinvrs = D2d_Mmatinvrs + D2d_Mmatinvrs';
            D2d_Mmatinvrs = D2d_Mmatinvrs - Mmatinv*D2d_Mmatrs*Mmatinv;
            D2d_fvecrs = (Mmatinv    *D2d_Amatrs' + ...
                          Dd_Mmatinvr*Dd_Amats'   + ...
                          Dd_Mmatinvs*Dd_Amatr'   + ...
                          D2d_Mmatinvrs*Amat'     )*yfine;
            D2d_fvec(:,r,s) = D2d_fvecrs;
            D2d_yhat(:,r,s) = Amat      *D2d_fvecrs  + ...
                              Dd_Amats  *Dd_fvec(:,r) + ...
                              Dd_Amatr  *Dd_fvec(:,s) + ...
                              D2d_Amatrs*fvec;
            D2d_fvec(:,s,r) = D2d_fvec(:,r,s);
            D2d_yhat(:,s,r) = D2d_yhat(:,r,s);
        end
    end

    %  compute Hessian according to the choice of criterion
        
    if crit == 1
        Fhess = 2.*(Dd_yhat'*Dd_yhat)./nfine;
    else
        Fhess = zeros(ndvec,ndvec);
        m = 0;
        for r=1:ndvec
            for s=1:r
                m = m + 1;
                 D2d_bbrs    =   sum(D2d_yhat(:,r,s).*yfine)/nfine;
                 D2d_ccrs    = 2*sum(D2d_yhat(:,r,s).*yhat + ...
                                     Dd_yhat(:,r).*Dd_yhat(:,s))/nfine;
                 Dd_ggrs     = 4*Dd_bb(r).*Dd_bb(s) + 4*bb*D2d_bbrs + ...
                               Dd_cc(r).*Dd_cc(s) - (aa - cc).*D2d_ccrs;
                 D2d_ddrs    = (Dd_ggrs*dd - gg(r)*Dd_dd(s))/dd^2;
%                  Fhess(r,s)  = sum(D2h(:,m));
%                  Fhess(r,s)  = sum(D2d_h(:,r,s));
%                  Fhess(r,s)  = sum(D2d_Amat(:,1,r,s));
%                  Fhess(r,s)  = sum(D2d_fvec(:,r,s));
%                  Fhess(r,s)  = sum(D2d_yhat(:,r,s));
                 Fhess(r,s)  = D2d_ccrs;
                 Fhess(r,s)  = Fhess(r,s) - D2d_ddrs;
                 Fhess(s,r)  = Fhess(r,s);
            end
        end
    end
    
    if ~isempty(Kmat) && Wlambda > 0
        Fhess = Fhess + 2.*Wlambda.*Kmat;
    end
    
    if ~periodic
        Fhess(1,:) = 0;
        Fhess(:,1) = 0;
        Fhess(1,1) = 1;
    end
    
end

%  ------------------------------------------------------------------------
%  compute the second cross-derivatives and H-update if required
%  ------------------------------------------------------------------------

if nargout > 3
    
    D_Phimat = eval_basis(hfine, Abasis, 1);

    %    compute gradient wrt cvec
    
    nAbasis = getnbasis(Abasis);
    Dc_Amat = zeros(nfine,K,K,nAbasis);
    for p=1:K
        for q=1:nAbasis
            Dc_Amatpq = zeros(nfine,K);
            Dc_Amatpq(:,p)   = Phimat(:,q);
            Dc_Amat(:,:,p,q) = Dc_Amatpq;
        end
    end
    
    %  compute the first c-derivative Dc_yhat wrt cvec
    
    Dc_Mmat    = zeros(K,K,K,nAbasis);
    Dc_Mmatinv = zeros(K,K,K,nAbasis);
    Dc_fvec = zeros(K,K,nAbasis);
    Dc_yhat = zeros(nfine,K*nAbasis);
    m = 0;
    for p=1:K
        for q=1:nAbasis
            m = m + 1;
            temp             = Amat'*Phimat(:,q);
            Dc_Mmatpq        = zeros(K);
            Dc_Mmatpq(p,:)   = temp;
            Dc_Mmatpq        = Dc_Mmatpq + Dc_Mmatpq';
            Dc_Mmat(:,:,p,q) = Dc_Mmatpq;
            Dc_Mmatinvpq     = -Mmatinv*Dc_Mmatpq*Mmatinv;
            Dc_Mmatinv(:,:,p,q) = Dc_Mmatinvpq;
            Dc_Amatpq        = squeeze(Dc_Amat(:,:,p,q));
            Dc_fvec(:,p,q)   = (Dc_Mmatinvpq *Amat' + ...
                               Mmatinv*Dc_Amatpq')*yfine;
            Dc_yhat(:,m)     = Amat*Dc_fvec(:,p,q) + ...
                               Dc_Amatpq*fvec;
        end
    end
    
    %  compute second cross-derivative of Dd_yhat wrt cvec
    
    Fcross = zeros(nWbasis,nAbasis*K);
%     DdDc_Mmat    = zeros(K,K,ndvec,K,nAbasis);
%     DdDc_Mmatinv = zeros(K,K,ndvec,K,nAbasis);
%     DdDc_fvec = zeros(K,ndvec,K,nAbasis);
%     DdDc_yhat = zeros(nfine,ndvec,K,nAbasis);
    for r=1:ndvec
        Dd_yhatr = Dd_yhat(:,r);
        Dd_bb    =   Dd_yhatr'*yfine/nfine;
        Dd_cc    = 2*Dd_yhatr'*yhat /nfine;
        Dd_Amatr = squeeze(Dd_Amat(:,:,r));
        Dd_Mmatr = squeeze(Dd_Mmat(:,:,r));
        Facr     = Mmatinv*Dd_Mmatr*Mmatinv;
        Dd_Mmatinvr = squeeze(Dd_Mmatinv(:,:,r));
        Dd_Phimatr  = D_Phimat.*(Dd_h(:,r)*ones(1,nAbasis));
        m = 0;
        for p=1:K
            for q=1:nAbasis
                m = m + 1;
                Dc_yhatpq    = squeeze(Dc_yhat(:,m));
                Dc_bb        =   Dc_yhatpq'*yfine/nfine;
                Dc_cc        = 2*Dc_yhatpq'*yhat/nfine;
                Dc_Mmatpq    = squeeze(Dc_Mmat(:,:,p,q));
                Dc_Mmatinvpq = squeeze(Dc_Mmatinv(:,:,p,q));
                temp = Amat'*Dd_Phimatr(:,q) + Dd_Amatr'*Phimat(:,q);
                DdDc_Mmatrpq = zeros(K);
                DdDc_Mmatrpq(p,:) = temp;
                DdDc_Mmatrpq = DdDc_Mmatrpq + DdDc_Mmatrpq';
%                 DdDc_Mmat(:,:,r,p,q) = DdDc_Mmatrpq;
                DdDc_Mmatinvrpq = Facr*Dc_Mmatpq*Mmatinv;
                DdDc_Mmatinvrpq = DdDc_Mmatinvrpq + DdDc_Mmatinvrpq';
                DdDc_Mmatinvrpq = DdDc_Mmatinvrpq - ...
                    Mmatinv*DdDc_Mmatrpq*Mmatinv;
%                 DdDc_Mmatinv(:,:,r,p,q) = DdDc_Mmatinvrpq;
                DdDc_fvecrpq = (DdDc_Mmatinvrpq  *Amat'            + ...
                                Dd_Mmatinvr(p,:)'*Phimat(:,q)'     + ...
                                Dc_Mmatinvpq     *Dd_Amatr'        + ...
                                Mmatinv(p,:)'    *Dd_Phimatr(:,q)')*yfine;
%                 DdDc_fvec(:,r,p,q) = DdDc_fvecrpq;
                DdDc_yhatrpq = Amat        *DdDc_fvecrpq   + ...
                               Dd_fvec(p,r)*Phimat(:,q)    + ...
                               Dd_Amatr    *Dc_fvec(:,p,q) + ...
                               fvec(p)     *Dd_Phimatr(:,q);
%                 DdDc_yhat(:,r,p,q) = DdDc_yhatrpq;
                DdDc_bb =    DdDc_yhatrpq'*yfine/nfine;
                DdDc_cc = 2*(DdDc_yhatrpq'*yhat + ...
                             Dd_yhatr'*Dc_yhatpq)/nfine;
                Dc_dd   = (4*bb*Dc_bb - (aa-cc)*Dc_cc)/dd;
                DdDc_dd = ((4*Dc_bb*Dd_bb + 4*bb*DdDc_bb + ...
                            Dc_cc*Dd_cc - (aa-cc)*DdDc_cc)*dd - ...
                           (4*bb*Dd_bb - (aa-cc)*Dd_cc)*Dc_dd)/dd^2;
                Fcross(r,m) = DdDc_cc;
                Fcross(r,m) = Fcross(r,m) - DdDc_dd;
            end
        end
    end
    
    if ~periodic
        Fcross(1,:) = 0;
    else
        Fcross(1,:) = 1;
    end
    
    %  compute the update of H and Dc_H
    
    reshat  = yfine - yhat;
    Hi      = trapz(tfine,reshat.^2);
    tempmat = Dc_yhat - (Fhess\Dd_yhat')'*Fcross;
    Dc_Hi   = -2*trapz(tfine,tempmat.*(reshat*ones(1,nAbasis*K)))';
    
end

