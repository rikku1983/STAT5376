function [Fval, Fgrad, Fhess, Fcross, fvec, yhat, Hi, Dc_Hi, Wmat] = ...
                   PCA_Reg_trapz(dvec, tfine, yfine, AfdPar, WfdPar, ...
                           PenStruct, periodic, crit)
%  Arguments:
%  DVEC    ... A vector of coefficients defining the function W_i(t) that
%              in turn defines the time-warping function h_i(t).
%  TFINE   ... Column vector of NFINE times of observation for the data 
%              values
%  YFINE   ... Column vector of NFINE data values corresponding to
%              observation tines in TFINE
%  AFDPAR  ... A functional parameter object for the K principal component 
%              functions.  The coefficient mat_rix is the transpose of
%              coefficient mat_rix CMAT defining these functions
%  WFDPAR  ... A functional parameter object for the function W_i(t) 
%              defining the initial warping functions h_i(t).
%  PENSTRUCT ... A struct object containing roughness penalty mat_rices
%              and their associated parameters
%  PERIODIC... 1 implies the data are periodic and a shift is estimated.
%              In this case the initial coefficient in DVEC is the
%              time-shift
%              0 implies data are not periodic.  In this case the initial
%              coefficient in DVEC is not used.
%  CRIT    ... 1, 2, or 3 defining the fitting criterion. 
%              2 implies the minimum eigenvalue of the cross-product
%              mat_rix
%
%  Last modified 25 July 2012 by Jim

if nargin < 8, crit = 2;     end
if nargin < 7, periodic = 0; end
if nargin < 6 || isempty(PenStruct)
    FRowMat = [];
    FRowPar = 0;
else
    FRowMat = PenStruct.FRowMat;
    FRowPar = PenStruct.FRowPar;
    %  check dimensions of roughness penalty mat_rices
    if ~isempty(FRowMat) && any(size(FRowMat) ~= K)  
        error('FRowMat is not of order K');  
    end
end

Afd     = getfd(AfdPar);
Abasis  = getbasis(Afd);
Cmat    = getcoef(Afd);

Wfd     = getfd(WfdPar);
Wbasis  = getbasis(Wfd);
Kmat    = getpenmat(WfdPar);
Wlambda = getlambda(WfdPar);

trng     = getbasisrange(Wbasis);
nfine    = length(tfine);
tlo      = trng(1);
thi      = trng(2);
width    = thi - tlo;

ndvec   = length(dvec);
K       = size(getcoef(Afd),2);

%  W function for input coefficients

nWbasis = getnbasis(Wbasis);
Wfd     = fd(dvec, Wbasis);

%  Warping function values defined by these coefficients

ffine        = monfn(tfine, Wfd);
fmax         = ffine(nfine);
hfine        = tlo + width.*ffine./fmax;
hfine(1)     = tlo;
hfine(nfine) = thi;

%  evaluate A_of_hfine ... the principal component values 

%  column penalties are not possible because PC_Reg is
%  called only with data from a single row.

Phi_of_hfine = eval_basis(hfine, Abasis);
A_of_hfine   = Phi_of_hfine*Cmat;

%  compute fvec ... the K factor scores for this case

nAbasis = getnbasis(Abasis);
Wmat = zeros(nAbasis,nAbasis);
for k=1:nAbasis
    for l=1:k
        Wmat(k,l) = trapz(tfine,Phi_of_hfine(:,k).*Phi_of_hfine(:,l));
        Wmat(l,k) = Wmat(k,l);
    end
end
CtWCmat = Cmat'*Wmat*Cmat;
Mmat    = CtWCmat;
if ~isempty(FRowMat) && FRowPar > 0
    Mmat = Mmat + FRowPar.*FRowMat;
end
Mmatinv = inv(Mmat);
ndvec = length(dvec); 
onecoef = ones(1,ndvec);
oneK    = ones(1,K);
fvec    = Mmatinv*trapz(tfine,A_of_hfine.*(yfine*oneK))';

%  compute predicted values 

Fmat = ones(nfine,1)*fvec';
yhat = sum(Fmat.*A_of_hfine,2);

%  compute criterion value according to the choice of criterion

aa = trapz(tfine,yfine.^2);
bb = trapz(tfine,yfine.*yhat);
cc = trapz(tfine,yhat.^2);

res = yfine - yhat;
if crit == 1
    %  least squares criterion
    Fval = aa - 2*bb + cc;
elseif crit == 2
    % minimum eigenvalue criterion
    ff   = aa - cc;
    dd   = sqrt(ff^2 + 4*bb^2);
%     Fval = sum(ffine);
%     Fval = sum(hfine);
%     Fval = sum(A_of_hfine(:,1));
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
    Dd_hfine = width.*(fmax.*Dd_ffine - ffine*Dd_fmax)./fmax^2;
    if periodic
        Dd_hfine(:,1) = 1;
    else
        Dd_hfine(:,1) = 0;
    end
    
    %  compute derivative of A_of_hfine wrt dvec
    
    Dd_A_of_hfine = zeros(nfine,K,ndvec);
    DA_of_hfine   = eval_fd(hfine, Afd, 1);
    for k=1:K
        Dd_A_of_hfine(:,k,:) = (DA_of_hfine(:,k)*onecoef).*Dd_hfine;
    end
    
    %  compute derivative of Mmat wrt dvec
    %  compute derivative of fvec and yhat wrt dvec
    
    Dd_Mmat    = zeros(K,K,ndvec);
    Dd_Mmatinv = Dd_Mmat;
    Dd_fvec    = zeros(K,ndvec);
    Dd_yhat    = zeros(nfine,ndvec);
    for r=1:ndvec
        Dd_A_of_hfine_r = squeeze(Dd_A_of_hfine(:,:,r));
        Dd_Mmat_r = zeros(K,K);
        for k1=1:K
            for k2=1:K
                Dd_Mmat_r(k1,k2) = ...
                    trapz(tfine,A_of_hfine(:,k1).*Dd_A_of_hfine_r(:,k2) + ...
                                Dd_A_of_hfine_r(:,k1).*A_of_hfine(:,k2));
            end
        end
        Dd_Mmatinv_r      = -Mmatinv*Dd_Mmat_r*Mmatinv;
        Dd_Mmat(:,:,r)    = Dd_Mmat_r;
        Dd_Mmatinv(:,:,r) = Dd_Mmatinv_r;
        tempmat = Dd_Mmatinv_r*A_of_hfine' + Mmatinv*Dd_A_of_hfine_r';
        Dd_fvec(:,r) = trapz(tfine,tempmat'.*(yfine*oneK))';
        Dd_yhat(:,r) = A_of_hfine*Dd_fvec(:,r) + Dd_A_of_hfine_r*fvec;
    end
    
    %  compute gradient according to the choice of criterion
    
    if crit == 1
        %  least squares criterion
        Fgrad = -2.*trapz(tfine,Dd_yhat.*(res*onecoef))';
    elseif crit == 2
        % minimum eigenvalue criterion
        Dd_bb   =    trapz(tfine,Dd_yhat.*(yfine*onecoef))';
        Dd_cc   = 2.*trapz(tfine,Dd_yhat.*(yhat *onecoef))';
        gg      = 4.*bb.*Dd_bb - (aa-cc).*Dd_cc;
        Dd_dd   = gg./dd;
%         Fgrad   = sum(Dd_ffine)';
%         Fgrad   = sum(Dd_hfine)';
%         Fgrad   = sum(squeeze(Dd_A_of_hfine(:,1,:)))';
%         Fgrad   = sum(Dd_fvec)';
%         Fgrad   = sum(Dd_yhat)';
        Fgrad   = Dd_cc - Dd_dd;
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
    D2hfine  = monhess_Reg(tfine, Wfd);
    D2d_fmax = D2hfine(nfine,:);
     m = 0;
     for r=1:ndvec
         for s=1:r
             m = m + 1;
             D2hfine(:,m) = (width./fmax.^3).* ...
                (2.*ffine.*Dd_fmax(r).*Dd_fmax(s) -    ...
                 fmax.*(Dd_ffine(:,r).*Dd_fmax(s)  + ...
                        Dd_ffine(:,s).*Dd_fmax(r)) + ...
                 fmax.^2.*D2hfine(:,m) - ffine.*fmax.*D2d_fmax(m));
            D2d_h(:,r,s) = D2hfine(:,m);
            D2d_h(:,s,r) = D2d_h(:,r,s);
        end
    end

    if ~periodic
        D2d_h(1,:) = 0;
        D2d_h(:,1) = 0;
        D2d_h(1,1) = 1;
    end

    %  compute second derivative of A_of_hfine wrt dvec
    
    D2d_A_of_hfine = zeros(nfine,K,ndvec,ndvec);
    D2A_of_hfine   = eval_fd(hfine, Afd, 2);
    for k=1:K
        for r=1:ndvec
            for s=1:r
                D2d_A_of_hfine(:,k,r,s) = ...
                    DA_of_hfine(:,k).*squeeze(D2d_h(:,r,s)) + ...
                    D2A_of_hfine(:,k).*Dd_hfine(:,r).*Dd_hfine(:,s);
                D2d_A_of_hfine(:,k,s,r) = D2d_A_of_hfine(:,k,r,s);
            end
        end
    end
    
    %  compute second derivative of yhat wrt dvec
    
    D2d_fvec = zeros(K,ndvec,ndvec);
    D2d_yhat = zeros(nfine,ndvec,ndvec);
    for r=1:ndvec
        Dd_A_of_hfine_r = squeeze(Dd_A_of_hfine(:,:,r));
        Dd_Mmat_r       = squeeze(Dd_Mmat(:,:,r));
        Dd_Mmatinv_r    = squeeze(Dd_Mmatinv(:,:,r));
        Facmat_r = Mmatinv*Dd_Mmat_r*Mmatinv;
        for s=1:r
            Dd_A_of_hfine_s   = squeeze(Dd_A_of_hfine(:,:,s));
            Dd_Mmat_s         = squeeze(Dd_Mmat(:,:,s));
            Dd_Mmatinv_s      = squeeze(Dd_Mmatinv(:,:,s));
            D2d_A_of_hfine_rs = squeeze(D2d_A_of_hfine(:,:,r,s));
            D2d_Mmat_rs = zeros(K,K);
            for k1=1:K
                for k2=1:k2
                    tempvec = ...
                        A_of_hfine(:,k1).*D2d_A_of_hfine_rs(:,k2) + ...
                        Dd_A_of_hfine_r(:,k1).*Dd_A_of_hfine_s(:,k2);
                    D2d_Mmat_rs(k1,k2) = trapz(tfine,tempvec);
                    D2d_Mmat_rs(k2,k1) = D2d_Mmat_rs(k1,k2);
                end
            end 
            D2d_Mmat_rs    = D2d_Mmat_rs + D2d_Mmat_rs';
            D2d_Mmatinv_rs = Mmatinv*Dd_Mmat_s*Facmat_r;
            D2d_Mmatinv_rs = D2d_Mmatinv_rs + D2d_Mmatinv_rs';
            D2d_Mmatinv_rs = D2d_Mmatinv_rs - Mmatinv*D2d_Mmat_rs*Mmatinv;
            tempmat = Mmatinv*D2d_A_of_hfine_rs' + ...
                      Dd_Mmatinv_r*Dd_A_of_hfine_s'   + ...
                      Dd_Mmatinv_s*Dd_A_of_hfine_r'   + ...
                      D2d_Mmatinv_rs*A_of_hfine';
            D2d_fvecrs = trapz(tfine, tempmat'.*(yfine*oneK))';
            D2d_fvec(:,r,s) = D2d_fvecrs;
            D2d_yhat(:,r,s) = A_of_hfine       *D2d_fvecrs  + ...
                              Dd_A_of_hfine_s  *Dd_fvec(:,r) + ...
                              Dd_A_of_hfine_r  *Dd_fvec(:,s) + ...
                              D2d_A_of_hfine_rs*fvec;
            D2d_fvec(:,s,r) = D2d_fvec(:,r,s);
            D2d_yhat(:,s,r) = D2d_yhat(:,r,s);
        end
    end

    %  compute Hessian according to the choice of criterion
        
    Fhess = zeros(ndvec,ndvec);
    if crit == 1
        for l1=1:ndvec
            for l2=1:l1                
                Fhess(l1,l2) = 2.*trapz(tfine,Dd_yhat(:,l1).*Dd_yhat(:,l2));
                Fhess(l2,l1) = Fhess(l1,l2); 
            end
        end
    else
        m = 0;
        for r=1:ndvec
            for s=1:r
                m = m + 1;
                 D2d_bbrs    =   trapz(tfine,D2d_yhat(:,r,s).*yfine);
                 D2d_ccrs    = 2*trapz(tfine,D2d_yhat(:,r,s).*yhat + ...
                                             Dd_yhat(:,r).*Dd_yhat(:,s));
                 Dd_ggrs     = 4*Dd_bb(r).*Dd_bb(s) + 4*bb*D2d_bbrs + ...
                               Dd_cc(r).*Dd_cc(s) - (aa - cc).*D2d_ccrs;
                 D2d_ddrs    = (Dd_ggrs*dd - gg(r)*Dd_dd(s))/dd^2;
%                  Fhess(r,s)  = sum(D2hfine(:,m));
%                  Fhess(r,s)  = sum(D2d_h(:,r,s));
%                  Fhess(r,s)  = sum(D2d_A_of_hfine(:,1,r,s));
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
    
    D_Phi_of_hfine = eval_basis(hfine, Abasis, 1);

    %    compute gradient wrt cvec
    
    nAbasis = getnbasis(Abasis);
    Dc_A_of_hfine = zeros(nfine,K,K,nAbasis);
    for p=1:K
        for q=1:nAbasis
            Dc_A_of_hfinepq = zeros(nfine,K);
            Dc_A_of_hfinepq(:,p)   = Phi_of_hfine(:,q);
            Dc_A_of_hfine(:,:,p,q) = Dc_A_of_hfinepq;
        end
    end
    
    %  compute the first c-derivative Dc_yhat wrt cvec
    
    Dc_Mmat    = zeros(K,K,K,nAbasis);
    Dc_Mmatinv = zeros(K,K,K,nAbasis);
    Dc_fvec    = zeros(K,K,nAbasis);
    Dc_yhat    = zeros(nfine,K*nAbasis);
    m = 0;
    for p=1:K
        for q=1:nAbasis
            m = m + 1;
            temp = trapz(tfine,A_of_hfine.*(Phi_of_hfine(:,q)*oneK))';
            Dc_Mmatpq        = zeros(K);
            Dc_Mmatpq(p,:)   = temp;
            Dc_Mmatpq        = Dc_Mmatpq + Dc_Mmatpq';
            Dc_Mmat(:,:,p,q) = Dc_Mmatpq;
            Dc_Mmatinvpq     = -Mmatinv*Dc_Mmatpq*Mmatinv;
            Dc_Mmatinv(:,:,p,q) = Dc_Mmatinvpq;
            Dc_A_of_hfinepq        = squeeze(Dc_A_of_hfine(:,:,p,q));
            tempmat = Dc_Mmatinvpq*A_of_hfine' + Mmatinv*Dc_A_of_hfinepq';
            Dc_fvec(:,p,q)   = trapz(tfine,tempmat'.*(yfine*oneK))';
            Dc_yhat(:,m)     = A_of_hfine*Dc_fvec(:,p,q) + ...
                               Dc_A_of_hfinepq*fvec;
        end
    end
    
    %  compute second cross-derivative of Dd_yhat wrt cvec
    
    Fcross = zeros(nWbasis,nAbasis*K);
%     DdDc_Mmat    = zeros(K,K,ndvec,K,nAbasis);
%     DdDc_Mmatinv = zeros(K,K,ndvec,K,nAbasis);
%     DdDc_fvec = zeros(K,ndvec,K,nAbasis);
%     DdDc_yhat = zeros(nfine,ndvec,K,nAbasis);
    for r=1:ndvec
        Dd_yhat_r = Dd_yhat(:,r);
        Dd_bb     =  trapz(tfine,Dd_yhat_r.*yfine);
        Dd_cc    = 2*trapz(tfine,Dd_yhat_r.*yhat);
        Dd_A_of_hfine_r = squeeze(Dd_A_of_hfine(:,:,r));
        Dd_Mmat_r = squeeze(Dd_Mmat(:,:,r));
        Facr     = Mmatinv*Dd_Mmat_r*Mmatinv;
        Dd_Mmatinv_r = squeeze(Dd_Mmatinv(:,:,r));
        Dd_Phi_of_hfine_r  = D_Phi_of_hfine.*(Dd_hfine(:,r)*ones(1,nAbasis));
        m = 0;
        for p=1:K
            for q=1:nAbasis
                m = m + 1;
                Dc_yhatpq    = squeeze(Dc_yhat(:,m));
                Dc_bb        =   trapz(tfine,Dc_yhatpq.*yfine);
                Dc_cc        = 2*trapz(tfine,Dc_yhatpq.*yhat);
                Dc_Mmatpq    = squeeze(Dc_Mmat(:,:,p,q));
                Dc_Mmatinvpq = squeeze(Dc_Mmatinv(:,:,p,q));
                tempmat = A_of_hfine.*...
                         (Dd_Phi_of_hfine_r(:,q)*oneK) + ...
                          Dd_A_of_hfine_r.*(Phi_of_hfine(:,q)*oneK);
                DdDc_Mmat_rpq = zeros(K);
                DdDc_Mmat_rpq(p,:) = trapz(tfine,tempmat);
                DdDc_Mmat_rpq = DdDc_Mmat_rpq + DdDc_Mmat_rpq';
%                 DdDc_Mmat(:,:,r,p,q) = DdDc_Mmat_rpq;
                DdDc_Mmatinv_rpq = Facr*Dc_Mmatpq*Mmatinv;
                DdDc_Mmatinv_rpq = DdDc_Mmatinv_rpq + DdDc_Mmatinv_rpq';
                DdDc_Mmatinv_rpq = DdDc_Mmatinv_rpq - ...
                    Mmatinv*DdDc_Mmat_rpq*Mmatinv;
%                 DdDc_Mmatinv(:,:,r,p,q) = DdDc_Mmatinv_rpq;
                tempmat = DdDc_Mmatinv_rpq  *A_of_hfine'            + ...
                          Dd_Mmatinv_r(p,:)'*Phi_of_hfine(:,q)'     + ...
                          Dc_Mmatinvpq     *Dd_A_of_hfine_r'        + ...
                          Mmatinv(p,:)'    *Dd_Phi_of_hfine_r(:,q)';
                DdDc_fvec_rpq = trapz(tfine,tempmat'.*(yfine*oneK))';
%                 DdDc_fvec(:,r,p,q) = DdDc_fvec_rpq;
                DdDc_yhat_rpq = A_of_hfine        *DdDc_fvec_rpq   + ...
                                Dd_fvec(p,r)*Phi_of_hfine(:,q)     + ...
                                Dd_A_of_hfine_r    *Dc_fvec(:,p,q) + ...
                                fvec(p)     *Dd_Phi_of_hfine_r(:,q);
%                 DdDc_yhat(:,r,p,q) = DdDc_yhat_rpq;
                DdDc_bb =   trapz(tfine,DdDc_yhat_rpq.*yfine);
                DdDc_cc = 2*trapz(tfine,(DdDc_yhat_rpq.*yhat + ...
                                         Dd_yhat_r.*Dc_yhatpq));
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
    
    Hi      = trapz(tfine,res.^2);
    tempmat = Dc_yhat - (Fhess\Dd_yhat')'*Fcross;
    Dc_Hi   = -2*trapz(tfine, tempmat.*(res*ones(1,nAbasis*K)))';
    
end

