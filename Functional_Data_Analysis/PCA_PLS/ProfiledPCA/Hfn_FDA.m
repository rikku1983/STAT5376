function [Hval, Hgrad, Fmat, Yhat, Bmat] = ...
    Hfn_FDA(parvec, Xfd, AfdPar, Ymat, gamval, PenStruct, ...
            radius, anglewrd, GrssMnWrd)
% function [Hval, Hgrad, Fmat] = Hfn_FDA(cvec, Xfd, Ymat, AfdPar, ...
%                                        gamval, PenStruct)
%  HFN_FDA evaluates the loss function and its gradient for a parameter
%  cascading approach to principal components analysis of Xfd which is
%  either a functional data object or a data matrix.
%  K factors are used.
%
%  Inner loss function for the functional data object is 
%
%  J(F|A,X) = (1-\gamval)*||Xfd - F*Afd||^2 + 
%                \gamval *||Ymat'[I-F(A)(F'(A)F(A))^{-1}F'(A)]||^2 +
%              FRowPar ||F  FRowMat F'||^2 + FColPar ||F' FColMat F ||^2
%
%  Outer loss function in the is
%
%  H(A|X) = (1-\gamval)*||Xfd - F(Afd)*Afd||^2 + 
%              \gamval *||Ymat'[I-F(A)(F'(A)F(A))^{-1}F'(A)]||^2 +
%                       ARowPar ||C' RowMat C||^2
%  where K-vector-valued function a has the expansion C \phi where 
%  \phi is a vector of L basis functions.
%
%  In this version, parvec has length one less than avec, and contains
%  the angles in the map of a vector into surface of a hypersphere.
%
%  Arguments:
%
%  PARVEC    ... Vectorized K by L coefficient matrix for principal 
%                component functions in Afd converted to angular measure
%  XFD       ... Either an N-vector of functional data objects or
%                a N by n matrix of multivariate observations
%  AFDPAR    ... A functional parameter object defining the basis system
%                and smoothing of the coefficient functions.  
%                Cannot be empty if XFD is a functional data object.
%                Must be empty if XFD is a matrix.
%  YMAT      ... An N by p matrix of dependent variable values.  Not used
%                if GAMVAL == 0.
%  GAMVAL    ... weight in [0,1] placed on fitting Ymat
%  PENSTRUCT ... A struct object with a field for each of four roughness
%                penalties and each corresponding roughness matrix.
%                If no penalty is to be applied, the parameter value
%                should be 0 and the matrix value [].
%  RADIUS    ... Radius for angular transformation.         Default 1.
%  ANGLEWRD  ... If nonzero, apply angular transformation.  Default 1.
%  GRSSMNWRD ... If nonzero, apply Grassmann map.           Default 1.
%
%  An alternative argument list for pure PCA:
%
%  PARVEC    ... Vectorized K by L coefficient matrix for principal 
%                component functions in Afd.
%  XFD       ... Either an N-vector of functional data objects or
%                a N by n matrix of multivariate observations
%  AFDPAR    ... A functional parameter object defining the basis system
%                and smoothing of the coefficient functions.  Not used if
%                XFD is a matrix.
%  PENSTRUCT ... A struct object with a field for each of four roughness
%                penalties and each corresponding roughness matrix.
%                If no penalty is to be applied, the parameter value
%                should be 0 and the matrix value [].
%  RADIUS    ... Radius for angular transformation
%  ANGLEWRD  ... If nonzero, apply angular transformation
%  GRSSMNWRD ... If nonzero, apply Grassmann map

%  Last modified 17 June 2014

%  Set default values for penalized least squares estimation of PC scores

if isstruct(Ymat)
    %  pure PCA arguments, no Ymat and gamval input
    if nargin < 6
        anglewrd = 1;   
    else
        anglewrd = PenStruct;
    end 
    if nargin < 5 
        radius = 1;   
    else
        radius = gamval;
    end 
    if nargin < 7 
        GrssMnWrd = 1;
    else
        GrssMnWrd = radius;
    end
    PenStruct = Ymat;  
    gamval    = 0;
    Ymat      = [];
else
    %  PCA/PLS arguments
    if nargin <  9, GrssMnWrd = 1;   end
    if nargin <  8, anglewrd  = 1;   end
    if nargin <  7, radius    = 1;   end
    if nargin <  6, PenStruct = [];  end
    if nargin <  5, gamval    = 0;   end
    if nargin <  4, Ymat      = [];  end
end

if gamval < 0 || gamval > 1
    error('GAMMAVAL is not within [0,1[.');
end

%  declare persistent variables

persistent GrssMnMap Umat Wmat XXmat

%  ------------------------------------------------------------------------
%            Xfd is a fd and AfdPar is not empty
%  ------------------------------------------------------------------------

if isa_fd(Xfd) && ~isempty(AfdPar)
    
    %  -----------------------------------------------------
    %              Xfd is a functional data object
    %  -----------------------------------------------------
    
    %  update AfdPar by replacing the coefficient matrix for Afd
    
    AfdPar       = fdParcheck(AfdPar);
    Afd          = getfd(AfdPar);
    Acoef        = getcoef(Afd);
    [nAbasis, K] = size(Acoef);
    nK           = nAbasis*K;
    
    %  Compute Grassman map, which in this case is 
    %    nAbasis*K by (nAbasis-K)*K
    
    if isempty(GrssMnMap)        
        GrssMnPts = linspace(0,1,nAbasis*K)';
        GrssBasis = create_fourier_basis([0,1],(nAbasis-K)*K);
        GrssMnMap = eval_basis(GrssMnPts,GrssBasis);
        GrssMnMap = GrssMnMap(:,1:(nAbasis-K)*K);
        [GrssMnMap,S,V] = svd(GrssMnMap,0);
    end
    
    %  convert PARVEC to CVEC, and CVEC to K by n matrix AMAT
    
    if anglewrd
        if GrssMnWrd
            [Gvec,DGvec] = sphere2cartes(parvec,radius);
            cvec  = GrssMnMap*Gvec;
            Dcvec = GrssMnMap*DGvec;
        else
            [cvec,Dcvec] = sphere2cartes(parvec,radius);
        end
    else
        if GrssMnWrd
            cvec  = GrssMnMap*parvec;
            Dcvec = GrssMnMap;
        else
            cvec  = parvec;
            Dcvec = eye(nK);
        end
    end
    
    %  set up Cmat
    
    Cmat = reshape(cvec, nAbasis, K);
    
    % define Afd and AfdPar with this value of CMAT
    
    Afd    = putcoef(Afd,Cmat);
    AfdPar = putfd(AfdPar, Afd);
    Abasis = getbasis(Afd);
    
    %  set up ARow penalty
    
    ARowPar      = getlambda(AfdPar);
    ARowMat      = getpenmat(AfdPar);
    if ARowPar > 0 && isempty(ARowMat)
        ALfd    = getLfd(AfdPar);
        ARowMat = eval_penalty(Abasis, ALfd);
    end
    
    %  get dimensions of data and convert AVEC to K by n matrix AMAT
    
    Xcoef  = getcoef(Xfd);
    Xbasis = getbasis(Xfd);
    N      = size(Xcoef,2);
    
    %  Retrieve roughness penalty parameters and matrices
    
    if nargin < 3,  PenStruct = [];  end
    
    if isempty(PenStruct)
        FRowMat = [];
        FRowPar = 0;
        FColMat = [];
        FColPar = 0;
        FVecMat = [];
        FVecPar = 0;
    else
        if isfield(PenStruct, 'FRowMat')
            FRowMat = PenStruct.FRowMat;
        else
            FRowMat = [];
        end
        if isfield(PenStruct, 'FRowPar')
            FRowPar = PenStruct.FRowPar;
        else
            FRowPar = 0;
        end
        if isfield(PenStruct, 'FColMat')
            FColMat = PenStruct.FColMat;
        else
            FColMat = [];
        end
        if isfield(PenStruct, 'FColPar')
            FColPar = PenStruct.FColPar;
        else
            FColPar = 0;
        end
        if isfield(PenStruct, 'FVecMat')
            FVecMat = PenStruct.FVecMat;
        else
            FVecMat = [];
        end
        if isfield(PenStruct, 'FVecPar')
            FVecPar = PenStruct.FVecPar;
        else
            FVecPar = 0;
        end        
    end
    
    if ~isempty(ARowMat) && any(size(ARowMat) ~= nAbasis)
        error('ARowMat is not of order nAbasis');
    end
    
    %  Compute persistent inner product matrices
    
    if isempty(Umat)
        Umat    = inprod_basis(Xbasis,Abasis);
        Vmat    = eval_penalty(Xbasis,0);
        Wmat    = eval_penalty(Abasis,0);
        XXmat   = Xcoef'*Vmat*Xcoef;
    end
    
    %  set up matrices depedent on Cmat
    
    XAmat   = Xcoef'*Umat*Cmat;
    CtWmat  = Cmat'*Wmat;
    CtWCmat = CtWmat*Cmat;
    
    %  Compute N by K matrix FMAT containing smoothed PC scores
    
    if FColPar > 0 && ~isempty(FColMat)
        % N*K equations in vec(F)
        I_N = eye(N);
        I_K = eye(K);
        Mmat = kron(I_N,CtWCmat);
        if FRowPar > 0 && ~isempty(FRowMat)
            Mmat = Mmat + kron(I_N,FRowPar*FRowMat);
        end
        if FColPar > 0 && ~isempty(FColMat)
            Mmat = Mmat + kron(FColPar*FColMat,I_K);
        end
        AXtvec = reshape(XAmat',N*K,1);
        Mmatinv = inv(Mmat);
        Fvec   = AXtvec'*Mmatinv;
        Fmat   = reshape(Fvec,K,N)';
        %     Fvec   = (Mmat\AXtvec)';
        %     Fmat   = reshape(Fvec,K,N)';
    else
        %  K equations in Fmat'
        Mmat = CtWCmat;
        if FRowPar > 0 && ~isempty(FRowMat)
            Mmat = Mmat + FRowPar*FRowMat;
        end
        Mmatinv = inv(Mmat);
        Fmat    = XAmat*Mmatinv;
        %     Fmat    = (Mmat\XAmat')';
    end
    
    Xhatcoef = Cmat*Fmat';
    
    %  Evaluate the outer loss function H(A) if gamval < 1
    
    if gamval < 1
        XAterm = XAmat*Fmat';
        Hval   = (1-gamval)*trace(XXmat + Xhatcoef'*Wmat*Xhatcoef - ...
                 (XAterm + XAterm'));
    else
        Hval = 0;
    end
    
    %  evaluate PLS fit to Ymat if it is supplied
    
    if ~isempty(Ymat)
        Bmat = Fmat\Ymat;
        Yhat = Fmat*Bmat;
    else
        Yhat = [];
        Bmat = [];
    end
    if gamval > 0 && ~isempty(Ymat)
        Emat = Ymat - Yhat;
        Hval = Hval + gamval*trace(Emat'*Emat);
    end
    
    %  add penalty to Hval if required
    
    if ARowPar > 0 && ~isempty(ARowMat)
        Hval = Hval + ARowPar*trace(Cmat'*ARowMat*Cmat);
        diPenRow = Cmat'*ARowMat;
    end
    
    %  Evaluate the gradient if required
    
    if nargout > 1
        CtWmat   = Cmat'*Wmat;
        CtWCmat  = CtWmat*Cmat;
        Hgrad    = zeros(nAbasis*K,1);
        %  loop through the gradient elements
        m = 0;
        for p=1:K
            for q=1:nAbasis
                m = m + 1;
                %  Evaluate the partial derivative of F(A) wrt Cqp
                if FColPar > 0 && ~isempty(FColMat)
                    %  N*K equations in vec(F)
                    %  Evaluate the derivative of A*A' wrt Cqp
                    diAAtdiCqp      = zeros(K);
                    diAAtdiCqp(:,p) = CtWmat(:,q);
                    diAAtdiCqp(p,:) = diAAtdiCqp(p,:) + CtWmat(:,q)';
                    %  Evaluate the derivative of M wrt Cqp
                    diMdiCqp        = kron(I_N,diAAtdiCqp);
                    %  Evaluate the derivative of A*X' wrt Cqp
                    diAXtdiCqp      = zeros(K,N);
                    diAXtdiCqp(p,:) = (Xcoef'*Umat(:,q))';
                    diAXtdiCqp      = reshape(diAXtdiCqp,1,N*K);
                    %  Evaluate the derivative of F wrt Cqp
                    diFdiCqp = (diAXtdiCqp - Fvec*diMdiCqp)*Mmatinv;
                    % diFdiCqp = (Mmat\(diAXtdiCqp - Fvec*diMdiCqp))';
                    diFdiCqp = reshape(diFdiCqp,K,N)';
                else
                    %  K equations in Fmat'
                    %  Evaluate the derivative of Mmat wrt Cqp
                    diMdiCqp      = zeros(K);
                    diMdiCqp(:,p) = CtWmat(:,q);
                    diMdiCqp(p,:) = diMdiCqp(p,:) + CtWmat(:,q)';
                    %  Evaluate the derivative of M(A)^{-1} wrt Cqp
                    %                 diMinvdiCqp   = -Mmat\diMdiCqp;
                    diMinvdiCqp   = -Mmatinv*diMdiCqp;
                    %  Evaluate the derivative of Fmat wrt Cmat(p,q)
                    diXAmatdiCqp = zeros(N,K);
                    diXAmatdiCqp(:,p) = Xcoef'*Umat(:,q);
                    diFdiCqp = (XAmat*diMinvdiCqp + diXAmatdiCqp)*Mmatinv;
                    % diFdiCqp = 
                    %       (Mmat\(XAmat*diMinvdiCqp + diXAmatdiCqp)')';
                end
                %  This has to be made simpler !!
                %  derivative of AA'
                diaatdiCqp = zeros(K,K);
                diaatdiCqp(p,:) = CtWmat(:,q)';
                %  derivative of AY'
                diaytdiCqp = zeros(K,N);
                diaytdiCqp(p,:) = (Umat(:,q)'*Xcoef)';
                %  Evaluate total derivative of yhat wrt Cqp
                Hgrad(m) = -2*(1-gamval)*trace( ...
                    XAmat'*diFdiCqp + ...
                    diaytdiCqp*Fmat         - ...
                    Fmat'*diFdiCqp*CtWCmat  - ...
                    Fmat'*Fmat*diaatdiCqp);
                if gamval > 0 && ~isempty(Ymat)
                    DYhat = diFdiCqp*Bmat;
                    DGmat = Emat'*DYhat;
                    Hgrad(m) = Hgrad(m) - 2*gamval*trace(DGmat);
                end
                if ARowPar > 0 && ~isempty(ARowMat)
                    Hgrad(m) = Hgrad(m) + 2*ARowPar*diPenRow(p,q);
                end
            end
        end
        
        Hgrad = Dcvec'*Hgrad;
        
    end
    
%  ------------------------------------------------------------------------
%            Xfd is a matrix and AfdPar is empty
%  ------------------------------------------------------------------------

elseif ismatrix(Xfd) && isempty(AfdPar)
    
    %  -----------------------------------------------------
    %              Xfd is an N by n matrix
    %  -----------------------------------------------------
    
    %  get dimensions of data
    
    Xmat  = Xfd;
    [N,n] = size(Xmat);

    if GrssMnWrd
        if anglewrd
            J  = length(parvec) + 1;
        else
            J  = length(parvec);
        end
        K  = n*(1-sqrt(1-4*J/n^2))/2;
        nK = n*K;
        if isempty(GrssMnMap)
            GrssMnPts = linspace(0,1,n*K)';
            GrssBasis = create_fourier_basis([0,1],(n-K)*K);
            GrssMnMap = eval_basis(GrssMnPts,GrssBasis);
            GrssMnMap = GrssMnMap(:,1:(n-K)*K);
        end

    else
        if anglewrd
            nK = length(parvec) + 1;
        else
            nK = length(parvec);
        end
        K  = nK/n;
    end
    
    %  convert PARVEC to AVEC,  and AVEC to K by n matrix AMAT
    
    if anglewrd
        if GrssMnWrd
            [Gvec,DGvec] = sphere2cartes(parvec,1);
            avec  = GrssMnMap*Gvec;
            Davec = GrssMnMap*DGvec;
        else
            [avec,Davec] = sphere2cartes(parvec,1);
        end
    else
        if GrssMnWrd
            avec  = GrssMnMap*parvec;
            Davec = GrssMnMap;
        else
            avec  = parvec;
            Davec = eye(nK);
        end
    end
    
    Amat = reshape(avec, K, n);
    
    %  Retrieve roughness penalty parameters and matrices
    
    if nargin < 3,  PenStruct = [];  end
    
    if isempty(PenStruct)
        FRowMat = [];
        FRowPar = 0;
        FColMat = [];
        FColPar = 0;
        FVecMat = [];
        FVecPar = 0;
        ARowMat = [];
        ARowPar = 0;
        AColMat = [];
        AColPar = 0;
        AVecMat = [];
        AVecPar = 0;
    else
        %  ------  F regularization ---------
        %  FRow values
        if isfield(PenStruct,'FRowMat')
            FRowMat = PenStruct.FRowMat;
        else
            FRowMat = [];
        end
        if isfield(PenStruct,'FrowPar')
            FRowPar = PenStruct.FRowPar;
        else
            FRowPar = 0;
        end
        %  FCol values        
        if isfield(PenStruct,'FColMat')
            FColMat = PenStruct.FColMat;
        else
            FColMat = [];
        end
        if isfield(PenStruct,'FColPar')
            FColPar = PenStruct.FColPar;
        else
            FColPar = 0;
        end
        %  FVec values        
        if isfield(PenStruct,'FVecMat')
            FVecMat = PenStruct.FVecMat;
        else
            FVecMat = [];
        end
        if isfield(PenStruct,'FVecPar')
            FVecPar = PenStruct.FVecPar;
        else
            FVecPar = 0;
        end
        %  ------  A regularization ---------
        %  ARow values
        if isfield(PenStruct,'ARowMat')
            ARowMat = PenStruct.ARowMat;
        else
            ARowMat = [];
        end
        if isfield(PenStruct,'ARowPar')
            ARowPar = PenStruct.ARowPar;
        else
            ARowPar = 0;
        end
        %  ACol values
        if isfield(PenStruct,'AColMat')
            AColMat = PenStruct.AColMat;
        else
            AColMat = [];
        end
        if isfield(PenStruct,'AColPar')
            AColPar = PenStruct.AColPar;
        else
            AColPar = 0;
        end
        %  AVec values
        if isfield(PenStruct,'AVecMat')
            AVecMat = PenStruct.AVecMat;
        else
            AVecMat = [];
        end
        if isfield(PenStruct,'AVecPar')
            AVecPar = PenStruct.AVecPar;
        else
            AVecPar = 0;
        end
        %  check dimensions of roughness penalty matrices
        %  F matrices
        if ~isempty(FRowMat) && any(size(FRowMat) ~= K)
            error('FRowMat is not of order K');
        end
        if ~isempty(FColMat) && any(size(FColMat) ~= N)
            error('FColMat is not of order N');
        end
        if ~isempty(FVecMat) && any(size(FVecMat) ~= N*K)
            error('FColMat is not of order N*K');
        end
        %  A matrices
        if ~isempty(ARowMat) && any(size(ARowMat) ~= K)
            error('ARowMat is not of order K');
        end
        if ~isempty(AColMat) && any(size(AColMat) ~= n)
            error('AColMat is not of order n');
        end
        if ~isempty(AVecMat) && any(size(AVecMat) ~= n*K)
            error('AColMat is not of order n*K');
        end
    end
    
    %  Compute N by K matrix FMAT containing smoothed PC scores
    
    if (FColPar > 0 || FVecPar > 0) && ~isempty(FColMat)
        % N*K equations in vec(F)
        I_N = eye(N);
        I_K = eye(K);
        Mmat = kron(I_N,Amat*Amat');
        if FRowPar > 0 && ~isempty(FRowMat)
            Mmat = Mmat + kron(I_N,FRowPar*FRowMat);
        end
        if FColPar > 0 && ~isempty(FColMat)
            Mmat = Mmat + kron(FColPar*FColMat,I_K);
        end
        if FVecPar > 0 && ~isempty(FVecMat)
            Mmat = Mmat + FVecPar*FVecMat;
        end
        Mmatinv = inv(Mmat);
        AXtvec = reshape(Amat*Xmat',N*K,1);
        Fmat   = reshape(AXtvec'*Mmatinv,K,N)';
    else
        %  K equations in Fmat'
        Mmat = Amat*Amat';
        if FRowPar > 0 && ~isempty(FRowMat)
            Mmat = Mmat + FRowPar*FRowMat;
        end
        Mmatinv = inv(Mmat);
        Fmat    = Xmat*Amat'*Mmatinv;
    end
    Xhat = Fmat * Amat;         %  fit to the data
    Hval = sum(sum((Xmat - Xhat).^2));   %  error sum of squares   SSE
    if gamval > 0 && ~isempty(Ymat)
        Bmat = Fmat\Ymat;
        Yhat = Fmat*Bmat;
        Emat = Ymat - Yhat;
        Hval = Hval + gamval*trace(Emat'*Emat);
    end
    if ARowPar > 0 && ~isempty(ARowMat)
        Hval = Hval + ARowPar*trace(Amat*ARowMat*Amat');
    end
    
    %  Evaluate the outer loss function H(A)
    
    if ARowPar > 0 && ~isempty(ARowMat)
        diPenRow = Amat*ARowMat;
    end
    
    %  Evaluate the gradient if required
    
    if nargout > 1
        Hgrad = zeros(nK,1);
        %  loop through the gradient elements
        m = 0;
        for q=1:n
            for p=1:K
                m = m + 1;
                %  Evaluate the partial derivative of F(A) wrt Apq
                if (FColPar > 0 || FVecPar > 0) && ~isempty(FColMat)
                    %  N*K equations in vec(F)
                    %  Evaluate the derivative of A*A' wrt Apq
                    diAAtdiApq      = zeros(K);
                    diAAtdiApq(:,p) = Amat(:,q);
                    diAAtdiApq(p,:) = diAAtdiApq(p,:) + Amat(:,q)';
                    %  Evaluate the derivative of M wrt Apq
                    diMdiApq        = kron(I_N,diAAtdiApq);
                    %  Evaluate the derivative of A*X' wrt Apq
                    diAXtdiApq      = zeros(K,N);
                    diAXtdiApq(p,:) = Xmat(:,q);
                    diAXtdiApq      = reshape(diAXtdiApq,1,N*K);
                    %  Evaluate the derivative of F wrt Apq
                    diFdiApq = ...
                        (diAXtdiApq - AXtvec'*Mmatinv*diMdiApq)*Mmatinv;
                    diFdiApq = reshape(diFdiApq,K,N)';
                else
                    %  K equations in Fmat'
                    %  Evaluate the derivative of Mmat wrt Apq
                    diMdiApq      = zeros(K);
                    diMdiApq(:,p) = Amat(:,q);
                    diMdiApq(p,:) = diMdiApq(p,:) + Amat(:,q)';
                    %  Evaluate the derivative of M(A)^{-1} wrt Apq
                    diMinvdiApq   = -Mmatinv*diMdiApq*Mmatinv;
                    %  Evaluate the derivative of Fmat wrt Apq
                    Dffacpq       = Amat'*diMinvdiApq;
                    Dffacpq(q,:)  = Dffacpq(q,:) + Mmatinv(p,:);
                    diFdiApq      = Xmat*Dffacpq;
                end
                %  Evaluate total derivative of Xhat wrt Apq
                diXhatdiApq = diFdiApq*Amat;
                diXhatdiApq(:,q) = diXhatdiApq(:,q) + Fmat(:,p);
                Emat = Xmat - Xhat;
                Hgrad(m) = -2*(1-gamval)*trace(Emat'*diXhatdiApq);
                if gamval > 0 && ~isempty(Ymat)
                    DYhat = diFdiAqp*Bmat;
                    Emat  = Ymat - Yhat;
                    Hgrad(m) = Hgrad(m) - 2*gamval*trace(Emat'*DYhat);
                end
                if ARowPar > 0 && ~isempty(ARowMat)
                    Hgrad(m) = Hgrad(m) + 2*ARowPar*diPenRow(p,q);
                end
            end
        end
        
        Hgrad = Davec'*Hgrad;
        
    end
    
else
    
    error(['Second argument is neither a functional data object', ...
        ' nor a matrix.']);
    
end


