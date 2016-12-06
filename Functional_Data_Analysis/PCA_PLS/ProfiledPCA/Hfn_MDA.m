function [Hval, Hgrad, fmat] = Hfn_MDA(parvec, ymat, PenStruct)
%  HFN evaluates the loss function and its gradient for a parameter
%  cascading approach to principal components analysis of 
%  multivariate N by n data YMAT using K << n factors.
%  Inner loss function is 
%
%  J(F|A,Y) = ||Y - FA||^2 +  FRowPar ||F  FRowMat F'||^2 + 
%                             FColPar ||F' FColMat F ||^2
%
%  Outer loss function is H(A|Y) = ||Y - F(A)A||^2
%
%  H(A|Y) = ||Y - F(A)A||^2 + ARowPar ||A ARowMat A'||^2
%
%  In this version, parvec has length one less than avec, and contains
%  the angles in the map of a vector in to surface of a hypersphere.
%  Arguments:
%
%  AVEC      ... Vectorized K by n principal component coefficient matrix
%  YMAT      ... N by n matrix of data
%  PENSTRUCT ... A struct object with a field for each of four roughness
%                penalties and each corresponding roughness matrix.
%                If no penalty is to be applied, the parameter value
%                should be 0 and the matrix value [].
%
%  Last modified 22 September 2012

%  set default values for penalized least squares estimation of PC scores

if nargin < 3, PenStruct = [];  end

%  get dimensions of data

[N,n] = size(ymat);
nK    = length(parvec) + 1;
K     = nK/n;

%  convert PARVEC to AVEC,  and AVEC to K by n matrix AMAT

[avec,Davec] = sphere2cartes(parvec,1);

amat  = reshape(avec, K, n);

%  Retrieve roughness penalty parameters and matrices

if nargin < 3,  PenStruct = [];  end

if isempty(PenStruct)
    FRowMat = [];
    FRowPar = 0;
    FColMat = [];
    FColPar = 0;
    ARowMat = [];
    ARowPar = 0;
    ANrmPar = 0;
else
    FRowMat = PenStruct.FRowMat;
    FRowPar = PenStruct.FRowPar;
    FColMat = PenStruct.FColMat;
    FColPar = PenStruct.FColPar;
    ARowMat = PenStruct.ARowMat;
    ARowPar = PenStruct.ARowPar;
    ANrmPar = PenStruct.ANrmPar;
    %  check dimensions of roughness penalty matrices
    if ~isempty(FRowMat) && any(size(FRowMat) ~= K)  
        error('FRowMat is not of order K');  
    end
    if ~isempty(FColMat) && any(size(FColMat) ~= N)  
        error('FColMat is not of order N');  
    end
    if ~isempty(ARowMat) && any(size(ARowMat) ~= n)  
        error('ARowMat is not of order K');  
    end
end

%  Compute N by K matrix FMAT containing smoothed PC scores

if FColPar > 0 && ~isempty(FColMat)
    % N*K equations in vec(F)
    I_N = eye(N);
    I_K = eye(K);
    Mmat = kron(I_N,amat*amat');
    if FRowPar > 0 && ~isempty(FRowMat)
        Mmat = Mmat + kron(I_N,FRowPar*FRowMat);
    end    
    if FColPar > 0 && ~isempty(FColMat)
        Mmat = Mmat + kron(FColPar*FColMat,I_K);
    end
    Mmatinv = inv(Mmat);
    AXtvec = reshape(amat*ymat',N*K,1);
    fmat   = reshape(AXtvec'*Mmatinv,K,N)';
else
    %  K equations in fmat'
    Mmat = amat*amat';
    if FRowPar > 0 && ~isempty(FRowMat)
        Mmat = Mmat + FRowPar*FRowMat;
    end
    Mmatinv = inv(Mmat);
    fmat    = ymat*amat'*Mmatinv;
end   
yhat = fmat * amat;         %  fit to the data
emat = ymat - yhat;         %  residuals
Hval = sum(sum(emat.^2));   %  error sum of squares   SSE
if ARowPar > 0 && ~isempty(ARowMat)
    Hval = Hval + ARowPar*trace(amat*ARowMat*amat');
end

%  Evaluate the outer loss function H(A)

diHdiA = 2*fmat'*emat;    % -1/2 partial derivative of H    wrt A
diHdiF = 2*amat *emat';   %  1/2 partial derivative of SSE  wrt F
if ARowPar > 0 && ~isempty(ARowMat)
    diPenRow = amat*ARowMat;
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
            if FColPar > 0 && ~isempty(FColMat)
                %  N*K equations in vec(F)
                %  Evaluate the derivative of A*A' wrt Apq
                diAAtdiApq      = zeros(K);
                diAAtdiApq(:,p) = amat(:,q);
                diAAtdiApq(p,:) = diAAtdiApq(p,:) + amat(:,q)'; 
                %  Evaluate the derivative of M wrt Apq
                diMdiApq        = kron(I_N,diAAtdiApq);
                %  Evaluate the derivative of A*X' wrt Apq
                diAXtdiApq      = zeros(K,N);
                diAXtdiApq(p,:) = ymat(:,q);
                diAXtdiApq      = reshape(diAXtdiApq,1,N*K);
                %  Evaluate the derivative of F wrt Apq
                diFdiApq = (diAXtdiApq - AXtvec'*Mmatinv*diMdiApq)*Mmatinv;
                diFdiApq = reshape(diFdiApq,K,N)';
            else
                %  K equations in fmat'
                %  Evaluate the derivative of Mmat wrt Apq
                diMdiApq      = zeros(K);
                diMdiApq(:,p) = amat(:,q);
                diMdiApq(p,:) = diMdiApq(p,:) + amat(:,q)';
                %  Evaluate the derivative of M(A)^{-1} wrt Apq
                diMinvdiApq   = -Mmatinv*diMdiApq*Mmatinv;
                %  Evaluate the derivative of Fmat wrt Apq
                Dffacpq       = amat'*diMinvdiApq;
                Dffacpq(q,:)  = Dffacpq(q,:) + Mmatinv(p,:);
                diFdiApq      = ymat*Dffacpq;
            end
            %  Evaluate total derivative of yhat wrt Apq
            diyhatdiApq = diFdiApq*amat;
            diyhatdiApq(:,q) = diyhatdiApq(:,q) + fmat(:,p);
            Hgrad(m) = -2*trace(emat'*diyhatdiApq);
            if ARowPar > 0 && ~isempty(ARowMat)
                Hgrad(m) = Hgrad(m) + 2*ARowPar*diPenRow(p,q);
            end
        end
    end
    
    Hgrad = Davec'*Hgrad;
    
end



