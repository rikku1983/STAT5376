function [T,P,Q,B,B0] = pls1(X, y, l)
[m,n] = size(X);
W = zeros(n,l); % W = W_0,...,W_{l-1}
P = zeros(n,l); % P = P_0,...,P_{l-1}
Q = zeros(l,1); % Q = q_0,...,q_{l-1}
T = zeros(m,l); % T = T_0,...,T_{l-1}

Xk = X;             % initialize Xk to X
Wk = X'*y/norm(y);  % Wk = normalized dot products of cols of X and y
Tk = X*Wk;          % Tk = weighted sum of cols of X
for k=0:l
    if k < l
        W(:,k+1) = Wk;
    end
    tk = Tk'*Tk; % squared norm of Tk
    Tk = Tk./tk; % Tk divided by its squared norm
    Pk = Xk'*Tk; % Pk dot products of cols of Xk and Tk
    qk =  y'*Tk; % qk dot products of y and Tkl
    if qk == 0
        l = k;
        break;
    else
        Xk = Xk - tk*Tk*Pk';  % subtract approximation of Xk from Xk
        Wk = Xk'*y;           % Wk = dot products of cols of Xk and y
        Tk = Xk*Wk;           % Tk = weighted sum of cols of X
    end
    if k < l
        P(:,k+1) = Pk;
        Q(k+1)   = qk;
        T(:,k+1) = Tk;
    end
end
B = W*((P'*W)\Q);   % regression coefficients for model Y = X*B + B0
B0 = Q(1) - P(:,1)'*B;  % intercept for model