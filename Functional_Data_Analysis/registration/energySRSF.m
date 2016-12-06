function cost=energySRSF(q1,q2,k,l,i,j)
%This is the function to calculate 2-norm distence between two fuction
%between path from (k,l) to (i,j) used in dynamic programming
n=length(q1);
m=length(q2);
slope=(j-l)/(i-k);
q2idx=round((l+((k+1:i)-k).*slope)/n*m);
%q2idx=min(q2idx,length(q2));
cost=norm(q1(k+1:i)-q2(q2idx)*sqrt(slope))^2;