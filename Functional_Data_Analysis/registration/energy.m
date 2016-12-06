function cost=energy(f,g,k,l,i,j)
%This is the function to calculate 2-norm distence between two fuction
%between path from (k,l) to (i,j) used in dynamic programming
n=length(f);
m=length(g);
slope=(j-l)/(i-k);
gidx=round((l+((k+1:i)-k).*slope)/n*m);
cost=norm(f(k+1:i)-g(gidx))^2/n;


