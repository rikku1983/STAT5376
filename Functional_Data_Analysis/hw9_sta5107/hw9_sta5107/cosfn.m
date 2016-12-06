function E=cosfn(f,g,k,l,i,j)

x=[k:1:i];
m=(j-l)/(i-k);
y=(x-k)*m+l;
idx=round(y);
vec=g(idx);

E=(norm(f(x)-vec))^2/length(f);