%Codes for Gibbs sampling
function [gib] = Gibbs(I,D,sigma1,sigma2)
[n1,n2]=size(I);
for j=1:n1;
   for k=1:n2; 
    mid=mean1(j,k,I); 
    mu = (mid/sigma1+D(j,k)/sigma2)*(1/((1/sigma1)^2 +(1/sigma2)^2)); 
    sd = sqrt(1/((1/sigma1)^2+(1/sigma2)^2)); 
    gib(j,k)=random('normal',mu,sd,1,1);
   end;
end;