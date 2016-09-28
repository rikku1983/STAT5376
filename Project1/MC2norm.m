function [avg,mcm,y,z]=MC2norm(m, sig1, d, sig2, step2)
%function of calculating mean of distribution of product of two normal
%density with mean m and d, variance sig1^2 and sig2^2.
y=normrnd(d,sig2,1, step2);
z=1/sqrt(2*pi*sig1^2)*exp(-1/2/(sig1^2).*((y-m).^2));
mcm=y.*(1/sqrt(2*pi*sig1^2)*exp(-1/2/(sig1^2).*((y-m).^2)));
%subplot(221);hist(y,200);title('y');
%subplot(222);hist(z,200);title('z');
%subplot(223);hist(mcm,200);title('mcm');
avg=mean(mcm)/mean(z);