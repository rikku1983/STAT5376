function [conm,x,r,conma] = MH(x0, m, sig1, d, sig2, step)
x(1)=x0;
%rej=0;
for n=1:step
    %Generate y from N(m, sigma1) 
    y=normrnd(m,sig1);
    rho=min(1,exp(((x(n)-d)^2-(y-d)^2)/(2*sig2^2)));
    r(n)=rho;
    if rand()>rho
        x(n+1)=x(n);
        %rej=rej+1;
    else
        x(n+1)=y;
    end
end
%rejr=rej/step;
conma=cumsum(x)./(1:step+1);
conm=conma(step+1);
%tm=MC2norm(m,sig1,d,sig2,10000);
%close all;
%subplot(221);plot(1:step+1,x);title('trace of x');hold on; plot(get(gca,'xlim'),[tm,tm]);
%subplot(222);plot(1:step+1,conm);title('convergence of mean');;hold on; plot(get(gca,'xlim'),[tm,tm]);
%subplot(223);plot(1:step,r);title('trace of acceptance ratrio');



