%1
x0=rand();
x(1)=sqrt(abs(x0))+normrnd(0,0.01);
y(1)=x(1)^2+normrnd(0,0.01);
for t=2:20
    x(t)=sqrt(abs(x(t-1)))+normrnd(0,0.01);
    y(t)=x(t)^2+normrnd(0,0.01);
end

clc;clear all ;close all;
n=5000;t=20;
x(1,:)=rand(1,n);
for i=1:t
    xtemp=sqrt(abs(x(i,:)))+normrnd(0,0.1,1,n);
    wtemp(i,:)=xtemp.^2+normrnd(0,0.1,1,n);
    wt(i,:)=wtemp(i,:)/(sum(wtemp(i,:)));
    for j=1:n
        u=rand;
        idx=min(find(u<cumsum(wt(i,:))));
        x(i+1,j)=xtemp(idx);
    end
end

fig=figure();
set(fig,'Position', [20 20 1000 3200]);
subplot(321);hist(x(2,:),100);title(['x for t=1']);
subplot(322);hist(x(6,:),100);title(['x for t=5']);
subplot(323);hist(x(11,:),100);title([' x for t=10']);
subplot(324);hist(x(16,:),100);title(['x for t=15']);
subplot(325);hist(x(21,:),100);title([' x for t=20']);

        
%2
clc;clear all ;close all;
n=1000;t=20;
x(1,:)=rand(1,n);
for i=1:t
    xtemp=x(i,:)/2+25*x(i,:)./(1+x(i,:).^2)+8*cos(1.2*(i-1))+normrnd(0,10,1,n);
    wtemp(i,:)=xtemp.^2/20+normrnd(0,1,1,n);
    wt(i,:)=wtemp(i,:)/(sum(wtemp(i,:)));
    for j=1:n
        u=rand;
        idx=min(find(u<cumsum(wt(i,:))));
        x(i+1,j)=xtemp(idx);
    end
    theta(i)=norm(x(i+1,:)-mean(x(i+1,:)));
end

fig=figure();
set(fig,'Position', [20 20 1000 3200]);
subplot(321);hist(x(2,:),30);title(['x for t=1']);
subplot(322);hist(x(6,:),30);title(['x for t=5']);
subplot(323);hist(x(11,:),30);title([' x for t=10']);
subplot(324);hist(x(16,:),30);title(['x for t=15']);
subplot(325);hist(x(21,:),30);title([' x for t=20']);
subplot(326);plot(theta);title('SD of Posterior');xlabel('t');