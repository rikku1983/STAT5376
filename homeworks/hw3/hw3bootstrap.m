clear all; close all;
load('lawschooldata.mat');
% Q1
a1=datasample(X,15);
mean(a1);
corr(a1);

B=[25 50 100 200 500 1000 2000];
for j=1:7
    for i=1:B(j)
        sample=datasample(a1,15);
        mx(i)=mean(sample(:,1));
        my(i)=mean(sample(:,2));
        rho(i)=corr(sample(:,1), sample(:,2));
    end
    mxse(j)=std(mx);
    myse(j)=std(my);
    rhose(j)=std(rho);
end
subplot(2,2,1),plot(B,mxse);title('SE of LSAT');
subplot(2,2,2),plot(B,myse);title('SE of GPA');
subplot(2,2,3),plot(B,rhose);title('SE of Correlation');

close all;
for i=1:2000
    sample=datasample(X,15);
    rhopop(i)=corr(sample(:,1), sample(:,2));
end
subplot(121),hist(rho, 60);title('histogram of BS Corr (B=2000)');xlim([0,1]);
subplot(122),hist(rhopop, 60);title('histogram of POP Corr (B=2000)');xlim([0,1]);
%rhotrue=corr(X(:,1),X(:,2));
%std(rho)
%std(rhopop)

%Q2
close all; clear all;
dat2=[1 2 3.5 4 7 7.3 8.6 12.4 13.8 18.1];
B=[25 100 200 500 1000 2000];
for j=1:6
    for i=1:B(j)
        sample=sort(datasample(dat2,10));
        tm(i)=mean(sample(3:8));
    end
    se(j)=std(tm);
end

close all; clear all;
dat2=[1 2 3.5 4 7 7.3 8.6 12.4 13.8 18.1];
B=[2,20,100,200,300,500,750,1000,1300,1600,2000];
for j=1:length(B)
    for r=1:20
        for i=1:B(j)
            sample=sort(datasample(dat2,10));
            tm(i,r)=mean(sample(3:8));
        end
    end
    finalse(:,j)=std(tm);
end

scatter(reshape(repmat(B, 20,1), length(B)*20,1),reshape(finalse,length(B)*20,1),15,'filled');
title('Bootstrap SE estimates');xlabel('Bootstrap repeats');ylabel('Estimated SE of trimed mean');

%Q3
clear all;close all;
n=10;
dat3=normrnd(0,1,n,1);
theta=mean(dat3);
%Calculate BS estimate
B=1000;
for i=1:B
    bsdat=datasample(dat3,n);
    bstheta(i)=median(bsdat);
end
bias=theta-mean(bstheta);
bias2=theta-median(dat3);

