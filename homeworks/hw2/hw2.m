%Q1, generate 4 sample paths with frequency of each states versus i up to
%1000
clear all;
close all;
pai=[0.2 0.2 0.1 0.5; 0.1 0.3 0.4 0.2; 0.3 0.2 0.3 0.2; 0.1 0.3 0.1 0.5];
%generate MC chain with 4 different starting states
x=MCsim(pai, 1000);
%Calculate frequency versus i
y=@(n,i)[mean(x(i,1:n)==1), mean(x(i,1:n)==2),mean(x(i,1:n)==3),mean(x(i,1:n)==4)];
steps=1000;
for n=1:steps
    xper1(n,1:4)=y(n,1);
    xper2(n,1:4)=y(n,2);
    xper3(n,1:4)=y(n,3);
    xper4(n,1:4)=y(n,4);
end
subplot(221);plot(1:steps,xper1);title('start at 1');xlim([1,steps+1]);legend('1','2','3','4');
subplot(222);plot(1:steps,xper2);title('start at 2');xlim([1,steps+1]);legend('1','2','3','4');
subplot(223);plot(1:steps,xper3);title('start at 3');xlim([1,steps+1]);legend('1','2','3','4');
subplot(224);plot(1:steps,xper4);title('start at 4');xlim([1,steps+1]);legend('1','2','3','4');

%Calculate dominant eigenvector
[V,D,W]=eig(pai)
W(:,1)/sum(W(:,1))
%Simulated frequency
xper1(1000,:)

%question 2, for another pi which is reducible
clear all;
close all;
pai=[0.5 0.5 0 0; 0.1 0.9 0 0; 0 0 0.3 0.7; 0 0 0.2 0.8];
x=MCsim(pai, 1000);
%Calculate frequency versus i
y=@(n,i)[mean(x(i,1:n)==1), mean(x(i,1:n)==2),mean(x(i,1:n)==3),mean(x(i,1:n)==4)];
steps=1000;
for n=1:steps
    xper1(n,1:4)=y(n,1);
    xper2(n,1:4)=y(n,2);
    xper3(n,1:4)=y(n,3);
    xper4(n,1:4)=y(n,4);
end
subplot(221);plot(1:steps,xper1);title('start at 1');xlim([1,steps+1]);legend('1','2','3','4');
subplot(222);plot(1:steps,xper2);title('start at 2');xlim([1,steps+1]);legend('1','2','3','4');
subplot(223);plot(1:steps,xper3);title('start at 3');xlim([1,steps+1]);legend('1','2','3','4');
subplot(224);plot(1:steps,xper4);title('start at 4');xlim([1,steps+1]);legend('1','2','3','4');

%Calculate dominant eigenvector
[V,D,W]=eig(pai)
W(:,2)/sum(W(:,2))
W(:,4)/sum(W(:,4))
%Simulated frequency
xper1(1000,:)
xper2(1000,:)
xper3(1000,:)
xper4(1000,:)

% For the pi whose periodicity is 2
clear all;
close all;
pai=[0 0.5 0 0.5;0.5 0 0.5 0; 0 0.5 0 0.5; 0.5 0 0.5 0];
x=MCsim(pai, 1000);
%Calculate frequency versus i
y=@(n,i)[mean(x(i,1:n)==1), mean(x(i,1:n)==2),mean(x(i,1:n)==3),mean(x(i,1:n)==4)];
steps=1000;
for n=1:steps
    xper1(n,1:4)=y(n,1);
    xper2(n,1:4)=y(n,2);
    xper3(n,1:4)=y(n,3);
    xper4(n,1:4)=y(n,4);
end
subplot(221);plot(1:steps,xper1);title('start at 1');xlim([1,steps+1]);legend('1','2','3','4');
subplot(222);plot(1:steps,xper2);title('start at 2');xlim([1,steps+1]);legend('1','2','3','4');
subplot(223);plot(1:steps,xper3);title('start at 3');xlim([1,steps+1]);legend('1','2','3','4');
subplot(224);plot(1:steps,xper4);title('start at 4');xlim([1,steps+1]);legend('1','2','3','4');

%Calculate dominant eigenvector
[V,D,W]=eig(pai)
W(:,4)/sum(W(:,4))
%Simulated frequency
xper1(1000,:)
xper2(1000,:)
xper3(1000,:)
xper4(1000,:)

%Q3
clear all;
close all;
pai=[0.1 0.3 0.4 0.2;0.2 0.1 0.3 0.4;0.4 0.2 0.1 0.3; 0.3 0.4 0.2 0.1];
[V,D,W]=eig(pai);
W(:,1)/sum(W(:,1));
%simulation
x=MCsim(pai, 1000, 1);
%Calculate frequency versus i
f=[2,1,2.5,-1];
steps=1000;
for n=1:steps
    xper1(n)=mean(f(x(1:n)));
end
plot(1:steps,xper1);title('start at 1');xlim([1,steps+1]);
