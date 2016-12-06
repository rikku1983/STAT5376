%Part I
clear all;
close all;
%define functions
n=100;
m=100;
xx1=(0:n)/n;
xx2=(0:m)/m;
f=normpdf(xx1,0.3,0.04)+normpdf(xx1,0.5,0.04);
g=normpdf(xx2,0.4,0.04)+normpdf(xx2,0.6,0.04);

%Smooth x1 and x2 by splines
smthpara=1;
fs=fit(xx1', f', 'smoothingspline', 'SmoothingParam', smthpara);
gs=fit(xx2', g', 'smoothingspline', 'SmoothingParam', smthpara);
%Generate q1 from f and q2 from g

%x2=0:1/m:1;
for i = 1:length(xx1)
    q1(i)=sign((fs(xx1(i)+0.0001)-fs(xx1(i)-0.0001))/(2*0.0001)).*sqrt(abs((fs(xx1(i)+0.0001)-fs(xx1(i)-0.0001))/(2*0.0001)));
end
for i=1:length(xx2)
    q2(i)=sign((gs(xx2(i)+0.0001)-gs(xx2(i)-0.0001))/(2*0.0001)).*sqrt(abs((gs(xx2(i)+0.0001)-gs(xx2(i)-0.0001))/(2*0.0001)));
end

[path, E]=sldpSRSF2(q1,q2);

close all;
fig=figure();
set(fig,'Position', [200 200 1600 450]);
subplot(131);plot(f);hold on;plot(g((1:n)*m/n));title('before registration');
set(gca,'XTickLabel',[],'YTickLabel',[]);
subplot(132);imagesc(E');colormap(gray);hold on;plot(path(:,1),path(:,2));axis xy;
axis equal;set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Warping function');
subplot(133);plot(f);hold on;
plot((1:n),g(round(interp1(path(:,1),path(:,2),1:n)*(m+1)/(n+1))));title('after registered');
set(gca,'XTickLabel',[],'YTickLabel',[]);

%Calculate da and dp
%da is E
da=E(n+1,n+1)
%plot(path(:,1),path(:,2));
pathd=slderi(path,1,1);
%plot(pathd(:,1),pathd(:,2));
dp=acos(sum(sqrt(pathd(:,2)))/n)


%%%Part2, Growth data.
close all, clear all;
dat=csvread('bgd.csv',1,0);
age=dat(:,1);
boy=dat(:,2:40);
girl=dat(:,41:length(dat));
%Visualization of data
fig1=figure();
set(fig1,'Position', [200 200 1100 450]);
subplot(121);plot(age,boy);title('Heights of Boys');xlabel('age');ylabel('Heights in cm');
subplot(122);plot(age,girl);title('Heights of Girls');xlabel('age');ylabel('Heights in cm');
%Derivatives(Growth rates)
%Smooth function and find growth rates.
%Pick two curve.
c1=girl(:,29);
c2=girl(:,28);
smthpara=1;
c1s=fit(age,c1, 'smoothingspline', 'SmoothingParam', smthpara);
c2s=fit(age,c2, 'smoothingspline', 'SmoothingParam', smthpara);
%Generate new functions gr1 and gr2 for derivative
x=1:0.2:18;
for i = 1:length(x)
    gr1(i)=(c1s(x(i)+0.02)-c1s(x(i)-0.02))/(2*0.02);
end
for i=1:length(x)
    gr2(i)=(c2s(x(i)+0.02)-c2s(x(i)-0.02))/(2*0.02);
end
%plot the functions of growth rate to be registered
close all;
plot(x,gr1);hold on;plot(x,gr2);title('Growth rates');xlabel('age');ylabel('Growth rate');

%Now lets calculate Q
%Smooth x1 and x2 by splines and scale the curve to between 0 and 1
smthpara=1;
fs=fit(((x-1)/max(x-1))', gr1', 'smoothingspline', 'SmoothingParam', smthpara);
gs=fit(((x-1)/max(x-1))', gr2', 'smoothingspline', 'SmoothingParam', smthpara);
close all;
fig1=figure();
set(fig1,'Position', [200 200 1100 450]);
subplot(121);plot(fs);hold on; plot(gs, 'b');title('Growth rate');
%Generate q1 from f and q2 from g
%x2=0:1/m:1;
n=100;
m=1000;
xx1=(0:n)/n;
xx2=(0:m)/m;

for i = 1:length(xx1)
    q1(i)=sign((fs(xx1(i)+0.0001)-fs(xx1(i)-0.0001))/(2*0.0001)).*sqrt(abs((fs(xx1(i)+0.0001)-fs(xx1(i)-0.0001))/(2*0.0001)));
end
for i=1:length(xx2)
    q2(i)=sign((gs(xx2(i)+0.0001)-gs(xx2(i)-0.0001))/(2*0.0001)).*sqrt(abs((gs(xx2(i)+0.0001)-gs(xx2(i)-0.0001))/(2*0.0001)));
end
subplot(122);plot(xx1,q1);hold on; plot(xx2,q2);title('Functions q1 and q2');

% Ready to register
[path, E]=sldpSRSF2(q1,q2);

close all;
f=fs(xx1);
g=gs(xx2);
fig=figure();
set(fig,'Position', [200 200 1600 450]);
subplot(131);plot(f);hold on;plot(g((1:n)*m/n));title('before registration');
set(gca,'XTickLabel',[],'YTickLabel',[]);
subplot(132);imagesc(E');colormap(gray);hold on;plot(path(:,1),path(:,2));axis xy;
axis equal;set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Warping function');
subplot(133);plot(f);hold on;
plot((1:n),g(round(interp1(path(:,1),path(:,2),1:n)*(m+1)/(n+1))));title('after registered');
set(gca,'XTickLabel',[],'YTickLabel',[]);

da=E(n+1,n+1)
pathd=slderi(path,1,1);
dp=acos(sum(sqrt(pathd(:,2)))/n)