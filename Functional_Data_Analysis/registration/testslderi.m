close all;
clear all;
x1=0:0.1:2*pi;
y1=sin(x1);
fig1=figure();
f=[x1',y1'];
f1=slderi(f,1,1);
f2=slderi(f,2,1);
f3=slderi(f,3,1);
f4=slderi(f,4,1);
f5=slderi(f,5,1);



set(fig1,'Position', [200 200 600 450]);
subplot(321);plot(x1,y1);
subplot(322);plot(f1(:,1),f1(:,2));
subplot(323);plot(f2(:,1),f2(:,2));
subplot(324);plot(f3(:,1),f3(:,2));
subplot(325);plot(f4(:,1),f4(:,2));
subplot(326);plot(f5(:,1),f5(:,2));