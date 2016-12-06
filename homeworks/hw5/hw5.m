clear all;
close all;
%define functions
n=100;
m=1000;
xx1=(0:n)/n;
xx2=(0:m)/m;
f=[xx1',normpdf(xx1,0.7,0.06)']*0.8;
g=[xx2',normpdf(xx2,0.3,0.05)'];

%order f
[x1,idx1]=sort(f(:,1));
ytemp=f(:,2);
y1=ytemp(idx1);
f=[x1,y1];
%order g
[x2,idx2]=sort(g(:,1));
ytemp=g(:,2);
y2=ytemp(idx2);
g=[x2,y2];


[path,E]=sldpSRSF(f,g,1);

close all;
fig=figure();
set(fig,'Position', [200 200 1600 450]);
subplot(131);plot(x1,y1);hold on;plot(x2,y2);title('before registration');
set(gca,'XTickLabel',[],'YTickLabel',[]);
subplot(132);imagesc(E');colormap(gray);hold on;plot(path(:,1),path(:,2));axis xy;
axis equal;set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Warping function');
subplot(133);plot(x1,y1);hold on;
plot(x2,y2(round(interp1(path(:,1)-1,path(:,2)-1,0:n/m:n)*m/n)+1));title('after registered');
set(gca,'XTickLabel',[],'YTickLabel',[]);