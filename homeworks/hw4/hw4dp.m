clear all;
close all;
%define functions
n=100;
m=1000;
xx1=(0:n)/n;
xx2=(0:m)/m;
f=normpdf(xx1,0.7,0.05)*0.8;
g=normpdf(xx2,0.3,0.05);


[path, E]=sldp(f,g);
close all;
fig=figure();
set(fig,'Position', [200 200 1600 450]);
subplot(131);plot(f);hold on;plot(g((1:n)*m/n));title('before registration');
set(gca,'XTickLabel',[],'YTickLabel',[]);
subplot(132);imagesc(E');colormap(gray);hold on;plot(path(:,1),path(:,2));axis xy;
axis equal;set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Warping function');
subplot(133);plot(f);hold on;
plot((1:n),g(round(interp1(path(:,1),path(:,2),1:n)*m/n)));title('after registered');
set(gca,'XTickLabel',[],'YTickLabel',[]);