clear all
clc;
close all;
n=100;
x=[0:n]/n;
f=normpdf(x,0.3,0.05)*0.4;
g=normpdf(x,0.5,0.05)*0.6;

figure(3);
plot(f);
hold on;
plot(g);
H=DPalog(f,g);


