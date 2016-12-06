%HW1
%Classical MC
clear all;
close all;

a=5;
i=0;
x=normrnd(0,1,[1,2000]);
for n=5:5:2000
    i=i+1;
    y=x(1:n);
    theta1(i)=mean(y>a);
end

i=0;
x=normrnd(5,1,[1,2000]);
for n=5:5:2000
    i=i+1;
    y=x(1:n);
    theta2(i)=mean((y>a).*exp(a^2/2-a*y));
end

real=normcdf(5, 'upper');
close all;
subplot(1,2,1);plot(5:5:2000, theta1);hold on; plot([5 2000], [real real]);legend('Simulated theta hat','Real Value');
ylabel('Classical MC estimator');xlabel('number of n simulated');axis([0 2000 0 0.000001]);
subplot(1,2,2);plot(5:5:2000, theta2);hold on; plot([5 2000], [real real]);legend('Simulated theta hat','Real Value');
ylabel('Tilt sampling estimator');axis([0 2000 0 0.000001]);xlabel('number of n simulated');