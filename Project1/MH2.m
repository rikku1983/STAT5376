%This is the function sample from density of 1/z*N(m,sig1^2)*N(d,30)
function x100 = MH2(x0, m, sig1, d)
sig2=30;
step=100;
x100=x0;
for n=1:step
    %Generate y from N(m, sigma1) 
    y=normrnd(m,sig1);
    rho=exp(((x100-d)^2-(y-d)^2)/(2*sig2^2));
    if rand()<rho
        x100=y;
    end
end

