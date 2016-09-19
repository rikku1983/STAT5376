function x = MCsim(pai, n, start)
nstat=size(pai,1);
steps=n;
    if nargin < 3
        start=1:nstat;
    end
nsim=length(start);
for i=1:nsim
    x(i,1)=start(i); %start status at t1;
    for t=1:steps
        u=rand;
        if(u<pai(x(i,t),1))
            x(i,t+1)=1;
        elseif(u<sum(pai(x(i,t),1:2)))
            x(i,t+1)=2;
        elseif(u<sum(pai(x(i,t),1:3)))
            x(i,t+1)=3;
        else
            x(i,t+1)=4;
        end
    end
end

