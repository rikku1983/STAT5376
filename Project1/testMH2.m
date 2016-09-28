for i=1:2000
    t(i)=MH2(300,10,30,15);
end
hist(t,100);
mean(t);