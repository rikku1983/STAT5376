load DataFile1.mat 
I=D1;
[n1,n2]=size(D1); 
sigma1=10; 
sigma2=30;
figure(1)
imagesc(I(:,:));
title('Initial Image'); 

for i=1:8;
    for j=1:n1;
        for k=1:n2;
            mid=mean1(j,k,I);
            if rand>=0.5;
            I(j,k) = random('normal',mid,sigma1,1,1);
            else
            I(j,k)=I(j,k);
            end;
        end;
     end;
W=random('normal',0,sigma2,n1,n2); 
D=I+W; 
          I2= Gibbs(I,D,sigma1,sigma2); 
          I=I2;
figure(2)
subplot(2,4,i);
imagesc(I)
figname = sprintf('Image of sweep %d',i);
title (figname);
end;




