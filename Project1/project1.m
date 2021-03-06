%Project 1
clear all;
close all;
load('DataFile1.mat');
load('DataFile2.mat');
load('DataFile3.mat');
load('DataFile4.mat');
load('DataFile5.mat');
imagesc(D1);
imagesc(D2);
imagesc(D3);
imagesc(D4);
imagesc(D5);
subplot(321);hist(reshape(D1, 192*128,1),200);title('D1');
subplot(322);hist(reshape(D2, 192*128,1),200);title('D1');
subplot(323);hist(reshape(D3, 192*128,1),200);title('D1');
subplot(324);hist(reshape(D4, 192*128,1),200);title('D1');
subplot(325);hist(reshape(D5, 192*128,1),200);title('D1');


%Validify MH algorithm for sigma1=10
m=-50:350;
d=-50:350;
for i=1:401;
    for j=1:401;
        tm(i,j)=MC2norm(i-51,10,j-51,30,10000);
        mcm(i,j)=MH(600,i-51,10,j-51,30,500);
    end
end
image(abs(mcm-tm));colorbar;xlabel('Mean of prior/proposal normal');ylabel('Mean of likelihood norm');title('Diff between MH mean and MC mean of Posterior');
m=-50:350;
d=-50:350;
for i=1:401;
    for j=1:401;
        tm(i,j)=MC2norm(i-51,10,j-51,30,1000);
        mcm(i,j)=MH(600,i-51,10,j-51,30,100);
    end
end
image(abs(mcm-tm));colorbar;xlabel('Mean of prior/proposal normal');ylabel('Mean of likelihood norm');title('Diff between MH mean and MC mean of Posterior');

%Validify MH algorithm for sigma1=20
m=-50:350;
d=-50:350;
for i=1:401;
    for j=1:401;
        tm(i,j)=MC2norm(i-51,20,j-51,30,1000);
        mcm(i,j)=MH(600,i-51,20,j-51,30,500);
    end
end
image(abs(mcm-tm));colorbar;xlabel('Mean of prior/proposal normal');ylabel('Mean of likelihood norm');title('Diff between MH mean and MC mean of Posterior'),
m=-50:350;
d=-50:350;
for i=1:401;
    for j=1:401;
        tm(i,j)=MC2norm(i-51,20,j-51,30,1000);
        mcm(i,j)=MH(600,i-51,20,j-51,30,100);
    end
end
image(abs(mcm-tm));colorbar;xlabel('Mean of prior/proposal normal');ylabel('Mean of likelihood norm');title('Diff between MH mean and MC mean of Posterior')       

%Validify MH algorithm for sigma1=100
m=-50:350;
d=-50:350;
for i=1:401;
    for j=1:401;
        tm(i,j)=MC2norm(i-51,100,j-51,30,1000);
        mcm(i,j)=MH(600,i-51,100,j-51,30,500);
    end
end
image(abs(mcm-tm));colorbar;xlabel('Mean of prior/proposal normal');ylabel('Mean of likelihood norm');title('Diff between MH mean and MC mean of Posterior'),
m=-50:350;
d=-50:350;
for i=1:401;
    for j=1:401;
        tm(i,j)=MC2norm(i-51,100,j-51,30,1000);
        mcm(i,j)=MH(600,i-51,100,j-51,30,100);
    end
end
image(abs(mcm-tm));colorbar;xlabel('Mean of prior/proposal normal');ylabel('Mean of likelihood norm');title('Diff between MH mean and MC mean of Posterior')       


%Gibbs Sampling
%Initialization
clear all;
close all;
sig1=10;
nr=128;
nc=192;
npic=5;
m=NaN(nr+2,nc+2);
load('DataFile5.mat')
m(2:nr+1,2:nc+1)=D5;
finalm(:,:,1)=m(2:nr+1,2:nc+1);
steps=20;
for n=1:steps
   for i=2:(nr+1)
       for j=2:(nc+1)
           %sample from density(full conditionals)
           im=mean([m(i-1,j),m(i+1,j),m(i,j-1),m(i,j+1)],'omitnan');
           m(i,j)=MH2(1/2*(im+D5(i-1,j-1)),im,sig1,D5(i-1,j-1));
       end
   end
   finalm(:,:,n+1)=m(2:nr+1,2:nc+1);
end
save('D1_1000m.mat', 'finalm');
%Draw figures
npic=6;
for i=1:npic
    intv=floor(20/(npic-1));
    subplot(3,2,i);imagesc(mean(finalm(:,:,1:max(1,intv*(i-1))),3));title(['iteration= ' num2str(intv*(i-1)) '.']);
end