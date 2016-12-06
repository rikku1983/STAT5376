function H=DPalog(f,g)
c=1;
n=length(f);
H=zeros(n,n);
H(1,:)=c;
H(:,1)=c;
H(1,1)=0;
v=[1 1;2 1; 1 3; 2 3; 3 2; 3 1];
for i=2:n;
    for j=2:n;
        for r=1:size(v,1);
            k=i-v(r,1);
            l=j-v(r,2);
            if (k>0 && l>0)
                CandE(r) = H(k,l) + cosfn(f,g,k,l,i,j);
            else
                CandE(r)=c;
            end
        end
        [H(i,j),idx] =min(CandE);
        path(i,j,1) = i-v(idx,1);
        path(i,j,2) = j-v(idx,2);
    end
end

x(1) = n;
y(1) = n;
cnt = 1;

while x(cnt)>1;
    x(cnt+1) = path(x(cnt),y(cnt),1);
    y(cnt+1) = path(x(cnt),y(cnt),2);
    cnt = cnt+1;
end

figure(1);
imagesc(H');
axis xy;
axis equal;
colormap(gray);
hold on
plot(x,y);
figure(2);
plot(f);hold on
plot((1:100),g(round(interp1(x,y,1:100))));ylim([0,8]);title('registered');

    
    