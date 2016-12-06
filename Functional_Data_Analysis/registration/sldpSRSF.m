%f,g: 2 functions in the form of n by 2 matrix, where f(:,1)=x, and
%g(:,2)=y, will be sorted according to x
%smthpara: from 0 to 1, 1 means completely interpolation.
%This function search neiborghood v
%n is the resolution of the grid, default is 100
function [path,E]=sldpSRSF(f,g,smthpara,n,m)
% Fill in unset optional values.
%define energy for impossible places
c=inf;
n=length(f);
m=length(g)
%Smooth x1 and x2 by splines and smthpara
fs=fit(f(:,1), f(:,2), 'smoothingspline', 'SmoothingParam', smthpara);
gs=fit(g(:,1), g(:,2), 'smoothingspline', 'SmoothingParam', smthpara);
%Generate q1 from f and q2 from g
x1=(0:n-1)/(n-1);
x2=0:1/((m-1)/(n-1)*n-1):1;
%x2=0:1/m:1;
for i = 1:length(x1)
    q1(i)=sign((fs(x1(i)+0.0001)-fs(x1(i)-0.0001))/(2*0.0001)).*sqrt(abs((fs(x1(i)+0.0001)-fs(x1(i)-0.0001))/(2*0.0001)));
end
for i=1:length(x2)
    q2(i)=sign((gs(x2(i)+0.0001)-gs(x2(i)-0.0001))/(2*0.0001)).*sqrt(abs((gs(x2(i)+0.0001)-gs(x2(i)-0.0001))/(2*0.0001)));
end
%Initialize E grid
E=zeros(n+1,n+1);
E(1,:)=c;
E(:,1)=c;
E(1,1)=0;
%define neigbors
v=[1,1;2,1;3,1;4,1;5,1;6,1;1,2;1,3;1,4;1,5;1,6;2,3;3,2;3,4;4,3;2,5;3,5;4,5;5,2;5,3;5,4;5,6;6,5;
    1,7;2,7;3,7;4,7;5,7;6,7;7,1;7,2;7,3;7,4;7,5;7,6;1,8;3,8;5,8;7,8;8,7;8,5;8,3;8,1];
for i=2:n;
    for j=2:n;
        for r=1:size(v,1);
            k=i-v(r,1);
            l=j-v(r,2);
            if (k>0 && l>0)
                CandE(r) = E(k,l) + energySRSF(q1,q2,k,l,i,j);
            else
                CandE(r)=c;
            end
        end
        [E(i,j),idx] =min(CandE);
        path(i,j,1) = i-v(idx,1);
        path(i,j,2) = j-v(idx,2);
    end
end
%reconstruct gamma
x(1) = n;
y(1) = n;
cnt = 1;
while x(cnt)>1;
    x(cnt+1) = path(x(cnt),y(cnt),1);
    y(cnt+1) = path(x(cnt),y(cnt),2);
    cnt = cnt+1;
end
path=[x',y'];
