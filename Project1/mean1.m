
%%%M-H

function [mid] = mean1(j,k,x)
[n1,n2]=size(x);
if (j==1) && (k==1)
    mid=(x(j,k+1)+x(j+1,k))/2; 
end;
if (j==1)&&(k==n2)
    mid=(x(1,k-1)+x(j+1,k))/2; 
end;
 
if (j==n1)&&(k==1)
mid=(x(j-1,k)+x(j,k+1))/2; 
end;
if (j==n1) && (k==n2)
    mid=(x(j,k-1)+x(j-1,k))/2; 
end;
if (j==1 && k~=1 && k~=n2)
    mid=(x(j+1,k)+x(j,k-1)+x(j,k+1))/3;
end;
if (j==n1 && k~=1 && k~=n2)
mid=(x(j-1,k)+x(j,k-1)+x(j,k+1))/3; 
end;
if (j~=1 && j~=n1 && k==1)
    mid=(x(j-1,k)+x(j+1,k)+x(j,k+1))/3; 
end;
if (j~=1 && j~=n1 && k==n2)
mid=(x(j-1,k)+x(j+1,k)+x(j,k-1))/3;
end;
if (j~=1 && j~=n1) && (k~=1 && k~=n2)
    mid=(x(j-1,k)+x(j+1,k)+x(j,k+1)+x(j,k+1))/4;
end;
