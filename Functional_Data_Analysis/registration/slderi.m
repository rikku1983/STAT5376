%function to calculate derivative based on inpute discrete data in f as a 2
%column matrix, x is col1 and y is col2. d is order of derivative from 0 and above. smthpara is smooth parameter from
%0 to 1. 
function q=slderi(f,d,smthpara)
if d==0
    q=f;
else
    ftemp=slderi(f,d-1,smthpara);
    x=ftemp(:,1);
	y=ftemp(:,2);
    tempx=(x-min(x))/range(x);
    q(:,1)=x;
    fs=fit(tempx, y, 'smoothingspline', 'SmoothingParam', smthpara);
    for i = 1:length(tempx)
        q(i,2)=(fs(tempx(i)+0.00001)-fs(tempx(i)-0.00001))/(0.00002)/range(x);
    end
end


