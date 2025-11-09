%copyright by J. E Hasbun and T. Datta
%rombergInt.m
function [result]=rombergInt(a,b,funcInt)
T=zeros(15,5);         %Initial T
error=9999;            %start with a large error
N = 1;
h  = (b-a)/N;
m = 1;
T(m,1) = 0.5*h*(funcInt(a)+funcInt(b));
%fprintf('m=%2.0f, T=%15.10f\n',m,T(m,1))
while ( (m < 15 & abs(error) > 1.e-8) | (m <= 3)) %require a min # of steps
    m = m + 1;
    N = 2*N;
    h  =(b-a)/N;
    extra = 0.0;
    for i=1:2:N-1
        extra = extra + funcInt(a+i*h);
    end
    T(m,1) = 0.5*T(m-1,1)+h*extra;
    T(m,2) = (4.0*T(m,1) - T(m-1,1) )/ 3.0;
    if(m>=3), T(m,3) = (16.0*T(m,2) - T(m-1,2)) / 15.0; end
    if(m>=4), T(m,4) = (64.0*T(m,3) - T(m-1,3)) / 63.0; end
    if(m>=5), T(m,5) =(256.0*T(m,4) - T(m-1,4)) /255.0; end
    %fprintf('m=%2.0f, T=%15.10f,%15.10f,%15.10f,%15.10f,%15.10f\n',...
    %    m,T(m,1),T(m,2),T(m,3),T(m,4),T(m,5))
    %Use error to check convergence...
    %add a small part ot the denominator to avoid division by zero
    error =(T(m,min(m,5))-T(m,min(m,5)-1))/sqrt(T(m,min(m,5))^2+1.e-12^2);
end
result = T(m,min(m,5));
