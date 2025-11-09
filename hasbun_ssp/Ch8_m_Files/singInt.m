%copyright by J. E Hasbun and T. Datta
%singInt.m
function y=singInt(cnum,cden,h)
%Estimates the integral of f(x)/g(x) dx on [-h,h]
%where h=(b-a)/(n-1) and n=number of function evaluations.
%The program sums all the contributions on interval [a,b]
%This is Roth's method of integrating singular functions
%[L. M. Roth Phys. Rev. B Vol. 7, p4321-4337 (1973)]
%as implemented by J. E. Hasbun.
%
n=length(cnum); %cden must be the same length as well
y=complex(0.0,0.0);
for k=1:2:n-2
    f0=cnum(k); f1=cnum(k+1); f2=cnum(k+2);
    g0=cden(k); g1=cden(k+1); g2=cden(k+2);
    a=abs(g0+g2-2.0*g1);
    %The amount added to g1 in the next line is slightly adjustable
    if(a < 1.e-7), g1=g1+1i*1.e-7; end
    a=(g0+g2-2.0*g1)/2.0;
    b=(g2-g0)/2.0;
    c=g1;
    q=sqrt(b*b-4.0*a*c);
    xp=(-b+q)/a/2.0;
    xm=(-b-q)/a/2.0;
    xpmxm=xp-xm;
    i1=log((1-xp)*(1+xm)/(1-xm)/(1+xp))/xpmxm/a;
    i2=(log(g2/g0)-b*i1)/a/2.0; %note: a+b+c=g2, a-b+c=g0
    i3=(2.0-b*i2-c*i1)/a;
    y=y+3.0*(f0*(i3-i2)+f2*(i3+i2))/2.0+3.0*f1*(i1-i3);
end
y=y*h/3.0;
