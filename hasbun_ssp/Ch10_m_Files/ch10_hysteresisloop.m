%copyright by J. E Hasbun and T. Datta
%ch10_hysteresisloop.m

clear;

%here x=phi, y=theta, z=h
fH =@(x,y,z) 0.5*sin(2*(x-y))+z*sin(x); % define fH(x,y,z)
dfH =@(x,y,z) cos(2*(x-y))+z*cos(x);    % define dfH(x,y,z) (the derivative)
y=pi/4;
z1=-2; z2=2; Nz=40; zs=(z2-z1)/(Nz-1);  %h range to work with
z=z1:zs:z2;
x1=-pi; x2=pi; Nx=20; xs=(x2-x1)/(Nx-1);%seek roots using this range of guesses
x=x1:xs:x2;
k=0;
for ic=1:length(z)
  %seek all the aeros of fH and collect all the one that make dfH >=0
  for ix=1:length(x)   %let's sweep through many guesses to find many roots
    xx=fzero(@(x) fH(x,y,z(ic)),x(ix));
    if(dfH(xx,y,z(ic)) >=0)  %collect all the roots for which dfH is >=0
      k=k+1;
      phi(k)=xx;
      h(k)=z(ic);
    end
  end
end
plot(h,cos(phi),'k*')
xlabel('h')
ylabel('Cos(\phi)')
