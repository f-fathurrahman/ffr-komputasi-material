%copyright by J. E Hasbun and T. Datta
%Green_1DdosCalc.m
%In addition to calculating the one-dimensional Green's function in
%the nearest neighbor tight binding model analytically, a comparison
%is made with the purely numerical calculation. The numerical
%calculation uses the Roth's integration scheme for singular
%functions; i.e., singint.m (separate listing in the Appendix).
function Green_1DdosCalc
clear
el=-2;                         %upper energy value
eu=2;                          %lower energy value
Ne=100;                        %number of energy points
es=(eu-el)/(Ne-1);             %energy step
gam=1/2;                       %band energy parameter (in Ha)
zim=complex(0.,1.0);           %imaginary number
delta=1.e-4;                   %small part for plottting, integrating
for i=1:Ne
  e0(i)=el+(i-1)*es;
  g00(i)=e0(i)/(sqrt((e0(i)+zim*delta)^2-4*gam^2)*abs(e0(i)));
  dosa(i)=-imag(g00(i))/pi;    %analytic density of states
end
%The plots are done below
%Numerical calculation follows
a=1;
kl=0;                          %upper limit
ku=2*pi/a;                     %lower limit
Nk=301;                        %number of k-points
dk=(ku-kl)/(Nk-1);
factor=a/2.0/pi;               %factor for the numerical integral
for j=1:Nk
  k(j)=kl+(j-1)*dk;
  Ek(j)=-2*gam*cos(k(j)*a);    %the band energy defined only once
  top(j)=1.0;                  %the numerator for integration
end
%Integration over k for each energy value
for i=1:Ne
  for j=1:Nk
    deno(j)=e0(i)-Ek(j)+zim*delta;    %the denominator
  end
  ge(i)=factor*singInt(top,deno,dk);  %ge is the numerical integral
  dosn(i)=-imag(ge(i))/pi;            %numerical density of states
end
%
%Plotting
subplot(1,2,1)
plot(e0,real(g00),'k');               %Real part of g00 - analytic
hold on
plot(e0,real(ge),'ko','MarkerSize',2) %numeric
legend ('analytic','numeric',0)
xlabel('E (Ha)')
ylabel('Real(G_{00}) (1/Ha)')
title('Real part of Green''s function (1D)')
%
subplot(1,2,2)
plot(e0,dosa,'k')                     %density of states - analytic
hold on
plot(e0,dosn,'ko','MarkerSize',2)     %numeric
xlabel('E (Ha)')
ylabel('D(E) (1/Ha)')
title('Density of States (1D)')
legend ('analytic','numeric',0)
