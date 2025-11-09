%copyright by J. E Hasbun and T. Datta
%F0_simple_cubic.m
%f0=(1/N)sum_over_k (1/(E0-Ek)) where
%Ek=-cos(kx*a)-cos(ky*a)-cos(kz*a). If we use the density of states
%corresponding to the simple cubic, then
%f0(E)=int(g(E') dE'/(E-E'+im*delta))

function F0_simple_cubic
clear, clc;
delta=1.5e-2;
im=complex(0.0,1.0);
e2=4.0;               %range for energy E
e1=-e2;
ntmax=101;
es=(e2-e1)/(ntmax-1); %energy step
a=1.0;                %lattice constant
tpa=2*pi/a;
VBZ=tpa^3;            %total SC BZ volume
for nt=1:ntmax
  e0(nt)=e1+(nt-1)*es;                  %energy E
  ge(nt)=jelittoScDosAnal(e0(nt))/VBZ;  %exact SC dos
end
plot(e0,ge,'k:','LineWidth',2), hold on %exact dos plot
xlabel('E (Hartrees)'), ylabel('D(E), F_0(E)')
str=cat(2,'D(E) and F_0 for the simple cubic (no spin)');
title(str, 'Fontsize',12)
%Repeat the e loop
e2p=3.0;                                %range for E'
e1p=-e2p;
ntpmax=201;
esp=(e2p-e1p)/(ntpmax-1);               %E' step
for nt=1:ntmax
  %Integrate over the e' loop
  for ntp=1:ntpmax
    if(nt==1)                      %calculate this part only once
      e0p(ntp)=e1p+(ntp-1)*esp;    %energy E'
      top(ntp)=jelittoScDosAnal(e0p(ntp))/VBZ; %use SC density of states
                                               %as the numerator
    end
    deno(ntp)=e0(nt)-e0p(ntp)+im*delta;        %the denominator
  end;
  f0(nt)=singInt(top,deno,esp);    %F0 at the energy E integration
end
plot(e0,real(f0),'k--')            %real part of F0 vesus E
hold on
plot(e0,imag(f0),'k-')             %imaginary part of f0 versus E
axis tight
legend('Jelitto-dos','real(F_0)','Im(F_0)',0)
