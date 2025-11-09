%copyright by J. E Hasbun and T. Datta
%Green_1D.m
%It calculates the real and imaginary parts of a one-dimensional
%Green's function for a nearest neighbor tight binding model
%where the band energy Ek=-2*gamma*cos(ka)
clear
el=-2;
eu=2;
Ne=100;
es=(eu-el)/(Ne-1);
gam=1/2;                       %band energy parameter (in Ha)
zim=complex(0.,1.0);
delta=1.e-3;                   %for the small imaginary part
for i=1:Ne
  e0(i)=el+(i-1)*es;
  g00(i)=e0(i)/(sqrt((e0(i)+zim*delta)^2-4*gam^2)*abs(e0(i)));
end
plot(e0,real(g00),'k-.')
hold on
plot(e0,imag(g00),'k:','LineWidth',2)
axis tight
legend('Real(G_{00})','Im(G_{00})',0)
xlabel('E (Ha)','FontSize',12)
ylabel('G_{00} (1/Ha)','FontSize',12)
title('One-Dimensional Green''s Function','FontSize',12)
