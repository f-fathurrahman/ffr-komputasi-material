%copyright by J. E Hasbun and T. Datta
%Green_1Totdos.m
%Uses the density of states from the analytic Green's function
%for the one-dimensional tight binding model and integrates it
%to obtain the total density of states N(E). The Romberg method
%is used to perform the integration over the density of states D(E).
%Romberg uses the interpolated D(E) to produce N(E).
function Green_1Totdos
clear
global e0 dosa                 %needed for interpolation
el=-2;                         %upper energy value
eu=2;                          %lower energy value
Ne=100;                        %number of energy points
es=(eu-el)/(Ne-1);             %energy step
gam=1/2;                       %band energy parameter (in Ha)
zim=complex(0.,1.0);           %imaginary number
delta=1.e-4;                   %small part for plottting, integrating
for i=1:Ne
  e0(i)=el+(i-1)*es;
  g00=e0(i)/(sqrt((e0(i)+zim*delta)^2-4*gam^2)*abs(e0(i)));
  dosa(i)=-imag(g00)/pi;       %analytic density of states
end
%Total density of states and plots. Integrate with Romberg's method.
x=0.1;           %energy step
eT=el:x:eu;      %energy range (using previous e1, eu limits)
for nt=1:length(eT)
  intdos(nt)=rombergInt(el,eT(nt),@fForRomb);  %integrate on [el,eT]
  fprintf('E0=%9.4f, integrated dos=%14.6e\n',eT(nt),intdos(nt));
end
%density of states (dosa)
plot(e0,dosa,'k')
hold on
%total dos multiplied by 2 to include electron spin
plot(eT,2*intdos,'k:','LineWidth',2)
legend ('D(E)','2*N(E)','Location','North')
xlabel('E (Ha)')
ylabel('D(E) (1/Ha) and 2*N(E)')
title('Density of States and 2*Total Density of States')

function y=fForRomb(p)
%Function used by rombergInt integration and which interpolates
%the density of states versus e0.
global e0 dosa            %variables passed from the main program
y=interpFunc(p,e0,dosa);  %the Langrange interpolator for Romberg
