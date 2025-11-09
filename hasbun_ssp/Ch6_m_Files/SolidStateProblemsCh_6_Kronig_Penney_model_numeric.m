%copyright by J. E Hasbun and T. Datta
%Kronig_Penney_model_numeric.m
%Solves for the energy eigenvalues versus k' in the kronig-Penney
%model of a one dimensional crystal

function Kronig_Penney_model_numeric
clear
global fL fR u0 a b
h=6.62606896e-34;           %Planck'constant (J.s)
hbar=h/2./pi;               %hbar (J.s)
me=9.10938215e-31;          %electron mass (kg)
e=1.602176487e-19;          %electronic charge
ab=1e-10;                   %1 angstom unit of distance
kb=1/ab;                    %wevevector unit
Eb=hbar^2*kb^2/2/me;        %energy unit in joules
Eb_eV=Eb/e;                 %energy unit in eV
u0=150;                     %potential height in Eb units
a=2; b=0.025;               %well, barrier widths in ab units
kp=(0:0.005:3.0)*pi/(a+b);  %vary k'
k=(-6:0.05:6)*pi/a;         %vary k (associated with energy)
%right hand side of energy equation (x=a, y=b)
fR=@(kv,x,y) cos(kv*(x+y));
%left hand side of energy equation (x=a, y=b)
fL=@(xE,V,x,y) ((V-2*xE)./(2*sqrt(u0-xE).*sqrt(xE))).*...
  sinh(sqrt(V-xE)*y).*sin(sqrt(xE)*x)+cosh(sqrt(V-xE)*y).* ...
  cos(sqrt(xE)*x);
%Evaluation and plotting (yR=RHS, yL=LHS)
yR=fR(kp,a,b);
yL=fL(k.^2,u0,a,b);  %Note: in dimensionless units energy=k^2
subplot(1,2,1)       %LHS and min(RHS), max(RHS) plotted versus k
line([min(k*a) max(k*a)],[max(yR), max(yR)],'Color','k',...
  'LineStyle','--')  %min(RHS)
hold on
plot(k*a,yL,'k')     %LHS
for i=1:length(k)
  if ((yL(i) <=  max(yR)) & (yL(i) >=  min(yR)))
    plot(k(i)*a,yL(i),'k.','MarkerSize',5)
  end
end
%legend('RHS','LHS','roots',1)
line([min(k*a) max(k*a)],[min(yR), min(yR)],'Color','k',...
  'LineStyle','--')  %max(RHS)
xlabel('ka','FontSize',14)
ylabel('LHS, RHS','FontSize',14)
str1=cat(2,'Kronig-Penney Model: u0=',num2str(u0,'%5.2f'),...
  ' \epsilon_b',', a=',num2str(a,'%5.3f'),' a_b',', b=',...
  num2str(b,'%5.3f'),' a_b, ');
title(['         ',str1])
axis tight
hold off
%
subplot(1,2,2)       %Energy eigenvalues versus k_prime
%find the lowest energy guess for the first kp point
eps_guess=searchguess(1.e-3,kp(1));
nr=0;               %root counter
test_old=0.0;       %variable to check when kp*pi/(a+b)=integer
for i=1:length(kp)
  yRi=fR(kp(i),a,b);
  if (yRi <= 1 & yRi >= -1)         %if solutions exist
    nr=nr+1;                        %root counter
    kr(nr)=kp(i);                   %store the related k'
    %find the energies for which FofE=fL-fR=0, for each k'
    eps(nr) = fzero(@(xE) FofE(xE,kp(i)),eps_guess);
    eps_guess=eps(nr);              %use this as next guess
    test_new=mod(kp(i)*(a+b)/pi,1); %to check if k'=integer
    %There is an energy gap when kp*pi/(a+b)=integer, so we need
    %to search for a higher energy guess at those points
    if (test_new < test_old)
      eps_guess=searchguess(eps_guess,kp(i));
      %fprintf(' ********** xguess=%8.3f\n',xguess);
    end
    test_old=test_new;
  end
end
%
%use the calculated values and apply symmetry to get the whole
%BZ plot for the energy bands
kp_a=kr*(a+b);         %k' variable where a root was found
ep_a=eps;              %corresponding found energy root
kp_b=kp_a-kp_a(end);   %get the reflection of k'
ep_b=ep_a(end:-1:1);   %use reflection on energy also
plot(kp_a,ep_a,'k.','MarkerSize',1)
hold on
plot(kp_b,ep_b,'k.','MarkerSize',1)
xlabel('k{''}(a+b)','FontSize',14)
ylabel('\epsilon(k'') (\epsilon_b)','FontSize',14)
str2=cat(2,'\epsilon_b=',num2str(Eb_eV,'%5.2f'),...
  ' eV',', a_b=',num2str(ab/1e-10,'%3.2f'),...
  ' Angs');
title(['            ',str2])

function [y]=searchguess(x,kk)
%This function does a simple search for a root by tracking
%a sign change in the FofE function
global u0 a
del=13.2e-3*u0/a+0.2*x;    %step size to search for a root
                           %notice it depends on u0, a - may tweek
xi=x;
xf=xi;
for i=1:50*del
  xf=xf+del;
  if(FofE(xi,kk)*FofE(xf,kk) <= 0.0)
    y=xf;         %root exits, so return as guess
    return
  end
  xi=xf;
end

function [y]=FofE(eps,kk)
global fL fR u0 a b
%This function is the difference of left and right functions
y=fL(eps,u0,a,b)-fR(kk,a,b);
