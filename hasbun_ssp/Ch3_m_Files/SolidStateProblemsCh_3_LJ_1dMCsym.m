%copyright by J. E Hasbun and T. Datta
%LJ_1dMCsym.m
%Program to do a simulation of the distance between two
%particles in 1 d using the Lennard-Jones potential,
%where U=4*epsilon*((sigma/R)^12-(sigma/R)^6)
clear
%We work with Argon
eps=0.0104;           %eV
sig=3.4;              %Angstroms
%Define the LJ potential as an anonymous function
U=@(eps,sig,r) 4*eps*((sig./r).^12-(sig./r).^6);
Rguess=2.0*sig;
Rold=Rguess;          %guess starting distance in terms of sigma
Uold=U(eps,sig,Rold); %starting energy
cmod=0.25;            %moderator for the random number used below
tol=1.e-3*eps;        %convergence tolerance
itermax=500;          %maximum iterations
diffU=10*tol;
plot(0,0,'ko','MarkerSize',5,'MarkerFaceColor','k')
hold on
axis([-1 1.5*Rguess -1 1])
plot(Rold,0,'bo','MarkerSize',5,'MarkerFaceColor','b')
iter=0;
while (diffU > tol & iter < itermax)
  iter=iter+1;
  cla
  h(1)=plot(0,0,'ko','MarkerSize',5,'MarkerFaceColor','k');
  rn=-1.0+2.0*rand(); %random number  -1 < r < 1
  Rnew=Rold*(1+cmod*rn);   %new R based on rand #, with moderation
  Unew=U(eps,sig,Rnew);    %new energy
  diffU=abs(Unew-Uold)/abs(Unew+Uold);
  if (Unew < Uold)
    Rold=Rnew;
    Uold=Unew;
    diffUold=diffU;
    h(2)=plot(Rold,0,'ko','MarkerSize',5);
    pause(0.05)
  end
end
h(2)=plot(Rold,0,'ko','MarkerSize',5);
h(3)=line([0,Rold],[0,0],'LineStyle','--','Color','k');
fprintf('diffU=%8.3e, Min U=%9.6f (eV), R=%9.6f (A), iter=%6i\n',...
    diffUold,Uold,Rold,iter)
Ns=500;
rl=0.1*Rold;
ru=3.0*Rold;
rs=(ru-rl)/(Ns-1);
r=rl:rs:ru;
Ur=U(eps,sig,r);
Umin=min(Ur);
h(4)=plot(r,Ur,'k-');
axis([-1 1.5*Rguess Umin abs(Umin)])
xlabel('r (\AA)','interpreter','latex')
ylabel('U (eV)')
str=cat(2,'Monte-Carlo Lennard-Jones U(r) Simulation:',' Umin=',...
  num2str(Umin,'%9.4g'),'eV, R=',num2str(Rold,'%9.4g'),' Angstroms');
title(str)
legend(h,'origin atom','other atom','common bond','L-J: Potential')
