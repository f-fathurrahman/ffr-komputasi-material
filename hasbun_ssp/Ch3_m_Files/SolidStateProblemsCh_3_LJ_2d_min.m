%copyright by J. E Hasbun and T. Datta
%LJ_2d_min.m
%This program finds the equilibrium distance between three particles
%in 2d. It does so by minimizing the Lennard-Jones potential
%U=4*epsilon* [Sum of i,j > i of ((sigma/Rij)^12-(sigma/Rij)^6)]
function LJ_2d_min
clear
global eps sig
e=1.602176487e-19;      %electronic charge
eps=5.0e-22;            %for Neon in Joules
eps=eps/e;              %for Neon in eV
sig=2.74;               %in Angstroms
np=3;                   %Number of particles
nd=2;                   %number of dimensions
r_guess(1,1)=0;         %x - origin particle is fixed
r_guess(1,2)=0;         %y
r_guess(2,1)=2;         %x - 2nd particle
r_guess(2,2)=0;         %y
r_guess(3,1)=2;         %x - 3rd particle
r_guess(3,2)=2;         %y
[rnew,Umin]=fminsearch(@LJ_funNd,r_guess,[],np,nd);
%show the final positions of the particles
hold on
axis([-1*max(abs(rnew(:,1))) 2*max(abs(rnew(:,1))) ...
  -2*max(abs(rnew(:,2))) 2*max(abs(rnew(:,2)))])
for i=1:np
  plot(rnew(i,1),rnew(i,2),'ko','MarkerSize',10,'MarkerFaceColor','k')
end
fprintf('Umin=%9.6g (eV)\n',Umin)
%Final particle distances
rave=0.0;
icord=0;
for i=1:np
  fprintf(' p#%2i x=%9.4f, y=%9.4f (A)\n',i,rnew(i,1),rnew(i,2))
end
for i=1:np-1
  for j=i+1:np
    icord=icord+1;
    rij=((rnew(i,1)-rnew(j,1))^2+(rnew(i,2)-rnew(j,2))^2)^(1/2);
    rave=rij+rave;
    fprintf(' r(%2i,%2i)=%9.4f (A)\n',i,j,rij)
    line([rnew(i,1),rnew(j,1)],[rnew(i,2),rnew(j,2)],'LineStyle','--',...
      'Color','k','LineWidth',3)
  end
end
axis tight
rave=rave/icord;
fprintf(' average r=%9.4f (A)\n',rave)
xlabel('x (\AA)','interpreter','latex')
ylabel('y (\AA)','interpreter','latex')
str=cat(2,'3 atom Molecule 3D',' Umin=',...
     num2str(Umin,'%9.4g'),'eV, r_{ave}=',num2str(rave,'%9.4g'),'A');
title(str)

function U_LJ=LJ_funNd(r,Np,Nd)
%r is a Np x Nd vector
%Np equal number of particles and
%Nd=number of dimensions
global eps sig
%The Lennard-Jones potential evaluated at vector r
U_LJ=0.0;
sig12=sig^12;
sig6=sig^6;
for i=1:Np-1
  for j=i+1:Np
    rij=0.0;
    for k=1:Nd
      rij=(r(i,k)-r(j,k))^2+rij;
    end
    U_LJ=(sig12*(1./rij).^6-sig6*(1./rij).^3)+U_LJ;
  end
end
U_LJ=4*eps*U_LJ;
