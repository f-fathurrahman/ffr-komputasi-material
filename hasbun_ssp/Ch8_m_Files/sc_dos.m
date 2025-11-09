%copyright by J. E Hasbun and T. Datta
%sc_dos.m
%Density of states for the single band, simple cubic system
%according to the ray approach for the tetrahedron method of
%An-Ban Chen, Phys Rev. B V16, 3291 (1977). We also use the
%singInt integration method for the Green's function.
%For the ray method, the 2D vectors on the faces of the thin
%tetrahedrons are obtained from the TTareas function. These vectors are
%average vectors.
function sc_dos
clear, clc;
global e0 dos
delta=1.5e-2;
im=complex(0.0,1.0);
x=0.1;             %energy step
e2=3.0;
e1=-e2;
e0=e1:x:e2;        %energy range
ntmax=length(e0);
a=1.0;             %lattice constant
tpa=2*pi/a;
VBZ=tpa^3;         %total SC BZ volume
%SC case tetrahedron (only one needed) in units of 2*pi/a;
X=[1/2,0,0]; R=[1/2,1/2,1/2]; M=[1/2,1/2,0]; %tetrahedron symmetry points
%Tetrahedron q vectors according to the method of An-Ban Chen & B. I. Reser.
q(1,:)=R*tpa; q(2,:)=(X-R)*tpa ; q(3,:)=(M-X)*tpa;
%total tetrahedron volume
VSC_t=abs(dot(X,cross(R,M))/6)*VBZ;
lim=9; lom=lim-1; %lim must be odd
%number of divisions along q2, and q3 => total number of TT's is n^2
n=25;
%VSC_t is the standard tetrahedron volume, and there are 48 of them in
%the SC cube's total BZ
factor=VSC_t*48;      %initital factor used in the full integral result
%TT=thin tetrahedron, T=tetrahedron
DD=1/n^2;   %ratio of TT area to T area (since all thin triangle are the same)
factor=3.0*DD*factor; %factor is modified further by 3*DD
par1=0.0;
par2=1.0;
dk=(par2-par1)/lom;
al=par1:dk:par2;   %alpha range to integrate (main TT axis)
dv=1/n;            %TT division size
be=0:dv:1;         %beta
ga=be;             %gamma
%The function TTareas produces the average vectors "va" on the faces of
%the TT's. Notation: va(ith TT vector,coordinate(x,y,z))
%Ntt=number of TT's, areas(ith TT area)
[Ntt,areas,va]=TTareas(q(2,:),q(3,:),n,be,ga);
terr=0.0;          %total error
dos=zeros(1,ntmax);
for nt=1:ntmax
  dos(nt)=0.0;
  for it=1:Ntt       %Ntt TT's
    for ko=1:lim   %loop over alpha
      for v=1:3  %build the k vector, kx=k(1), ky=k(2), kz=k(3)
        %The va vectors are average vectors on the TT faces
        k(v)=al(ko)*(q(1,v)+va(it,v)); %va vectors used here
      end
      top(ko)=al(ko)*al(ko);
      deno(ko)=e0(nt)+cos(k(1))+cos(k(2))+cos(k(3))+im*delta;
    end
    gy=singInt(top,deno,dk);
    dos(nt)=-imag(gy)/pi+dos(nt);  %single TT dos, add all n^2 TT's
  end
  dos(nt)=factor*dos(nt)/VBZ;      %numeric dos
  ge(nt)=jelittoScDos(e0(nt))/VBZ; %exact dos
  err=abs(dos(nt)-ge(nt));         %the error per energy
  terr=terr+err;                   %cumulativs error
  fprintf('E0=%9.4f, dos=%9.4f, ge=%9.4f, err=%14.6e\n',...
      e0(nt),dos(nt),ge(nt),err)
end
terr=terr/ntmax;
fprintf('Total error=%14.6e\n',terr)
%Total integrated density of states, and plot
fprintf('Integrated Density of States for the simple cubic')
intdos=zeros(1,ntmax);
for nt=1:ntmax
  intdos(nt)=rombergInt(e1,e0(nt),@fForRomb);  %integrate on [e1,e0]
  fprintf('E0=%9.4f, integrated dos=%14.6e\n',e0(nt),intdos(nt));
end
plot(e0,ge,'k'), hold on         %exact dos
plot(e0,dos,'ko','MarkerSize',5) %numeric dos
%Next we do the total density of states without the factor of 2 for spin.
plot(e0,intdos,'k:','LineWidth',2);
xlabel('E (Hartrees)'), ylabel('D(E) (states/energy), N(E) (states)')
str=cat(2,'D(E) and N(E) for the simple cubic (no spin)');
title(str, 'Fontsize',12)
legend('Exact D(E)','numeric D(E)','N(E)',0)

function ge=jelittoScDos(ee)
%Jelitto's DOS for the simple cubic (exact)
eaa=abs(ee);
if (eaa <= 3.) & (eaa >=  1.)
  a1=3.-eaa;
  a2=a1^2;
  a=sqrt(a1);
  b=80.3702-16.3846*a1;
  d=0.78978*(a2);
  f=-44.2639+3.66394*a1;
  h=-0.17248*(a2);
  ge=a*((b+d)+(f+h)*sqrt(eaa-1.));
else
  if(eaa < 1.)
    ge=70.7801+1.0053*(ee^2);
  else
    ge=0.;
  end
end

function y=fForRomb(p)
%function used by romberg integration and which interpolates
%tdos versus e0
global e0 dos
y=interpFunc(p,e0,dos);
