%copyright by J. E Hasbun and T. Datta
%fcc_dos.m
%Density of states for the single band, face centered cubic system
%according to the ray approach for the tetrahedron method of
%An-Ban Chen, Phys Rev. B V16, 3291 (1977). We also use the
%singInt integration method for the Green's function.
%For the ray method, the 2D vectors on the faces of the thin
%tetrahedrons are obtained from the TTareas function. These vectors are
%average vectors.

function fcc_dos
clear, clc;
global e0 dos
delta=1.5e-3;
im=complex(0.0,1.0);
x=0.05;            %energy step
e2=-3.0;
e1=1.0-delta;      %not to go too close to 1.0
e0=e2:x:e1;        %energy range
ntmax=length(e0);
a=1.0;             %lattice constant
tpa=2*pi/a;
VBZ=4*tpa^3;       %total FCC BZ volume
%FCC case tetrahedrons (three needed) in units of 2*pi/a;
%points used in units of 2*pi/a
L=[1/2,1/2,1/2]; K=[3/4,3/4,0]; U=[1,1/4,1/4]; W=[1,1/2,0]; X=[1,0,0];
%1st tetrahedron (Vectors in units of 2*pi/a)
A(1,:)=L; B(1,:)=K; C(1,:)=W;  %L, K, W points
%tetrahedron volume = 1/32 of the total BZ vol, so use in corresp. integral
VFCC_t1=abs(dot(A(1,:),cross(B(1,:),C(1,:)))/6);
factor(1)=VFCC_t1;
%2nd tetrahedron (Vectors in units of 2*pi/a)
A(2,:)=L; B(2,:)=U; C(2,:)=W;  %L, U, W points
%tetrahedron volume = 1/32 of the total BZ vol, so use in corresp. integral
VFCC_t2=abs(dot(A(2,:),cross(B(2,:),C(2,:)))/6);
factor(2)=VFCC_t2;
%3rd tetrahedron (Vectors in units of 2*pi/a)
A(3,:)=X; B(3,:)=U; C(3,:)=W;  %X, U, W points
%tetrahedron volume = 1/48 of the total BZ vol, so use in corresp. integral
VFCC_t3=abs(dot(A(3,:),cross(B(3,:),C(3,:)))/6);
factor(3)=VFCC_t3;
lim=21; lom=lim-1; %lim must be odd
%number of divisions along q2, and q3 => total number of TT's is n^2
n=10;
%VFCC_t=(VFCC_t1+VFCC_t2+VFCC_t3) %total volume=sum of 3 tetrahedrons
%VFCC_t is the standard tetrahedron volume part, and there are 12 VFCC_t's
%in the face center cube's total BZ (1/32+1/32+1/48=1/12)
factor_t=12*VBZ;         %to be used in the full integral result
%TT=thin tetrahedron, T=tetrahedron
DD=1/n^2;   %ratio of TT area to T area (since all thin triangles are the same)
factor_t=3.0*DD*factor_t; %factor is modified further by 3*DD
par1=0.0;
par2=1.0;
dk=(par2-par1)/lom;
al=par1:dk:par2;   %alpha range to integrate (main TT axis)
dv=1/n;            %TT division size
be=0:dv:1;         %beta
ga=be;             %gamma
%Tetrahedron q vectors according to the method of An-Ban Chen & B. I. Reser.
%Notation: q(tetrahedron(1,2,3),coordinate(x,y,z),vector(1,2,3))
q(:,:,1)=A*tpa; q(:,:,2)=(B-A)*tpa; q(:,:,3)=(C-B)*tpa;
%The function TTareas produces the average vectors "va" on the faces of the TT's
%Notation: va(ith TT vector,coordinate(x,y,z),tetrahedron(1,2,3))
%Ntt=number of TT's, areas(tetrahedron(1,2,3),ith TT area)
for tet=1:3      %the fcc has three tetrahedrons
  [Ntt(tet),areas(tet,:),va(:,:,tet)]=TTareas(q(tet,:,2),q(tet,:,3),n,be,ga);
end
terr=0.0;                  %total error
dos=zeros(1,ntmax);
for nt=1:ntmax
  dos(nt)=0.0;
  for tet=1:3            %the fcc has three tetrahedrons
    dos_tet(tet)=0.0;
    for it=1:Ntt(tet)  %Ntt TT's for the tet tetrahedron
      for ko=1:lim   %loop over gamma
        for v=1:3  %build the k vector, kx=k(1), ky=k(2), kz=k(3)
          %The va vectors are average vectors on the TT faces
          k(v)=al(ko)*(q(tet,v,1)+va(it,v,tet)); %va vectors used here
        end
        top(ko)=al(ko)*al(ko);
        deno(ko)=e0(nt)+cos(k(1)/2)*cos(k(2)/2)...
          +cos(k(1)/2)*cos(k(3)/2)...
          +cos(k(2)/2)*cos(k(3)/2)+im*delta;
      end
      gy=singInt(top,deno,dk);
      dos_tet(tet)=-imag(gy)/pi+dos_tet(tet);  %single TT dos, add all n^2 TT's
    end
    %weigh the corresponding tetrahedron contribution by its factor
    dos(nt)=factor(tet)*dos_tet(tet)+dos(nt);  %dos=sum over 3 tetrahedrons
  end
  %finally, multiply by the total three-tetrahedon contribution factor
  dos(nt)=factor_t*dos(nt)/VBZ;
  ge(nt)=jelittoFccDos(e0(nt))/VBZ;
  err=abs(dos(nt)-ge(nt));
  terr=terr+err;
  fprintf('E0=%9.4f, dos=%9.4f, ge=%9.4f, err=%14.6e\n',...
    e0(nt),dos(nt),ge(nt),err)
end
terr=terr/ntmax;
fprintf('Total error=%14.6e\n',terr)
%Total integrated density of states, and plot
fprintf('Integrated Density of States for the face centered cubic')
intdos=zeros(1,ntmax);
for nt=1:ntmax
  intdos(nt)=rombergInt(e2,e0(nt),@fForRomb);  %integrate on [eL,e0]
  fprintf('E0=%9.4f, integrated dos=%14.6e\n',e0(nt),intdos(nt));
end
plot(e0,ge,'k.'), hold on        %exact dos
plot(e0,dos,'ko','MarkerSize',5) %numeric dos
%Next we do the total density of states without the factor of 2 for spin.
plot(e0,intdos,'k:','LineWidth',2);
xlabel('E (Hartrees)'), ylabel('D(E) (states/energy), N(E) (states)')
str=cat(2,'D(E) and N(E) for the face centered cubic (no spin)');
title(str, 'Fontsize',12)

%legend('Exact D(E)','numeric D(E)','N(E)',0)

function ge=jelittoFccDos(ee)
%Jelitto's DOS for the face centered cubic (exact)
if((ee > -3.0)) & (ee < 0.0)
    eaa=abs(ee);
    a1=3.0+ee;
    a2=a1^2;
    a=sqrt(a1);
    b=-85.9325+101.103*a1;
    d=-16.2885*(a2);
    f=56.8683-47.1215*a1;
    h=2.9045*(a2);
    ge=4.*a*((b+d)+(f+h)*sqrt(eaa));
else
    if((ee >= 0.0) & (ee <= 1.0))
        b11=122.595-19.4100*ee+1.76011*ee^2;
        b22=(-44.8100+7.18628*ee)*log(1.-ee);
        ge=4.0*(b11+b22);
    else
        ge=0.;
    end
end

function y=fForRomb(p)
%function used by romberg integration and which interpolates
%tdos versus e0
global e0 dos
y=interpFunc(p,e0,dos);
