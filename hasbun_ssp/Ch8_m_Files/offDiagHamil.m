%copyright by J. E Hasbun and T. Datta
%offDiagHamil.m
function offDiagHamil(bx,by,bz,dasb,ess,esp,exx,exy)
%c=cation, a=anion, es=s-energy, ep=p-energy
%builds the off-diagonal part of the 8x8 Harrison hamiltonian
%bx,by,bz are the wavevector directions
global H zim NB
do=dasb/4.;   %scale for anion positions is lattice contant/4
cx=cos(bx*do);
cy=cos(by*do);
cz=cos(bz*do);
sx=sin(bx*do);
sy=sin(by*do);
sz=sin(bz*do);
go=4.*(cx*cy*cz-zim*(sx*sy*sz));
g1=4.*(-cx*sy*sz+zim*(sx*cy*cz));
g2=4.*(-sx*cy*sz+zim*(cx*sy*cz));
g3=4.*(-sx*sy*cz+zim*(cx*cy*sz));
H(1,2)=ess*go;
H(1,6)=esp*g1;
H(1,7)=esp*g2;
H(1,8)=esp*g3;
H(2,3)=-esp*conj(g1);
H(2,4)=-esp*conj(g2);
H(2,5)=-esp*conj(g3);
H(3,6)=exx*go;
H(3,7)=exy*g3;
H(3,8)=exy*g2;
H(4,6)=exy*g3;
H(4,7)=exx*go;
H(4,8)=exy*g1;
H(5,6)=exy*g2;
H(5,7)=exy*g1;
H(5,8)=exx*go;
for i=1:NB-1
  for j=i+1:NB
    H(j,i)=conj(H(i,j)); %hermitian matrix
  end
end
