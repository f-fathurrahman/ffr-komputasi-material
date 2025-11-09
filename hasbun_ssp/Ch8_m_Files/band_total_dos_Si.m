%copyright by J. E Hasbun and T. Datta
%band_total_dos_Si.m
%Here, we also find the total integrated density of states
%versus energy. This program uses Harrison's s-p parametrized
%scaling approach for 3-5 semiconductors. Off diagonal elements
%scale 1/(bond length)^2. This is a nearest neighbor tight binding
%method.
function band_total_dos_Si
global H zim NB
global e0 tdos
%*********** constants ***********
NB=8;          %number of bands (hamiltonian dimension =NBxNB also)
H=zeros(NB,NB);%initilize the hamiltonian
delta=1.5e-2;
zim=complex(0.0,1.0);
%a=lattice constant
%Initialize energies, compound: c=cation, a=anion, es=s-energy, ep=p-energy
%The compound is "system"
[a,esc,esa,epc,epa,ess,esp,exx,exy,system]=initialize_Si();
%x=0.25;            %energy step
%eL=-23.0; eU=5.0;  %lowest and highest energy from the band structure
disp('Return to accept the default values within brackets')
disp('Once the calculation begins, wait for the results')
eL =input('  low energy limit [-23] = ');
if(isempty(eL)), eL=-23.0; end
eU =input('  high energy limit [5] = ');
if(isempty(eU)), eU=5.0; end
x =input('  energy step [0.25] = ');
if(isempty(x)), x=0.25; end
e0=eL:x:eU;        %energy range
ntmax=length(e0);
diagHamil(esc,esa,epc,epa) %diagonal elements of H
tpa=2*pi/a;
VBZ=4*tpa^3;    %total FCC BZ volume
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
%Initialize arrays
dos=zeros(NB,ntmax);
dos_tet=zeros(NB,3);
top=zeros(1,lim);
deno=zeros(NB,lim);
gy=zeros(1,NB);
tdos=zeros(1,ntmax);
%Main Loop
fprintf('system: %s\n',system)
for nt=1:ntmax
  for nB=1:NB            %initialize each band dos at each energy
    dos(nB,nt)=0.0;
  end
  for tet=1:3            %the fcc has three tetrahedrons
    for nB=1:NB        %initialize each tetrahedron dos for each band
      dos_tet(nB,tet)=0.0;
    end
    for it=1:Ntt(tet)  %Ntt TT's for the tet tetrahedron
      for ko=1:lim   %loop over gamma
        for v=1:3  %build the k vector, kx=k(1), ky=k(2), kz=k(3)
          %The va vectors are average vectors on the TT faces
          k(v)=al(ko)*(q(tet,v,1)+va(it,v,tet)); %va vectors used here
        end
        top(ko)=al(ko)*al(ko);
        %Tight binding interations - nearest neighbors
        offDiagHamil(k(1),k(2),k(3),a,ess,esp,exx,exy) %off-diagonal els of H
        [z,w]=eig(H); %z=eigenvectors, w=eigenvalues
        %ww contains the sorted eigenvalues and zz  the sorted eigenvectors
        [ww,zz]=sorter(diag(w),z);
        deno(:,ko)=e0(nt)-ww(:)+zim*delta; %w-Ek for the NB energy bands
      end
      for nB=1:NB
        gy(nB)=singInt(top,deno(nB,:),dk);    %integrate each band over k
        dos_tet(nB,tet)=-imag(gy(nB))/pi+dos_tet(nB,tet);%single TT dos,  TT's
      end
    end
    %for each band weigh the corresponding tetrahedron contribution by its factor
    for nB=1:NB
      dos(nB,nt)=factor(tet)*dos_tet(nB,tet)+dos(nB,nt);  %dos=sum over tetraheds
    end
  end
  %finally, for each band multiply by the three-tetrahedon contribution factor
  for nB=1:NB
    dos(nB,nt)=factor_t*dos(nB,nt);
    %fprintf('E0=%9.4f, dos=%9.4f\n',e0(nt),dos(nB,nt)
  end
end
%
%total dos sum over all the bands
fprintf('Density of States for the system: %s\n',system)
for nt=1:ntmax
  tdos(nt)=0.0;
  for nB=1:NB
    tdos(nt)=tdos(nt)+dos(nB,nt);
  end
  %Total dos=integral of d^3k over the BZ / BZ volume
  tdos(nt)=tdos(nt)/VBZ;
  fprintf('E0=%9.4f, tdos=%14.6e\n',e0(nt),tdos(nt));
end
%
%For the Harrison model, we know the band edges, and gap
Ev1=(epc+epa)/2-sqrt(((epc-epa)/2)^2+(4*exx)^2); %Top of VB
Ev2=(epc+epa)/2+sqrt(((epc-epa)/2)^2+(4*exx)^2); %Top of VB
Ec1=(esc+esa)/2+sqrt(((esc-esa)/2)^2+(4*ess)^2); %Bottom of CB
Ec2=(esc+esa)/2-sqrt(((esc-esa)/2)^2+(4*ess)^2); %Bottom of CB
disp('Energies in eV ')
fprintf('Ev1,Ev2,Ec1,Ec2=%8.4f %8.4f %8.4f %8.4f\n',Ev1,Ev2,Ec1,Ec2)
Ev=max(Ev1,Ec2);  %highest of these two roots
Ec=min(Ec1,Ev2);  %lowest of these two roots
EgH=Ec-Ev; %gap - based on Harrison's formulas at the gamma (k=0) point
fprintf('Ev=%9.4f, Ec=%9.4f, Eg=%9.4f\n',Ev,Ec,EgH);
indEv=1+floor((Ev-eL)/x);
indEc=1+floor((Ec-eL)/x);
strv=cat(2,'Val. band edge: Ev=',num2str(Ev,'%9.4f'),' eV');
strc=cat(2,'Cond. band edge: Ec=',num2str(Ec,'%9.4f'),' eV');
strg=cat(2,'Eg=',num2str(EgH,'%9.4f'),' eV');
mxd=max(tdos);
%
%plot the total density of states versus energy
figure(1)
p0(1)=plot(e0,tdos,'k-');
hold on

%
% %Back to figure 1 and add the new information
% figure(1)

%place a line at Ev and Ec
p0(2)=line([Ev Ev],[0 mxd*0.1],'LineStyle','--','color','k','LineWidth',2);
p0(3)=line([Ec Ec],[0 mxd*0.1],'LineStyle','-.','color','k','LineWidth',2);
pL=legend(p0,'DOS',strv,strc,0);
set(pL,'FontSize',9);
text(Ev-abs(Ev)*0.1,mxd*0.5,strg,'FontSize',9)
xlabel('E(eV)'), ylabel('DOS (eV^{-1})')
str=cat(2,'Density of States for the ',system,' system');
title(str)
hold off
pause(1)
%
%Total integrated density of states, and plot
fprintf('Integrated Density of States for the system: %s\n',system)
intdos=zeros(1,ntmax);
for nt=1:ntmax
  intdos(nt)=rombergInt(eL,e0(nt),@fForRomb);  %integrate on [eL,e0]
  fprintf('E0=%9.4f, integrated dos=%14.6e\n',e0(nt),intdos(nt));
end
%
figure(2)
hold on
p1(1)=plot(e0,intdos,'k-');
%place a line at Ev and Ec of respective heights tdos(indEv) & tdos(indEc)
p1(2)=line([Ev Ev],[0 intdos(indEv)],'LineStyle','--','color','k','LineWidth',2);
p1(3)=line([Ec Ec],[0 intdos(indEc)],'LineStyle','-.','color','k','LineWidth',2);
pL=legend(p1,'integated DOS',strv,strc,0);
set(pL,'FontSize',9);
text(Ec+0.1*abs(Ec),intdos(indEc)*(1-0.3),strg,'FontSize',9)
xlabel('E (eV)'), ylabel('Integrated DOS (number of states)')
str=cat(2,'Integrated Density of States for the ',system,' system');
title(str)
hold off

function y=fForRomb(p)
%function used by romberg integration and which interpolates
%tdos versus e0
global e0 tdos
y=interpFunc(p,e0,tdos);
