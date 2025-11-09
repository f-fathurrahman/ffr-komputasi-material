%copyright by J. E Hasbun and T. Datta
%band_structure_Si.m
%This program uses Harrison's parametrized scaling approach.
%Off diagonal elements scale 1/(bond length)^2
%This is a nearest neighbor tight binding method.
%It's an s-p parametrized approach for 3-5 semiconductors.
%
function band_structure_Si
global H zim NB
%*********** constants ***********
NB=8;          %number of bands (hamiltonian dimension =NBxNB also)
Ns=12;         %number of k steps to do in a given direction
zim=1i;        %the complex number i by itself
H=zeros(NB,NB);%initilize the hamiltonian
%Initialize energies, compound: c=cation, a=anion, es=s-energy, ep=p-energy
%The compound is "system"
[a,esc,esa,epc,epa,ess,esp,exx,exy,system]=initialize_Si();
diagHamil(esc,esa,epc,epa) %diagonal elements of H
%Note: symmetry directions (points): Lambda (L),delta (X), Sigma (K)
sympt=['L';'X';'K'];                 %symmetry points
ba0=2*pi/a;                          %reciprocal lattice vector magnitude
bkpt=ba0*[1/2,1,3./4.];              %corresponding maximum k point
% *******************************************************
% *************** four different directions --loop *******
% *******************************************************
direction=['Lamb';'Delt';'Sigm';'XtoK']; %chosen directions for calculations
for idir=1:4
  %reset maximum k value when the direction is changed
  if strcmp(direction(idir,:),'Delt')       %Delta direction
    a1=1.;
    b1=0.;
    c1=0.;
    ba=ba0;                                 %X point
  else
    if strcmp(direction(idir,:),'Sigm')     %Sigma direction
      a1=1.;
      b1=1.;
      c1=0.;
      ba=ba0*(3./4.);                       %K point
    else
      if strcmp(direction(idir,:),'Lamb')   %Lambda direction
        a1=1.;
        b1=1.;
        c1=1.;
        ba=ba0/2.;                          %L point
      else
        %the vector direction from X to K=[(1,1,0)(3/4)-(1,0,0)](2*pi/a)=[-1/4,3/4,0](2pi/a)
        %we will need to start at X, and move in this direction to end up at the U point
        if strcmp(direction(idir,:),'XtoK') %direction toward K point
          a1=-1./4.;
          b1=3/4;
          c1=0;
          ba=ba0;
        else
          break
        end
      end
    end
  end
  st=ba/Ns;
  fprintf('system: %s\n',system)
  fprintf('k-direction a1=%6.3f, b1=%6.3f, c1=%6.3f, (%s)\n',a1,b1,c1,...
    direction(idir,:))
  fprintf('upper limit ba=%5.2f, step size st =%5.2f\n',ba,st)
  % *******************************************************
  % *************** k--loop ****************
  % *******************************************************
  %the first direction set of bands are done from high k value to zero
  %the second direction set of bands are done from zero to high k value
  if strcmp(direction(idir,:),'Lamb')
    bkL=ba;
    bkS=-st;
    bkU=0.0;
  else
    bkL=0.0;
    bkS=st;
    bkU=ba;
  end
  kc=0;               %reset counter for each direction
  for bk=bkL:bkS:bkU
    kc=kc+1;
    %Notice bx,by,bz includes a factor of 2*pi/a here (see definition of ba above)
    if strcmp(direction(idir,:),'XtoK') %in the XtoK direction, start at X (2pi/a)[1,0,0]
      bx=ba0+bk*a1;%and continue in the XtoK direction (2pi/a)[-1/4,3/4,0] in steps of bk
      by=0+bk*b1;
      bz=0+bk*c1;
    else
      bx=bk*a1;
      by=bk*b1;
      bz=bk*c1;
    end
    %Tight binding interations - nearest neighbors
    offDiagHamil(bx,by,bz,a,ess,esp,exx,exy) %off-diagonal elements of H
    [z,w]=eig(H); %z=eigenvectors, w=eigenvalues
    %let ww contain the sorted eigenvalues and zz contain the sorted eigenvectors
    [ww,zz]=sorter(diag(w),z);
    kval(kc,idir)=bk;
    bandEnergy(:,kc,idir)=ww;
    fprintf('eig.val. at %6.2f %6.2f %6.2f with direction: %s\n',bx,by,bz,direction(idir,:))
    fprintf('%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',ww(1:NB))
  end
end
[tvb,kiv]=max(bandEnergy(NB/2,:,2));    %top of the valence band
[bcb,kic]=min(bandEnergy(NB/2+1,:,2));  %bottom of the conduction band
Eg=bcb-tvb;                             %band gap
fprintf('Top Val. Band tvb =%8.4f %s, and occurs at k=%i\n',tvb,'eV',kval(kiv,2))
fprintf('Bot. Cond. Band bcb =%8.4f %s, and occurs at k=%i\n',bcb,'eV',kval(kic,2))
fprintf('Band gap Eg=bcb-tvb =%8.4f %s\n',Eg,'eV')
%Next, check with Harrison's formulas for the top of the VB and bottom of the CB
%Most of the time, Ev=Ev1 and Ec=Ec1 for 3-5 semiconductors.
%To make Carbon work, Ev2 can be lower than Ec1 so we fixed that here as
%follows:
Ev1=(epc+epa)/2-sqrt(((epc-epa)/2)^2+(4*exx)^2); %Top of VB
Ev2=(epc+epa)/2+sqrt(((epc-epa)/2)^2+(4*exx)^2); %Top of VB
Ec1=(esc+esa)/2+sqrt(((esc-esa)/2)^2+(4*ess)^2); %Bottom of CB
Ec2=(esc+esa)/2-sqrt(((esc-esa)/2)^2+(4*ess)^2); %Bottom of CB
fprintf('Ev1,Ev2,Ec1,Ec2=%8.4f %8.4f %8.4f %8.4f\n',Ev1,Ev2,Ec1,Ec2)
Ev=max(Ev1,Ec2);  %highest of these two roots
Ec=min(Ec1,Ev2);  %lowest of these two roots
EgH=Ec-Ev; %based on Harrison's formulas at the gamma (k=0) point
fprintf('Harrison formulas of the top of the VB and Bottom of the CB, and gap:\n')
fprintf('Ev=%8.4f %s, at Gamma\n',Ev,'eV')
fprintf('Ec%8.4f %s,  at Gamma\n',Ec,'eV')
fprintf('Eg=Ec-Ev =%8.4f %s\n',EgH,'eV')
%Below we shift k by reciprocal lattice vector mag (ba0) to plot correctly
plot(ba0-kval(:,1),bandEnergy(:,:,1),'k')           %Lamda direction
hold on
plot(ba0+kval(:,2),bandEnergy(:,:,2),'k')           %Delta direction
plot(3.5*ba0-kval(:,3),bandEnergy(:,:,3),'k')       %Sigma direction
%max k value at the K point is (3/4)*(2pi/a)+, and shift by recip. lat.
%to plot correctly
plot(2*ba0+(3./4.)*kval(:,4),bandEnergy(:,:,4),'k') %XtoK direction
str1=cat(2,'Energy Band structure versus wavevector k for the ',system,' system');
title(str1,'Fontsize',12)
ylabel('Energy (eV)','Fontsize',12)
%xlabel('(k-space)','Fontsize',8)
emin1=min(bandEnergy(1,:,1));
emin2=min(bandEnergy(1,:,2));
emin3=min(bandEnergy(1,:,3));
emin4=min(bandEnergy(1,:,4));
emax1=max(bandEnergy(NB,:,1));
emax2=max(bandEnergy(NB,:,2));
emax3=max(bandEnergy(NB,:,3));
emax4=max(bandEnergy(NB,:,4));
xmin1=min(ba0-kval(:,1));
xmin2=min(ba0+kval(:,2));
xmin3=min(3.5*ba0-kval(:,3));
xmin4=min(2*ba0+(3./4.)*kval(:,4));
xmax1=max(ba0-kval(:,1));
xmax2=max(ba0+kval(:,2));
xmax3=max(3.5*ba0-kval(:,3));
xmax4=max(2*ba0+(3./4.)*kval(:,4));
emin=min([emin1,emin2,emin3,emin4]); emax=max([emax1,emax2,emax3,emax4]);
xmin=min([xmin1,xmin2,xmin3,xmin4]); xmax=max([xmax1,xmax2,xmax3,xmax4]);
axis ([xmin xmax emin emax])
str2=cat(2,'Eg=',num2str(Eg),' eV');
text(ba0+0.05,tvb*(1-0.05),str2)
text(ba0,min(emin1,emin2)*(1+0.075),'\Gamma')
text(ba0-bkpt(1),emin*(1+0.075),sympt(1,:))               %L point
text(ba0+bkpt(2),emin*(1+0.075),sympt(2,:))               %X point
text(3.5*ba0-bkpt(3),emin*(1+0.075),sympt(3,:))           %K point
text(3.5*ba0,emin*(1+0.075),'\Gamma')
line ([ba0 ba0],[emin emax],'LineStyle','-','Color','k')  %Gamma line
line ([ba0+bkpt(2) ba0+bkpt(2)],[emin emax],...           %X line
  'LineStyle','-','Color','k')
line ([3.5*ba0-bkpt(3) 3.5*ba0-bkpt(3)],[emin emax],...   %K line
  'LineStyle','-','Color','k')
