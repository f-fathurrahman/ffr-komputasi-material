%copyright by J. E Hasbun and T. Datta
%CPA.m
%We calculate the coherent potential approximation (CPA) sigma.
%The virtual crystal approximation (VCA) is used as guess
%The resulting sigma is used to obtain the density of states
%for the alloy.
function CPA
clear; clc;
delta=1.e-3;
im=complex(0.0,1.0);
Ea=-2.5;           %alloy system species 1 diagonal energy
Eb=2.5;            %alloy system species 2 diagonal energy
scbw=3.0;          %simple cubic known band width
x=1;
y=1-x;             %y=concentration of Eb species
e2=Ea-scbw;        %energy range to work with
e1=-e2;
ntmax=201;
es=(e2-e1)/(ntmax-1);
a=1.0;             %lattice constant
tpa=2*pi/a;
VBZ=tpa^3;         %total SC BZ volume
e2p=3.0;
e1p=-e2p;
ntpmax=251;
esp=(e2p-e1p)/(ntpmax-1);
Evca=x*Ea+y*Eb;                     %VCA energy
sig_guess=Evca;                     %initial sigma guess - VCA energy
cv=0.85;                            %new guess helper parameter
ncmax=25;                           %maximum iterations
tol=1.e-4;
str1=cat(2,'nt=%4i, e0=%4.2f, nc=%2i, sig=%5.4g, drs=%5.4g');
str2=cat(2,', dis=%5.4g, dos=%5.4f\n');
str=cat(2,str1,str2);
for nt=1:ntmax                      %Main energy loop
  nc=0;
  converge=0;
  drs=10*tol;
  dis=10*tol;
  e0(nt)=e1+(nt-1)*es;              %the energy E
  sig(nt)=sig_guess;                %use former sigma guess
  while (converge==0 & nc < ncmax)
    nc=nc+1;
    if(drs <  tol & dis < tol)
      converge=1;            %use the converged sigma as guess
      sig_guess=sig(nt);     %for next energy
    end
    %Integrate over the E' loop
    for ntp=1:ntpmax
      if(nt==1)                   %calculate this part only once
        e0p(ntp)=e1p+(ntp-1)*esp;
        top(ntp)=jelittoScDosAnal(e0p(ntp))/VBZ;
      end
      deno(ntp)=e0(nt)-e0p(ntp)-sig(nt)+im*delta;
    end;
    f0(nt)=singInt(top,deno,esp);
    %CPA rule for sigma's next guess. The 2nd form converges better.
    %signew=Evca-(Ea-sig(nt))*f0(nt)*(Eb-sig(nt));
    signew=(Evca-(Ea*Eb-Evca*sig(nt))*f0(nt))/...
      (1-(y*Ea+x*Eb-sig(nt))*f0(nt));
    %The actual guess we make is a mixture between new and old
    sig(nt)=cv*signew+(1-cv)*sig(nt);
    drs=abs(real(sig(nt)-signew));    %use for convergence criteria
    dis=abs(imag(sig(nt)-signew));
  end
  dos(nt)=-imag(f0(nt))/pi;
  fprintf(str,nt,e0(nt),nc,sig(nt),drs,dis,dos(nt));
end
str=cat(2,'Coherent Potential Approximation (CPA), \epsilon_A=',...
  num2str(Ea,'%5.2f'),'H_a, \epsilon_B=',num2str(Eb,'%5.2f'),...
  'H_a, x=',num2str(x,'%4.2f'));
title(str)
plot(e0,dos,'k')
xlabel('E (H_a)'), ylabel('D(E) (states/energy)')
title(str)
