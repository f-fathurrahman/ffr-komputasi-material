clear all
close all

%Constants (all MKS, except energy which is in eV)
hbar=1.055e-34;q=1.602e-19;a=3e-10;m=9.110e-31;

%Luttinger-Kohn parameters
g1=6.85;g2=2.1;g3=2.9;%GaAs
w1=(hbar^2)*g1/(2*m*q*(a^2));
w2=(hbar^2)*g2/(2*m*q*(a^2));
w3=(hbar^2)*g3/(2*m*q*(a^2));

g1=3.45;g2=0.68;g3=1.29;%AlAs
a1=(hbar^2)*g1/(2*m*q*(a^2));b1=(.7*w1)+(.3*a1);
a2=(hbar^2)*g2/(2*m*q*(a^2));b2=(.7*w2)+(.3*a2);
a3=(hbar^2)*g3/(2*m*q*(a^2));b3=(.7*w3)+(.3*a3);
Ev=0;Evb=(0.7*0)+(0.3*0.75),kx=0*pi;ky=0*pi;k2=(kx^2)+(ky^2);

for nk=1:20
    Nw = nk+10;
    Nb = Nw;
    Np = Nb + Nw + Nb;
    W(nk) = (Nw-1)*a*1e9;
    Z = zeros(Np);
    X(nk) = Nw-1;
    t=[b1*ones(1,Nb) w1*ones(1,Nw-1) b1*ones(1,Nb)];tt=[0 t]+[t 0];
    Ebk=Evb+(b1*k2);Ewk=(w1*k2);Ebwk=(Ebk+Ewk)/2;
    U=Ev+[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    P=-diag(t,1)-diag(t,-1)+diag(tt)+diag(U);

    t=-2*[b2*ones(1,Nb) w2*ones(1,Nw-1) b2*ones(1,Nb)];tt=[0 t]+[t 0];
    Ebk=b2*k2;Ewk=w2*k2;Ebwk=(Ebk+Ewk)/2;
    U=[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    Q=-diag(t,1)-diag(t,-1)+diag(tt)+diag(U);

    Ebk=-(sqrt(3)*b2*((kx^2)-(ky^2)))+(i*2*b3*sqrt(3)*kx*ky);
    Ewk=-(sqrt(3)*w2*((kx^2)-(ky^2)))+(i*2*w3*sqrt(3)*kx*ky);
    Ebwk=(Ebk+Ewk)/2;
    U=[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    R=diag(U);

    t=2*i*sqrt(3)*(kx-(i*ky))*[b3*ones(1,Nb) w3*ones(1,Nw-1) b3*ones(1,Nb)];
    S=diag(t,1)-diag(t,-1);

    H=[P+Q Z;Z P+Q];HL=[P-Q Z;Z P-Q];
    HC=[-S R;R' S'];    
    H=-[H HC;HC' HL];

    [V,D]=eig(H);D=diag(D);D=-(sort(real(-D)))';
    E1(nk)=D(1);E2(nk)=D(2);E3(nk)=D(3);E4(nk)=D(4);
    E5(nk)=D(5);E6(nk)=D(6);E7(nk)=D(7);E8(nk)=D(8);
end

%Analytical results for infinite well
Ean1=-(w1-(2*w2))*(pi^2)./(X.^2);
Ean2=-(w1+(2*w2))*(pi^2)./(X.^2);

hold on
%h=plot(W,Ean1,'b');
%h=plot(W,Ean2,'b');
h = plot(W,E1,'b');
set(h,'linewidth',[2.0])
%h=plot(W,E2,'bx');
h = plot(W,E3,'b');
set(h,'linewidth',[2.0])
%h=plot(W,E4,'b+');
h = plot(W,E5,'b');
set(h,'linewidth',[2.0])
%h=plot(W,E6,'bx');
h = plot(W,E7,'b');
%h=plot(W,E8,'b+');
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[24])
xlabel('W (nm)')
ylabel('Energy (eV)')
axis([3 9 -0.25 0])
grid on
