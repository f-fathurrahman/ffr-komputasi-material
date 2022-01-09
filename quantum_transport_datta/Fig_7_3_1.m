clear all
close all

%Constants (all MKS, except energy which is in eV)
hbar = 1.06e-34;
q = 1.6e-19;
eps0 = 8.85E-12;
epsr = 4;
m = 0.25*9.1e-31;
mu = 0;
kT = .025;
n0 = m*kT*q/(2*pi*(hbar^2));

%inputs
a = 3e-10;
t0 = (hbar^2)/(2*m*(a^2)*q);
e0 = q*a/eps0;
Nox = 7;
Nc = 10;%use Nc=10,30 for 3,9nm channel respectively
Np = Nox + Nc + Nox;
XX = a*1e9*[1:1:Np];
Ec = [3*ones(Nox,1); 0*ones(Nc,1); 3*ones(Nox,1)];

%Hamiltonian matrix
T = (2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));

%dielectric constant matrix
D2 = epsr*((2*diag(ones(1,Np)))-(diag(ones(1,Np-1),1))-(diag(ones(1,Np-1),-1)));
iD2 = inv(D2);

Vg = linspace(-0.25, 0.25, 26);
for kg = 1:26
    Ubdy = -epsr*[Vg(kg); zeros(Np-2,1); Vg(kg)];
    %Ubdy=-epsr*[0;zeros(Np-2,1);Vg(kg)];;%for asymmetric bias

    U0 = iD2*Ubdy;

    %self-consistent calculation
    U1 = 1e-6*ones(Np,1);
    UU = U1;
    change = 1;
    while change > 1e-3
        U1 = U1 + (0.1*(UU-U1));
        [P,D] = eig(T + diag(Ec) + diag(U1));
        D = diag(D);
        rho=log(1+exp((mu-D)./kT));rho=P*diag(rho)*P';
        n=2*n0*diag(rho);
        for kp=1:Np
            ncl(kp) = a*2*(n0^1.5)*Fhalf((mu-Ec(kp)-U1(kp))/kT);
        end
        %n=ncl';%use for semiclassical calculation
        UU = U0 + (iD2*e0*n);
        change = max(max((abs(UU-U1))));
        U = Ec + U1; %self-consistent band profile
    end
    %electron density in channel per cm2    
    ns(kg)=1e-4*sum(sum(n.*[zeros(Nox,1);ones(Nc,1);zeros(Nox,1)]));
    nn(:,kg)=1e-6*n./a;%electron density per cm3
    Fn(:,kg)=mu*ones(Nc+Nox+Nox,1);
end

C = q*(ns(26)-ns(25))/(Vg(26)-Vg(25));
d = 1e-4*epsr*eps0*2/C;

d
C

%ns=log10(ns)

hold on
h = plot(Vg,ns,'b');
set(h, 'linewidth',[2.0])
set(gca, 'Fontsize',[24])
xlabel('Vg (V)')
ylabel('ns (/cm2)')
%axis([0 .3 0 3.5e12])
grid on