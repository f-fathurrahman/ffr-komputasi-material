clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.06e-34;
q = 1.6e-19;
eps0 = 8.85E-12;
epsr = 4;
m = 0.25*9.1e-31;
mu = 0;
kT = 0.025;
n0 = m*kT*q/(2*pi*(hbar^2));

n0

% inputs
a = 3e-10;
t0 = (hbar^2)/(2*m*(a^2)*q);
e0 = q*a/eps0;
Nox = 7;
Nc = 10; %use Nc=10,30 for 3,9nm channel respectively
Np = Nox + Nc + Nox;
XX = a*1e9*[1:1:Np];
Ec = [ 3*ones(Nox,1); 0*ones(Nc,1); 3*ones(Nox,1) ];

% Hamiltonian matrix
T = (2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));

% dielectric constant matrix
D2 = epsr*((2*diag(ones(1,Np)))-(diag(ones(1,Np-1),1))-(diag(ones(1,Np-1),-1)));
iD2 = inv(D2);
Vg = 0.25;
Ubdy = -epsr*[Vg; zeros(Np-2,1); Vg];
%Ubdy=-epsr*[0;zeros(Np-2,1);Vg];;%for asymmetric bias
U0 = iD2*Ubdy;

%self-consistent calculation
U1 = 1e-6*ones(Np,1);
UU = U1;
change = 1;
while change > 1e-3
    U1 = U1 + (0.1*(UU - U1));
    [P,D] = eig(T + diag(Ec) + diag(U1));
    D = diag(D);
    rho = log(1 + exp((mu-D)./kT));
    rho = P*diag(rho)*P';
    n = 2*n0*diag(rho);
    for kp = 1:Np
        ncl(kp) = a*2*(n0^1.5)*Fhalf((mu-Ec(kp)-U1(kp))/kT);
    end
    %n=ncl';%use for semiclassical calculation
    UU = U0 + (iD2*e0*n);
    change = max(max((abs(UU-U1))));
    U = Ec + U1;%self-consistent band profile
end

% electron density in channel per cm2
ns = 1e-4*sum( sum( n.*[zeros(Nox,1); ones(Nc,1); zeros(Nox,1)] ) );
Vg, ns
nn = 1e-6*n./a; %electron density per cm3
Fn = mu*ones(Nc+Nox+Nox,1);
hold on
h = plot(XX,nn,'g');
%h=plot(XX,Ec,'g--');
%h=plot(XX,Ec+U1,'g');
%h=plot(XX,Fn,'g:');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [24])
xlabel('z (nm)')
%ylabel('E (eV)')
ylabel('n (/cm3)')
%axis([0 8 -.5 3])
axis([0 8 0 15e18])
grid on
