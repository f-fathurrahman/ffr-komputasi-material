clear all
close all

%Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
q = 1.602e-19;
mu = 0.25;
kT = 0.025;% 0.025 for Part (c),(e) and 0.0025 for Part (d),(f)

% Lattice
Np = 100;
a = 2e-10;
X = a*[1:1:Np];t0=(hbar^2)/(2*m*(a^2))/q;
U = linspace(0,0,Np);

T = (2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));
T(1,Np) = -t0;
T(Np,1) = -t0; %Periodic boundary conditions for Parts (d), (f)
U(Np/2) = U(Np/2)+10; %Impurity potential with Parts (d), (f)

[V,D] = eig(T+diag(U));
E = sort(diag(D)');

D = diag(D) - mu;
rho = 1./(1+exp(D./kT));
rho = V*diag(rho)*V';
rho = diag(rho)/a;

hold on
grid on
%h=plot(E,'b');h=plot(mu*ones(Np/2,1),'b');% Part (b)
h = plot(X,rho,'b'); % Part (c)-(f)
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
grid on
%xlabel('Eigenvalues number');% Part (b)
%ylabel('Energy (eV)');% Part (b)
xlabel('x (m)');% Part (c)-(f)
ylabel('Electron density (/m^3)');% Part (c)-(f)
%axis([0 100 0 4]);% Part (b)
axis([0 2e-8 0 1e9]);% Part (c)-(f)
