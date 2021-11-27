clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
epsil = 8.854e-12;
q = 1.602e-19;
a0 = 4*pi*epsil*hbar*hbar/(m*q*q);
E0 = q/(8*pi*epsil*a0);
F = linspace(0,1e9,11);
A = (a0*128*sqrt(2)/243)*F;
B = (-3*a0)*F;

for kF = 1:11
    M = [-E0 0 A(kF);
         0 -E0/4 B(kF);
         A(kF) B(kF) -E0/4];
    [V,D] = eig(M);
    D = diag(D);
    [DD,ind] = sort(D);
    E1(kF) = D(ind(1));
    E2(kF) = D(ind(2));
    E3(kF) = D(ind(3));
end

% perturbation theory results
E1s = -E0-((A.^2)/(3*E0/4));
E2s = -(E0/4)+B;
E2p = -(E0/4)-B;

hold on

h = plot(F, E1, 'b');% Fig.3.4.1
h = plot(F, E1s, 'bx');% Fig.3.4.1

set(h,'linewidth', [2.0])
set(gca,'Fontsize', [25])
grid on
xlabel('Field (V/m)');
ylabel('Energy (eV)');
