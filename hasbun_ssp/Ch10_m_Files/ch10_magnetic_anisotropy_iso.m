%copyright by J. E Hasbun and T. Datta
%ch10_magnetic_anisotropy_iso.m

% ch10_magnetic_anisotropy.m - by J. Hasbun and T. Datta
% code to simulate the magnetic anisotropy energy surface
% for an isotropic cubic crystal - with and without magnetic field

clear;

% declaring the anisotropy constants
% isotropic case
K0 = 1; K1 = 0;

% definining the anisotropy energy surface projections along x,y, and z

ecubic = @(x,y,K0,K1,H) (K0+ K1*(((sin(x).*cos(y)).^2).*...
         ((sin(x).*sin(y)).^2 + ((cos(x)).^2).*((sin(x).*sin(y)).^2)+...
         ((cos(x)).^2).*((sin(x).*cos(y)).^2)))+ ...
         (H/sqrt(3))*(sin(x).*cos(y)+sin(x).*sin(y)+cos(x)));

fcubicx = @(x,y,K0,K1,H) (sin(x).*cos(y)).*ecubic(x,y,K0,K1,H);
fcubicy = @(x,y,K0,K1,H) (sin(x).*sin(y)).*ecubic(x,y,K0,K1,H);
fcubicz = @(x,y,K0,K1,H) (cos(x)).*ecubic(x,y,K0,K1,H);

rotate3d on;
colormap(copper);
subplot(1,2,1)
ycubic = ezsurf(@(x,y) fcubicx(x,y,K0,K1,0),@(x,y) ...
   fcubicy(x,y,K0,K1,0),@(x,y) fcubicz(x,y,K0,K1,0),[0 pi],[0 2*pi]);
title('H = 0 (Zero field)');
% view(45,45)

% plot with magnetic field for the simple cubic lattice

rotate3d on;
colormap(copper);
subplot(1,2,2)
ycubich = ezsurf(@(x,y) fcubicx(x,y,K0,K1,1.5),@(x,y) ...
   fcubicy(x,y,K0,K1,1.5),@(x,y) fcubicz(x,y,K0,K1,1.5),[0 pi],[0 2*pi]);
title('H = 1.5 (Strong field)')
% view(45,45)
