%copyright by J. E Hasbun and T. Datta
%hybrid_one_bond.m
%hybrid wavefunctions plotted in units of sqrt(1/[32*pi*a0^3])
%versus distance in units of a0. We work with the 2s, 2p hybridization
%
clear;
N=128;                    %grid points
ul=1.0; us=2.*ul/N;       %range and step size
[x,y,z]=meshgrid(-ul:us:ul,-ul:us:ul,-ul:us:ul);
r=sqrt(x.^2+y.^2+z.^2);   %distance
fr=0.5*exp(-r/2.0);       %fp(r)/r/2 - needed throughout
% One of the hybrids plot follows
figure
ff1=(2.0-r+x+y+z).*fr;    %one hybrid
s=0.95;                   %isosurface value plotted
isosurface(x,y,z,ff1,s);
colormap(gray)
lighting gouraud          %lighting control type for smooth looks
az=110; el=14;
camlight (az,el)          %camera light from azimuth, elevation
view(az,el)
box on
axis ([-0.25 0.25 -0.25 0.25 -0.25 0.25])
xlabel('x')
ylabel('y')
zlabel('z')
title('sp^3 Hybrid Orbitals')
