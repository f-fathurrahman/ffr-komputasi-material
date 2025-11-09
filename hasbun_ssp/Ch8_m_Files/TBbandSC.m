%copyright by J. E Hasbun and T. Datta
%TBbandSC.m
%Tight binding model band for the simple cubic system
clear, clc;
a=1.0;
q=2*pi/a;    %range of k
kz=0.0;      %Keep kz fixed
[kx,ky]=meshgrid(-q:0.1:q,-q:0.1:q);  %variables for the 2D plot
E=-cos(kx*a)-cos(ky*a)-cos(kz*a);     %SC energy band
figure;
surfl(kx,ky,E);                       %2D plot
colormap(gray)
shading interp
xlabel('kx'), ylabel('ky'), zlabel('E (kx,ky,kz)')
str=cat(2,'SC Energy Band E_k=-cos(kx*a)-cos(ky*a)-cos(kz*a) vs kx, ky at kz=',...
    num2str(kz,'%4.3f'));
title(str)
axis([-q q -q q min(min(E)) max(max(E))])
view(115,70)                          %viewpoint azimuth, elevation
