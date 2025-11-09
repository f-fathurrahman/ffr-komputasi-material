%copyright by J. E Hasbun and T. Datta
%TBfermiIsoSurfBzSc.m, Simple Cubic
%To find the full Fermi surface, we use MATLAB's isosurface
%command for a given Fermi energy value.
%The obtained Fermi Surface is superimposed with the Brillouin Zone.
%Given kx, ky, and kz, surfaces of contant energy E=Ef are sought.
function TBfermiIsoSurfBzSc
clear,clc
a=1.0;
q=pi/a;
qq=(7./5.)*q;
Ef=0.0;                             %half filled band
[kx,ky,kz]=meshgrid(-qq:0.1:qq,-qq:0.1:qq,-qq:0.1:qq);
E=-cos(kx*a)-cos(ky*a)-cos(kz*a);   %SC energy band
hf=figure;
%isosurface connects points that have the specified value
%much the way contour lines connect points of equal elevation.
isosurface(kx,ky,kz,E,Ef);    %the isosurface to get the Fermi surface
colormap(gray)
xlabel('kx','FontSize',14), ylabel('ky','FontSize',14)
str=cat(2,'Fermi Energy Surface and 1st BZ - Simple Cubic Ef=',...
    num2str(Ef,'%4.2f'));
zlabel('kz','FontSize',14), title(str,'FontSize',14)

%Let's draw the first BZ
%SC BZ - with volume: (2*pi/a)^3
set(hf,'Position',[90 78 560 420])
view(130,20)                  %viewpoint azimuth, elevation
grid on
axis([-qq qq -qq qq -qq qq])

%corners
c1=[1/2,1/2,1/2]; c2=[-1/2,1/2,1/2]; c3=[-1/2,-1/2,1/2]; c4=[1/2,-1/2,1/2];
c5=[1/2,1/2,-1/2]; c6=[-1/2,1/2,-1/2]; c7=[-1/2,-1/2,-1/2]; c8=[1/2,-1/2,-1/2];
c1=c1*2*q; c2=c2*2*q; c3=c3*2*q; c4=c4*2*q; c5=c5*2*q; c6=c6*2*q;
c7=c7*2*q; c8=c8*2*q;
%c1 top corner connectors
liner(c1,c2,'-','k',1.0)
liner(c1,c4,'-','k',1.0)
liner(c1,c5,'-','k',1.0)
%c2 top corner connectors
liner(c2,c6,'-','k',1.0)
liner(c2,c3,'-','k',1.0)
%c3 top corner connectors
liner(c3,c4,'-','k',1.0)
liner(c3,c7,'-','k',1.0)
%c5 bottom corner connectors
liner(c5,c6,'-','k',1.0)
liner(c5,c8,'-','k',1.0)
%c7 bottom corner connectors
liner(c7,c6,'-','k',1.0)
liner(c7,c8,'-','k',1.0)
%c8 bottom corner connectors
liner(c8,c4,'-','k',1.0);

function liner(v1,v2,lin_style_txt,lin_color_txt,lin_width_num)
%Draws a line given initial vector v1 and final vector v2
%lin_style_txt: line style text format
%lin_color_txt: line color text format
%lin_width_num: line width number format
line([v1(1),v2(1)],[v1(2),v2(2)],[v1(3),v2(3)],...
    'LineStyle',lin_style_txt,'color',...
    lin_color_txt,'linewidth',lin_width_num)
