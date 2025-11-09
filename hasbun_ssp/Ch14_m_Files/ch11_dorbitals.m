%copyright by J. E Hasbun and T. Datta
% ch11_dorbitals.m

% Boundary (isosurfaces) of atomic d orbitals
% A substantial portion (say 90%) of the total electron density of the orbital
% is enclosed within the boundary surface

% Wavefunctions defined below

% Normalized wave functions of hydrogen atom
% NOTE: Z =1, a0 = 1, 0<= l < n %-l<= m <=l

% d orbital
% ylm choice for d orbitals
 ylm_20 = @(theta,phi) (3*cos(theta).^2 - 1.0)*sqrt(5/pi)/4;
 ylm_21 = @(theta,phi) sqrt(2).*sin(theta).*cos(theta).*cos(phi)...
          *sqrt(15/pi)/2;
 ylm_22 = @(theta,phi) (sin(theta).^2).*(cos(2*phi))*sqrt(15/pi)/4;

% 3d
 psi_3dz2 = @(r,theta,phi) exp(-(1/3)*r).*(r.^2).*(3*cos(theta).^2 - 1.0)...
           /(81*sqrt(6*pi));
 psi_3dxz = @(r,theta,phi) exp(-(1/3)*r).*(r.^2).*sqrt(2).*sin(theta)...
           .*cos(theta).*cos(phi)/(81*sqrt(pi));
 psi_3dyz = @(r,theta,phi) exp(-(1/3)*r).*(r.^2).*sqrt(2).*sin(theta)...
            .*cos(theta).*sin(phi)/(81*sqrt(pi));
 psi_3dx2y2 = @(r,theta,phi) exp(-(1/3)*r).*(r.^2).*(sin(theta).^2)...
           .*(cos(2*phi))/(81*sqrt(2*pi));
 psi_3dxy = @(r,theta,phi) exp(-(1/3)*r).*(r.^2).*(sin(theta).^2)...
            .*(sin(2*phi))/(81*sqrt(2*pi));

%configuring the range
[x, y , z] = meshgrid(-30:0.5:30,-30:0.5:30,-30:0.5:30);

% Cartesian to spherical coordinates conversion
R=sqrt(x.^2+y.^2+z.^2);
THETA=acos(z./R);
PHI=atan2(y,x);

% Plotting orbtials

figure

subplot(1,5,1);
colors = ylm_22(THETA,PHI); psi = psi_3dxy(R,THETA,PHI); psisq = psi.^2;
set(gcf,'color',[1 1 1]);
daspect([1 1 1]); axis off; view(3);
camlight('left'); camzoom(0.75); lighting phong;
axis vis3d; colormap jet; rotate3d on; brighten(1);
isosurface(psisq,1E-5,colors);
title('(a) 3d_{xy}','FontName','Times','FontSize',12)

subplot(1,5,2)
colors = ylm_21(THETA,PHI); psi = psi_3dyz(R,THETA,PHI); psisq = psi.^2;
set(gcf,'color',[1 1 1]);
daspect([1 1 1]); axis off; view(3);
camlight('left'); camzoom(0.75); lighting phong;
axis vis3d; colormap jet; rotate3d on; brighten(1);
isosurface(psisq,1E-5,colors);
title('(b) 3d_{yz}','FontName','Times','FontSize',12)

subplot(1,5,3);
colors = ylm_21(THETA,PHI); psi = psi_3dxz(R,THETA,PHI); psisq = psi.^2;
set(gcf,'color',[1 1 1]);
daspect([1 1 1]); axis off; view(3);
camlight('left'); camzoom(0.75); lighting phong;
axis vis3d; colormap jet; rotate3d on; brighten(1);
isosurface(psisq,1E-5,colors);
title('(c) 3d_{xz}','FontName','Times','FontSize',12)

subplot(1,5,4)
colors = ylm_22(THETA,PHI); psi = psi_3dx2y2(R,THETA,PHI); psisq = psi.^2;
set(gcf,'color',[1 1 1]);
daspect([1 1 1]); axis off; view(3);
camlight('left'); camzoom(0.75); lighting phong;
axis vis3d; colormap jet; rotate3d on; brighten(1);
isosurface(psisq,1E-5,colors);
title('(d) 3d_{x^2 - y^2}','FontName','Times','FontSize',12)

subplot(1,5,5)
colors = ylm_20(THETA,PHI); psi = psi_3dz2(R,THETA,PHI); psisq = psi.^2;
set(gcf,'color',[1 1 1]);
daspect([1 1 1]); axis off; view(3);
camlight('left'); camzoom(0.75); lighting phong;
axis vis3d; colormap jet; rotate3d on; brighten(1);
isosurface(psisq,1E-5,colors);
title('(e) 3d_{3z^2 - r^{2}}','FontName','Times','FontSize',12)
