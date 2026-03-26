% harmonic osscilator 3D

function y = V_osc_3D(x,y,z)
r = (x^2+y^2+z^2)^(1/2);
y = r^2;
%y = r^(1/2);
return;