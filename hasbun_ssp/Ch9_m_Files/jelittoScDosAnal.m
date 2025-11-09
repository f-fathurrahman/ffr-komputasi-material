%copyright by J. E Hasbun and T. Datta
%jelittoScDosAnal.m
function ge=jelittoScDosAnal(ee)
%Jelitto's DOS for the simple cubic (exact)
eaa=abs(ee);
if (eaa <= 3.) & (eaa >=  1.)
  a1=3.-eaa;
  a2=a1^2;
  a=sqrt(a1);
  b=80.3702-16.3846*a1;
  d=0.78978*(a2);
  f=-44.2639+3.66394*a1;
  h=-0.17248*(a2);
  ge=a*((b+d)+(f+h)*sqrt(eaa-1.));
else
  if(eaa < 1.)
    ge=70.7801+1.0053*(ee^2);
  else
    ge=0.;
  end
end
