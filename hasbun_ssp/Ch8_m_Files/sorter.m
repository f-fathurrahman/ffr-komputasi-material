%copyright by J. E Hasbun and T. Datta
%sorter.m
function [ww,zz]=sorter(ww,zz)
%sorts eigenvalues and eigenvectors
%On input ww=unsorted eigenvalues, zz=unsorted eigenvectors
%On output ww=sorted eigenvalues, zz=sorted eigenvectors
NB=length(ww);
for io=1:NB
  for jo=io:NB-1
    if  ww(jo+1) < ww(io)
      wo=ww(jo+1);
      ww(jo+1)=ww(io);
      ww(io)=wo;
      z0=zz(:,jo+1);
      zz(:,jo+1)=zz(:,io);
      zz(:,io)=z0;
    end
  end
end
