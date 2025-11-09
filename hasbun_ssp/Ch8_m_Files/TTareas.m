%copyright by J. E Hasbun and T. Datta
%TTareas.m
function [Ntt,areas,va]=TTareas(q2,q3,n,be,ga)
%Trangles based on the ray method of An-Ban Chen, Phys Rev. B V16, 3291 (1977)
%Divides a large triangular area into smaller triangles and finds their
%twodimensional vector directions and the smaller triangles ares
%finds the areas of the thin tetrahedron trangles according to how it's
%described in the Chen paper
%Here is an example of how the vertices 1,2,3,4,5,6 of a triangle is
%constructed. Example, we divide a large triangle in 2^2 smaller triangles
%1=(beta1*q2,gamma1*q3): (1,1)         /1\
%2=(beta2*q2,gamma1*q3): (2,1)        /   \
%3=(beta3*q2,gamma1*q3): (3,1)       /2---4\
%4=(beta2*q2,gamma2*q3): (2,2)      / \ * / \
%5=(beta2*q2,gamma2*q3): (2,3)     /   \ /   \
%6=(beta3*q2,gamma3*q3): (3,3)    3-----5-----6
%On the first loop pass, v1=1, v2=2, v3=4. On the next pass, v1=2
%v2=3, v3=5, v4=4 (extra triangle). On the last pass, we only have
%one triangle v1=4, v2=5, and v3=6.
%The triangle with the * is referred here to as an extra triangle.
%Once the vertices are found, the areas of each triangle can be found
%The numbers in parenthesis are the indices of the beta and gamma
%coefficients
%Inputs
%3D vectors q2, q3,
%arrays beta and gamma made according to be=0:dv:1 and ga=be
%number of subtriangles desired n.
%Output
%Nt: number of triangles
%areas: their areas
%va: average vector positions of these triangles based on the
%    gamma and beta as in Chen's paper
Ntt=0;
for io=1:n           %gamma loop
  for jo=io:n        %beta loop
    Ntt=Ntt+1;
    v1=be(jo)*q2+ga(io)*q3;     %small triangle vertices 1,2,3
    v2=be(jo+1)*q2+ga(io)*q3;
    v3=be(jo+1)*q2+ga(io+1)*q3;
    va(Ntt,:)=(v1+v2+v3)/3;     %average vertix vector (triangle centroid)
    vv1=v2-v1; vv2=v3-v2;       %small triangle vectors
    areas(Ntt)=norm(cross(vv1,vv2)/2); %area of little triangle
    if(jo > io)                 %extra triangle, with indices of v1, v3, v4
      %[jo,io]
      %[jo+1,io+1]
      %[jo,io+1]
      Ntt=Ntt+1;
      v4=be(jo)*q2+ga(io+1)*q3;
      va(Ntt,:)=(v1+v3+v4)/3;   %average vertix vector (triangle centroid)
      vv1=v3-v1; vv2=v4-v3;     %its associated vectors
      areas(Ntt)=norm(cross(vv1,vv2))/2; %area of little triangle
    end
  end
end
