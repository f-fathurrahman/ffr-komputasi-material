%copyright by J. E Hasbun and T. Datta
%interpFunc.m
function y=interpFunc(rho,x,fx)
%Interpolates the function fx using Lagrange interpolation
%at point rho given the function values fx at x in array form
%search the index of x where rho lies
no=3;     %interpolator type, 1=linear, 2=quadratic, 3=cubic, etc.
xlen=length(x);
xmin=min(x); xmax=max(x);
dx=(xmax-xmin)/xlen;
i0=1+floor((rho-xmin)/dx);
if(i0 >= xlen-no),i0=xlen-no; end
y = 0.;
for j = i0:i0+no
  %Evaluate the j-th coefficient
  Lj = 1.0;
  for k =i0:i0+no
    if(j ~= k)
      Lj = Lj * (rho-x(k) )/( x(j)-x(k) );
    end
  end
  %Add contribution of j-th term to the polynomial
  y = y + Lj * fx(j);
end
