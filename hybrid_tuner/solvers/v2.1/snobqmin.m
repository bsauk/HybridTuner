%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobqmin.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = snobqmin(a,b,xl,xu)
% minimization of the quadratic polynomial p(x) = a*x^2+b*x over [xl,xu]
%
% Input:
% a, b    coefficients of the polynomial
% xl,xu   bounds (xl < xu)
% 
% Output:
% x       minimizer in [xl,xu]
%
function x = snobqmin(a,b,xl,xu)
if a > 0
  x = -0.5*b/a;
  x = min(max(xl,x),xu);
else
  fl = a*xl^2+b*xl;
  fu = a*xu^2+b*xu;
  if fu <= fl
    x = xu;
  else
    x = xl;
  end
end
