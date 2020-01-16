%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobround.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = snobround(x,u,v,dx)
% a point x is projected into the interior of [u,v] and x(i) is
% rounded to the nearest integer multiple of dx(i)
%
% Input:
% x	vector of length n
% u,v	vectors of length n such that u < v
% dx    vector of length n
%
% Output:
% x	projected and rounded version of x
%
function x = snobround(x,u,v,dx)
x = min(max(x,u),v);
x = round(x./dx).*dx;
i1 = find(x<u);
if ~isempty(i1)
  x(i1) = x(i1) + dx(i1);
end
i2 = find(x>v);
if ~isempty(i2)
  x(i2) = x(i2) - dx(i2);
end
