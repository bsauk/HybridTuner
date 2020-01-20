%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobpoint.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = snobpoint(x,xl,xu,f0,g,sigma,u,v,dx)
% for a box [xl,xu] containing a point x, a point y in the intersection 
% of [xl,xu] and [u,v] is constructed such that it is both not close to 
% x and to the boundary of [xl,xu] and its function value is estimated
% from a local quadratic model around x
%
% Input:
% x	     point contained in [xl,xu]
% xl,xu	     box bounds
% u,v        the point is to be generated in [u,v]
% f0         f0(1) is the function value at x, f0(2) is its uncertainty
% g,G,sigma  the local quadratic model around x is given by
%            q(y)=f0(1)+g*(y-x)'+sigma*((y-x)*diag(D)*(y-x)'+f0(2))
%            for a row vector y, where D = f0(2)./dx.^2
% dx         resolution vector
%
% Output:
% y	     point in the intersection of [xl,xu] and [u,v]
% f          corresponding estimated function value
%
function [y,f] = snobpoint(x,xl,xu,f0,g,sigma,u,v,dx)
n = length(x);
for i=1:n
  if x(i) - xl(i) > xu(i) - x(i)
    y(i) = 0.5*(xl(i)+x(i));
  else
    y(i) = 0.5*(x(i)+xu(i));
  end
end
y = min(max(y,u),v);
y = snobround(y,max(xl,u),min(xu,v),dx);
D = f0(2)./dx.^2;
f = f0(1) + g*(y-x)' +sigma*((y-x)*diag(D)*(y-x)'+f0(2));
