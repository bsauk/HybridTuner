%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobnewb.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xl,xu,vol,u,v] = snobnewb(xnew,xl,xu,small,u,v,u1,v1)
% if xnew contains points outside the box [u,v], the box is made larger
% such that it just contains all new points; moreover, if necessary, 
% [u,v] is made larger to contain [u1,v1]
% the box bounds xl and xu of the boxes on the boundary and the volume 
% measures small of these boxes are updated if necessary
%
% Input:
% xnew	the rows of xnew contain the new points
% xl	xl(j,:) is the lower bound of box j
% xu	xu(j,:) is the upper bound of box j
% small	small(j) is the volume measure of box j
% u,v	old box bounds
% u1,v1 box in which the points are to be generated
%
% Output:
% xl	updated version of xl
% xu 	updated version of xu
% small	updated version of small
% u,v	updated box bounds (the old box is contained in the new one)
%
function [xl,xu,small,u,v] = snobnewb(xnew,xl,xu,small,u,v,u1,v1)
nx = size(xl,1);
n = size(xl,2);
uold = u;
vold = v;
u = min([xnew;u],[],1);
v = max([xnew;v],[],1);
u = min(u,u1);
v = max(v,v1);
i1 = find(u<uold);
i2 = find(v>vold);
ind = [];
for j=1:length(i1)
  j1 = find(xl(:,i1(j))==uold(i1(j)));
  ind = [ind; j1];
  xl(j1,i1(j)) = u(i1(j));
end
for j=1:length(i2)
  j2 = find(xu(:,i2(j))==vold(i2(j)));
  ind = [ind; j2];
  xu(j2,i2(j)) = v(i2(j));
end
if length(i1)+length(i2)  % at least one of the bounds was changed
  small = -sum(round(log2((xu-xl)./(ones(nx,1)*(v-u)))),2);
end


