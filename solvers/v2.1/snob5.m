%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snob5.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x1 = snob5(x,u,v,dx,nreq)
% generates nreq points of class 5 in [u,v]
%
% Input:
% x       the rows are the points already chosen by Snobfit (can be 
%         empty)
% u,v     bounds of the box in which the points are to be generated
% dx      resolution vector, i.e. the ith coordinate of a point to be
%         generated is an integer-valued multiple of dx(i) 
% nreq    number of points to be generated
%
% Output:
% x1      the rows are the requested points
%         x1 is of dimension nreq x n (n = dimension of the problem)
%
function x1 = snob5(x,u,v,dx,nreq)
n = length(u);   % dimension of the problem
nx = size(x,1);
nx1 = 100*nreq;
xnew = ones(nx1,1)*u+rand(nx1,n).*(ones(nx1,1)*(v-u));
if nx
  for j=1:nx1
    xnew(j,:) = snobround(xnew(j,:),u,v,dx);
    d(j) = min(sum((x-ones(nx,1)*xnew(j,:)).^2,2));
  end
  ind = find(~d);
  xnew(ind,:) = [];
  d(ind) = [];
  x1 = [];
  nx1 = size(xnew,1);
  if size(d,2) > 1, d = d'; end
else
  x1 = xnew(1,:);
  xnew(1,:) = [];
  nx1 = nx1-1;
  d = sum((xnew-ones(nx1,1)*x1).^2,2);
  nreq = nreq -1;
end
for j = 1:nreq
  [dmax,i] = max(d);
  y = xnew(i,:);
  x1 = [x1; y];
  xnew(i,:) = [];
  d(i) = [];
  nx1 = nx1 - 1;
  d1 = sum((xnew-ones(nx1,1)*y).^2,2);
  d = min(d,d1);
end
