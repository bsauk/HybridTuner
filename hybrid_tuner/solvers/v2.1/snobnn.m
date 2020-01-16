%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobnn.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [near,d] = snobnn(x0,x,m,dx)
% computes a safeguarded set of m nearest neighbors to a point x0
% for each coordinate, the nearest point differing from x0 in the ith
% coordinate by at least 1 is chosen
% 
% Input:
% x0	point for which the neighbors are to be found
% x	the rows are the points from which the neighbors are to be
%	chosen (possibly x0 is among these points)
% m	number of neighbors to be found
% dx    resolution vector
%
% Output:
% near	vector pointing to the nearest neighbors (i.e. to the rows
%	of x)
% d	maximal distance between x0 and a point of the set of m
% 	nearest neighbors
% 
function [near,d] = snobnn(x0,x,m,dx)
n = length(x0);  % dimension of the problem
d = sqrt(sum((x-ones(size(x,1),1)*x0).^2,2));
[d1,ind] = sort(d);
if ~d1(1), ind(1) = []; end	% eliminate x0 if it is in the set
near = [];
for i=1:n
  j = min(find(abs(x(ind,i)-x0(i))-dx(i)>=0));
  near = [near ind(j)];
  ind(j) = [];
end
j = 1;
while length(near) < m
  if isempty(find(near==ind(j))) & max(abs(x0-x(ind(j),:))-dx)>= 0
    near = [near ind(j)]; 
  end
  j = j + 1;
end
[d,ind]=sort(d(near));
near = near(ind);
d = max(d);
