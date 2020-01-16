%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobinput.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,np,t] = snobinput(x,f)
% checks whether there are any duplicates among the points given by the
% rows of x, throws away the duplicates and computes their average
% function values and an estimated uncertainty
%
% Input:
% x	the rows of x are a set of points
% f	f(j,1) is the function value of x(j,:) and f(j,2) is its
%	uncertainty
%
% Output:
% x	updated version of x (possibly some points have been deleted)
% f	updated version of f (f(j,1) is the average function value and
%	f(j,2) the estimated uncertainty pertaining to x(j,:))  
% np	np(j) is the number of times the row x(j,:) appeared in the
%	input version of x
% t	t(j) is np(j) times the variance of the function values measured
%	for point x(j,:)
%
function [x,f,np,t] = snobinput(x,f)
sx = size(x,1);
n = size(x,2);
i = 1;
if isempty(x)
  np = [];
  t = [];
end
while i <= sx
  j = i + 1;
  ind = [];
  while j <= sx
    if sum(x(i,:)==x(j,:)) == n
      ind = [ind j];
    end
    j = j + 1;
  end
  if ~isempty(ind)
    ind = [i ind];
    ind1 = find(isnan(f(ind,1)));
    if length(ind1) < length(ind)
      ind(ind1) = [];
      np(i) = length(ind);
      fbar = sum(f(ind,1))/np(i);
      t(i) = sum((f(ind,1)-fbar).^2);
      f(i,1) = fbar;
      f(i,2) = sqrt((sum(f(i,2).^2)+t(i))/np(i));
    else
      np(i) = 1;
      t(i) = 0;
    end
    x(ind(2:end),:) = [];
    f(ind(2:end),:) = [];
    sx = size(x,1);
  else
    np(i) = 1;
    t(i) = 0;
  end 
  i = i + 1;
end






