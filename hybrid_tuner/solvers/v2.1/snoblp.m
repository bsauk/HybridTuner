%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snoblp.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [local,nlocal] = snoblp(f,near,ind)
% computes a pointer to all `local' points (i.e. points whose neighbors
% have `significantly larger' function values)
%
% Input:
% f	f(j) is the function value of point j
% near	near(j,:) contains the indices of the size(near,2) neighbors of
%	point j
% ind	pointer to the boxes to be considered (optional, default
%       1:length(f))
%
% Output:
% local		vector containing the indices of all local points
% nlocal	vector containing the indices of all nonlocal points
%
function [local,nlocal] = snoblp(f,near,ind)
if nargin < 3, ind = 1:length(f); end
local = [];
nlocal = ind';
jj = [];
for j = 1:length(ind)
  fmi = min(f(near(ind(j),:)));
  fma = max(f(near(ind(j),:)));
  if f(ind(j)) < fmi - 0.2*(fma-fmi)
    local = [local ind(j)];
    jj = [jj j];
  end
end
if ~isempty(jj)
  nlocal(jj) = [];
end
nloc = length(local);
