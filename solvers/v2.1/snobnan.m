%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobnan.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = snobnan(fnan,f,near,inew)
% replaces the function values NaN of a set of points by a value 
% determined by their nearest neighbors with finite function values, 
% with a safeguard for the case that all neighbors have function value
% NaN
%
% Input:
% fnan	  	vector containing the pointers to the points where the
%               function value could not be obtained
% f	  	f(:,1) set of available function values
%		f(:,2) their uncertainty/variation
% near(j,:)	vector pointing to the nearest neighbors of point j
% inew          vector pointing to the new boxes and boxes whose nearest
%               neighbors have changed
%
% Output:
% f		updated version of f
%
function f = snobnan(fnan,f,near,inew)
lnn = size(near,2);
notnan = 1:size(f,1);
notnan(fnan) = [];
[fmx,imax] = max(f(notnan,1));
fmn = min(f(notnan,1));
dfmax = f(imax,2);
for j=1:length(fnan)
  l = fnan(j);
  if ~isempty(find(inew==l))
% a substitute function value is only computed for new points and for 
% points whose function values have changed
    ind = near(l,:); 
% eliminate neighbors with function value NaN
    ind1 = []; 
    for i=1:length(ind)
      if ~isempty(find(fnan==ind(i))),ind1 = [ind1 i]; end
    end
    if ~isempty(ind1)
      ind(ind1) = [];
    end
    if isempty(ind)
      f(l,1) = fmx + 1.e-3*(fmx-fmn);
      f(l,2) = dfmax;
    else
      [fmax1,k] = max(f(ind,1));
      fmin1 = min(f(ind,1));
      f(l,1) = fmax1 + 1.e-3*(fmax1-fmin1);
      f(l,2) = f(ind(k),2);
    end
  end
end

