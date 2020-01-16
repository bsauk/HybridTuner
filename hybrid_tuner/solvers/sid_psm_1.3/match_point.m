function [match,x,f] = match_point(x,xnorm,X,F_values,X_norms,tol_match);
%
% Purpose:
%
%    Function match_point scans a list of previously evaluated points
%    to try to match a point provided by the optimizer. When matching 
%    is successful, the function returns the corresponding objective
%    function value.
%
% Input: 
%
%         x (point to be checked).
%
%         xnorm (1-norm of the point to be checked).
%
%         X (matrix of points to be used in point comparisons,
%           storing the points columnwise).
%
%         F_values (vector storing the objective function values at the
%                  points stored in X).
%
%         X_norms (vector storing the 1-norm of the points stored in X).
%
%         tol_match (tolerance value within which two points are
%                   considered as equal).
%
% Output: 
%
%         match (0-1 variable: 1 if x was previously evaluated; 0 otherwise).
%
%         x (vector storing the matched point coordinates).
%
%         f (objective function value at the matched point).
%
% Copyright (C) 2005 A. L. Custodio and L. N. Vicente.
%
% http://www.mat.uc.pt/sid-psm
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
%
%
% Prune from X the points that do not satisfy a 1-norm criterion.
%
f     = +inf;
index = find(abs(X_norms - xnorm) <= tol_match);
if ~isempty(index)
   X        = X(:,index);
   F_values = F_values(index);
end
%
% Finish search.
%
nX    = size(X,2);
index = find(max(abs(X-repmat(x,1,nX)),[],1) <= tol_match);
match = ~isempty(index);
%
% Retrieve the point coordinates and the objective function value. 
%
if match
   x = X(:,index(1));
   f = F_values(index(1));
end
%
% End of match_point.