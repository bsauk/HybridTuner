function [deriv] = simplex_deriv(Y,F_values,sd_order,min_norm,deriv_old);
%
% Purpose:
%
%    Given a set of points, simplex_deriv computes the simplex
%    gradient and, when specified, it also computes a diagonal 
%    approximation to the simplex Hessian matrix. The
%    determined case, as well as the underdetermined and
%    the overdetermined cases are considered.
%
% Input:  
%
%         Y (matrix storing the points columnwise, being the first
%           column the simplex center).
%
%         F_values (vector storing the objective function values 
%                  at the points).
%
%         sd_order (1-2 variable: represents the order of the simplex
%                  derivatives to be computed; 1 if only a simplex
%                  gradient is required; 2 if a diagonal approximation
%                  to the simplex Hessian matrix should be computed,
%                  together with a simplex gradient).
%
%         min_norm (0-1 variable: 1 if a minimum norm solution is computed 
%                  when solving underdetermined systems corresponding to 
%                  simplex derivatives calculations; 0 if in this situation
%                  the closest solution to the previously calculated simplex
%                  derivatives is computed).
%
%         deriv_old (last simplex derivatives computed; only used when
%                   min_norm is set to 0).
%
% Output: 
%
%         deriv (vector storing the components of the simplex derivatives).
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
% Translate the points to a ball centered at Y(:,1).
%
[n,nsimp] = size(Y);
S         = Y(:,2:nsimp) - diag(Y(:,1))*ones(n,nsimp-1);
if sd_order == 2
   S = [S; 1/2*S.^2];
end
delta = F_values(2:nsimp) - F_values(1);
%
% Check if in presence of determined simplex derivatives.
%
determined  = 0;
overdeterm  = 0;
underdeterm = 0;
if sd_order == 2
   if nsimp == (2*n+1)
      determined = 1;
   else
       if nsimp > (2*n+1)
           overdeterm = 1;
       else
           underdeterm = 1;
       end
   end
else
   if nsimp == (n+1)
      determined = 1;
   else
       if nsimp > (n+1)
           overdeterm = 1;
       else
           underdeterm = 1;
       end
   end
end
%
% Compute the simplex derivatives.
%
if determined
   deriv = S' \ delta;
else
    if underdeterm & ~min_norm
       delta = delta - S'*deriv_old;  
    end
    deriv = pinv(S') * delta;
    if underdeterm & ~min_norm
        deriv = deriv_old + deriv;
    end
end
%
% End of simplex_deriv.