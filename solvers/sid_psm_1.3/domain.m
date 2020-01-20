function [feasible,c_const] = domain(x);
%
% Purpose:
%
%    Function domain checks the feasibility of a given point.
%
% Input: 
%
%         x (point to be checked).
%
% Output: 
%
%         feasible (0-1 variable: 1 if x is feasible; 0 if x is not feasible).
%
%         c_const  (vector storing the values of the functions defining the 
%                  constraints at the given point).
%
% Functions called: func_const (application, user provided).
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
c_const = func_const(x);
m       = size(c_const,1);
%
if sum(c_const <= 0) == m
   feasible = 1;
else
   feasible = 0;
end
%
% End of domain.