function [x,feasible] = proj_ort(x);
%
% Purpose:
%
%    Function proj_ort computes the orthogonal projection of a given point
%    in a feasible region defined only by imposing bounds in the problem 
%    variables.
%
% Input: 
%
%         x (point to be projected).
%
% Output: 
%
%         x        (projection).
%
%         feasible (0-1 variable: 1 if x is feasible; 0 if x is not feasible).
%
% Functions called: func_const, grad_const (application, user provided).
%
% Copyright (C) 2008 A. L. Custodio and L. N. Vicente.
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
index   = logical(c_const > 0);
if sum(index) == 0
   feasible = 1;
else
   n        = size(x);
   grad_c   = grad_const(x);
   new_x    = - func_const(zeros(n))'./sum(grad_c,1);
   new_var  = grad_c(:,index);
   new_x    = new_x(index);
   m        = size(new_var,2);
   for i = 1:m
       x(logical(abs(new_var(:,i)))) = new_x(i);
   end
   feasible = 1;
end
%
% End of proj_ort.