function [r,A] = grad_act(x,c_const,epsilon);
%
% Purpose:
%
%    Function grad_act computes the approximated active constraints
%    at a point provided by the optimizer.
%
% Input:  
%
%         x (point given by the optimizer).
%
%         c_const (vector storing the values of the functions
%                 defining the constraints at the given point).
%
%         epsilon (tolerance for a constraint be considered as
%                 approximated active).
%
% Output: 
%
%         r (number of approximated active constraints).
%
%         A (n-times-r matrix storing columnwise the gradients of the 
%           functions defining the approximated active constraints).
%
% Functions called: grad_const (application, user provided).
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
Acomplete = grad_const(x);
%
A         = Acomplete(:,find(c_const >= -epsilon));
r         = size(A,2);
%
% End of grad_act.