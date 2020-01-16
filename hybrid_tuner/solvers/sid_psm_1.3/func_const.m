function [c_const] = func_const(x);
%
% Purpose:
%
%    Function func_const is an user provided function which
%    computes the values of the functions defining the 
%    constraints c_i(x) <= 0, i = 1,...,m, at a point provided
%    by the optimizer.
%
% Input:  
%
%         x (point given by the optimizer).
%
% Output: 
%
%         c_const (vector storing the values of the functions c_i
%                 defining the constraints at the given point).
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
c_const = [ ];
%
% -----------------------------------------------------------------------------
%  Block to be user modified.
% -----------------------------------------------------------------------------
c_const(1) = x(1);
c_const(2) = -x(1)-2;
c_const(3) = x(2)-1;
% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
c_const = c_const';
%
% End of func_const.