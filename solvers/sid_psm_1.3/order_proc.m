function [D,di_cosines] = order_proc(D,di);
%
% Purpose:
%
%    Function order_proc reorders the columns of D according to the 
%    angles that they make with the vector di (by decreasing order 
%    of the corresponding cosines).
%
% Input:  
%
%         D (matrix whose columns are to be reordered).
%
%         di (vector with as many components as rows in D).
%
% Output: 
%
%         D (matrix whose columns are the columns of D,
%           after reordering).
%
%         di_cosines (vector storing the cosines of the angles
%                    between the columns of D and the descent
%                    indicator di).
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
di_cosines         = ( di'*D ) ./ ( sqrt(sum(D.^2)) * norm(di) );
[di_cosines,index] = sort(-di_cosines);
D                  = D(:,index);
%
% End of order_proc.