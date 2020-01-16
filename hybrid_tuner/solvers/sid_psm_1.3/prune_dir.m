function [D] = prune_dir(D,di_cosines,pruning);
%
% Purpose:
%
%    Function prune_dir prunes the columns of D according to the 
%    option specified by pruning.
%
% Input:  
%
%         D (matrix whose columns will be pruned).
%
%         di_cosines (vector storing the cosines of the angles
%                    between the columns of D and the descent
%                    indicator).
%
%         pruning (0-2 variable: 0 if no pruning is required;
%                 1 if only the most promising direction will
%                 be selected; 2 if all the descent directions,
%                 according to the descent indicator, will be
%                 considered).
%
% Output: 
%
%         D (pruned matrix).
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
if pruning == 0
    return;
end
if pruning == 1
    D = D(:,1);
else
    D = D(:,find(di_cosines > 0));
end
%
% End of prune_dir.