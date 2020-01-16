% driver_const.m script file
%
% Purpose:
%
% File driver_const applies the SID-PSM algorithmic framework to 
% solve the constrained optimization problem:
%
%                      min (x_2 - x_1^2)^2 
%
%                      s.t.  -2 <= x_1 <= 0
%                                  x_2 <= 1
%
% The initial point considered is [-1.2 1]'. The optimizer uses the
% default options specified in the file parameters.m. A complete 
% output report is produced, both at the screen and in the text file
% sid_psm_report.txt (stored at the directory sid_psm_1.1).
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
sid_psm([-1.2 1]',1,2);
%
% End of driver_const.