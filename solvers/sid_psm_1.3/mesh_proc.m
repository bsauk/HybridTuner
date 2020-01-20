function [alfa] = mesh_proc(alfa,mesh_option,success,phi,theta,match_dir,...
                           poised,x,xtemp,f,ftemp,g,H);
%
% Purpose:
%
%    Function mesh_proc updates the mesh size parameter, according to the
%    strategy specified by mesh_option.
%
% Input: 
%
%         alfa (current mesh size parameter).
%       
%         mesh_option (0-3 variable: 0 if mesh size is increased in every 
%                     successful iteration; 1 if the mesh update for 
%                     successful iterations is based on a sufficient
%                     decrease condition, allowing mesh contractions;
%                     2 if the mesh update for successful iterations is 
%                     based on a sufficient decrease condition, but 
%                     contractions are not allowed; 3 if the mesh size is
%                     maintained in successful iterations, except when two
%                     consecutive successful iterations are found using 
%                     the same direction, where the mesh size is increased).
%
%        success (0-1 variable: 1 if the current iterate is successful; 0 otherwise).
%
%        phi (coefficient for mesh expansion).
%
%        theta (coefficient for mesh contraction).
%
%        match_dir (0-1 variable used only if mesh_option is set to 3: 1 if
%                  the same direction lead to a success in two consecutive
%                  iterations; 0 otherwise).
%
%        poised (0-2 variable used only if mesh_proc is based on a sufficient
%               decrease condition: 1 if a linear \Lambda-poised set was identified;
%               2 if quadratic \Lambda-poisedness was achieved; 0 otherwise).
%
%        x (current iterate, used only if mesh_proc is based on a sufficient
%          decrease condition).
%
%        xtemp (new iterate, used only if mesh_proc is based on a sufficient
%              decrease condition).
%
%        f (function value at the current iterate, used only if mesh_proc is
%          based on a sufficient decrease condition).
%
%        ftemp (function value at the new iterate, used only if mesh_proc is
%              based on a sufficient decrease condition).
%
%        g (simplex gradient at the current iterate, used only if mesh_proc is
%          based on a sufficient decrease condition).
%
%        H (simplex Hessian at the current iterate, used only if mesh_proc is
%          based on a sufficient decrease condition).
%
%
% Output: 
%
%        alfa (new value for the mesh size parameter).
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set mesh update tolerances and constants (not to be user-modified).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
if (mesh_option == 1 | mesh_option == 2)
   gamma1  = 0.25;      % Bad predicted reduction when the mesh size update
                        % is based on a sufficient decrease condition.
   gamma2  = 0.75;      % Good predicted reduction when the mesh size update
                        % is based on a sufficient decrease condition.
   tol_rho = 10^-8;     % Tolerance for degenerated \rho value.
end
%
if success
   if (mesh_option == 1 | mesh_option == 2) & (poised == 1 | poised == 2)
%
% Check for any degenerated situation.
%
      deg_rho = 0;
      if poised == 2
         if abs(g'*(xtemp-x) + 1/2*(xtemp-x)'*H*(xtemp-x)) <= tol_rho
            deg_rho = 1;
         end
      else
         if abs(g'*(xtemp-x)) <= tol_rho 
            deg_rho = 1;
         end
      end
%         
% Compute achieved reduction vs. predicted reduction ratio.
%
      if ~deg_rho
         if poised == 2
            rho = (ftemp-f)/(g'*(xtemp-x)+1/2*(xtemp-x)'*H*(xtemp-x));
         else
            rho = (ftemp-f)/(g'*(xtemp-x));
         end
         if rho > gamma2
            alfa = phi*alfa;
         end 
         if mesh_option == 1 & rho <= gamma1
            alfa = theta*alfa;
         end
      else
         if mesh_option == 1
            alfa = theta*alfa;
         end
      end
   else
      if (mesh_option ~= 3  | (mesh_option == 3 & match_dir))
         alfa = phi*alfa;
      end
   end 
%          
else
   alfa = theta*alfa;
end
%
% End of mesh_proc.