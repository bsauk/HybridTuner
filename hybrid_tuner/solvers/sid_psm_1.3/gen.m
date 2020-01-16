function [D,deg] = gen(pss,n,const,r,A,tol_degconst);
%
% Purpose:
%
%    Function gen:
%
%    (i) computes a positive spanning set for R^n (if the problem to be solved
%        is unconstrained, if the only constraints are variable bounds or 
%        if there are no approximated active constraints).
%
%    (ii) computes a set of positive generators for the tangent cone 
%         defined by the approximated active constraints gradients.
%
%    (iii) detects any degenerated situation (numerical linear dependence
%          of constraints gradients). 
%
% Input:  
%
%         pss (0-4 variable: 0 for the minimal positive basis 
%             [-ones(n,1) eye(n)]; 1 for the maximal positive basis
%             [eye(n) -eye(n)]; 2 for the positive spanning set
%             [ones(n,1) -ones(n,1) eye(n) -eye(n)]; 3 for a positive basis
%             with angles of uniform amplitude among vectors;
%             4 for an asymptotic dense set of directions in the unit sphere).
%
%         n (number of variables).
%
%         const (0-2 variable: 0 if the problem is unconstrained; 
%               1 if the problem has general constraints; 2 if the
%               constraints are only bounds).
%
%         r (number of approximated active constraints).
%
%         A (n-times-r matrix storing columnwise the gradients of the 
%           functions defining the approximated active constraints).
%
%         tol_degconst (tolerance for detecting any degenerated situation).
%
% Output: 
%
%         D (matrix storing columnwise a positive basis for R^n, in the
%           unconstrained case or when there are no approximated active
%           constraints, or a positive generator set for the tangent cone 
%           defined by the approximated active constraints gradients).
%
%         deg (0-1 variable: 1 if the approximated active constraints
%             gradients are numerically linearly dependent; 0 otherwise).
%
% Copyright (C) 2009 A. L. Custodio and L. N. Vicente.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unconstrained or constrained optimization, when is required a dense set
% of directions in the unit sphere.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (pss == 4)
   deg   = 0;
   v     = 2*rand(n,1)-1;
   [Q,R] = qr(v);
   if ( R(1) > 1 )
      D = Q * [ eye(n) -eye(n) ];
   else
      D = Q * [ -eye(n) eye(n) ];
   end
else    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unconstrained case, bound constrained optimization or general constrained
% optimization with no approximated active constraints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (const == 0) | (const == 2) | (r == 0)   
      deg = 0;
      if pss == 0
         D = [-ones(n,1) eye(n)];
      else
         if pss == 1
            D = [eye(n) -eye(n)];
         else
            if pss == 2
               D = [ones(n,1) -ones(n,1) eye(n) -eye(n)]; 
            else
%               
% Generation of a positive basis with angles of uniform amplitude between
% vectors.
%
               D = eye(n)*(1/n+1)+(zeros(n)-1/n);
               D = chol(D);  
               D = [sum(D')' -D];  
            end
         end
      end
   else
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained case with approximated active constraints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if r > n
         deg = 1;
         D   = [];
      else
%  
%  Case of r (<= n) approximated active constraints.
%
         [U,S,V] = svd(A);
         S1      = S(1:r,1:r);
         if min(diag(S1)) < tol_degconst
            deg = 1;
            D   = [];
         else
            deg = 0;
            U1  = U(1:n,1:r); 
            W   = U1 * diag(1./diag(S1)) * V';
            if r < n      
               U2 = U(1:n,r+1:n);  
               d  = sum(U2',1)'; 
               D  = [ -U2 d -W ];
%
% Other alternatives:
%
%           D = [ U2 -d -W ]; 
%
            else 
               D = -W;
            end
         end
      end  
   end
%   
end
%
% End of gen.