function [quad,H,g] = quad_Frob(X,F_values,always,regopt,H_calc);
%
% Purpose:
%
%    Function quad_Frob computes a quadratic interpolation model by 
%    minimizing the Frobenius norm of the Hessian matrix. When the number
%    of points in the interpolation set exceeds the total number of points
%    required for computing a complete quadratic model, two variants are 
%    implemented: to compute a regression model or to discard
%    some of the points and build a complete determined model.
%
% Input:  
%
%        X (matrix storing the points columnwise).
%
%        F_values (vector storing the objective function values at the
%        points).
%
%        always (0-1 variable: 1 if a quadratic model should always be used, 
%               once it has been computed; 0 otherwise).
%
%        regopt (0-1 variable: 1 if a regression model is considered when 
%               the number of points in the interpolation set exceeds the
%               number of points required for complete quadratic 
%               interpolation; 0 if some of the interpolation points are
%               discarded, in order to only compute determined quadratic 
%               interpolation models).
%
%        H_calc (number of matrix Hessians computed).
%
%
% Output: 
%
%        quad (0-2 variable: 1 if a quadratic model is computed; 2 if a
%             previous computed model should be used; 0 otherwise).
%
%        H (matrix storing the Hessian of the model).
%
%        g (vector storing the model gradient).
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
% Initialization step.
%
tol_svd = eps;   % Minimum value accepted for a singular value.
%
if (regopt == 0)
   perc = 0.8;   % Percentage of nearest points to the current iterate, 
                 % which should be used in model computation.
end
%
[n,m] = size(X);
H     = zeros(n);
g     = zeros(n,1);
%
% Check if there are enough points in X.
%
if (m <= n+1)
    quad = 0;
    if (always & (H_calc ~= 0))
        quad = 2;
    end
else
%
    quad     = 1;
    if (m <= (n+1)*(n+2)/2) | (regopt == 1)
       Y        = X;
       Y_values = F_values;
    else
       Y        = [X(:,1:floor(perc*(n+1)*(n+2)/2)),X(:,(m-(n+1)*(n+2)/2+floor(perc*(n+1)*(n+2)/2)+1):m)];
       Y_values = [F_values(:,1:floor(perc*(n+1)*(n+2)/2)),F_values(:,(m-(n+1)*(n+2)/2+floor(perc*(n+1)*(n+2)/2)+1):m)];
    end
    m_Y = size(Y,2);    
    Y   = Y - diag(X(:,1))*ones(n,m_Y);
%
% Compute a quadratic model by minimizing the Frobenius norm of the Hessian.
%    
    if (m_Y <= (n+1)*(n+2)/2) | (regopt == 0)
       b = [Y_values zeros(1,n+1)]';
       A = ((Y'*Y).^2)/2;
       W = [A ones(m_Y,1) Y';[ones(1,m_Y);Y] zeros(n+1,n+1)];
%
% Compute the model coefficients.
%
       [U,S,V]        = svd(W);
       Sdiag          = diag(S);
       indeces        = find(Sdiag < tol_svd);      
       Sdiag(indeces) = tol_svd;
       Sinv           = diag(1./Sdiag);
       lambda         = V * Sinv *U'* b;
%
% Retrieve the model coefficients.
%                
       g = lambda(m_Y+2:m_Y+n+1);
       H = zeros(n,n);
       for j = 1:m_Y
          H = H + lambda(j)*(Y(:,j)*Y(:,j)');
       end
    else
%
% Compute a complete quadratic model by solving a minimum least squares problem.
%         
       phi_Q = [ ];
       for i = 1:m_Y
           y      = Y(:,i);
           aux_H  = y*y'-1/2*diag(y.^2);
           aux    = [ ];
           for j = 1:n
              aux = [aux;aux_H(j:n,j)];
           end
           phi_Q  = [phi_Q; aux'];
       end
       W = [ones(m_Y,1) Y' phi_Q];
       b = Y_values';
%
% Compute the coefficients of the model.
%
       [U,S,V]        = svd(W,0);
       Sdiag          = diag(S);
       indeces        = find(Sdiag < tol_svd);      
       Sdiag(indeces) = tol_svd;
       Sinv           = diag(1./Sdiag);    
       lambda         = V * Sinv' *U'* b;
%
% Retrieve the model coefficients.
%        
       g    = lambda(2:n+1);
       H    = zeros(n,n);
       cont = n+1;
       for j = 1:n
          H(j:n,j) = lambda(cont+1:cont+n-(j-1));
          cont     = cont + n - (j-1);
       end
       H = H + H' - diag(diag(H));
    end        
end
%
% End of quad_Frob.