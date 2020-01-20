function [poised,Y,Y_values] = lambda_poised(X,F_values,s_min,s_max,Delta,...
                                       lambda,tol_degset,shessian,economic);
%
% Purpose:
%
%    Function lambda_poised tries to identify a \Lambda-poised subset in X.
%
% Input:  
%
%        X (matrix storing the points columnwise).
%
%        F_values (vector storing the objective function values at the points).
%
%        s_min (minimum number of points required in the \Lambda-poised set).
%
%        s_max (maximum number of points allowed in the \Lambda-poised set).
%
%        Delta (trust-region radius).
%
%        lambda (\Lambda-poisedness constant).
%
%        tol_degset (tolerance value for considering a \Lambda-poised set
%                   as degenerated).
%
%        shessian (0-2 variable: 0 if linear \Lambda-poisedness should be 
%                 checked; 1 if quadratic \Lambda-poisedness is required;
%                 2: if linear \Lambda-poisedness should be tested when
%                 quadratic \Lambda-poisedness can not be achieved).
%
%        economic (0-2 variable: 0 if a sequential process is used for 
%                 \Lambda-poised set computation; 1 if block 
%                 \Lambda-poisedness is tested before starting the
%                 sequential process; 2 if block checking is allowed 
%                 and additionally the \Lambda-poisedness condition is 
%                 tested using the QR decomposition).
%
% Output: 
%
%        poised (0-2 variable: 1 if a linear \Lambda-poised set is identified;
%               2 if quadratic \Lambda-poisedness is achieved; 0 otherwise).
%
%        Y (matrix storing columnwise the points in the \Lambda-poised set).
%
%        Y_values (vector storing the objective function values at
%                 the points stored in Y).
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
[n,m] = size(X);
%
% Check if there are enough points in X.
%
if  ((shessian ~= 2) & m < s_min) | ((shessian == 2) & (m < floor((s_min+1)/2)))
   poised   =  0;
   Y        = zeros(n,m);
   Y_values = zeros(1,n);
else
%
% Translate the points to a ball of radius one centered at X(:,1).
%
   Sinitial = X(:,2:m) - diag(X(:,1))*ones(n,m-1);
   Sinitial = Sinitial/Delta;
%
% Check whether the points are or not in the trust region.
%
  index_initial = find(sqrt(sum(Sinitial.^2,1)) <= 1);
  test_poised   = 0;
%
  while test_poised ~= 2
%
% Check if there are enough points in the trust region.
%
   index = index_initial;
   S     = Sinitial(:,index);
   m_S   = size(S,2);
   if ((shessian ~= 2) & m_S < (s_min-1)) | ((shessian == 2) & ...
                                   (m_S < (floor((s_min+1)/2)-1)))
      poised      =  0;
      Y           = zeros(n,m);
      Y_values    = zeros(1,n);
      test_poised = 2;
   else
       poised = 0;
%
% Check for block \Lambda-poisedness.
%
         if economic ~= 0
            test_eco = 0;
            while ~poised & test_eco <= 1
               if test_eco == 0
                  m_block = min(m_S,s_max-1);
               else
                  m_block = s_min-1;
               end      
               if shessian ~= 0
                 Y = [S(:,1:m_block);1/2*S(:,1:m_block).^2]';
               else
                 Y = S(:,1:m_block)';
               end
%    
               if economic == 2
                  [Q,R]  = qr(Y,0);
                  if min(diag(R)) >= tol_degset
                     Y_lambda = norm(1./diag(R));
                  else
                     Y_lambda = Inf;
                  end
               else
                     Y_lambda = cond(Y)/norm(Y);
               end
%        
               if Y_lambda <= lambda
                  poised = 1;
                  index  = [1,index(1:m_block)+1];
               end
%     
               test_eco = test_eco + 1;
            end
         end
%
% Start the pointwise \Lambda-poised set computation.
%
         if ~poised
%
% Set counters and variables.
%
            cont_S = 1;
            cont_Y = 0;
            Y      = [ ];
%
            while (cont_S <= m_S) & (cont_Y < s_max-1)
%
% Check whether the point is well-poised or not.
%
               if shessian ~= 0 
                 aux = [S(:,cont_S);1/2*S(:,cont_S).^2];      
               else
                 aux = S(:,cont_S);  
               end
               Y_aux = [Y,aux]';
%
               if economic == 2
                  [Q,R]  = qr(Y_aux,0);
                  if min(diag(R)) >= tol_degset
                     Y_lambda = norm(1./diag(R));
                  else
                     Y_lambda = Inf;
                  end
               else
                     Y_lambda = cond(Y_aux)/norm(Y_aux);
               end
 %
               if Y_lambda <= lambda
                  poised      = 1;
                  cont_Y      = cont_Y+1;
                  Y(:,cont_Y) = aux;
               else
                  index(cont_S) = 0;
               end
               cont_S = cont_S + 1;            
            end 
%
% Check if there are enough well-poised points in Y.
%
            if cont_Y < (s_min-1)
              poised   = 0;
              Y        = zeros(n,m);
              Y_values = zeros(1,n);
            else    
%
% Adjust values at the centering point.
%
              index = [1,index(find(index(1:cont_S-1)))+1];
            end
%       
         end   
%
% When shessian = 2 and quadratic \Lambda-poisedness fails,
% set variables for testing linear \Lambda-poisedness.
%
         if (shessian == 2) & (~poised) & (test_poised == 0)
            shessian    = 0;
            s_min       = floor((s_min+1)/2);
            s_max       = floor((s_max+1)/2);
            test_poised = 1;
         else
%
% Indicate that a set verifying quadratic \Lambda-poisedness was identified.
%
            if (shessian ~= 0) & poised & (test_poised == 0)
                poised = 2;
            end
            test_poised = 2;
         end
      end
   end
%
% Retrieve the original coordinates.
%
   if poised ~= 0  
      Y        = X(:,index);
      Y_values = F_values(index);
   end
end
%
% End of lambda_poised.