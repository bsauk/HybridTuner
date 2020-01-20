function [x_final,f_final,histout] = sid_psm(x_initial,const,output);
%
% Purpose:
%
% Function sid_psm applies a pattern search method to the problem:
%
%          min f(x)  s.t.  c_i(x) <= 0, i = 1,...,m,
%
% where x is a real vector of dimension n. The derivatives of the function
% f are unknown. Only function values are provided for f. Both function and 
% gradient values are provided for c_i, i = 1,...,m. 
%
% The user must provide: func_f (for f function values),
%
%                        and, if const = 1,
%
%                        func_const (for c_i, i = 1,...,m, function values),
%                        grad_const (for c_i, i = 1,...,m, gradient values).
%
% Input:  
%
%         x_initial (the initial point to start the optimizer).
%
%         const (0-2 variable: 0 if the problem is unconstrained; 
%               1 if the problem has general constraints; 2 if the
%               constraints are only bounds).
%
%         output (0-2 variable: 0 - no output; 1 - verbose; 
%                2 - very verbose).
%
%
% Output:
%
%         sid_psm_report.txt (output text file).
%
% Functions called: func_f (application, user provided)
%
%                   domain, gen, grad_act, lambda_poised, match_point,
%                   mesh_proc, order_proc, proj_ort, prune_dir, quad_Frob, 
%                   simplex_deriv (provided by the optimizer).
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
time = clock;
warning off all;
format compact;
format long;
fprintf('\n\n');
rand('state', 0);
n = size(x_initial,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set tolerances, parameters and constants (not to be user-modified).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
parameters;            % Accept the parameters defined by the user or
                       % load the default values.
%
tol_degconst = 10^-3;  % Tolerance for degenerated constraints.
tol_hess     = 10^-3;  % Tolerance for degenerated Hessian element.
tol_Delta    = 10^-5;  % Minimum trust-region radius used in model 
                       % computation.
%
% For \Lambda-poised set computation.
%
lambda     = 100;       % \Lambda-poised constant.
sigma      = 2;         % Coefficient used in trust-region definition.
tol_degset = sqrt(eps); % Tolerance value for considering a\Lambda-poised
                        % set as degenerated.
%                      
% For simplex derivatives computation.
% 
if shessian
    if store_all
       s_min = 2*n+1;
       s_max = 2*n+1;
       p_max = 8*(n+1);
    else
       s_min = n;   
       s_max = 2*n+1;
       p_max = 4*(n+1);
    end 
else   
    if store_all
       s_min = n+1;
       s_max = n+1;
       p_max = (n+1)*(n+2);
    else
       s_min = floor((n+1)/2);   
       s_max = n+1;
       p_max = 2*(n+1);
    end 
end
%
% For cache implementation.
%
if cache == 1
   n_cache = 50*(n+1);
else
   n_cache = p_max;
end
%
tol_match = 10^-2*tol_alfa; % Tolerance used in point comparisons.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set counters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
func_eval = 0; % Function evaluation counter.
iter      = 0; % Iteration counter.
iter_suc  = 0; % Successful iteration counter.
iter_uns  = 0; % Unsuccessful iteration counter.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute a set of positive generators at the initial point.
%
% For the unconstrained case.
%
if ~const
  [D,deg] = gen(pss,n,const);
  r       = -1;
else
%
% For the constrained case.   
%   
% Check whether the initial point is feasible or not.
%
   [feasible,c_const] = domain(x_initial);
   if ~feasible
      fprintf('Initial point provided is not feasible. \n\n');
      fprintf('Please call function sid_psm with a feasible point. \n\n\n');
      return
   end
%
   epsilon = epsilon_ini;
   [r,A]   = grad_act(x_initial,c_const,epsilon);
   if (const == 2) & (pss == 0 | pss == 3)
      fprintf('For bound constrained optimization, the positive generating set ')
      fprintf('considered should include the coordinate directions. \n\n');
      fprintf('Please call function sid_psm with a different positive spanning set. \n\n\n');
      return
   end
   [D,deg] = gen(pss,n,const,r,A,tol_degconst);
   if deg
      print_format = ['Degenerated situation encountered. r = %2d. \n\n\n'];
      fprintf(print_format, r);
      return
   end
end
%
% Determine the number of positive generators.
%
nD    = size(D,2);
max_D = max(sqrt(sum(D.^2)));
%
% Compute f at the initial point.
%
x         = x_initial;
f         = func_f(x);
func_eval = func_eval+1;
%
% Record the initial point for printing.
%
histout(iter+1,1) = 1;
histout(iter+1,2) = f;
if (output < 2)
   REPORT(1,1) = 0;
   REPORT(2,1) = f;
   REPORT(3,1) = alfa;
else
   REPORT(1,1)     = 0;
   REPORT(2,1)     = -1;
   REPORT(3,1)     = -1;
   REPORT(4:n+3,1) = x;
   REPORT(n+4,1)   = f;
   REPORT(n+5,1)   = alfa;
   REPORT(n+6,1)   = r;
   REPORT(n+7,1)   = -1;
   REPORT(n+8,1)   = -1;
end
%
% Initialize the cache and/or the list used in simplex derivatives computation.
%
if search_option | (order_option == 1) | (order_option == 5) | stop_grad |...
        (mesh_option == 1) | (mesh_option == 2) | (cache ~= 0)
   X(:,1)      = x;
   F_values(1) = f;
end
if cache ~= 0
   X_norms(1) = norm(x,1);
end
if cache == 1
   label(1) = 1;
end
%
% Print the iteration report header.
%
if (output > 0)    
   fprintf('Iteration Report: \n\n');
   if (output == 1)
      fprintf('| iter  |     f_value      |      alpha       |\n');
      print_format = ['| %5d | %+13.8e | %+13.8e |\n'];
      fprintf(print_format, iter, f, alfa); 
   else
      fprintf('| iter  | success | #fevals |     f_value      |');
      fprintf('      alpha       | active | search | poised |\n');
      print_format = ['| %5d |    %2s   |    %2s   | %+13.8e | %+13.8e |'];
      if r == -1
         print_format = strcat(print_format,'   %2s   |   %2s   |   %2s   |\n');
         fprintf(print_format, iter, '--', '--', f, alfa, '--', '--', '--');
      else
         print_format = strcat(print_format,'   %2d   |   %2s   |   %2s   |\n');
         fprintf(print_format, iter, '--', '--', f, alfa, r, '--', '--');
      end
   end 
end
%
% Initialize some auxiliary variables.
%
halt        = 0;
poised      = -1;
search_step = -1;
quad        = 0;
deriv_old   = zeros(n,1);
if order_option == 2
   first_flag_colD = 0;
end
if mesh_option == 3
   dir_past = zeros(n,1);
end
%
% Initialize some auxiliary variables used in model computation. 
%
if search_option & always
   H_old  = zeros(n); 
   g_old  = zeros(n,1); 
   x_old  = zeros(n,1);
   H_calc = 0;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START THE PATTERN SEARCH METHOD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
while (~halt)
   success          = 0;
   func_iter        = 0;
   func_iter_finite = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEARCH STEP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           
   match = 0;
   if search_option
     search_step = 0;
%       
% Compute quadratic model.
%
     if iter == 0
        max_D_old = max_D;
     end
     Delta = alfa * sigma * max_D_old;
     if Delta < tol_Delta
        Delta = tol_Delta;
     end
     if (Delta > tol_match)
         if cache == 1
            X_aux        = X(:,find(label == 1));
            X_aux        = X_aux(:,1:min(size(X_aux,2),p_max));
            F_values_aux = F_values(find(label == 1));
            F_values_aux = F_values_aux(:,1:size(X_aux,2));
         else
            X_aux        = X;
            F_values_aux = F_values;
         end
        [quad,H,g] = quad_Frob(X_aux,F_values_aux,always,regopt,H_calc);
%
% If a quadratic model exists, compute a trust region trial point.
%   
        if (quad ~= 0)
           if quad == 1
              H_old  = H; 
              g_old  = g; 
              x_old  = X(:,1); 
              H_calc = H_calc + 1;
           else
              H = H_old; 
              g = g_old; 
           end
           if norm(g) <= eps
               info       = 1;
               [U,e]      = eig(H);
               [e,indexe] = sort(diag(e));
               xtemp      = Delta* U(:,indexe(1))/norm(U(:,indexe(1)));
           else
              if trs_solver
                 [fmodel,xtemp,iter_model,info] = soltr(H,g,Delta);
              else
                 info    = 1;
                 [xtemp] = trust(g,H,Delta);
              end
           end
           if (info >= 0)
              xtemp = xtemp + x;
%
% Check whether the trial point is feasible or not.
%
              if const
      	         [feasible,c_const_temp] = domain(xtemp);
              end
%
% In case of bound constrained optimization, project the trial point into 
% the feasible region.
%              
              if (const == 2) & ~feasible 
                 [xtemp,feasible] = proj_ort(xtemp);
              end
%      
% Compute the function value at the point or find a match.
%
              if (~const) | feasible 
                 if cache ~= 0
                    xtempnorm           = norm(xtemp,1);
                    [match,xtemp,ftemp] = match_point(xtemp,xtempnorm,X,...
                                          F_values,X_norms,tol_match);
                 end
                 if ~match
                    ftemp     = func_f(xtemp);
                    func_eval = func_eval + 1;
                    func_iter = func_iter + 1;
                 end
%
% Test for a better function value.
%
                 if ftemp < f 
                    success     = 1;
                    search_step = 1;      
                 end
%
% Update the cache and/or the list used in simplex derivatives
% computation (store_all version).      
% 
                 if isfinite(ftemp)& ~match & ((store_all & (search_option | (order_option == 1) |...
                 (order_option == 5) | stop_grad | (mesh_option == 1) |...
                 (mesh_option == 2) | (cache == 2))) | (cache == 1)) 
                    colX = size(X,2);
                    if cache == 1
                       if store_all
                          label = [1,label(:,1:min(n_cache-1,colX))];
                       else 
                          label = [0,label(:,1:min(n_cache-1,colX))];
                       end
                    end
                    func_iter_finite = func_iter_finite + 1;
                    X        = [xtemp,X(:,1:min(n_cache-1,colX))];
                    F_values = [ftemp,F_values(:,1:min(n_cache-1,colX))];
                    if cache ~= 0
                       X_norms = [xtempnorm,X_norms(:,1:min(n_cache-1,colX))];
                    end
                 end
              end 
           end
        end
     end
   end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POLL STEP.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   count_col = 0;
   if ((order_option ~= 4) & (order_option ~= 5) & (order_option ~= 8) &...
      (order_option ~= 9)) | const | (iter == 0)
      col_index = 1;
   else
      if (col_index == nD+1)
         col_index = 1;    
      end
   end
   if const & (order_option == 2)
      max_nD = nD+1;
   else
      max_nD = nD;
   end
   if const & (order_option == 2) & ~first_flag_colD
       count_col = count_col +1;
   end
%
% Permutes the set of positive generators.
%  
   if order_option == 3
       D = D(:,randperm(nD));
   end
   while ~success & (count_col < max_nD)
      match  = 0;
%
% Test the last successful direction (dynamic polling option).
%
      if (order_option == 2) & first_flag_colD
         if count_col == 0 
            colD      = flag_colD;	
            col_index = 0; 
         else
            colD = D(1:n,col_index);
            if colD == flag_colD
               if const 
                  count_col = count_col+1;
                  if col_index == nD
                      break
                  end
               end
               col_index = col_index + 1;
               colD      = D(1:n,col_index);   
            end
         end
      else
         colD = D(1:n,col_index);
      end
      xtemp = x + alfa*colD;
%
% Check whether the trial point is feasible or not.
%
      if const
      	[feasible,c_const_temp] = domain(xtemp);
      end
%      
% Compute the function value at the point or find a match.
%
      if (~const) | feasible
         if cache ~= 0
            xtempnorm           = norm(xtemp,1);
            [match,xtemp,ftemp] = match_point(xtemp,xtempnorm,X,...
                                             F_values,X_norms,tol_match);
         end
         if ~match
            ftemp     = func_f(xtemp);
            func_eval = func_eval + 1;
            func_iter = func_iter + 1;
         end
%
% Test for a better function value.
%
         if ftemp < f 
            success = 1;
            if (order_option == 2) & (col_index ~= 0)
               first_flag_colD = 1;
               flag_colD       = colD;
            end
         end   
%
% Update the cache and/or the list of points used in simplex 
% derivatives computation (store_all version).
%
         if isfinite(ftemp) & ~match & ((store_all & (search_option | (order_option == 1) |...
         (order_option == 5) | stop_grad | (mesh_option == 1) |...
         (mesh_option == 2) | (cache == 2))) | (cache == 1)) 
            colX = size(X,2);
            if cache == 1
               if store_all
                  label = [1,label(:,1:min(n_cache-1,colX))];
               else 
                  label = [0,label(:,1:min(n_cache-1,colX))];
               end
            end
            func_iter_finite = func_iter_finite + 1;
            X          = [xtemp,X(:,1:min(n_cache-1,colX))];
            F_values   = [ftemp,F_values(:,1:min(n_cache-1,colX))];
            if cache ~= 0
               X_norms = [xtempnorm,X_norms(:,1:min(n_cache-1,colX))];
            end            
         end
      end
      col_index = col_index + 1;
      if (col_index == nD+1)
         col_index = 1;    
      end
      count_col = count_col + 1;
   end
   iter = iter + 1;
%
% Update the mesh size parameter.
%
   match_dir = 0;
   if mesh_option == 3
       dir_new = (xtemp-x) / alfa;
       if success & (dir_past == dir_new)
           match_dir = 1;
       end
       dir_past = dir_new;
   end
   if (mesh_option == 0) | (mesh_option == 3)
       alfa = mesh_proc(alfa,mesh_option,success,phi,theta,match_dir);
   else
       if (poised == 1) | (poised == 2)
           if poised == 2
              alfa = mesh_proc(alfa,mesh_option,success,phi,theta,match_dir,...
                              poised,x,xtemp,f,ftemp,g,H);
           else
              alfa = mesh_proc(alfa,mesh_option,success,phi,theta,match_dir,...
                              poised,x,xtemp,f,ftemp,g); 
           end
       else
            alfa = mesh_proc(alfa,mesh_option,success,phi,theta,match_dir,poised); 
       end
   end
%
   if success
      iter_suc = iter_suc + 1;
%
% Update the cache and/or the list of points used in simplex 
% derivatives computation (store-successful version).
%     
      if (search_option | (order_option == 1) | (order_option == 5) |...
              stop_grad | (mesh_option == 1) | (mesh_option == 2) |...
              (cache ~= 0)) & (~store_all) & (~match)
         if cache == 1 
            label(1) = 1;
         else
            colX     = size(X,2);
            X        = [xtemp,X(:,1:min(n_cache-1,colX))];
            F_values = [ftemp,F_values(:,1:min(n_cache-1,colX))];
            if cache ~= 0
               X_norms = [xtempnorm,X_norms(:,1:min(n_cache-1,colX))];
            end            
         end
      end   
      f = ftemp;
      x = xtemp;
      if const
         c_const = c_const_temp;
      end
   else
      iter_uns = iter_uns + 1;
   end
%
% Record the iteration results for printing.
%
   histout(iter+1,1) = histout(iter,1) + func_iter;
   histout(iter+1,2) = f;
   if (output < 2)
      REPORT(1,iter+1)     = iter;
      REPORT(2,iter+1)     = f;
      REPORT(3,iter+1)     = alfa;
   else
      REPORT(1,iter+1)     = iter;
      REPORT(2,iter+1)     = success;
      REPORT(3,iter+1)     = func_iter;
      REPORT(4:n+3,iter+1) = x;
      REPORT(n+4,iter+1)   = f;
      REPORT(n+5,iter+1)   = alfa;
      REPORT(n+6,iter+1)   = r;
      REPORT(n+7,iter+1)   = search_step;
      REPORT(n+8,iter+1)   = -1;
   end
   max_D_old = max_D;
%   
% If in the constrained case, compute a set of positive 
% generators at the new point.
%
   if const
      epsilon  = min(epsilon_ini,10*alfa);
      [r,A]    = grad_act(x,c_const,epsilon);
      [D,deg]  = gen(pss,n,const,r,A,tol_degconst);
      nD       = size(D,2); 
      max_D    = max(sqrt(sum(D.^2)));
      if deg
         print_format = ['Degenerated situation encountered. r = %2d. \n\n\n'];
         fprintf(print_format, r);
         return
      end
   end
%
% If in the unconstrained case and pruning was selected, compute a set of  
% positive generators at the new point.
%
   if (~const) & (pruning ~= 0)
      [D,deg]  = gen(pss,n,const);
      nD       = size(D,2); 
      max_D    = max(sqrt(sum(D.^2)));
   end   
%      
   if search_option | (order_option == 1) | (order_option == 5) | stop_grad |...
           (mesh_option == 1) | (mesh_option == 2)
%
% Move the last successful iterate for the first column of
% the matrix used in simplex derivatives computation.
%
      if store_all & (~success)
          aux_x                          = X(:,func_iter_finite+1);
          X(:,2:func_iter_finite+1)      = X(:,1:func_iter_finite);
          X(:,1)                         = aux_x;  
          aux_f                          = F_values(func_iter_finite+1);
          F_values(2:func_iter_finite+1) = F_values(1:func_iter_finite);
          F_values(1)                    = aux_f;
          if cache ~= 0
             aux_n                         = X_norms(func_iter_finite+1);
             X_norms(2:func_iter_finite+1) = X_norms(1:func_iter_finite);
             X_norms(1)                    = aux_n;
          end
          if cache == 1
             aux_l                       = label(func_iter_finite+1);
             label(2:func_iter_finite+1) = label(1:func_iter_finite);
             label(1)                    = aux_l;
          end
      end  
   end
   if (order_option == 1) | (order_option == 5) | stop_grad |...
           (mesh_option == 1) | (mesh_option == 2)
%
% Try to build a \Lambda-poised set.
%
      if cache == 1
         X_aux        = X(:,find(label == 1));
         X_aux        = X_aux(:,1:min(size(X_aux,2),p_max));
         F_values_aux = F_values(find(label == 1));
         F_values_aux = F_values_aux(:,1:size(X_aux,2));
      else
         X_aux        = X;
         F_values_aux = F_values;
      end
      Delta               = sigma * alfa * max_D_old;
      [poised,Y,Y_values] = lambda_poised(X_aux,F_values_aux,s_min,...
                            s_max,Delta,lambda,tol_degset,shessian,economic);       
      if (output == 2)
         REPORT(n+8,iter+1) = poised;
      end
%
% Compute the simplex derivatives.
%         
      if poised ~= 0 
          deriv = simplex_deriv(Y,Y_values',poised,min_norm,deriv_old);
          if ~min_norm
             deriv_old = deriv;
          end
          g = deriv(1:n);
          if poised == 2
             aux_hess        = deriv(n+1:2*n);
             logic           = logical(abs(aux_hess)<tol_hess);
             aux_hess(logic) = aux_hess(logic)+sign(aux_hess(logic))*tol_hess;
             logic           = logical(aux_hess==0);
             aux_hess(logic) = aux_hess(logic) + tol_hess;
             H               = diag(aux_hess);
          end
%
% Compute the descent indicator to be used in the search step and/or for 
% reordering the polling vectors.
%
          if (order_option == 1) | (order_option == 5) 
             if poised == 2
                di = -diag(1./diag(H))*g;
             else
                di = -g;
             end
%         
% Reorder the vectors in the positive generator set using the
% descent indicator and prunes the positive generator set.
%            
             [D,di_cosines] = order_proc(D,di);
             col_index      = 1;
             if pruning ~= 0
                 D  = prune_dir(D,di_cosines,pruning);
                 nD = size(D,2);
             end
          end
%
% Build a stopping indicator based on simplex gradients.
%
          if stop_grad == 2
             nsimp        = size(Y,2);
             dir          = Y(:,2:nsimp) - diag(Y(:,1))*ones(n,nsimp-1);
             sg_dir_deriv = dir'*g
             if sg_dir_deriv >= -tol_grad
                 halt = 1;
             end
          else 
             if (stop_grad == 1) & (max(abs(Delta*g)) <= tol_grad)
                 halt = 1;
             end
          end       
      end    
   end
   if ((order_option == 6) | (order_option == 8)) & (quad ~= 0)
      if shessian
         di = -H_old*g_old;
      else
         di = -g_old;
      end
      [D]       = order_proc(D,di);
      col_index = 1;
   end
   
%         
% Reorder the vectors in the positive generator set according to the
% model values.
%
  if ((order_option == 7) | (order_option == 9)) & (quad ~= 0)
     xvec = diag(x)*ones(n,nD) + alfa*D;
     for j = 1:nD
        fvec(1,j) = g_old'*xvec(:,j) + xvec(:,j)'*H_old*xvec(:,j)/2;
     end
     [fval,ind] = sort(fvec);
     D          = D(:,ind);
     col_index  = 1;
  end
%
% Test for ending.
%
   if (stop_alfa & (alfa < tol_alfa))
      halt = 1;
   end
   if (stop_fevals & (func_eval >= fevals_max))
      halt = 1;
   end
   if (stop_iter & (iter >= iter_max))
      halt = 1;
   end
%
% Print iteration report.
% 
   if (output > 0)
      if (output == 1)
         print_format = ['| %5d | %+13.8e | %+13.8e |\n'];
         fprintf(print_format, iter, f, alfa); 
      else
         print_format = ['| %5d |    %2d   |    %2d   | %+13.8e | %+13.8e |'];         
         if r == -1
            if poised == -1
               if search_step == -1
                  print_format = strcat(print_format,'   %2s   |   %2s   |   %2s   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, '--', '--', '--');
               else
                  print_format = strcat(print_format,'   %2s   |   %2d   |   %2s   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, '--', search_step, '--'); 
               end
            else
               if search_step == -1
                  print_format = strcat(print_format,'   %2s   |   %2s   |   %2d   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, '--' ,'--', poised);
               else
                  print_format = strcat(print_format,'   %2s   |   %2d   |   %2d   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, '--' ,search_step, poised);
               end
            end
         else
           if poised == -1
               if search_step == -1
                  print_format = strcat(print_format,'   %2d   |   %2s   |   %2s   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, r, '--', '--');
               else
                  print_format = strcat(print_format,'   %2d   |   %2d   |   %2s   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, r, search_step, '--'); 
               end
            else
               if search_step == -1
                  print_format = strcat(print_format,'   %2d   |   %2s   |   %2d   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, r ,'--', poised);
               else
                  print_format = strcat(print_format,'   %2d   |   %2d   |   %2d   |\n');
                  fprintf(print_format, iter, success, func_iter, f, alfa, r ,search_step, poised);
               end
           end     
         end
      end
   end
%
end   % End of while.
%
time    = etime(clock,time);
x_final = x;
f_final = f;
%
% Print final report.
%
if (output~=0)
   fprintf('\n\nFinal Report: \n\n');
   print_format = 'Elapsed Time = %10.3e \n\n';
   fprintf(print_format,time);
   fprintf('| #iter | #isuc | #fevals |  final f_value   |   final alpha    |\n');
   print_format = ['| %5d | %5d |  %5d  | %+13.8e | %+13.8e |\n\n'];
   fprintf(print_format, iter, iter_suc, func_eval, f, alfa);
   fprintf('Minimum Point:\n');
   for i = 1:n
      fprintf('%13.8e \n',x(i));
   end
%
% Print report file.
%
   fresult = fopen('sid_psm_report.txt','w');
   fprintf(fresult,'Report: \n\n');
   if (output == 1)
      fprintf(fresult,'| iter  |     f_value      |      alpha       |\n');
      print_format  = ['| %5d | %+13.8e | %+13.8e |\n'];
      fprintf(fresult,print_format,REPORT(:,1:iter+1));
   else
      print_format   = '| iter  | success | #fevals |';
      string_counter = 1;
      for i = 1:n
         if fix(i/10^string_counter) == 1
            string_counter = string_counter + 1;  
         end
         print_format = strcat(print_format,'       x(',int2str(i),')');
         for j = 1:(8-string_counter)
            print_format = strcat(print_format,char(160));
         end
         print_format = strcat(print_format,'|');
      end
      print_format = strcat(print_format,'     f_value      |      alpha       |');
      print_format = strcat(print_format,' active | search | poised |\n');
      fprintf(fresult,print_format);
      print_format  = ['| %5d |    %2d   |    %2d   |'];
      print_format1 = ['| %5d |    %2s   |    %2s   |'];
      for i = 1:n+2
         print_format  = strcat(print_format,' %+13.8e |');
         print_format1 = strcat(print_format1,' %+13.8e |');
      end
      if r == -1
         print_format1 = strcat(print_format1,'   %2s   |   %2s   |   %2s   |\n');
         if poised == -1
             if search_step == -1
               print_format  = strcat(print_format,'   %2s   |   %2s   |   %2s   |\n');
             else
               print_format  = strcat(print_format,'   %2s   |   %2d   |   %2s   |\n');
             end
         else
             if search_step == -1
               print_format  = strcat(print_format,'   %2s   |   %2s   |   %2d   |\n');
             else
               print_format  = strcat(print_format,'   %2s   |   %2d   |   %2d   |\n');
             end
         end
      else
         print_format1 = strcat(print_format1,'   %2d   |   %2s   |   %2s   |\n'); 
         if poised == -1
             if search_step == -1
               print_format  = strcat(print_format,'   %2d   |   %2s   |   %2s   |\n');
             else
               print_format  = strcat(print_format,'   %2d   |   %2d   |   %2s   |\n');
             end
         else
             if search_step == -1
               print_format  = strcat(print_format,'   %2d   |   %2s   |   %2d   |\n');
             else
               print_format  = strcat(print_format,'   %2d   |   %2d   |   %2d   |\n');
             end
         end
      end
      if r == -1 
         fprintf(fresult,print_format1,REPORT(1,1),'--','--',...
                REPORT(4:n+5,1),'--','--','--'); 
         if poised == -1 
            if search_step == -1
               for l = 2:iter+1
                  fprintf(fresult,print_format,REPORT(1:n+5,l),'--','--','--');
               end
            else
               for l = 2:iter+1
                  fprintf(fresult,print_format,REPORT(1:n+5,l),'--',REPORT(n+7,l),'--');
               end
            end
         else
            if search_step == -1
               for l = 2:iter+1
                  fprintf(fresult,print_format,REPORT(1:n+5,l),'--','--',REPORT(n+8,l));
               end
            else
               for l = 2:iter+1
                  fprintf(fresult,print_format,REPORT(1:n+5,l),'--',REPORT(n+7,l),REPORT(n+8,l));
               end
            end
         end
      else
         fprintf(fresult,print_format1,REPORT(1,1),'--','--',...
                REPORT(4:n+6,1),'--','--'); 
         if poised == -1
            if search_step == -1
               for l = 2:iter+1
                  fprintf(fresult,print_format,REPORT(1:n+6,l),'--','--');
               end
            else
               for l = 2:iter+1
                  fprintf(fresult,print_format,REPORT(1:n+7,l),'--');
               end
            end
         else
            if search_step == -1
               for l = 2:iter+1
                  fprintf(fresult,print_format,REPORT(1:n+6,l),'--',REPORT(n+8,l));
               end
            else
               fprintf(fresult,print_format,REPORT(:,2:iter+1));
            end
         end
      end
   end
   fprintf(fresult,'\n\nFinal Report: \n\n');
   print_format = 'Elapsed Time = %10.3e \n\n';
   fprintf(fresult,print_format,time);
   fprintf(fresult,'| #iter | #isuc | #fevals |  final f_value   |');
   fprintf(fresult,'   final alpha    |\n');
   print_format = ['| %5d | %5d |  %5d  | %+13.8e | %+13.8e |\n\n'];
   fprintf(fresult,print_format, iter, iter_suc, func_eval, f, alfa);
   fprintf(fresult,'Minimum Point:\n');
   for i = 1:n
      fprintf(fresult,'%13.8e \n',x(i));
   end
   fclose(fresult);
end
%
% End of sid_psm.