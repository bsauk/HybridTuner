% parameters.m script file
%
% Purpose:
%
% File parameters sets the algorithmic strategies, parameters, constants
% and tolerances values to be used by the function sid_psm, that can be 
% user modified.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithmic strategies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
always        = 1;  % 0-1 variable: 1 if a quadratic model should always be 
                    % considered at the search step of the algorithm, once  
                    % it has been computed; 0 otherwise.
cache         = 0;  % 0-2 variable: 0 if no point comparisons are performed   
                    % prior to function evaluation; 1 if all the points   
                    % where the objective function has been evaluated are used  
                    % in point comparisons; 2 if the only points used for 
                    % comparisons are the ones stored for simplex derivatives
                    % computation.
economic      = 0;  % 0-2 variable: 0 if a sequential process is used when 
                    % computing a \Lambda-poised set; 1 if block
                    % \Lambda-poisedness is checked before starting the 
                    % sequential process; 2 if block checking is allowed and  
                    % additionally the test for the \Lambda-poisedness
                    % condition is based on the QR decomposition. 
mesh_option   = 0;  % 0-3 variable: 0 if mesh size is increased in every 
                    % successful iteration; 1 if the mesh update for 
                    % successful iterations is based on a sufficient decrease
                    % condition, allowing mesh contractions; 2 if the mesh
                    % update for successful iterations is based on a
                    % sufficient decrease condition, but contractions are
                    % not allowed; 3 if the mesh size is maintained in 
                    % successful iterations, except when two consecutive
                    % successful iterations are found using the same
                    % direction, where the mesh size is increased.
min_norm      = 1;  % 0-1 variable: 1 if a minimum norm solution is computed
                    % when solving underdetermined systems corresponding to 
                    % simplex derivatives calculations; 0 if the closest 
                    % solution to the previously calculated simplex 
                    % derivatives is computed in this situation.
order_option  = 5;  % 0-9 variable: 0 if the vectors in the poll set are
                    % tested at each iteration in a consecutive order 
                    % of storage; 1 if the reordering of the poll
                    % directions is based on a simplex descent indicator;
                    % 2 if dynamic polling is considered, testing first
                    % the last successful poll direction; 3 if a random
                    % order is considered for the poll set; 4 if a
                    % consecutive, cycling order of storage is
                    % considered through all the iterations (only 
                    % available for unconstrained optimization); 
                    % 5 ordering strategy identical to 4, except when 
                    % $\Lambda$-poisedness is achieved. In this case the
                    % poll directions are ordered according to the 
                    % simplex descent indicator; 6 if the reordering of
                    % the poll directions is based on a model descent 
                    % indicator; 7 if the reordering of the poll 
                    % directions is based on the model values; 8 ordering
                    % strategy identical to 4, except when a model is 
                    % available for use. In this case the poll directions
                    % are ordered according to the model descent indicator;
                    % 9 ordering strategy identical to 4, except when a 
                    % model is available for use. In this case the
                    % poll directions are ordered according to the
                    % model values.
pruning       = 0;  % 0-2 variable: 0 if no pruning is required; 1 if only the 
                    % most promising direction is considered for polling; 2 if 
                    % all the potential descent directions, according to the
                    % descent indicator, are used for polling.
pss           = 2;  % 0-4 variable:  0 for the minimal positive basis 
                    % [-ones(n,1) eye(n)]; 1 for the maximal positive
                    % basis [eye(n) -eye(n)]; 2 for the positive spanning
                    % set [ones(n,1) -ones(n,1) eye(n) -eye(n)]; 3 for 
                    % a positive basis with angles of uniform amplitude
                    % among vectors; 4 for an asymptotic dense set of 
                    % directions in the unit sphere.
regopt        = 1;  % 0-1 variable: 1 if, at the search step, a regression 
                    % model is considered when the number of points in the
                    % interpolation set exceeds the number of points 
                    % required for complete quadratic interpolation; 0 if 
                    % some of the interpolation points are discarded, in 
                    % order to only compute determined quadratic 
                    % interpolation models.
search_option = 1;  % 0-1 variable: 0 for no search step; 1 for a search step 
                    % based on a minimum Frobenius norm model.
shessian      = 0;  % 0-2 variable: 0 if the simplex derivatives to be computed
                    % correspond to a simplex gradient; 1 if the simplex 
                    % derivatives to be computed correspond to a simplex gradient
                    % and a diagonal simplex Hessian; 2 if one first
                    % attempts shessian = 1 and, in case of failure, then 
                    % tries shessian = 0.
store_all     = 1;  % 0-1 variable: 1 if all function evaluations are stored
                    % for simplex derivatives and model computations; 0 if
                    % the only points stored with this purpose are the ones
                    % corresponding to successful iterations.
trs_solver    = 0;  % 0-1 variable: 1 if trs routine, based on dgqt solver 
                    % of MINPACK2 is used for solving the trust-region
                    % subproblems; 0 if the solution of the trust-region 
                    % subproblems is computed using trust.m Matlab function.                        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stopping criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
%
stop_alfa   = 1;      % 0-1 variable: 1 if the stopping criterion is based on
                      % the mesh size parameter; 0 otherwise.
tol_alfa    = 10^-5;  % Lowest value allowed for the mesh size parameter.
                      %                        
stop_fevals = 0;      % 0-1 variable: 1 if the stopping criterion is based on a 
                      % maximum number of function evaluations; 0 otherwise.
fevals_max  = 1500;   % Maximum number of function evaluations allowed.
                      %
stop_grad   = 0;      % 0-2 variable: 0 if the stopping criterion does not
                      % involve any simplex gradients; 1 if the stopping 
                      % criterion is related to the norm of a scaled simplex 
                      % gradient; 2 if the stopping criterion is based on the 
                      % signs of the approximated directional derivatives 
                      % (computed using simplex gradients).
tol_grad    = 10^-5;  % Convergence tolerance when the stopping criterion
                      % is based on simplex gradients.
                      % 
stop_iter   = 0;      % 0-1 variable: 1 if the stopping criterion is based 
                      % on a maximum number of iterations; 0 otherwise.
iter_max    = 1500;   % Maximum number of iterations allowed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh size parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
alfa  = max(1,norm(x_initial,inf));  % Initial mesh parameter value 
                                     % (set as in Moré and Wild).
phi   = 1;                           % Coefficient for mesh expansion.
theta = 0.5;                         % Coefficient for mesh contraction.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
epsilon_ini = 10^-1;  % Initial value for considering a constraint as
                      % approximated active.
%
% End of parameters.