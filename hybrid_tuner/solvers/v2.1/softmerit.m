%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% softmerit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fm = softmerit(f,F,F1,F2,f0,Delta,sigma)
% merit function of the soft optimality theorem
%
% Input:
% f        objective function value
% F        vector containing the values of the constraint functions
%          (m-vector)
% F1       m-vector of lower bounds of the constraints
% F2       m-vector of upper bounds of the constraints
% f0       scalar parameter in the merit function
% Delta    scalar, positive parameter in the merit function
% sigma    positive m-vector, where sigma(i) is the permitted violation 
%          of constraint i
%
function fm = softmerit(f,F,F1,F2,f0,Delta,sigma)
if ~isfinite(f) | any(~isfinite(F)) 
% if the objective function or one of the constraint functions is 
% infinite or NaN, set the merit function value to 3
  fm = 3;
  return
end
m = length(F);
delta = 0;
for i=1:m
  if F(i) < F1(i)
    delta = delta + (F1(i)-F(i))^2/sigma(i)^2;
  elseif F(i) > F2(i)
    delta = delta + (F(i)-F2(i))^2/sigma(i)^2;
  end
end
fm = (f-f0)/(Delta+abs(f-f0)) + 2*delta/(1+delta);
