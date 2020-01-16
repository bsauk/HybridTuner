%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snoblocf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [y,f,g,sigma] = snoblocf(j,x,f,near,dx,u,v)
% computes a local fit around the point x0 = x(j,:) and minimizes it
% on a trust region
%
% Input:
% j         index of the point around which the fit is to be computed 
% x         the rows contain the points where the function has been
%           evaluated
% f         the corresponding function values and their uncertainties,
%           i.e., f(j,1) = f(x(j,:)), f(j,2) = df(x(j,:))
% near      near(j,:) is a vector containing the indices of the nearest
%           neighbors of the point x(j,:)
% dx        resolution vector, i.e. the ith coordinate of a point to be
%           generated is an integer-valued multiple of dx(i) 
% u,v       bounds of the box where the points should be generated
%
% Output:
% y         estimated minimizer in the trust region
%           
% f1        its estimated function value
% g         estimated gradient for the fit
% sigma     sigma = norm(A*g-b)/sqrt(K-n), where A and b are the 
%           coefficients resp. right hand side of the fit, n is the
%           dimension and K the number of nearest neighbors considered
%           (estimated standard deviation of the model errors)
%
function [y,f1,g,sigma] = snoblocf(j,x,f,near,dx,u,v)
n = length(u);   % dimension of the problem
x0 = x(j,:);
f0 = f(j,1);
df0 = f(j,2);
D = df0./dx.^2;
x1 = x(near(j,:),:);
K = size(x1,1);
S = x1-ones(K,1)*x0;
d = 0.5*max(abs(S),[],1);
d = max(d,dx);
for i=1:K
  sc(i) = S(i,:)*diag(D)*S(i,:)'+f(near(j,i),2);
end
A = S./(sc'*ones(1,n));
b = (f(near(j,:),1)-f0)./sc';
[U,Sigma,V] = svd(A,0);
Sigma = diag(Sigma);
Sigma = max(Sigma,1.e-4*Sigma(1));
g = V*diag(1./Sigma)*U'*b;
sigma = sqrt(sum((A*g-b).^2)/(K-n));
pl = max(-d,u-x0);
pu = min(d,v-x0);
for i=1:n
  p(i) = snobqmin(sigma*D(i),g(i),pl(i),pu(i));
end
y = snobround(x0+p,u,v,dx);
nc = 0;
while min(max(abs(x-ones(size(x,1),1)*y)-ones(size(x,1),1)*dx,[],2))<0 & nc < 5
  p = pl+(pu-pl).*rand(1,n);
  y = snobround(x0+p,u,v,dx);
  nc = nc + 1;
end
p = y-x0;
err=p*diag(D)*p'+df0;
f1 = f0+p*g+sigma*err;
