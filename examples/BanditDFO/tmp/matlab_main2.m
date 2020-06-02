path(path, '/home/bsauk/Documents/research/5thyear/HybridTuner/solvers/v2.1');
file = 'test';
fcn = 'func_f';
fac = 0;
ncall = 21;
u = bl;
v = bu;
n = 4;
npoint = 1;
nreq = n + 6;
x = x_0;
dx = (v-u)'*1.e-5;
for i=1:n
   if(dx(i) > 1.e-5)
      dx(i) = 1.e-5
   end
end
p = 0.5;
prt = 0;
for j=1:npoint
   f(j,:) = [feval(fcn, x(j,:))+fac*randn max(sqrt(eps), 3*fac)];
end
ncall0 = npoint;
params = struct('bounds', {u,v}, 'nreq', nreq, 'p', p);
% repeated calls to Snobfit
while ncall0 < ncall %repeat till ncall function values are reached
   if ncall0 == npoint % intial call
      [request, xbest, fbest] = snobfit(file, x, f, params, dx);
      ncall0, xbest, fbest
   else
      [request, xbest, fbest] = snobfit(file, x, f, params);
   end
   if prt>0, request, end
   clear x
   clear f
   for j=1:size(request,1)
      x(j,:) = request(j,1:n);
      f(j,:) = [feval(fcn, x(j,:))+fac*randn max(sqrt(eps), 3*fac)];
   end
   ncall0 = ncall0 + size(f,1);
   [fbestn jbest] = min(f(:,1));
   if fbestn < fbest
      fbest = fbestn
;      xbest = x(jbest, :);
      ncall0, xbest, fbest % Display Current number of function values, best point and function values if fbest has changed
   end
end
ncall0, xbest, fbest %show number of function values, best point and function value
