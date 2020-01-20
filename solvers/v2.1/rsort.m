%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rsort.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,w,cdfx,dof]=rsort(x,w);
% sort x in increasing order
% remove multiple entries, and adapt weights w accordingly
% x and w must both be a row or a column
% default input weights are w=1
% 
% if nargout>2, the weighted empirical cdf is computed at x
% dof = length(x) at input is the original number of degrees of freedom
%
function [x,w,cdfx,dof]=rsort(x,w);

if nargin==1, w=ones(size(x)); end;

debug=0;
if debug>0, 
  save rsort.debug x w
elseif debug<0,
  load rsort.debug; 
  disp('debug mode');
end;


% turn into rows and sort
col=(size(x,2)==1);
if col, x=x';w=w'; end;
[x,ind]=sort(x);
w=w(ind);

% remove repetitions
n=length(x);
ind=find([x(2:n),inf]~=x);
nn=length(ind);
x=x(ind);
w(1)=sum(w(1:ind(1)));
for i=2:nn,
  w(i)=sum(w(ind(i-1)+1:ind(i)));
end;
w(nn+1:n)=[];

% restore original shape
if col, x=x';w=w'; end;

if nargout<3, return; end; 

% get cumulative sum of weights
cdfx=w;
for i=2:nn,
  cdfx(i)=cdfx(i-1)+w(i);
end;

% adjust for jumps and normalize
cdfx=(cdfx-0.5*w)/cdfx(nn);
dof=n;



return 

%%%%%%%% test documentation %%%%%%%%%
% test rsort.m (outside of rsort.m)

n=200;

u=[1 1 1 5 9 13 27 66 78 85 88 93]';

figure(1);clf
while 1,
  x=u(fix(10*rand(n,1))+1);
  [x,w,cdfx,dof]=rsort(x);
  [x,w]
  plot(x,cdfx)
  input('next>');
end;
