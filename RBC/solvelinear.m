% solves func(x)=0, assuming that func is linear in x
% Inputs:
%   func: handle or name of function, R^n->R^n
%   n: dimension of domain and range of func
%   varargin: will be passed to func
function x = solvelinear(func,n,varargin)

    z = zeros(n,1);
  f0 = feval(func,z,varargin{:});
  % J = zeros(n,n);
  J = 0*repmat(f0,1,n);
  for i=1:n
    x = z;
    x(i) = 1;
    J(:,i) = feval(func,x,varargin{:})-f0;
  end;

  if(n>100 & nnz(J)<0.1*n*n)
    x = -sparse(J)\f0;
  else
    x = -J\f0;
    f1 = feval(func,x,varargin{:});
    x = x - J\f1;
  end
