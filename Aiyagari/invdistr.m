% calculations invariant distribution from probability transition matrix Pi
% input is tranpose of Pi, where Pi(i,j) is probability to go to i from j
% uses eigs command to find the eigenvector of the matrix 
%   corresponding to the largest eigenvalue (which equals 1)
% Michael Reiter, IHS Vienna
function D = invdistr(Pi)
  assert(all(abs(sum(Pi)-1)<1e-10));
  opts.disp=0;
  % compute 6 largest eigenvalues, with eigenvectors:
  [x,eval] = eigs(Pi,[],6,1+1e-10,opts);
  % find the unit eigenvalue:
  ii = find(abs(diag(eval)-1)<1e-10);
  if(isempty(ii))
    error('eigs not properly converged; no unit eigenvalue');
  elseif(length(ii)>1)
    error('multiple unit eigenvalues (stationary states)');
  else
    D = x(:,ii);
    % normalize eigenvector to have sum=1:
    D = D/sum(D);
    % ensure non-negativity:
    assert(min(D)>-1e-12);
    D = max(D,0);
    D = D/sum(D);
  end
end
