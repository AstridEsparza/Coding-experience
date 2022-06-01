% Computes transition matrices and invariant cross-sectional distribution
%   in a model with one endogenous and one exogenous state
% Solution is captured in the structure Sol, which is input
function [TEndog,TExog,TT] = transmat(Sol)
  [nEndog,nExog] = size(Sol.V);
  nn = nEndog*nExog;
  % exogenous transition is the same independent of the endogenous state:
  TExog = kron(sparse(Sol.TransExog'),speye(nEndog));
  % for each exogenous state there is a different endogenous transition:
  TEndog = sparse(nn,nn);
  ee = (1:nEndog)';
  for j=1:nExog
    ii = (j-1)*nEndog+1 : j*nEndog; % index of endogenous grid points
    if(isfield(Sol,'pLowOpt'))  % endogenous transition as it arises from linear interpolation
      jj = Sol.IndxOpt(:,j);
      pLow = Sol.pLowOpt(:,j);
      TEndog(ii,ii) = sparse([jj;jj+1],[ee;ee],[pLow;1-pLow],nEndog,nEndog);
    else  % endogenous transition is deterministic
      TEndog(ii,ii) = sparse(Sol.IndxOpt(:,j),1:nEndog,ones(nEndog,1),nEndog,nEndog);
    end
  end
  assert(all(abs(sum(TEndog)-1)<1e-10))
  TT = TExog*TEndog;
end
