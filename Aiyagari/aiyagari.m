% solve for the steady state of the model with heterogeneous households
% labor supply is exogenous
% solve HH infinite horizon saving problem on finite grid, using acceleration steps
% household productivity follows Markov process
% Michael Reiter, IHS, February 2017
% last modified November 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A SIMPLE HETEROGENEOUS AGENT MODEL (a modification of Aiyagari 1994)
% 
% Economy with very many agents
% Agents (infinitely lived) are identical ex-ante:
% - have the same preferences
% - facing the same productivity process 
%   they offer one unit of labor inelastically
%   labor income = wage * productivity
%   productivity follows a Markov process (approximation to an AR1)
% Agents differ ex post:
%   different realizations of the productivity process
% 
% Agents decide on consumption/saving
% Borrowing limit (typically: you can save, but not borrow)
% 
% MACROECOOMY:
% 
% stationary economy: constant r and wage (stationary: no aggregate shocks)
% To do:
%   find equilibrium r and wage
%   compute cross-sectional distribution of asset holdings
%   compute aggregate variables (such as capital, output, etc.)
% 
% 
% FIRMS:
%   representative firm (like in RBC model)
%   constant returns to scale
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r2,sol2] = aiyagari(sigma,nK,nZ,Kmax_factor)
  if nargin<2
    nK = 1000;
  end
  if nargin<3
    nZ = 11;
  end
  if nargin<4
    Kmax_factor = 100;
  end
  params.sigma = sigma;  % standard deviation of shock
  params.nK = nK;
  params.nZ = nZ;
  params.rho = 0.9;
  params.beta = 0.97;
  params.Kmax_factor = Kmax_factor;
  [r2,sol2] = do_aiya(params);
end

function [rEqu,Sol] = do_aiya(params)
  rMax = 1/params.beta - 1;

  f = @(x) findR(x,params);
  options = optimset('TolX',1e-10);
  rEqu = fzero(f,[rMax-0.01 rMax-1e-7],options);
  [resid,Sol] = findR(rEqu,params);
  % get invariant distribution:
  [TEndog,Sol.invD] = stationaryDistribIter(Sol);
end

function [resid,Sol] = findR(r,params)
  persistent Vstart;
  if(isempty(Vstart))
    Vstart = zeros(params.nK,params.nZ);
  end
  % 1) prepare parameters:
  params.r = r;
  % production function:
  % use y = K^alpha L^(1-alpha)
  % r + delta = alpha (K/L)^(alpha-1)
  % wage = (1-alpha) (K/L)^alpha
  % parameters:
  delta = 0.1;
  alpha = 0.36;
  % from input r:
  KLratio = ((r+delta)/alpha)^(1/(alpha-1));
  params.wage = (1-alpha)* KLratio^alpha;
  % set Kmax relative to wage:
  Kmax = params.Kmax_factor*params.wage; 
  % 2) compute HH solution
  Sol = consinterp(Vstart,params,Kmax,100);
  Vstart = Sol.V; % take current solution as starting point for iteration next time
  % 3) compute stationary distribution of capital and productivity
  [TEndog,invD] = stationaryDistribIter(Sol);
  % 4) compute aggregate capital
  K = sum(Sol.gridEndog' * invD);
  % marginal productivity of capital
  Laggr = dot(Sol.gridExog,Sol.invdExog);
  % 5) compute interest rate from the marginal productivity of capital
  r2  = alpha * (K/Laggr)^(alpha-1) - delta;

  resid = r2 - r;
  fprintf(1,'Interest rate = %f; capital = %f; resid = %e\n',r,K,resid);
end


% Computes transition matrices and invariant cross-sectional distribution
%   in a model with one endogenous and one exogenous state
% Solution is captured in the structure Sol, which is input
function [TEndog,invD] = stationaryDistribIter(Sol)
  [nEndog,nExog] = size(Sol.V);
  nn = nEndog*nExog;
  % for each exogenous state there is a different endogenous transition:
  global TEndog TExog
  TExog = Sol.TransExog;
  TEndog = sparse(nn,nn);
  ee = (1:nEndog)';
  for j=1:nExog
    ii = (j-1)*nEndog+1 : j*nEndog; % index of endogenous grid points
    if(isfield(Sol,'pHiOpt'))  % endogenous transition as it arises from linear interpolation
      jj = Sol.IndxOpt(:,j);
      pHi = Sol.pHiOpt(:,j);
      TEndog(ii,ii) = sparse([jj;jj+1],[ee;ee],[1-pHi;pHi],nEndog,nEndog);
    else  % endogenous transition is deterministic
      TEndog(ii,ii) = sparse(Sol.IndxOpt(:,j),1:nEndog,ones(nEndog,1),nEndog,nEndog);
    end
  end
  assert(all(abs(sum(TEndog)-1)<1e-12))

  % start iteration with equal distribution:
  D = ones(nn,1) / nn;
  for iter = 1:10000; % maximum of 10000 iterations
    Dmid = TEndog*D;
    D2 = vec(reshape(Dmid,nEndog,nExog) * TExog);
    sumD = sum(D2);
    assert(abs(sumD-1)<1e-12);
    D2 = D2 / sumD;
    dist = max(abs(D-D2));
    D = D2;
    if(dist<1e-10)
      invD = reshape(D,nEndog,nExog);
      return;
    end
  end
  fprintf(1,'distance = %e',dist)
  warning('maximum iterations exceeded')
  invD = reshape(D,nEndog,nExog);
end

