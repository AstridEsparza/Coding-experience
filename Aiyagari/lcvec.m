% solve deterministic HH lifecycle problem on a FINITE GRID of capital, simplest version
% Problem: 
%    max sum_{t=1}^T beta^t u(c_t)
%      subject to
%    c = (1+r)*Kbeg + y(t) - Kend;
%    Kend >= 0
% Michael Reiter, IHS, February 2017; this version November 2020
% Input arguments:
%   Kmax: maximal capital level, must be checked that it is big enough (not binding)
%   nK:   size of capital grid; higher nK gives higher precision
%
%   INCREASE SPEED BY VECTORIZATION!
function [ V ,indxOptSaving,gridK] = lcvec(Kmax,nK)

  r = 0.02;
  beta = 0.97;
  gridK = linspace(0,Kmax,nK)';
  y = [ones(40,1);zeros(20,1)];
  T = length(y);

  % value function: sum of all discounted value achievable in the future, given level of wealth
  V = NaN(nK,T+1);
  indxOptSaving = zeros(nK,T);
  vTry = NaN(nK,1);

  V(:,T+1) = 0;
  for t=T:-1:1 % solve backward in time
    fprintf(1,'now solving period %d\n',t);
    for i=1:nK  % solve separately for each state
      Kbeg = gridK(i);
      % eliminate innermost loop by a vector operation:  VECTORIZATION increases speed 
      Kend = gridK;  % try all possible values simultaneously
      c = (1+r)*Kbeg + y(t) - Kend;
      vTry = util(c) + beta*V(:,t+1);

      [V(i,t),indxOptSaving(i,t)] = max(vTry);
    end
  end

end


function u = util(c)
  u = -Inf(size(c));  % default: - infinity
  indxPos = c>0;
  u(indxPos) = log(c(indxPos));
end

