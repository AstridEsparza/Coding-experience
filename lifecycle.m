% solve deterministic HH lifecycle problem on a FINITE GRID of capital, simplest version
% Problem: 
%    max sum_{t=1}^T beta^t u(c_t)
%      subject to
%    c = (1+r)*Kbeg + y(t) - Kend;
%    Kend >= 0
% Michael Reiter, IHS, February 2017; this version February 2020
% Input arguments:
%   Kmax: maximal capital level, must be checked that it is big enough (not binding)
%   nK:   size of capital grid; higher nK gives higher precision
function [ V,IndxOptSaving,gridK ] = lifecycle(Kmax, nK)

  Kmin = 0;  % borrowing constraint; will bind only at end of life
  gridK = linspace(Kmin,Kmax,nK)'; % finite grid of capital

  r = 0.02;
  beta = 0.97;
  % 40 years of work, 20 years of retirement:
  y = [ones(40,1);zeros(20,1)];
  T = length(y);

  V = NaN(nK,T+1);
  V(:,T+1) = 0;
  Vtry = NaN(nK,1);
  for t=T:-1:1  % go backwards in time
    fprintf(1,'now solving period %d\n',t);
    for i=1:nK % solve problem on each grid point (state)
      Kbeg = gridK(i);
      for j=1:nK  % try all possible levels of saving (decision)
        Kend = gridK(j);
        c = (1+r)*Kbeg + y(t) - Kend;
        Vtry(j) = util(c) + beta*V(j,t+1);
      end
      [V(i,t),IndxOptSaving(i,t) ] = max(Vtry);
    end
  end
end


function u = util(c)
  if(c<=0)
    u = -Inf;
  else
    u = log(c);
  end
end
