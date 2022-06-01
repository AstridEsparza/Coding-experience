% solve HH infinite horizon saving problem on finite grid
% labor income follows Markov process
% Michael Reiter, IHS, February 2017
% last update: February 2018
function Sol = consinfin( nK, nZ, Kmax)

  Kmin = 0;
  gridK = linspace(Kmin,Kmax,nK)'; % finite grid of capital

  r = 0.02;
  beta = 0.97;
  % finite state Markov process for productivity:
  m = 3.;   %width of grid in terms of standard deviations
  rho = 0.9;  % autocorrelation of AR(1)
  sigma = 0.05;  % standard deviation of shock
  [Pi,gridlogZ,invariantDistribution]=markovappr(rho,sigma,m,nZ);
  gridZ = exp(gridlogZ);  % level of productivity (wage)
  PiTransp = Pi';

  Vcont = zeros(nK,nZ); % value function of next period (as function of states)
  VcontExpected = zeros(nK,nZ);
  V = zeros(nK,nZ);
  Vtry = NaN(nK,1);
  IndxOptSaving = zeros(nK,nZ);
  for iter = 1:10000  % maximum number of steps = 10000  % think of it as iterating backward in time
    % expected continuation value, conditional on end-of-period states:
    % EXP( Vcont(i,znext) | zcurr = z_j) = sum_m Pi(j,m) Vcont(i,m)
    %                                       = Vcont(i,:)*Pi'  ; 
    VcontExpected = beta*Vcont*PiTransp; % discounted expected continuation value:
    for j=1:nZ % solve problem on each exogenous grid point 
      for i=1:nK % solve problem on each endogenous grid point 
        Kbeg = gridK(i);
        % try all possible levels of saving at once:
        Kend = gridK; % vector of all possible choices
        c = (1+r)*Kbeg + gridZ(j) - Kend; %vector of consumption levels
        Vtry = util(c) + VcontExpected(:,j);   % objective
        [V(i,j),kOpt] = max(Vtry);
        IndxOptSaving(i,j) = kOpt;
      end
      % Bellman operator maps Vcont into V;  V = Bellman(Vcont)
    end
    dist = max(max(abs(V-Vcont)));
    fprintf(1,'iter = %d; dist = %e\n',iter,dist);
    if(dist<1e-8) % criterion for convergence
      break;
    end
    Vcont = V; % step back one period; current value function becomes continuation value
  end
  Sol = struct( 'V',V,'IndxOpt',IndxOptSaving,'gridEndog',gridK,'gridExog',gridZ,'TransExog',Pi);
end


function u = util(c)
  u = -Inf(size(c));
  indx = c>0;
  u(indx) = log(c(indx));
end
