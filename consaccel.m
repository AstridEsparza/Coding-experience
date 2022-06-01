% solve HH infinite horizon saving problem on finite grid, using acceleration steps
% labor income follows Markov process
% Michael Reiter, IHS, February 2017
% last update: February 2018
function Sol = consaccel( nK, nZ, Kmax,nAccel)

  Kmin = 0;
  gridK = linspace(Kmin,Kmax,nK)'; % finite grid of capital

  r = 0.02;
  beta = 0.97;
  % finite state Markov process for productivity:
  m = 2.5;   %width of grid in terms of standard deviations
  rho = 0.9;  % autocorrelation of AR(1)
  sigma = 0.05;  % standard deviation of shock
  [TransZ,gridlogZ,invariantDistribution]=markovappr(rho,sigma,m,nZ);
  gridZ = exp(gridlogZ);  % level of productivity (wage)
  TransZTransp = TransZ';

  Vcont = zeros(nK,nZ);
  Vend = zeros(nK,nZ);
  Vtry = NaN(nK,1);
  IndxOptSaving = zeros(nK,nZ);
  for iter = 0:10000  % maximum number of steps = 10000;start with iter=0 to enforce optimization
    Vend = beta*Vcont*TransZTransp;
    if(mod(iter,nAccel)==0)
      for j=1:nZ % solve problem on each exogenous grid point 
        for i=1:nK % solve problem on each endogenous grid point 
          Kbeg = gridK(i);
          % try all possible levels of saving at once:
          Kend = gridK; % vector of all possible choices
          c = (1+r)*Kbeg + gridZ(j) - Kend; %vector of consumption levels
          Vtry = utilCheck(c) + Vend(:,j);
          [V(i,j),jOpt] = max(Vtry);
          IndxOptSaving(i,j) = jOpt;
        end
      end
      dist = max(max(abs(V-Vcont)));
      fprintf(1,'iter = %d; dist = %e\n',iter,dist);
      if(dist<1e-8)
        break;
      end
    else
      for j=1:nZ % solve problem on each grid point (state)
        for i=1:nK % solve problem on each grid point (state)
          Kbeg = gridK(i);
          % use best policy of earlier iteration:
          iOpt = IndxOptSaving(i,j);
          Kend = gridK(iOpt);
          c = (1+r)*Kbeg + gridZ(j) - Kend; %vector of consumption levels
          V(i,j) = util(c) + Vend(iOpt,j); % positive c is guaranteed
        end
      end
    end
    Vcont = V;
  end
  Sol = struct( 'V',V,'IndxOpt',IndxOptSaving,'gridEndog',gridK,'gridExog',gridZ,'TransExog',TransZ);
end


function u = util(c)
  u = log(c);
end
function u = utilCheck(c)
  u = -Inf(size(c));
  indx = c>0;
  u(indx) = util(c(indx));
end
