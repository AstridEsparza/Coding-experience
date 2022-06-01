% Solve model of Gomes, AER 2001, "Financing Investment", by DP
% Michael Reiter, this version February 2020
function Sol = fininv( nK, nZ, Kmax,nAccel)

  global A alphak alphal delta rho sigma lambda0 lambda1 betta f w;
  A = 1;
  Kmin = 0;
  gridK = logspaceshift(Kmin,Kmax,nK,0.01);
  % some variables
  f=0.001;
  beta = 1/1.065;
  lambda0=0.001; % 0.08;
  lambda1=0.028;
  alphak=0.3;
  alphal=0.65;  
  delta=0.145;
  w=1;

  % finite state Markov process for productivity:
  m = 2.5;   %width of grid in terms of standard deviations
  rho = 0.762;  % autocorrelation of AR(1)
  sigma = 0.0352;  % standard deviation of shock

  % TransZ is the transition matrix of the Markov chain
  % gridlogZ is the discretized state space of z_t
  [TransZ,gridlogZ,invariantDistribution]=markovappr(rho,sigma,m,nZ);
  gridZ = exp(gridlogZ);  % level of productivity (z, vector), 
  TransZTransp = TransZ';

  Vcont = zeros(nK,nZ);
  Vend = zeros(nK,nZ);
  Vtry = NaN(nK,1);
  IndxOptSaving = zeros(nK,nZ);
  optCF = zeros(nK,nZ);
  for iter = 0:10000  % maximum number of steps = 10000;start with iter=0 to enforce optimization
    Vend = beta*Vcont*TransZTransp;
    if(mod(iter,nAccel)==0)
      for j=1:nZ % solve problem on each exogenous grid point 
        for i=1:nK % solve problem on each endogenous grid point 
          Kbeg = gridK(i);
          % try all possible levels of saving at once:
          Kend = gridK; % vector of all possible choices

          z = gridZ(j);
          pr = profits(Kbeg,z);
          ic = invcosts(Kbeg,Kend);
          lam = lambdafunc(Kbeg,Kend,z);
          % cash flow:
          CFtry = pr - ic - lam;

          Vtry = CFtry + Vend(:,j);

          [V(i,j),jOpt] = max(Vtry);
          IndxOptSaving(i,j) = jOpt;
          optCF(i,j) = CFtry(jOpt);
        end
      end
      dist = max(max(abs(V-Vcont)));
      fprintf(1,'iter = %d; dist = %e\n',iter,dist);
      if(dist<1e-8)
        break;
      end
    else
      for j=1:nZ 
        Kbeg = gridK;
        iOpt = IndxOptSaving(:,j);
        Kend = gridK(iOpt);
        V(:,j) = optCF(:,j) + Vend(iOpt,j); % positive c is guaranteed
      end
    end
    Vcont = V;
  end
  Sol = struct( 'V',V,'IndxOpt',IndxOptSaving,'gridEndog',gridK,'gridExog',gridZ,'TransExog',TransZ,'mpk',margprodk(gridK,gridZ(4)));
end

function pr =  profits(k,z)
  global A alphak alphal w f;
  Ak = A*z.*(k.^alphak);
  l = (alphal.*Ak/w).^(1/(1-alphal)); % optimal labor
  pr = Ak.*l.^alphal - w*l - f; % f: fixed cost of production
  pr = max(pr,0);  % option of not producing; no fixed cost paid
end
% marginal product of capital:
function mpk =  margprodk(k,z)
  global A alphak alphal w f;
  Ak = A*z.*(k.^alphak);
  l = (alphal.*Ak/w).^(1/(1-alphal)); % optimal labor
  mpk = alphak*Ak .* l.^alphal ./ k; 
end

function ic = invcosts(k, knext)
  global delta;
  ic = knext - (1-delta)*k;
  indxNeg = knext<(1-delta)*k;
  salesPrice = 0.5;
  ic(indxNeg) = salesPrice*ic(indxNeg);
end
function lam =  lambdafunc(k,knext,z) % financing cost
  global lambda0 lambda1;
  pr = profits(k,z);
  ExternFunds = invcosts(k,knext) - pr;
  indx =  ExternFunds>0;  % indicator whether I need external financing
  lam = indx.*(lambda0 + lambda1*ExternFunds);
end
