% solve HH infinite horizon saving problem on finite grid, using acceleration steps
% and linear interpolation of the value function
% labor income follows Markov process
% Michael Reiter, IHS, February 2017
% last update: March 2018
function Sol = consinterp(Vcont,params,Kmax,nAccel)

  [nK,nZ] = size(Vcont);
  Kmin = 0;
  gridK = linspace(Kmin,Kmax,nK)'; % finite grid of capital

  r = params.r;
  m = 2.5;   %width of grid in terms of standard deviations
  [TransZ,gridlogZ,invariantDistrib]=markovappr(params.rho,params.sigma,m,nZ);
  gridZ = exp(gridlogZ);
  TransZTransp = TransZ';

  VendDiscounted = zeros(nK,nZ);
  Vtry = NaN(nK,1);
  OptCons = zeros(nK,nZ);
  IndxOptSaving = zeros(nK,nZ);
  pHiOptSaving = zeros(nK,nZ);
  for iter = 0:10000  % maximum number of steps = 10000
    VendDiscounted = params.beta*Vcont*TransZTransp;
    if(mod(iter,nAccel)==0)
      for j=1:nZ % solve problem on each exogenous grid point 
        for i=1:nK % solve problem on each endogenous grid point 
          Kbeg = gridK(i);
          % try all possible levels of saving at once:
          [cMin,cMax] = cbounds(Kbeg,gridZ(j),gridZ,r,Kmin,Kmax);
          OptCons(i,j) = golden(@obj2max,cMin,cMax,Kbeg,j,gridK,gridZ,VendDiscounted,r);
          [V(i,j),IndxOptSaving(i,j), pHiOptSaving(i,j)] = ...
            obj2max(OptCons(i,j),Kbeg,j,gridK,gridZ,VendDiscounted,r);
        end
      end
      dist = max(max(abs(V-Vcont)));
      fprintf(1,'iter = %d; dist = %e\n',iter,dist);
      if(dist<1e-10)
        break;
      end
    else
      for j=1:nZ % solve problem on each grid point (state)
        for i=1:nK % solve problem on each grid point (state)
          Kbeg = gridK(i);
          % use best policy of earlier iteration:
          iOpt = IndxOptSaving(i,j);
          pHi = pHiOptSaving(i,j);
          Kend = gridK(iOpt);
          c = OptCons(i,j);
          V(i,j) = util(c) + (1-pHi).*VendDiscounted(iOpt,j) + pHi.*VendDiscounted(iOpt+1,j);
        end
      end
    end
    Vcont = V;
  end
  Sol = struct( 'V',V,'OptCons',OptCons,'IndxOpt',IndxOptSaving,'pHiOpt',pHiOptSaving,...
    'gridEndog',gridK,'gridExog',gridZ,'TransExog',TransZ,'invdExog',invariantDistrib);
end

% objective function
function [val,indx,pHi] = obj2max(c,Kbeg,j,gridK,gridZ,VendDiscounted,r)
  y = gridZ(j);
  Kend = (1+r)*Kbeg + y - c; 
  [indx,pHi] = lookup_equidist(gridK,Kend);
  % linear interpolation of VendDiscounted at Kend:
  val = util(c) + (1-pHi).*VendDiscounted(indx,j) + pHi.*VendDiscounted(indx+1,j);
end

function u = util(c)
  u = log(c);
end
function u = utilCheck(c)
  u = -Inf(size(c));
  indx = c>0;
  u(indx) = util(c(indx));
end
function [cMin,cMax] = cbounds(Kbeg,y,gridY,r,Kmin,Kmax)
  cMin = max(1e-100,(1+r)*Kbeg + y - Kmax); %cannot consume less
  cMax = (1+r)*Kbeg + y - Kmin;
end

% table lookup: find position of X in grid gridX;
% gridX must be EQUIDISTANT !!
% output: 
%     iPos: index of grid point left to X
%     weightHi: (X-gridX(iPos))./ (gridX(iPos+1) - gridX(iPos))
function [iPos weightHi] = lookup_equidist(gridX,X)
  n = length(gridX);
  stepx = gridX(2)-gridX(1);
  % position of largest grid point that is smaller than X:
  iPos = max(1,floor( (X-gridX(1))/stepx ) + 1);
  iPos = min(iPos,n-1); %iPos must not be n
  % probability to go to high end of interval:
  weightHi = (X-gridX(iPos))./stepx;
  % test whether X are between gridX(iPos) and gridX(iPos+1):
  % this tests whether X must be within grid
  assert(all(weightHi>=-1e-10 & weightHi<=1+1e-10));
end

% maximize function f by golden search method, within bounds (a,b)
function [x1,f1] = golden(f,a,b,varargin)

  tol = 1e-8; % tolerance level

  alpha1 = (3-sqrt(5))/2;
  alpha2 = (sqrt(5)-1)/2;
  d  = b-a;
  x1 = a+alpha1*d;
  x2 = a+alpha2*d;
  f1 = feval(f,x1,varargin{:});
  f2 = feval(f,x2,varargin{:});

  d = alpha1*alpha2*d;
  while d>tol
    d = d*alpha2;
    if f2<f1 % x2 is new upper bound
      x2 = x1; x1 = x1-d; 
      f2 = f1; f1 = feval(f,x1,varargin{:});
    else     % x1 is new lower bound
      x1 = x2; x2 = x2+d; 
      f1 = f2; f2 = feval(f,x2,varargin{:});
    end
  end

  % Return the larger of the two
  if f2>f1, x1 = x2; f1 = f2; end
end
