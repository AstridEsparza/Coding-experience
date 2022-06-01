% Changes randomly the state in a model with discrete state space
% Inputs:
%   icurr: index of current state (1 to n)
%   randu: uniform random number (between 0 and 1)
%   transprob: transition probability matrix (row=from; col=to)
% Michael Reiter, October 2005
function inext = changestaterandomly(icurr,randu, transprob)
  if(randu<0 && randu>1)
    error('invalid random number');
  end
  sum = 0;
  n = size(transprob,2);
  inext=n;
  for(i=1:n-1)
    sum = sum + transprob(icurr,i);  
    if(randu<=sum)
      inext = i;
      return;
    end;
  end;

