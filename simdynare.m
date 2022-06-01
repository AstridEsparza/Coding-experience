% file to simulate model solved by dynare, as stored in dynare output oo_
% simulated for T periods, with stdeviation stdev given as column vector
% Michael Reiter, IHS, Feb. 2017
function simseries = simdynare(oo_,T,stdev)
  [A,B] =  dynare2ab(oo_);
  nz = size(B,2);
  eps = randn(nz,T) .* repmat(stdev,1,T);
  simseries = zeros(size(A,1),T);
  for t=1:T
    simseries(:,t+1) = A*simseries(:,t) + B*eps(:,t);
  end
  simseries = simseries';
end

