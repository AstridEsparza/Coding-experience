% return from Dynare output the effect of the variance on the variables 
% in the order they are declared in the mod-file
% Michael Reiter, IHS, Feb. 2017
function dx =  dxdsigma(oo_)
  d = oo_.dr;
  iiV = d.inv_order_var;
  dx = 0.5 * d.ghs2(iiV);  % ghs2 stores 0.5 times the effect of sigma^2
end
