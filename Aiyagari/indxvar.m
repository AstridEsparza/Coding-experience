function indx = indxvar(vname)
  global M_
  n = length(vname);
  nEndo = length(M_.endo_names(1,:));  % have all the same length
  for i=1:length(M_.endo_names)
    vn2 = M_.endo_names(i,:);
    if all(vn2(1:n)==vname) && (n==nEndo || vn2(n+1)==' ')
      indx = i;
      return;
    end
  end
  error(sprintf('model does not have a variable %s',vname))
end

