% -------------------------------------------------------------
% functions defining the jacobian of the BVP for limitcycles
% -------------------------------------------------------------

% integral contraint
function ic = BVP_PD2_jac_ic()
global lds
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  p = lds.dt(j)*(lds.vpoldp(:,range1).*lds.pwi);
  ic(range2) = ic(range2)+p(lds.cols_p1_coords);
  range1 = range1 + lds.ncol;
  range2 = range2 + lds.ncol_coord;
end

