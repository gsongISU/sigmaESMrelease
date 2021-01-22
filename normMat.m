function amp = normMat(x, dim)
% compute the norm of matrix x along dimension dim.
amp = sqrt(sum(x.^2,dim));
