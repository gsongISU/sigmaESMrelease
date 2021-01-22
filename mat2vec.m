function vec = mat2vec(mat)
% convert a matrix to a vector, extracting data *row-wise*.
vec = reshape(mat',1,prod(size(mat)));
% reshape extracts data column-wise, I would like it to 
% be row-wise. so mat' instead of mat.

%[row,col] = size(mat);
%for i=1:row
%    k = (i-1)*col;
%    vec(k+1:k+col) = mat(i,:);
%end

