function Y=rapidunsquare(mat)

%Converts a square distance matrix into a vector of distances (pdist
%format), containing only the non repeated distance values.

%Francisco Rodrigues Pinto, Oeiras, 2003

b=tril(ones(size(mat)),-1);
Y=mat(find(b))';