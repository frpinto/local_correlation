function mat=turncolrank(mat)

%this function transforms a matrix of values into a matrix of corresponding
%ranks inside each collumn.

%Francisco Rodrigues Pinto, Oeiras, 2003

for i=1:size(mat,2)
    mat(:,i)=turnrank(mat(:,i));
end
    
function rnk=turnrank(vec)

vec=[vec (1:length(vec))'];
vec=sortrows(vec,1);
vec=[vec (1:length(vec))'];
vec=sortrows(vec,2);
rnk=vec(:,3);