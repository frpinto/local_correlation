function mat=squareform(vec)

%this function transforms a distance vector (pdist format) into a square
%distance matrix.

%Francisco Rodrigues Pinto, Oeiras, 2003

bign=size(vec,2);

n=(1+(1+8*bign)^(0.5))/2;

mat=zeros(n,n);

for i=1:(n-1)
    mat(i,(i+1):end)=vec(1:(n-i));
    vec(1:(n-i))=[];
end

mat=mat+mat';
