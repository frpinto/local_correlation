function [lc]=lc_ratio(d1,d2)

%Calculates local ratio correlation coefficient for a set of objects
%characterized by distance vectors (pdist format) d1 and d2

%Francisco Rodrigues Pinto, Oeiras, 2003

d1=checkrow(d1);
d2=checkrow(d2);

mat1=squareform(d1);
mat2=squareform(d2);

n=length(mat1);
ratio=(n-1)*sum((mat1.*mat2),1)./(sum(mat1,1).*sum(mat2,1));

lc=ratio';