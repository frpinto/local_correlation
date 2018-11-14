function [lc]=lc_pearson(d1,d2)

%Calculates local pearson correlation coefficient for a set of objects
%characterized by distance vectors (pdist format) d1 and d2

%Francisco Rodrigues Pinto, Oeiras, 2003

d1=checkrow(d1);
d2=checkrow(d2);

mat1=squareform(d1);
mat2=squareform(d2);
n=length(mat1);

meanmat1=repmat(mean(mat1),n,1);
meanmat2=repmat(mean(mat2),n,1);
mat1=mat1-meanmat1;
mat2=mat2-meanmat2;
corr=sum((mat1.*mat2),1)./((sum((mat1.^2),1).*sum((mat2.^2),1)).^0.5);
lc=corr';
   