function lc=lc_score(d1,d2)

%Calculates local score correlation coefficient for a set of objects
%characterized by distance vectors (pdist format) d1 and d2

%Francisco Rodrigues Pinto, Oeiras, December 2003

d1=checkrow(d1);
d2=checkrow(d2);

mat1=squareform(d1);
mat2=squareform(d2);

n=length(mat1);

[ordmat1,index]=sort(mat1);
udindex=flipud(index);
linindex=index+repmat((0:n:(n*(n-1))),n,1);
linudindex=udindex+repmat((0:n:(n*(n-1))),n,1);
ordmat2=reshape(mat2(linindex(:)),n,n);
bestmat=turncolrank(ordmat2);
ordmat2=reshape(mat2(linudindex(:)),n,n);
worstmat=turncolrank(ordmat2);
worstmat=flipud(worstmat);
refmat=repmat((1:n)',1,n);
bestmat=bestmat.*((bestmat-refmat)>0)+refmat.*((bestmat-refmat)<=0);
worstmat=worstmat.*((worstmat-refmat)>0)+refmat.*((worstmat-refmat)<=0);
refmat=repmat((1:n),n,1);
bestmat=full(sparse(bestmat(:),refmat(:),1,n,n));
worstmat=full(sparse(worstmat(:),refmat(:),1,n,n));
bestmat=cumsum(bestmat);
worstmat=cumsum(worstmat);
bestscore=sum(bestmat,1);
worstscore=sum(bestmat,1);
commonneibsum=(bestscore'+worstscore')/2;


if n/2==round(n/2)
    minsum=sum(2:2:n);
else
    minsum=sum(2:2:(n-1))+n;
end

lc=(commonneibsum-minsum)/(sum(1:n)-minsum);