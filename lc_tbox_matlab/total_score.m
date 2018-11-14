function r=total_score(d1,d2)

%Calculates the total score correlation coefficient between distance matrices
%(in pdist vector format) d1 and d2 

%Francisco Rodrigues Pinto, Oeiras, 2003


d1=checkrow(d1)';
d2=checkrow(d2)';

n=length(d1);

[a,index]=sort(d1);
maxord=turncolrank(d2(index));
minord=flipud(turncolrank(d2(flipud(index))));

maxord=maxord.*((maxord-(1:n)')>0)+(1:n)'.*((maxord-(1:n)')<=0);
minord=minord.*((minord-(1:n)')>0)+(1:n)'.*((minord-(1:n)')<=0);

maxord=full(sparse(maxord,1,1,n,1));
minord=full(sparse(minord,1,1,n,1));

maxord=cumsum(maxord);
minord=cumsum(minord);
maxord=sum(maxord);
minord=sum(minord);
commonneibsum=(maxord+minord)/2;

if n/2==round(n/2)
    minsum=sum(2:2:n);
else
    minsum=sum(2:2:(n-1))+n;
end

r=(commonneibsum-minsum)/(sum(1:n)-minsum);
