function r=total_pearson(d1,d2)

%Calculates the total pearson correlation coefficient between distance matrices
%(in pdist vector format) d1 and d2 

%Francisco Rodrigues Pinto, Oeiras, 2003

d1=checkrow(d1)';
d2=checkrow(d2)';

n=length(d1);

m1=mean(d1);
m2=mean(d2);

r=sum((d1-m1).*(d2-m2))/((sum((d1-m1).^2)*sum((d2-m2).^2))^0.5);
