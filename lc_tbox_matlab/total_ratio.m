function r=total_ratio(d1,d2)

%Calculates the total ratio correlation coefficient between distance matrices
%(in pdist vector format) d1 and d2 

%Francisco Rodrigues Pinto, Oeiras, 2003

d1=checkrow(d1)';
d2=checkrow(d2)';

n=length(d1);

r=(n-1)*sum(d1.*d2)/(sum(d1)*sum(d2));