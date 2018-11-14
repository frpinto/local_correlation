function [lc,p,tc]=lc_calc(d1,d2,lcfun)

%Calculates local correlation coefficients for the set of objects
%characterized by the distance matrices d1 and d2, using coefficient lcfun

%d1 and d2 can be symetric distance matrices or vectors (pdist format)

%lcfun should be one of the following strings: 'ratio','pearson' or 'score'

%lc is the list of local correlation coefficients

%p is the list of p values associated with local correlation coefficients

%tc is the total correlation between d1 and d2

%Francisco Rodrigues Pinto, Oeiras 2003

size1=size(d1);
if size1(1)==size1(2)
    d1=rapidunsquare(d1);
else
    d1=checkrow(d1);
end

size2=size(d2);
if size2(1)==size2(2)
    d2=rapidunsquare(d2);
else
    d2=checkrow(d2);
end

bign=length(d1);
n=(1+(1+8*bign)^(0.5))/2;
bootsize=ceil(20000/n);

if strcmpi(lcfun,'ratio')
    tc=total_ratio(d1,d2);
    lc=lc_ratio(d1,d2);
    for i=1:bootsize
        bootsamp=ceil(bign*rand(bign,1));
        empdist(:,i)=lc_ratio(d1(bootsamp),d2(bootsamp));
    end
    empdist=empdist(:);
    samplesize=length(empdist);
    for i=1:n
        p(i,1)=sum(empdist>lc(i))/samplesize;
    end
elseif strcmpi(lcfun,'pearson')
    tc=total_pearson(d1,d2);
    lc=lc_pearson(d1,d2);
    for i=1:bootsize
        bootsamp=ceil(bign*rand(bign,1));
        empdist(:,i)=lc_pearson(d1(bootsamp),d2(bootsamp));
    end
    empdist=empdist(:);
    samplesize=length(empdist);
    for i=1:n
        p(i,1)=sum(empdist>lc(i))/samplesize;
    end
elseif strcmpi(lcfun,'score')
    tc=total_score(d1,d2);
    lc=lc_score(d1,d2);
    for i=1:bootsize
        bootsamp=ceil(bign*rand(bign,1));
        empdist(:,i)=lc_score(d1(bootsamp),d2(bootsamp));
    end
    empdist=empdist(:);
    samplesize=length(empdist);
    for i=1:n
        p(i,1)=sum(empdist>lc(i))/samplesize;
    end
end