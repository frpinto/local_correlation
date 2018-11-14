function [semantic_dist, annot]=lc_semantic(genelist,genome,annot)

%Calculates semantic distance vector (pdist format) between genes in
%genelist, from yeast (genome=1), mouse (genome=2), arabidopsis
%(genome=3) or else (genome=4). The user may supply his own annotation 
%in the annot input variable (see the manual for the format of this variable).

%Francisco Rodrigues Pinto, Oeiras, 2003

if genome==1
    load scchip;
    genomechip=scchip;
    clear scchip;
elseif genome==2
    load mmchip;
    genomechip=mmchip;
    clear mmchip;
elseif genome==3
    load atchip;
    genomechip=atchip;
    clear atchip;
elseif genome==4
    load genomechip;
end

n=length(genelist);

if nargin==2
    index_list=convert_list(genelist,genomechip);
    annot=search_annot(index_list,genomechip);
end

if genome==4
    ngenes=1;
else
    ngenes=length(genomechip.affyref);
end

load function_mat;
uniterms=unique(annot.f(:,2));
notvalid=find(ismember(annot.f(:,2),setdiff(uniterms,function_mat.id)));
annot.f(notvalid,:)=[];
jiangdist.f=lc_jiang(unique(annot.f(:,2)),function_mat,genomechip.all.f,ngenes);
clear function_mat;
load process_mat;
uniterms=unique(annot.p(:,2));
notvalid=find(ismember(annot.p(:,2),setdiff(uniterms,process_mat.id)));
annot.p(notvalid,:)=[];
jiangdist.p=lc_jiang(unique(annot.p(:,2)),process_mat,genomechip.all.p,ngenes);
clear process_mat;
load component_mat;
uniterms=unique(annot.c(:,2));
notvalid=find(ismember(annot.c(:,2),setdiff(uniterms,component_mat.id)));
annot.c(notvalid,:)=[];
jiangdist.c=lc_jiang(unique(annot.c(:,2)),component_mat,genomechip.all.c,ngenes);
clear component_mat;

dist.f=mingodistmat(squareform(jiangdist.f),annot.f,n);
dist.p=mingodistmat(squareform(jiangdist.p),annot.p,n);
dist.c=mingodistmat(squareform(jiangdist.c),annot.c,n);

semantic_dist.f=rapidunsquare(dist.f);
semantic_dist.p=rapidunsquare(dist.p);
semantic_dist.c=rapidunsquare(dist.c);

function jiangdist=lc_jiang(uniterms,gomat,freq,ngenes)

n=length(uniterms);
freq=freq/ngenes;
set=find(ismember(gomat.id,uniterms));
setfreq=freq(set);
setfreq=setfreq+(setfreq==0).*(1/ngenes);
distmat=zeros(n,n);

for i=1:(n-1)
    mat=gomat.children(:,set((i+1):end)).*repmat(gomat.children(:,set(i)),1,n-i)...
        .*repmat(freq,1,n-i);
    mat=mat+(mat==0);
    parents=min(mat);
    distmat(i,(i+1):end)=2*log(parents)-log(setfreq((i+1):end)')-log(setfreq(i));
end
jiangdist=rapidunsquare(distmat+distmat');

    

