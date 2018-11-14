function [lc,p,tc]=lc_expgo(genelist,genome,expdata,disttype,lcfun,annot)

%Calculates local correlation coefficients between expression data and GO
%annotation for a set of genes in genelist (cell array of strings)
%belonging to yeast (genome=1), mouse (genome=2), arabidopsis (genome=3) or
%else (genome=4).
%Expdata should be a gene expression data matrix, with the same number of
%rows as genelist, being each row correspondent to one gene, and as many
%collumns as the number of conditions under which gene expression was
%measured. Disttype has to be one of two strings: 'euc', for euclidean
%distances between expression profiles, and 'corr', for correlation
%distances between expression profiles. Lcfun has to be one of three
%strings: 'ratio','pearson' or 'score', identifying the desired local
%correlation coefficient to be used.
%The user may supply his own annotation in the annot input variable (see the 
%manual for the format of this variable).

%Francisco Rodrigues Pinto, Oeiras, 2003

d1=lc_dist(expdata,disttype);

n=size(genelist,1);

[semantic_dist,annot]=lc_semantic(genelist,genome,annot);

expmat=squareform(d1);

withtermset=unique(annot.f(:,1));
mat2=squareform(semantic_dist.f);
mat2=mat2(withtermset,withtermset);
mat1=expmat(withtermset,withtermset);
[minlc,minp,tc.f]=lc_calc(mat1,mat2,lcfun);
lc.f=zeros(n,1);
lc.f(withtermset,1)=minlc;
p.f=0.5*ones(n,1);
p.f(withtermset,1)=minp;


withtermset=unique(annot.p(:,1));
mat2=squareform(semantic_dist.p);
mat2=mat2(withtermset,withtermset);
mat1=expmat(withtermset,withtermset);
[minlc,minp,tc.p]=lc_calc(mat1,mat2,lcfun);
lc.p=zeros(n,1);
lc.p(withtermset,1)=minlc;
p.p=0.5*ones(n,1);
p.p(withtermset,1)=minp;

withtermset=unique(annot.c(:,1));
mat2=squareform(semantic_dist.c);
mat2=mat2(withtermset,withtermset);
mat1=expmat(withtermset,withtermset);
[minlc,minp,tc.c]=lc_calc(mat1,mat2,lcfun);
lc.c=zeros(n,1);
lc.c(withtermset,1)=minlc;
p.c=0.5*ones(n,1);
p.c(withtermset,1)=minp;
