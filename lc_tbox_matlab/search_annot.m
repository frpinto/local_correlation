function annot=search_annot(id_list,chipvar)

%this function searches the chipvar variable for the annotations of the
%genes in id_list (list of integers, each one is the index of the
%corresponding gene in the chipvar variable) and builds a stuctured
%variable, annot, with the annotations of those genes. annot has three
%fields: annot.f has the molecular function annotations; annot.p has the
%biological process annotations; annot.c has the cellular component
%annotations. Each of these three fields is a 2 collumn matrix. The first
%collumn contains indexes of genes (indexes for the id_list, not for the
%chipvar) and the second collumn contains corresponding GO codes with which
%the genes are annotated.

%Francisco Rodrigues Pinto, Oeiras, 2003


annot.p=[0 0];
annot.f=[0 0];
annot.c=[0 0];
for i=1:length(id_list)
    process_pos=find(chipvar.annot.p(:,1)==id_list(i));
    if isempty(process_pos)==0
        annot.p=[annot.p; i*ones(size(process_pos)) chipvar.annot.p(process_pos,2)];
    end
    
    function_pos=find(chipvar.annot.f(:,1)==id_list(i));
    if isempty(function_pos)==0
        annot.f=[annot.f; i*ones(size(function_pos)) chipvar.annot.f(function_pos,2)];
    end
    
    component_pos=find(chipvar.annot.c(:,1)==id_list(i));
    if isempty(component_pos)==0
        annot.c=[annot.c; i*ones(size(component_pos)) chipvar.annot.c(component_pos,2)];
    end
end
annot.p(1,:)=[];
annot.f(1,:)=[];
annot.c(1,:)=[];
