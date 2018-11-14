function index_list=convert_list(id_list,genomechip)

%this function converts a list of gene identifiers (id_list in cell array
%of strings format) into a list of corresponding indexes of those genes in
%the existing genomechip variable.

%Francisco Rodrigues Pinto, Oeiras, 2003


n=length(id_list);
index_mat=zeros(n,5);


    if isfield(genomechip,'affyref')
        searchlist=lower(genomechip.affyref);
        for i=1:n
            index=strmatch(lower(id_list{i}),searchlist,'exact');
            if isempty(index)==0
                index_mat(i,1)=index(1);
            end
        end
    end

    if isfield(genomechip,'sgd')
        searchlist=lower(genomechip.sgd);
        for i=1:n
            index=strmatch(lower(id_list{i}),searchlist,'exact');
            if isempty(index)==0
                index_mat(i,2)=index(1);
            end
        end
    end

    if isfield(genomechip,'name')
        searchlist=lower(genomechip.name);
        for i=1:n
            index=strmatch(lower(id_list{i}),searchlist,'exact');
            if isempty(index)==0
                index_mat(i,3)=index(1);
            end
        end
    end

    if isfield(genomechip,'locus')
        searchlist=lower(genomechip.locus);
        for i=1:n
            index=strmatch(lower(id_list{i}),searchlist,'exact');
            if isempty(index)==0
                index_mat(i,4)=index(1);
            end
        end
    end

    if isfield(genomechip,'ref')
        searchlist=lower(genomechip.ref);
        for i=1:n
            index=strmatch(lower(id_list{i}),searchlist,'exact');
            if isempty(index)==0
                index_mat(i,5)=index(1);
            end
        end
    end
    
    ngenes=length(searchlist);
    index_mat=index_mat+(index_mat==0)*2*ngenes;
    
    index_list=min(index_mat')';
    