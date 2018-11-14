function godistsq=mingodistmat(unigodistsq,annot,n)

%this function converts the semantic distance square matrix that was built
%by lc_semantic.m, that has distance values between GO terms, into a square
%semantic distance matrix with distance values between genes.

%Francisco Rodrigues Pinto, Oeiras, 2003

uniterms=unique(annot(:,2));
godistsq=zeros(n);
for i=1:n
    seti=annot(find(annot(:,1)==i),2);
    for j=(i+1):n
        setj=annot(find(annot(:,1)==j),2);
        goij=unigodistsq(find(ismember(uniterms,seti)),find(ismember(uniterms,setj)));
        if isempty(goij)
            godistsq(i,j)=666;
        else
            godistsq(i,j)=min(goij(:));
        end
    end
end

godistsq=godistsq+godistsq';