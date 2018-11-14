function dist=lc_dist(source,dist_type);

%generates vector of distances (pdist format) between rows of the data
%matrix: source. dist_type can be 'euc' or 'corr', for euclidean or
%correlation distances rescpectively.

%Francisco Rodrigues Pinto, Oeiras, 2003

n=size(source,1);

if strcmpi(dist_type,'euc')
    
    dist=(sum(((source(2:end,:)-repmat(source(1,:),(n-1),1)).^2),2)).^(0.5);
    for i=2:(n-1)
        dist=[dist;...
      (sum(((source((i+1):end,:)-repmat(source(i,:),(n-i),1)).^2),2)).^(0.5)];
    end
    dist=dist';
elseif strcmpi(dist_type,'corr')
    mean_source=mean(source')';
    dev_source=(source-repmat(mean_source,1,size(source,2)));
    sum_dev_source=sum(dev_source.^2,2);
    dist=sum((dev_source(2:end,:).*repmat(dev_source(1,:),(n-1),1)),2)./...
        ((sum_dev_source(2:end)*sum_dev_source(1)).^(0.5));
    for i=2:(n-1)
        dist=[dist; sum((dev_source((i+1):end,:).*repmat(dev_source(i,:),(n-i),1)),2)./...
        ((sum_dev_source((i+1):end)*sum_dev_source(i)).^(0.5))];
    end
    dist=1-dist';
end


