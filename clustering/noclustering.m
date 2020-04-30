function [idx,w,iterations] = noclustering(x,init_centroids,varargin)
%NOCLUSTERING assignes datapoints to their closer centroids 

    dmetric = 'squaredeuclidean';
    for i = 1:length(varargin)
        if isequal(varargin,'metric')
            dmetric = varargin{i+1};
        end
    end
    
    dist = pdist2(x,init_centroids,dmetric);
    [~,idx] = min(dist,[],2);
    w = ones(1,size(init_centroids,2));
    iterations = 1;      
end

