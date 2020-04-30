function [idx,centroids,iterations,iflag] = kmeans_lloyd(x,k,init_centroids,ITER,varargin)
%KMEANS_LLOYD carries out the K-means algorithm by Stuart P. Lloyd

%INPUT:
% x: datapoins; rows = observations, columns = attributes
% k: number of target clusters
% init_centroids: initial centroids positions; if empty then random positions
% iterations: number of iterations

%OUTPUT
% idx: vector; the cluster of each datapoint.
% centroids: final centroids locations.
% iterations: number of iterations before converged.
% iflag: 2 = no converged, 1 = empty cluster(s), 0 = converged

    % Assert initial centroids
    s = size(init_centroids);
    assert(s(1)==k);
    assert(s(2)==size(x,2));  
    
    dmetric = 'squaredeuclidean';
    for i = 1:length(varargin)
        if isequal(varargin,'metric')
            dmetric = varargin{i+1};
        end
    end
    
    % Initialise
    iflag = 0;
    ndata = size(x,1);
    ndim = size(x,2);
    idx = zeros(ndata,1);
    iterations = 1;
    
    % Create initial clusters
    % Assign datapoints to the nearest centroids
    dist = pdist2(x,init_centroids,dmetric);
    [~,idx] = min(dist,[],2);
    
    % Compute initial clusters
    centroids = nan(k,ndim);
    for ii = 1:k
        elements = find(idx == ii);
        if length(elements) > 1
            centroids(ii,:) = mean(x(elements,:));
        elseif length(elements) == 1 %if only 1 element
            centroids(ii,:) = x(elements,:);             
        elseif isempty(elements) %if no elements
            centroids(ii,:) = centroids(ii,:); 
        end                
    end     
    if length(unique(idx)) ~= k
        iflag = 1;
        return
    end        

    
    % Lloyd's main loop
    old_centroids = centroids;
    
    for T = 1:ITER + 1
        % Assign datapoints to clusters
        dist = pdist2(x,centroids,dmetric);
        [~,idx] = min(dist,[],2);        
        
        % Recompute the centroids
        un = unique(idx);
        if length(un) ~=k 
            iflag = 1;
            disp('Empty cluster(s) detected. Terminate!')
            return            
        end
        for i = 1:k
            centroids(i,:) = mean(x(idx == un(i),:),1);
        end      
        
        % Number of iterations
        iterations = iterations + 1;
        
        % Check for convergence in centroids
        if isequal(old_centroids,centroids)
            return
        else
            old_centroids = centroids;
        end
        
        % Check for convergence in iterations
        if T == ITER+1
            iflag = 2;
            iterations = iterations-1;
        end       
    end
end