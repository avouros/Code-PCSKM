function [idx,centroids,w,niter,C,iflag] = pcskmeans(x,k,s,constr,init_centers,varargin)
%SEMISUPERVISED_SPARSE_KMEANS implements a semi-supervised version of the 
%sparse k-means algorithm [1].
% Original R code can be found here: https://bit.ly/2PUFKim
% The MATLAB implementation uses the method of [2] to obtain the initial
% centroids in a deterministic way.

% References:

% [1] Witten, Daniela M., and Robert Tibshirani. "A framework for feature 
%     selection in clustering." Journal of the American Statistical 
%     Association 105.490 (2010): 713-726.

% [2] N. Nidheesh, K.A. Abdul Nazeer, P.M. Ameer, 
%     An enhanced deterministic K-Means clustering 
%     algorithm for cancer subtype prediction from 
%     gene expression data, Computers in Biology
%     and Medicine 91C (2017) pp. 213-221.
%     https://doi.org/10.1016/j.compbiomed.2017.10.014


% Input:
% - x : a matrix where rows are observations and columns are attributes.
% - k : number of target clusters.
% - s : sparsity parameter.
% - labels: element indexes and their respective clusters
% - init_centers: initial centroids, if it is empty then DKM++ init is used
% - varargin :

% Output:
% - idx       : vector specifying the cluster of each element.
% - centroids : final location of the centroids.
% - w         : vector specifying the final weight of each feature.
% - niter     : number of iterations until converge.
% - C         : initial centroids
% - iflag     : 0 = n/a, 1 = empty cluster(s), 2 = negative/zero obj func,
%                        3 = negative/nan weights
%               OK: -1 = weights, -2 = no re-assignments, -3 = max iters                        

    
    iflag = 0;
    % Algorithm properties
    metric = 'squaredeuclidean';  
    iters = 25; %sparse k-means iterations
    iterk = 25; %k-means iterations
    % Display
    DISPLAY = 0;
    
    % Custom parameters
    for i = 1:length(varargin)
        if isequal(varargin,'metric')
            metric = varargin{i+1};
            if ~isequal(metric,'squaredeuclidean') || ~isequal(metric,'euclidean')
                error('Euclidean or Euclidean^2 is allowed as metric')
            end
        elseif isequal(varargin{i},'iters')    
            iters = varargin{i+1};        
        elseif isequal(varargin{i},'iterk')    
            iterk = varargin{i+1};                        
        % Display    
        elseif isequal(varargin{i},'DISPLAY')   
            DISPLAY = varargin{i+1};        
        end
    end
    
    % Initialize
    [n,p] = size(x);
    niter = 0; %number of iterations
    
    % Initialize feature weights 
    w = (1/sqrt(p)) * ones(1,p);
	%w = ones(1,p);
    
    % Generate initial centroids if needed
    if isempty(init_centers)
        % Use DKM++ to init the centroids
        C = dkmpp_init(x,k);
        init_centers = x(C,:);
    else
        C = init_centers;
    end    
        
    % Execute PCK-Means (K-Means with constraints)
    if DISPLAY
        fprintf('\nK-Means iter %d',1);
    end    
    [idx,centroids,~,~,~] = mpckmeans(x,k,constr,...
        'centers',init_centers,'iterations',iterk,'metric_learning',0,...
        'transitiveML',0,'transitiveCL',0,'gap_iterations',1);   
    
    % Check if we have empty cluster(s)
    if length(unique(idx)) ~= k
        iflag = 1;
        return
    end
    
    % LOOP UNTIL CONVERGED
    idx_old = idx;
    while niter < iters
        niter = niter + 1;
        w_old = w;
        
        %% Update the clusters (after the first iteration)
        if niter > 1
            % Scale each feature by w
            wx = x .* repmat(sqrt(w),n,1);
            % Compute the new centroids
            for ii = 1:k
                elements = find(idx == ii);
                if ~isempty(elements)
                    centroids(ii,:) = mean(wx(elements,:),1);
                elseif length(elements) == 1 %if only 1 element
                    centroids(ii,:) = wx(elements,:); 
                elseif isempty(elements) %if no elements    
                    %Empty cluster(s)
                    iflag = 1;
                    return
                end                
            end
            % Do PCK-Means clustering using the new computed centroids as init centroids 
            if DISPLAY
                fprintf('\nK-Means iter %d',niter);
            end              
            [idx,centroids,~,~,~] = mpckmeans(wx,k,constr,...
                'centers',centroids,'iterations',iterk,'metric_learning',0,...
                'transitiveML',0,'transitiveCL',0,'gap_iterations',1);               
            % Check if we have empty cluster(s)
            if length(unique(idx)) ~= k
                iflag = 1;
                return
            end         
        end
        
        
        %% Update the weights
        % Compute within-cluster distance (per feature)
        WCSSp = zeros(1,p);
        for ii = 1:k
            switch metric
                case 'squaredeuclidean'
                    WCSSp = WCSSp + sum(normalizations(x(idx == ii,:),'mean').^2);
                case 'euclidean'
                    WCSSp = WCSSp + sum(sqrt(normalizations(x(idx == ii,:),'mean').^2) );
                otherwise
                    error('Undefined metric')
            end
        end
        % Compute global within-cluster distance (per feature)
        switch metric
            case 'squaredeuclidean'
                TSSp = sum(normalizations(x,'mean').^2);
            case 'euclidean'
                TSSp = sum( sqrt(normalizations(x,'mean').^2) );
            otherwise
                error('Undefined metric')
        end  
        % Check if objective function is negative
        if -WCSSp+TSSp <= 0 
            iflag = 2;
        end
        if DISPLAY
            fprintf('\nObjective function: %d',-WCSSp+TSSp);
        end        
        % Find Delta using binary search
        delta = BinarySearchDelta(-WCSSp+TSSp , s);
        % Compute the new weights
        w_tmp = sign(-WCSSp+TSSp) .* max(0,abs(-WCSSp+TSSp)-delta);     
        w = w_tmp/norm(w_tmp,2);
        % Check if weights are negative or nan
        if ~isempty(find(w<0 | isnan(w) | isinf(w)))
            iflag = 3;
            return
        end         
        
        
        %% Convergence
        if (sum(abs(w-w_old)) / sum(abs(w_old))) <= 10^-4
            if DISPLAY
                disp('Converged (weights)');
            end
            iflag = -1;
            return
        end
%         if isequal(idx_old,idx)
%             if DISPLAY
%                 disp('Converged (no re-assignments)');
%             end
%             iflag = -2;
%             return            
%         end
    end
    if DISPLAY
        disp('Maximum number of iterations reached');
    end    
    iflag = -3;
end

