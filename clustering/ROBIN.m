function [C,lof] = ROBIN(x,k,nn,varargin)
%The ROBIN clustering initialisation method based on the work of [1].
%On each iteration the original algorithm described in [1] suggests to pick
%points with LOF "around" 1. In [2] this is specified as points with 
%"LOF < 1.05". This MATLAB implementation uses instead the bounds: 
%"1-critRobin <= LOF <= 1+critRobin" with "critRobin < 1". Incase no such
%points exist then "LOF <= 1+critRobin" is used.

%References:
%[1] Al Hasan, Mohammad, et al. "Robust partitional clustering by outlier 
%    and density insensitive seeding." Pattern Recognition Letters 30.11 
%    (2009): 994-1002.
%[2] Brodinova, Sarka, et al. "Robust and sparse k-means clustering for 
%    high-dimensional data." Advances in Data Analysis and Classification 
%    13.4 (2019): 905-932.


% Input:
% - x : NxP matrix, ows are observations and columns are attributes.
% - k : number of target clusters.
% - nn: number of nearest neighborrs. [1] and [2] are using 10 as default.
% - varargin:
%       - 'LOF', Nx1 or 1xN vector: in case the LOF scores have been
%       computed they  can be given as input arguments.
%       - 'critRobin', scalar: LOF threshold, [0,1); default is 0.05 [2].
%       - 'DETERMINISTIC', scalar:
%         0: Based on [2] the reference point is chosen at random (default)
%         1: Based on [1] the reference point is chosen at origin.
%         2: The reference point is the point with LOF closest 1.

% Output:
% - C : vector, indeces of x (datapoints) to be used as initial
%       centroids. 


    DETERMINISTIC = 0;
	critRobin = 0.05;
    lof = [];
    C = [];
        
    for iii = 1:length(varargin)
        if isequal(varargin{iii},'DETERMINISTIC')
            DETERMINISTIC = varargin{iii+1};
        elseif isequal(varargin{iii},'critRobin')
            critRobin = varargin{iii+1};
        elseif isequal(varargin{iii},'LOF')
            lof = varargin{iii+1};
        end
    end  
    
    if critRobin >= 1 || critRobin < 0
        error('Error: ROBIN threshold should be in the interval [0,1).')
    end
    
    % Compute distance matrix
    dists = squareform(pdist(x));
    
    % Compute LOF either based on the code or the paper
    if isempty(lof)
        lof = lof_paper(x,nn);   
    end
    
    % Select reference point
    switch DETERMINISTIC
        case 0
            %Non deterministic: (pick a random point as reference
            %Brodinova implementation
            n = size(x,1);
            r = randsample(n,1);
        case 1
            %Deterministic: pick origin as reference
            %Default method in Al Hasan et al. (2009)
            r = mean(x,1);
        case 2
            %Deterministic: pick point with LOF closest to the threshold
            r = abs(1-lof);
            [~,r] = min(r);
        otherwise
            error('Error: Wrong reference point option.');
    end
    

    % Find centroids
    while length(C) < k
        if length(C) < 1
            if DETERMINISTIC ~= 1
                [~,sorted] = sort(dists(r,:),'descend');
            else
                %Default method in Al Hasan et al. (2009)
                [~,sorted] = sort(pdist2(r,x),'descend');
            end
        else
            [~,sorted] = sort(min(dists(C,:),[],1),'descend');
        end
        sorted_lof = lof(sorted);
        %Pick points with LOF around 1 (1 +/- critRobin)
		id = find( (1-critRobin < sorted_lof) & (sorted_lof < 1+critRobin) );
		if isempty(id)
			warning('ROBIN: no valid id point, try 1.');
            %Use Brodinova's implementation [2] and pick any point with
            %LOF < critRobin_1
			id = find((sorted_lof < 1+critRobin) == 1);    
			if isempty(id)
				error('ROBIN: cannot find valid id point.')
			end
		end
        id = id(1);
        r = sorted(id);
        C = union(C,r);
    end

end