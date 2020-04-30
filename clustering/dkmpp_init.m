function [C,optimal,rho,eps,eps_neighbors_] = dkmpp_init(x,k,L)
%DKMPP_INIT implements the N. Nidheesh et al. cluster initialization method
%Original code in R can be found here: https://github.com/nidheesh-n/dkmpp

% Matlab implementation:
% Avgoustinos Vouros <avouros1@sheffield.ac.uk>, <av.vouros@gmail.com>

% Reference:
%   N. Nidheesh, K.A. Abdul Nazeer, P.M. Ameer, 
%   An enhanced deterministic K-Means clustering 
%   algorithm for cancer subtype prediction from 
%   gene expression data, Computers in Biology
%   and Medicine 91C (2017) pp. 213-221.
%   https://doi.org/10.1016/j.compbiomed.2017.10.014

% Input:
% - x : a matrix where rows are observations and columns are attributes.
% - k : number of target clusters.
% - L : (optional) Minimum Spanning Tree computing a priori.

% Output:
% - C : vector of row indeces of x (datapoints) to be used as initial
%       centroids. 

%%
    optimal = 1;
    n = size(x,1); %number of observations
    
    % Compute the distance matrix of x
    %(distance between every pair of points)
    M = squareform(pdist(x));

    % Normalization of distances (min-max according to the paper)
    Mnorm = M./max(max(M)); 
    %Mnorm = normalizations(M,'scale');

    % Using M as an adjacency matrix of a completely connected weighted 
    %undirected graph, construct a Minimum Spanning Tree.
    if ~exist('L','var')
        L = genMinimumSpanningTree(Mnorm);
    end

    % Assign best value for the radius of the hypersphere (epsilon)
    eps_tmp = 3 * iqr(L) + prctile(L,75); %extreme outliers formula
    eps = min(max(L),eps_tmp);

    % Compute the local density of each datapoint
    rho = zeros(1,size(x,1));
    eps_neighbors_ = cell(1,size(x,1));
    for i = 1:n
        %epsilon neghbors of data point i
        eps_neighbors = setdiff(find(Mnorm(:,i) <= eps),i); 
        eps_neighbors_{i} = eps_neighbors;
        %local density of datapoint i
        for j = 1:length(eps_neighbors)
            rho(i) = rho(i) + exp(-1 * Mnorm(i,eps_neighbors(j)) / eps);
        end
    end

    % min-max normalization of the rho vector
    rho_n = normalizations(rho,'scale');

    % Assign datapoints as the initial centroids
    [~,max_rho] = max(rho_n); %index of the highest density point
    prospect = zeros(1,n); %prospectiveness of each datapoint to be selected
    mind2cent = inf*ones(1,n);
    C = max_rho;
    while length(C) < k
        for i = 1:n
            mind2cent(i) = min(Mnorm(i,max_rho), mind2cent(i));
            prospect(i) = rho_n(i) * mind2cent(i);        
        end
        %[~,max_rho] = max(prospect); %point with maximum prospect = a centroid
        [~,sorti] = sort(prospect,'descend');
        for z = 1:length(sorti)
            if ismember(sorti(z),C)
                % Extreme case where C already contains max_rho(i): throw
                % warning and try next datapoint
                %warning('dkmpp: Cannot find more optimal centroids.');
                optimal = 0;
            else
                max_rho = sorti(z); %point with maximum prospect = a centroid
                C = union(C,sorti(z),'stable'); %add centroid to C
                break;
            end
        end
    end
    
    
    
    function LL = genMinimumSpanningTree(Mnorm)
        G = graph(Mnorm); %generate graph from M
        mst = minspantree(G); %generate Minimum Spanning Tree
        tmp = table2array(mst.Edges); %table to matrix
        LL = tmp(:,3); %weights of edges of the Minimum Spanning Tree        
    end
end

