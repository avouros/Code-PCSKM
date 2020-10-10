function [LOF,ard,density,neighbors] = lof_paper(data,nn)
%LOF_PAPER computes the Local Outlier Factor (LOF) based on the work of
%[1]. This code is based on the paper formulas.

% REFERENCES:
% [1] Al Hasan, Mohammad, et al. "Robust partitional clustering by 
%     outlier and density insensitive seeding." Pattern Recognition 
%     Letters 30.11 (2009): 994-1002.

%INPUT: 
% nn: number of neighbors
% data: dataset, rows = observations (N), columns = attributes (p)
% n2: (optional) computes the LOF only for the first n2 data points

%OUTPUT:
% LOF: LOF score of each element (Nx1 vector).
% ard: average reachability distance (Nx1 vector).
% density: density of each datapoint (Nx1 vector).
% neighbors: a cell array (Nx1) where each cell contains the indexes of
%            the nn closest datapoints of each datapoint. 

% Author: Avgoustinos Vouros

    n = size(data,1);
    
    % Compute the density of each datapoint
    density = zeros(n,1);
    neighbors = cell(n,1);
    parfor i = 1:n
        % Find the distance between the i-th datapoint and all the other
        % elements
        tmp = pdist2(data(i,:),data);
        % Sort the distances using ascending order
        dtmp = sort(tmp,'ascend');
        % Take the nn closer neighbors of the i-th datapoint
        r = dtmp(nn+1); %+1 because we have the distance between the i-th and i-th datapoint
        tmp = find(tmp <=r );
        tmp(tmp==i) = []; %remove the i-th datapoint
        neighbors{i} = tmp;
        % Compute the density of the i-th datapoint
        density(i) = length(neighbors{i}) / sum(pdist2(data(i,:),data(neighbors{i},:)));
    end

    % Compute the average relative density of each datapoint
    ard = zeros(n,1);
    parfor i = 1:n
        ard(i) = density(i) / (sum(density(neighbors{i})) / length(neighbors{i}));
    end    
    
    % Compute the Local Outlier Factor score
    LOF = nan(n,1);
    parfor i = 1:n
        LOF(i) = 1 / ard(i);
    end    
end

