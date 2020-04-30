function [idx,centers,A,constraints,iter,init_centers,iflag,neighbours] = mpckmeans(x,k,constr,varargin)
%% IMPORTNAT NOTE:
%WORK IN PROGRESS. DO NOT USE FOR MPCK-MEANS AND METRIC K-MEANS. 
%USE IT ONLY FOR: Seeding, Lloyd's, PCK-Means

%REFERENCES:
%[1]
%Bilenko, Mikhail, Sugato Basu, and Raymond J. Mooney. 
%"Integrating constraints and metric learning in semi-supervised 
%clustering." Proceedings of the twenty-first international conference 
%on Machine learning. 2004.
%[2]
%Java implementation
%WekaUT: http://www.cs.utexas.edu/users/ml/risc/code/
%[3]
%R implementation:
%Tran Khanh Hiep and Nguyen Minh Duc (2016). 
%conclust: Pairwise Constraints Clustering. R package version 1.1
%https://rdrr.io/cran/conclust/



% - iflag     : 0 = ok, 1 = empty cluster(s), 2 = negative obj func,
%                       3 = negative/nan weights

    iterations = 25;          %number of iterations until exit
    gap_iterations = 10;      %number of iterations with no element changes
    ObjFunConvergenceDifference = -inf; %diff between old and new obj func; ignore if -inf
    metric_learning = 1;      %activate metric learning, if 0 = PCK-Means
    random_order = 0;         %process the elements in randon order
    %Infer extra constraints: 
    % - If 0 then no extra constraints are generated. If initial centers
    %   have not been provided then MUST-LINK constraints are overloaded
    %   only for the initialization process (MPCK-Means default if no 
    %   overloading is de-activated).
    % - If 1 then extra constraints are generated using the selected
    %   neighborhoods only if random centroids have not been generated.
    %   (MPCK-Means default).
    % - If 2 then extra constraints are generated using the selected
    %   neighborhoods.
    % - If > 2 then extra constraints are generated using all the
    %   neighborhoods (computational overkill).
    transitive_closure_ml = 0;%overload MUST-LINK constraints
    transitive_closure_cl = 0;%overload CANT-LINK constraints
    
    if isempty(constr)
        % K-Means
        iterations = 100;
        metric_learning = 0;
    end
    
    %The following options are active only with metric learning
    fix_normalizer = 1; %obj func is going negative use the Java method
    aML = 1; %weight of MUST-LINK constraints
    aCL = 1; %weight of CANNOT-LINK constraints
    
    DISPLAY = 0; %for debugging
    iflag = 0;
    init_only = 0; %return only the initial centroids
    
    [n,p] = size(x);    
    idx = [];             %cluster indexes
    centers = [];         %centroids
    constraints = [];     %constraints if extra are added from the transitive closure
    iter = 0;             %number of iterations
    A = diag(ones(1,p));  %distance matrix 
    prev_ObjFun = inf;    %previous objective function
    
    for i = 1:length(varargin)
        if isequal(varargin{i},'centers')
            centers = varargin{i+1};
            if size(centers,1) ~= k
                error('Number of given centers must be equal to k');
            end
            if size(centers,2) ~= size(x,2)
                error('Dimensionality of given centers must be equal to the dimensionality of the dataset x');
            end
        elseif isequal(varargin{i},'iterations')
            iterations = varargin{i+1};
        elseif isequal(varargin{i},'transitiveML')
            transitive_closure_ml = varargin{i+1};
        elseif isequal(varargin{i},'transitiveCL')
            transitive_closure_cl = varargin{i+1};  
        elseif isequal(varargin{i},'random_order')
            random_order = varargin{i+1};  
        elseif isequal(varargin{i},'metric_learning')
            metric_learning = varargin{i+1};        
        elseif isequal(varargin{i},'gap_iterations')
            gap_iterations = varargin{i+1};               
        elseif isequal(varargin{i},'objDiff')
            ObjFunConvergenceDifference = varargin{i+1};    
        elseif isequal(varargin{i},'init_only')
            init_only = varargin{i+1};
        elseif isequal(varargin{i},'MLweight')
            aML = varargin{i+1};
        elseif isequal(varargin{i},'CLweight')
            aCL = varargin{i+1};
        end
    end

    % Make the constraints
    %CANNOT-LINK a<b
    %MUST-LINK a>b
    cantLink = [];
    mustLink = [];
    for i = 1:size(constr,1)
        if constr(i,1) < constr(i,2)
            cantLink = [cantLink;constr(i,:)];
        elseif constr(i,1) > constr(i,2)
            mustLink = [mustLink;constr(i,:)];
        else
            error('Wrong constraints');
        end
    end

    % Adjacency matrix for MUST-LINK constraints
    nml = size(mustLink,1);
    M = zeros(n,n);
    for i = 1:nml
        u = mustLink(i,1);
        v = mustLink(i,2);
        M(u,v) = 1;
        M(v,u) = 1;
    end       
    
    % Adjacency matrix for CANNOT-LINK constraints
    ncl = size(cantLink,1);
    C = zeros(n,n);
    for i = 1:ncl
        u = cantLink(i,1);
        v = cantLink(i,2);
        C(u,v) = 1;
        C(v,u) = 1;
    end       
    
    % Transitive closure of the MUST-LINK constraints using DFS
    if isempty(centers)
        % Transitive closure of the MUST-LINK constraints using DFS
        [neighbours,idx] = depth_first_search(n,M);
        % Shrink neighbours variable
        sn = cellfun(@(x) length(x), neighbours);
        neighbours(sn==0) = [];
        sn(sn==0) = [];
        nn = length(neighbours);        
    end
    
    
    %% MPCK-Means centroids initialization
    %Note: exact Java results cannot be reproduced. This is (possibly) due 
    %to the numeric accuracy between the two languages. It was observed
    %that during the computation of weighted farthest first traversal there
    %were significant differences between Java's and MATLAB's MPCK-Means
    %implementations.
    if isempty(centers)
        % Neighbourhood centroids
        ncenters = nan(nn,p);
        for i = 1:nn
            ncenters(i,:) = mean(x(neighbours{i},:),1);
        end
        % Pick centroids
        if nn >= k
            %Pick the first centroid based on the largest neighbourhood
            [~,seln] = max(sn);  
            %Pick the rest based on weighted farthest first traversal
            for i = 2:k
                bestPoints = [];   
                maxDistanceSoFar = -inf;
                for j = 1:nn
                    if ~isempty(find(seln==j))
                        continue;
                    end 
                    minDistanceFromSet = +inf;
                    for jj = 1:length(seln)
                        % Weighted distance based on the cluster sizes
                        d = sqrt(sum((ncenters(j,:) - ncenters(seln(jj),:)) .^ 2)) * sqrt(sn(j) * sn(seln(jj)));
                        if d < minDistanceFromSet
                            minDistanceFromSet = d;
                        end                    
                    end
                    if minDistanceFromSet == maxDistanceSoFar
                        minDistanceFromSet = maxDistanceSoFar;
                        bestPoints = [bestPoints,j];
                    elseif minDistanceFromSet > maxDistanceSoFar
                        maxDistanceSoFar = minDistanceFromSet;
                        bestPoints = j;            
                    end                 
                end    
                %Pick the best point
                %If multiple pick one at random
                if length(bestPoints) == 1
                    seln = [seln,bestPoints];   
                else
                    seln = [seln,randsample(bestPoints,1)];     
                end
            end   
            %Pick the selected neighbourhoods
            centers = ncenters(seln,:);
        else
            %Pick centroids of neighbourhoods and rest at random from the data
            seln = 1:size(ncenters,1);
            centers = [ncenters ; x(randsample(1:n,k-nn),:)];
        end
    else
        seln = [];
    end
    init_centers = centers;
    if init_only
        return
    end
    
    % Assign points to selected neighbourhoods
    idx = zeros(n,1);
    dists = zeros(1,k);
    for i = 1:length(seln)
        elems = neighbours{seln(i)};
        for ii = 1:length(elems)
            for j = 1:k
                dists(j) = distFun(x(elems(ii),:), centers(j,:), A);
            end
            [~,idx(elems(ii))] = min(dists);   
        end
    end
    
    
    %% Collect the constraints and put datapoints in their clusters
    
    % Add inferred MUST-LINK constraints
    %MUST-LINK constraints are added among all the pairs of points in the 
    %same neighbourhood.
    if transitive_closure_ml
        switch transitive_closure_ml
            case 1
                %Infer constraints using the selected neighbourhoods only.
                %This happens if random centroids haven't been
                %genertaed (MPCK-Means default).  
                if nn >= k
                    nn_tmp = length(seln);
                else
                    nn_tmp = 0;
                end
            case 2
                %Infer constraints using the selected neighbourhood only.
                %Regardless of random centroids generation this process
                %will be executed.
                nn_tmp = length(seln);
            otherwise
                %Infer constraints using all the neighbourhoods
                nn_tmp = nn;
        end
        for i = 1:nn_tmp
            for j = 1:length(neighbours{seln(i)})
                for jj = 1:length(neighbours{seln(i)})
                    if neighbours{seln(i)}(j) < neighbours{seln(i)}(jj)
                        M(j,jj) = 1;
                        M(jj,j) = 1;
                    end
                end
            end
        end
        mustLink = [];
        for i = 1:n
            for j = i+1:n
                if M(i,j)==1
                    mustLink = [mustLink;[i,j]];
                end
            end
        end
        nml = size(mustLink,1);         
    end

    % Add inferred CANNOT-LINK constraints
    %CANT-LINK constraints are added among all the pairs of points in 
    %different neighbourhoods which have at least one CANNOT-LINK 
    %constraint
    if transitive_closure_cl
        switch transitive_closure_cl
            case 1
                %Infer constraints using the selected neighbourhoods only.
                %This happens if random centroids haven't been
                %genertaed (MPCK-Means default).  
                if nn >= k
                    nn_tmp = length(seln);
                else
                    nn_tmp = 0;
                end
            case 2
                %Infer constraints using the selected neighbourhood only
                %Regardless of random centroids generation this process
                %will be executed.
                nn_tmp = length(seln);
            otherwise
                %Infer constraints using all the neighbourhoods
                nn_tmp = nn;
        end  
        for i = 1:nn_tmp
            for j = i+1:nn_tmp
                if ~isempty(find(C(neighbours{seln(i)},:) == 1))
                    %There is at least one CANNOT-LINK constraint
                    %Generate CANNOT-LINK constraint between all the pairs
                    for ii = 1:length(neighbours{seln(i)})
                        for jj = 1:length(neighbours{seln(j)})
                            C(neighbours{seln(i)}(ii),neighbours{seln(j)}(jj)) = 1;
                            C(neighbours{seln(j)}(jj),neighbours{seln(i)}(ii)) = 1;
                        end
                    end
                end
            end
        end
        cantLink = [];
        for i = 1:n
            for j = i+1:n
                if C(i,j)==1
                    cantLink = [cantLink;[i,j]];
                end
            end
        end
        ncl = size(cantLink,1);        
    end  
    % Return also the new constraints
    constraints = {mustLink,cantLink};
    
    % Check if we have empty cluster(s)
    if length(find(unique(idx)~=0)) ~= k && ~isempty(seln)
        iflag = 1;
        return
    end    
    
    idx_ = nan(length(idx),1); %Tracking for elements that change cluster
    ccs = [mustLink;cantLink]; %Put all the constraints together
    
    
    %% Main loop
    for iter = 1:iterations
        
        % Maximally separated pair of points (enclosing hypercube)
        mins = min(min(x,[],1),zeros(1,size(x,2)));
        maxs = max(x,[],1);
        farthest = distFun(mins,maxs,A);
        farthest_vals = (maxs-mins).^2;
        
        % Iterate through elements and place them to clusters based on
        % penalty minimization
        movedp = 0;
        ObjFun = 0;
        if random_order
            %Process the elements in random order
            ord = datasample(1:n,n,'Replace',false);
        else
            ord = 1:n;
        end   
        for ii = 1:n
            % For each element...
            i = ord(ii);
            D = 1e15*ones(1,k);
            for j = 1:k
                % For each cluster...
                %Penalty without constraints
                if metric_learning
                    if fix_normalizer == 0
                        %R code
                        D(j) = distFun(x(i,:), centers(j,:), A) - ( log(prod(diag(A))) / log(2) );
                    else
                        %This is more towards the original java code
                        D(j) = distFun(x(i,:), centers(j,:), A) - 0.01*log(prod(diag(A)));
                    end
                else
                    %No metric learning
                    D(j) = distFun(x(i,:), centers(j,:), A);
                end
                %Penalty with constraints
                if ~isempty(ccs)
                    con = sort(find(ccs(:,1)==i | ccs(:,2)==i),'ascend');
                    for u = 1:length(con)
                        c1 = ccs(con(u),1); 
                        c2 = ccs(con(u),2);
                        if c1 == i
                            other = idx(c2);
                        else
                            other = idx(c1);
                        end
                        if other > 0 && other <= k
                            if other ~= j && c1 > c2
                                %MUST-LINK violation
                                D(j) = D(j) + aML * distFun(x(c1,:), x(c2,:), A);
                            elseif other == j && c1 < c2
                                %CANT-LINK violation
                                tmp = farthest - aCL * distFun(x(c1,:), x(c2,:), A);
                                if tmp < 0
                                    iflag = 2;
                                    return
                                end
                                D(j) = D(j) + tmp;                            
                            end
                        end
                    end
                end
                if DISPLAY
                    fprintf('Final penalty of element %d to cluster %d is: %8.9f\n',ii,j,D(j));
                end                
                
                %Find the cluster with the lowest penalty for this element
                [~,tmp] = min(D(1:j));
                if tmp ~= idx(i)
                    idx(i) = j;
                end
            end
            %Assign element to cluster
            if idx_(i) ~= idx(i)
                movedp = movedp+1;
            end
            if DISPLAY
                if ~isnan(idx_(i)) && idx_(i) ~= idx(i)
                    fprintf('Element %d is moved from cluster %d to cluster %d\n\n',ii,idx_(ii),idx(ii));
                else
                    fprintf('Element %d is assigned to cluster %d\n\n',ii,idx(ii));
                end
            end
            %Compute objective function
            ObjFun = ObjFun + D(idx(i));
        end
        
        if DISPLAY
            if iter == 1
                fprintf('\nNumber of elements assigned to clusters: %d\n',movedp);
            else
                fprintf('\nNumber of elements moved in this iteration (iter=%d): %d\n',iter,movedp);
            end
            fprintf('\nObjective function (old) = %3.5f\n',prev_ObjFun);
            fprintf('Objective function (new): %3.5f\n',ObjFun);
            fprintf('Difference: %3.5f\n',abs(prev_ObjFun - ObjFun));
        end
        
        
        %% Update centroids locations
        centers = zeros(k,p);
        un = unique(idx);
        if length(un) ~=k 
            iflag = 1;
            return            
        end
        for i = 1:k
            a = find(idx == un(i));
            if length(a) >= 1
                centers(i,:) = mean(x(a,:),1);
            else
                %error('Error in MPCKMeans M-step: empty cluster detected!');
                iflag = 1;
                return
            end
        end

        % Check if max iterations have been reached
        if iter == iterations
            if DISPLAY
                fprintf('\nMax number of iterations reached\n');
            end             
            break;
        end
        % Check if gap iterations have been reached
        if isequal(idx_,idx)
            gap_iterations = gap_iterations - 1;
            if gap_iterations == 0
                if DISPLAY
                    fprintf('\nConverged: Centroids\n');
                end     
                break;
            end
        end
        idx_ = idx;
        % Check if there is minor change in the objective function
        if abs(prev_ObjFun - ObjFun) <= ObjFunConvergenceDifference && ~isinf(ObjFunConvergenceDifference)
            if DISPLAY
                fprintf('\nConverged: ObjFunc\n');
            end            
            break;
        else
            prev_ObjFun = ObjFun;
        end
        
        
        %% Update the metric
        if metric_learning
            cA = A;
            W = 0;
            for i = 1:n
                %Penalty without constraints
                idiff = (x(i,:) - centers(idx(i),:)).^2;
                %Penalty with constraints
                cdiff = 0;
                con = sort(find(ccs(:,1)==i | ccs(:,2)==i),'ascend');
                for u = 1:length(con)
                    c1 = ccs(con(u),1); 
                    c2 = ccs(con(u),2);
                    if c1 == i
                        other = idx(c2);
                    else
                        other = idx(c1);
                    end
                    if other > 0
                        if other ~= idx(i) && c1 > c2
                            %MUST-LINK violation
                            cdiff = cdiff + aML/2 .* ((x(c1,:) - x(c2,:)).^2);
                        elseif other == idx(i) && c1 < c2
                            %CANT-LINK violation
                            tmp = (aCL/2 .* farthest_vals) - (aCL/2 .* ((x(c1,:) - x(c2,:)).^2));
                            if tmp < 0
                                error('Negative penalty')
                            end
                            cdiff = cdiff + tmp;                            
                        end
                    end
                end  
                weights = idiff + cdiff;
                W = W+weights;
            end
            %Check the weights: all must be >= 0
            for i = 1:p
                if W(i) > 0
                    A(i,i) = 0.01 .* (n / W(i));
                end
            end
        end
        
    end
   
    
    
    
    %% Function to compute the weighted distances
    function wdist = distFun(x,y,A)
        differ = x - y;
        wdiffer = differ * A; 
        wdist = wdiffer * differ'; 
    end

    %% Function for transitive closure with depth-first search
    function [neighbours,idx] = depth_first_search(n,M)
        %n = number of elements
        %M = adjacency matrix
        neighbours = cell(n,1); %neighbourhoods
        idx = zeros(n,1);       %neighbourhood of each element
        visited = zeros(n,1);   %mark visited nodes
        id = 1;                 %current neighbourhood
        for ui = 1:n
            aa = find(M(ui,:)>0);
            if visited(ui)==0 && ~isempty(aa) > 0
                [idx,neighbours,visited] = depth_first_search_visit(ui,visited,M,idx,neighbours,id);
                id = id + 1;
            end
        end
    end
    function [idx,neighbours,visited] = depth_first_search_visit(ui,visited,M,idx,neighbours,id)
        visited(ui) = 1;
        aa = find(M(ui,:)==1);
        for uj = 1:length(aa)
            if visited(aa(uj))==0
                %the vertex is still undiscovered
                [idx,neighbours,visited] = depth_first_search_visit(aa(uj),visited,M,idx,neighbours,id);
            end
        end
        %update stats for ui
        idx(ui) = id;
        neighbours{id} = [neighbours{id},ui];
        visited(ui) = 2;
    end    
end

