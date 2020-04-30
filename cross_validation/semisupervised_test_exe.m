function res_folds = semisupervised_test_exe(x,x_labs,constraints,ML,CL,method_clustering,method_centers,indices,nconstr,...
                                       citer,sstep,maxIter,DETERM,DISPLAY,fileID,CONSTR_PERC)
% Testing procedure for each fold:
% 1. Take a random (equal) sample of the MUST-LINK and CANNOT-LINK
%    constraints belonging to the training set until we have a total of 
%    'nconstr' constraints.
% 2. Shuffle the taken constraints.
% 3. Execute the algorithmic solution given 'maxIter' to converge.
% 4. Repeat for 'iter' times.
% An iteration x folds cell array will hold the results. Similar procedure
% was followed by [1] but with 50 repetitions of 5 folds.

% Constraints are belonging to the training set if their relevant labels
% are belonging to the training set. In the case of Pairwise Constrained
% Sparse K-Means algorithm different values of the sparsity parameter are
% tested starting for 1.1 with step 'sstep' up to sqrt(dims) as in [2].

%References:
%
%[1] Bilenko, M., Basu, S., & Mooney, R. J. (2004, July). Integrating 
%    constraints and metric learning in semi-supervised clustering. In 
%    Proceedings of the twenty-first international conference on Machine 
%    learning (p. 11). ACM.
%
%[2] Brodinová, Š., Filzmoser, P., Ortner, T., Breiteneder, C., & Rohm, M. 
%    (2017). Robust and sparse k-means clustering for high-dimensional 
%    data. Advances in Data Analysis and Classification, 1-28.


% Input:
% - x : a matrix where rows are observations and columns are attributes.
% - x_labs      : data labels
% - constraints : data contraints (*)
% - ML : (flag) use/don't use MUST-LINK constraints
% - CL : (flag) use/don't use CANNOT-LINK constraints
% - method_clustering : clustering algorithm
% - method_centers    : clustering initialization
% - indices           : indicates training and test sets
% - nconstr           : total number of constraints to be used
% - citer             : iterations per number of constraints
% - sstep             : sparsity parameter values to be tested
%                       [1.1 : sstep : sqrt(numFeatures)];
% - maxIter           : iterations for algorithm to reach convergence
% - DETERM            : (flag) start with fixed random seed

% (*) List of constraints is in Nx2 matrix format:
%     MUST-LINK = low -> high index
%     CANNOT-LINK = high -> low index
%     Original MPCK-Means:
%     MUST-LINK = high -> low index
%     CANNOT-LINK = low -> high index 

    if DETERM
        rng(47);
    end

    folds = unique(indices);
    k = length(unique(x_labs));
    res_folds = cell(citer,length(folds)); %iterations x folds
    
    if isequal(method_clustering,'Lloyd') || isequal(method_clustering,'SK-Means')
        if ~isequal(method_centers,'MPCK-Means')
            citer = 1;
        end
    end
    
    for i = 1:length(folds)
        % Split into training and test sets
        tmp_constraints = constraints;
        % Remove all the constraints generated using the test labels
        test = find(indices == i); 
        a = find(ismember(constraints(:,1),test));
        b = find(ismember(constraints(:,2),test));
        c = unique([a;b]);
        tmp_constraints(c,:) = [];     
        x_tmp = x;
        % Split to ML and CL
        a = find(tmp_constraints(:,1) > tmp_constraints(:,2));
        tmp_constraints_ML = tmp_constraints(a,:);
        a = find(tmp_constraints(:,1) < tmp_constraints(:,2));
        tmp_constraints_CL = tmp_constraints(a,:);
        if CONSTR_PERC && i==1
            nconstr = round(size(tmp_constraints,1)*nconstr/100);
            if rem(nconstr,2) == 1
                nconstr = nconstr + 1;
            end
        end
        % Take equal number of ML and CL
        for t = 1:citer
            if ML == -1 && CL == -1
                % Take random constraints without considering ML and CL
                cML = randsample(1:length([tmp_constraints_ML;tmp_constraints_CL]),nconstr/2);
                cCL = randsample(1:length([tmp_constraints_ML;tmp_constraints_CL]),nconstr/2);
            else            
                if ML && CL
                    % Take half and half
                    cML = randsample(1:length(tmp_constraints_ML),nconstr/2);
                    cCL = randsample(1:length(tmp_constraints_CL),nconstr/2);
                elseif ~ML && ~CL
                    % Take none
                    cML = [];
                    cCL = [];
                else
                    % Take full
                    if ML
                        cML = randsample(1:length(tmp_constraints_ML),nconstr);
                    else
                        cML = [];
                    end
                    if CL
                        cCL = randsample(1:length(tmp_constraints_CL),nconstr);
                    else
                        cCL = [];
                    end
                end
            end
            % Constraints to be used
            if ML == -1 && CL == -1
                constr = [tmp_constraints_ML;tmp_constraints_CL];
                constr = constr([cML,cCL],:);
            else
                constr = [tmp_constraints_ML(cML,:);tmp_constraints_CL(cCL,:)];
            end
            % Shuffle them
            constr = constr(randperm(length(constr)),:);   
            % Execute clustering solution
            if isequal(method_clustering,'PCSK-Means') || isequal(method_clustering,'SK-Means')
                S = 1.1:sstep:sqrt(size(x,2));
                if length(S) > 30
                    sstep = 1;
                    while length(S) > 30
                        S = 1.1:sstep:sqrt(size(x,2));
                        sstep = sstep + 1;
                    end   
                    sstep = S(2)-S(1);
                end
            else
                S = 0;
            end
            res_stats = struct('idx',[],'centroids',[],'w',[],'initCenters',[]);
            res_stats = repmat(res_stats,1,length(S));
            for ss = 1:length(S) 
                s = S(ss);
                %Clustering
                if isequal(method_clustering,'MPCK-Means')
                    [idx,centroids,w,iflag,init_centers] = semisupervised_test_exe_solution(method_centers,method_clustering,...
                        constr,x_tmp,k,s,maxIter,'RunFromJava');  
                else
                    [idx,centroids,w,iflag,init_centers] = semisupervised_test_exe_solution(method_centers,method_clustering,...
                        constr,x_tmp,k,s,maxIter);  
                end
                %Stats
                if iflag <= 0
                    [f_score,accuracy,recall,specificity,precision] = cl_FmeasureCL(x_labs(test),idx(test));
                else
                    f_score = NaN;
                    accuracy = NaN;
                    recall = NaN;
                    specificity = NaN;
                    precision = NaN;
                    if DISPLAY
                        str = sprintf('f%di%d: %s,%s,c=%d,k=%d,s=%0.2f,f=%d\n',...
                            i,t,method_centers,method_clustering,length(constr),k,s,iflag);
                        fprintf(str)
                    end
                    if fileID
                        str = sprintf('f%di%d: %s,%s,c=%d,k=%d,s=%0.2f,f=%d\n',...
                            i,t,method_centers,method_clustering,length(constr),k,s,iflag);                   
                        fprintf(fileID,'%s\n',str); 
                    end
                end
                %Collect
                res_stats(1,ss).idx = idx;
                res_stats(1,ss).centroids = centroids;
                res_stats(1,ss).w = w;
                res_stats(1,ss).flag = iflag;
                res_stats(1,ss).initCenters = init_centers;
                
                res_stats(1,ss).fscore = f_score;
                res_stats(1,ss).accuracy = accuracy;
                res_stats(1,ss).recall = recall;
                res_stats(1,ss).specificity = specificity;  
                res_stats(1,ss).precision = precision;
            end
            res_folds{t,i} = res_stats;
        end       
    end
end
