function [idx,centroids,w,iflag,init_centers] = semisupervised_test_exe_solution(method_centers,method_clustering,...
    constr,x,k,s,maxIter,varargin)
%SEMISUPERVISED_TEST_EXE_SOLUTION 

% Input:
% - method_centers   : method for centroids initialization
% - method_clustering: clustering algorithm
% - constr           : list of constraints in Nx2 matrix format:
%                      MUST-LINK = low -> high index
%                      CANNOT-LINK = high -> low index
%                      (opposite for original MPCK-Means)

% - x : a matrix where rows are observations and columns are attributes.
% - k : number of target clusters.
% - s : sparsity parameter.

% - maxIter     : number of iterations until algorithm converges.    
     
    JFlag = 0;
    if ~isempty(varargin)
        if isequal(varargin{1},'RunFromJava')
            %RUn MPCK-Means algorithm using Java
            JFlag = 1;
        end
    end
    
    
    % Generate centroids
    if isnumeric(method_centers)
        [t1,t2] = size(method_centers);
        assert(t1==k,'Wrong number of centroids (manual)');
        assert(t2==size(x,2),'Wrong dimensions of centroids (manual)');
        init_centers = method_centers;
    else
        switch method_centers
            case 'Density K-Means++'
                init_centers = x(dkmpp_init(x,k),:);
            case 'MPCK-Means'
                [~,~,~,~,~,init_centers,~] = mpckmeans(x,k,constr,'init_only',1);
            otherwise
                error('Unknown init method');
        end
    end
    if JFlag
        if exist(fullfile(pwd,'mycentroids.txt'),'file')
            delete(fullfile(pwd,'mycentroids.txt'));
        end        
        if ~isequal(method_centers,'MPCK-Means')
            try
                writematrix(init_centers,'mycentroids.txt','Delimiter',',');
            catch
                %MATLAB 2017
                dlmwrite('mycentroids.txt',init_centers,'Delimiter',',');
            end
        end
        method_clustering = 'JMPCK-Means';
    end
    
    
    % Run clustering algorithm
    switch method_clustering
        case 'PCK-Means'
            [idx,centroids,w,~,~,~,iflag] = mpckmeans(x,k,constr,...
                'centers',init_centers,'iterations',maxIter,'metric_learning',0,...
                'transitiveML',0,'transitiveCL',0,...
                'gap_iterations',1);
            w = diag(w);
            
%         case 'MPCK-Means' %Java here
%             [idx,centroids,w,~,~,~,iflag] = mpckmeans(x,k,constr,...
%                 'centers',init_centers,'iterations',maxIter,'metric_learning',1,...
%                 'transitiveML',0,'transitiveCL',0,...
%                 'gap_iterations',1);
%             w = diag(w);  
            
        case 'PCSK-Means'
            [idx,centroids,w,~,~,iflag] = pcskmeans(x,k,s,constr,init_centers,...
                'iters',maxIter,'iterk',maxIter,'transitiveML',0,'transitiveCL',0);
            
        case 'Lloyd'
            [idx,centroids,w,~,~,~,iflag] = mpckmeans(x,k,[],...
                'centers',init_centers,'iterations',maxIter,'metric_learning',0,...
                'transitiveML',0,'transitiveCL',0,...
                'gap_iterations',1);            
			w = diag(w);
            
		case 'SK-Means'
            [idx,centroids,w,~,~,iflag] = pcskmeans(x,k,s,[],init_centers,...
                'iters',maxIter,'iterk',maxIter,'transitiveML',0,'transitiveCL',0); 
            
        % SPECIAl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        case 'JMPCK-Means'
            %Run the MPCKM using the hacked version of weka
            constr = [constr(:,2),constr(:,1)];
            save('lastrun.mat','x','constr','k');
            [idx,centroids] = Jmpckmeans(x, constr, k);
            idx = double(idx')+1;
            iflag = -inf;
            fileID1 = fopen('initCents.txt');
            tline = fgetl(fileID1);
            init_centers = nan(k,size(x,2));
            kk = 1;
            while ischar(tline)
                if ischar(tline)
                    tmp = strsplit(tline,',');
                end                    
                for z = 1:length(tmp)
                    init_centers(kk,z) = str2double(tmp{z});
                end
                kk = kk+1;
                tline = fgetl(fileID1);
            end                
            fclose(fileID1);
            %final weights
            fileID1 = fopen('WeightsFinal.txt');
            tline = fgetl(fileID1);
            w = nan(1,size(x,2));
            kk = 1;
            while ischar(tline)
                w(1,kk) = str2double(tline);
                tline = fgetl(fileID1);
                kk = kk+1;
            end                
            fclose(fileID1);
            %delete(fullfile(pwd,'*.txt'));        
            
%         % TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                     
%         case 'MPCK-Means-test' %just run weka
%             constr = [constr(:,2),constr(:,1)];
%             [idx,centroids] = mpckmeans_(x, constr, k);
%             idx = double(idx');
%             centroids = centroids';
%             w = ones(1,size(x,2));   
%             iflag = -inf;
%             
%         case 'Lloyd-test' %kmeans_lloyd <-> mpckmeans
%             [idx,centroids,~,~] = kmeans_lloyd(x,k,init_centers,100);
%             [tidx,tcentroids,w,~,~,~,~] = mpckmeans(x,k,[],...
%                 'centers',init_centers,'iterations',maxIter,'metric_learning',0,...
%                 'transitiveML',0,'transitiveCL',0,...
%                 'gap_iterations',1);
%             w = diag(w);
%             idx = idx-tidx;
%             if ~isnan(sum(idx,'all'))
%                 centroids = centroids-tcentroids;
%                 assert(sum(idx,'all')==0,'Error: Lloyd-test, idx not equal');
%                 assert(sum(centroids,'all')==0,'Error: Lloyd-test, centroids not equal');
%                 iflag = -inf;
%             end
%             
%         case 'SK-Means-test' %sparse_kmeans <-> pcskmeans
%             [idx,centroids,w,~,~,~] = sparse_kmeans(x,k,s,'Start',init_centers,...
%                 'iters',maxIter,'iterk',maxIter);
%             [tidx,tcentroids,tw,~,~,~] = pcskmeans(x,k,s,[],init_centers,...
%                 'iters',maxIter,'iterk',maxIter,'transitiveML',0,'transitiveCL',0);  
%             idx = idx-tidx;
%             centroids = centroids-tcentroids;
%             w = w-tw;
%             if ~isnan(sum(idx,'all'))
%                 assert(sum(idx,'all')==0,'Error: Lloyd-test, idx not equal');
%                 assert(sum(centroids,'all')==0,'Error: Lloyd-test, centroids not equal');       
%                 assert(sum(w,'all')==0,'Error: SK-Means-test, weights not equal'); 
%             end
%             iflag = -inf;
 		otherwise
             error('Unknown clustering method');
    end

end

