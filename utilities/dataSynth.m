function [x,y,lb1,lb2] = dataSynth(varargin)
%DAYASYNTH. This function generates synthetic datasets for the purpose of 
%data clustering.

% Authors:
% Original function written in R by Sarka Brodinova (package: wrsk)
% <sarka.brodinova@tuwien.ac.at>
% MATLAB version written by Avgoustinos Vouros

% References:
% [1] S. Brodinova, P. Filzmoser, T. Ortner, C. Breiteneder, M. Zaharieva. 
%     Robust and sparse k-means clustering for high-dimensional data. 
%     Submitted for publication, 2017. 
%     Available at http://arxiv.org/abs/1709.10012
%
% [2] Campello, Ricardo JGB, et al. "Hierarchical density estimates for 
%     data clustering, visualization, and outlier detection." ACM 
%     Transactions on Knowledge Discovery from Data (TKDD) 10.1 (2015): 5.


%% Default values 
group_sizes = [40,40,40]; %datapoints and groups (3 groups of 40 points each)
p_info = 50;   %number of informative variables  
p_noise = 750; %number of uninformative variables   
groups_mu_range = [-6,-3;3,6]; %informative variables mu range; each row: [min,max]
groups_rho_range = [0.1,0.9];  %informative variables rho range; each row: [min,max]
p_noise_m = 0;  %uninformative variables mean
p_noise_s = 1;  %uninformative variables variance

MINDIST = 0; %minimum allowed distance between pairs in different clusters
REPEAT = 1;
END_REPEAT = 100000;

pn_info = 50;   %number of informative variables contaminated with outliers    
pp_info = 0.1;  %proportion of observations to be contaminated in the informative variables 
add_noise_info = 'first';
pn_noise = 75;  %number of uninformative variables contaminated with outliers    
pp_noise = 0.1; %proportion of observations to be contaminated in the uninformative variables.
add_noise_uninfo = 'last';

% Ouliers in informative variables:
scatter_out = 1; 
%1: scattered outliers are generated with the characteristics specified in out_scatter_range
out_scatter_range = [3,10];   
%0: uniformly distributed outliers are produced with the specification defined in out_unif_range
out_noise_range = [-12,6;6,12]; 

% Ouliers in uninformative variables:
out_unif_range = [-12,6;6,12]; %uninformative variables mu range; each row: [min,max]

%% Custom values
for i = 1:length(varargin)
    switch char(varargin{i})
        case 'dataset'
            group_sizes = varargin{i+1};
        case 'p_info'
            p_info = varargin{i+1};
        case 'p_noise'
            p_noise = varargin{i+1};
        case 'pn_info'
            pn_info = varargin{i+1};
        case 'pp_info'
            pp_info = varargin{i+1};
        case 'pn_noise'
            pn_noise = varargin{i+1};
        case 'pp_noise'
            pp_noise = varargin{i+1};
        case 'groups_mu_range'
            groups_mu_range = varargin{i+1};
        case 'groups_rho_range'
            groups_rho_range = varargin{i+1};
        case 'scatter_out'
            scatter_out = varargin{i+1};
        case 'out_scatter_range'
            out_scatter_range = varargin{i+1};
        case 'out_noise_range'
            out_noise_range = varargin{i+1};
        case 'out_unif_range'
            out_unif_range = varargin{i+1};
        case 'deterministic'
            rng(1,'twister');
        case 'MINDIST'
            MINDIST = varargin{i+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    ng = length(group_sizes); %number of clusters
    x = [];   %dataset
    y = [];   %group membership before contamination
    lb1 = []; %group and outlier (informative) membership
    lb2 = []; %group and outlier (uninformative) membership  
    
    %For the covarience generate ng numbers from continuous uniform 
    %random numbers within range specified in out_scatter_range
    ros = unifrnd(out_scatter_range(1,1),out_scatter_range(1,2),1,ng);
    %Number of elements in the informative variables to be contaminated
    n_out = round(group_sizes.*pp_info);   
    %Number of elements in the uninformative variables to be contaminated
    nu_out = round(group_sizes .* pp_noise);
    
    
    %% Generate the informative variables mu and sigma
    if p_info > 0
        %Generate ng continuous uniform random numbers from two 
        %distributions with lower and upper endpoints
        r1 = unifrnd(groups_mu_range(1,1),groups_mu_range(1,2),1,ng);
        r2 = unifrnd(groups_mu_range(2,1),groups_mu_range(2,2),1,ng);
        %Pick one number and clone it ng times
        groups_mu = repmat(randsample([r1,r2],1),1,ng);     
        %Construct the mu matix where each column is a group and each row an
        %informative feature. The means per informative variable of the groups 
        %will follow the format described in [1].
        info_mu_matrix = diag(groups_mu);
        while size(info_mu_matrix,1) < p_info
            info_mu_matrix = [info_mu_matrix;diag(groups_mu)];
        end
        info_mu_matrix = info_mu_matrix(1:p_info,:);

        %Generate ng ontinuous uniform random numbers from two distributions
        %with lower and upper endpoints
        groups_rho = unifrnd(groups_rho_range(1,1),groups_rho_range(1,2),1,ng);
        %Covaruence matrix formula ?t = Q * M * Q', where M is a matrix with 
        %ones in the diagonal and all the other elements are random numbers 
        %from groups_rho_range [2].
        %Construct ng such matrices, one for each group.
        diagonal = 1;
        groups_diag = diagonal * ones(1,ng);
        info_si_matrix = cell(1,ng); %store the covariance matrices of the groups
        for i = 1:ng
            gs = repmat(groups_rho(i),p_info,p_info);
            gs(1:p_info+1:end) = groups_diag(i); %replace the diagonal with groups_diag
            %Generate rotation matrix (based on R documentation, rRotationMatrix)
            if p_info == 1
                rm = 1;
            elseif p_info == 2
                p21 = rand;
                p11 = sqrt(1-p21^2);
                p22 = p11;
                p12 = -p21;
                rm = [p11,p12;p21,p22];
            else
                while 1
                    %We need the rank of A to be equal to the number of
                    %variables p_info
                    A = rand(p_info,p_info);
                    if rank(A) == p_info
                        break;
                    end
                end
                %QR decomposition of A
                [Q,~,~] = qr(A);
                rm = Q;
            end
            % Final covariance matrix Q x gs x Q', see [2]
            info_si_matrix{i} = rm * gs' * rm';
        end    
    end

    
    %% Generate the uninformative variables mu and sigma
    %nothing to do here...
 
    
    %% Generate the dataset 
    while REPEAT <= END_REPEAT
        % Terminate only if the distance between any pair of points
        % in two different clusters are more than MINDIST
        Xorig = [];
        Yorig = [];
        x = [];   %dataset
        y = [];   %group membership before contamination
        lb1 = []; %group and outlier (informative) membership
        lb2 = []; %group and outlier (uninformative) membership          
        for i = 1:ng
            lb = i*ones(1,group_sizes(i));
            lb_ = lb;
            y = [y,i*ones(1,group_sizes(i))];
            X = [];
            if p_info > 0
                %Generate the informative variables of the dataset using 
                %multivariate normal random numbers. The means per informative 
                %variable of the groups will follow the format described in [1].            
                X = mvnrnd(info_mu_matrix(:,i),info_si_matrix{i},group_sizes(i));
                Xorig = [Xorig;X];
                Yorig = y;
                %Check if we need to add outliers
                if pn_info > 0 && n_out(i) > 0
                    %If we have requested to contaminate informative variables
                    %and if the number of elements to be contaminated is enough...
                    if pn_info > p_info
                        pn_info = p_info;
                    end      
                    p_out = 1:pn_info;
                    %Create the outliers
                    if scatter_out
                        %Add them on the same location, same mean different
                        %variance. Contaminate first pn_info variables and 
                        %n_out(i) elements per group
                        covar = diag(ros(i)*ones(1,pn_info));
                        switch add_noise_info
                            case 'first'
                                %Add the outliers at the start
                                X(1:n_out(i),p_out) = mvnrnd(info_mu_matrix(p_out,i),covar,n_out(i)); 
                            case 'last'
                                %Add the outliers at the end
                                X(1:n_out(i)+1,p_out) = mvnrnd(info_mu_matrix(p_out,i),covar,n_out(i)); 
                        end                    
                    else
                        %Non-scatter outliers (outliers with differnt location than groups)
                        rtmp1 = unifrnd(out_noise_range(1,1),out_noise_range(1,2),1,n_out(i)*pn_info);
                        rtmp2 = unifrnd(out_noise_range(2,1),out_noise_range(2,2),1,n_out(i)*pn_info);
                        rtmp = randsample([rtmp1,rtmp2],n_out(i)*pn_info);  
                        switch add_noise_info
                            case 'first'
                                %Add the outliers at the start
                                X(1:n_out(i),p_out) = reshape(rtmp,n_out(i),pn_info);  
                            case 'last'
                                %Add the outliers at the end
                                X(end-n_out(i)+1:end,p_out) = reshape(rtmp,n_out(i),pn_info);
                        end
                    end
                    switch add_noise_info
                        case 'first'
                            %Add the outliers at the start
                            lb(1:n_out(i)) = 0;
                        case 'last'
                            %Add the outliers at the end
                            lb(end-n_out(i)+1:n_out(i)) = 0;
                    end                
                end
            end

            if p_noise > 0
                %Generare the uniformative variables of the dataset using 
                %normal random numbers
                rnorm = normrnd(p_noise_m,p_noise_s,group_sizes(i),p_noise); 
                %Check if we need to add outliers
                if p_noise > 0 && nu_out(i) > 0
                    %If we have requested to contaminate uninformative variables
                    %and if the number of elements to be contaminated is enough...
                    if pn_noise > p_noise
                        pn_noise = p_noise;
                    end   
                    p_out = 1:pn_noise;
                    %Contaminate first pn_noise variables and nu_out(i) 
                    %elements per group                
                    rtmp1_ = unifrnd(out_unif_range(1,1),out_unif_range(1,2),1,nu_out(i)*pn_noise);
                    rtmp2_ = unifrnd(out_unif_range(2,1),out_unif_range(2,2),1,nu_out(i)*pn_noise);
                    rtmp_ = randsample([rtmp1_,rtmp2_],nu_out(i)*pn_noise);  
                    switch add_noise_uninfo
                        case 'first'
                            %Add the outliers at the start
                            rnorm(1:nu_out(i),p_out) = reshape(rtmp_,nu_out(i),pn_noise); 
                            lb_(1:nu_out(i)) = 0;
                        case 'last'
                            %Add the outliers at the end
                            rnorm(end-nu_out(i)+1:end,p_out) = reshape(rtmp_,nu_out(i),pn_noise);
                            lb_(end-nu_out(i)+1:end) = 0;
                    end                
                end
                X = [X,rnorm];
            end
            x = [x;X]; 
            lb1 = [lb1,lb]; 
            lb2 = [lb2,lb_];
        end
    
        if MINDIST > 0 
            % Check if the distance between any pair of points
            % in two different clusters are more than MINDIST
            flag = 1;
            for k1 = 1:ng
                for k2 = k1+1:ng
                    dists = pdist2(Xorig(Yorig==k1,:),Xorig(Yorig==k2,:),'euclidean');
                    a = min(min(dists));
                    if a <= MINDIST
                        flag = 0;
                    end
                end
            end
            if flag
                y = y'; %Nx1 
                lb1 = lb1'; %Nx1
                lb2 = lb2'; %Nx1                
                return
            end
            REPEAT = REPEAT + 1;
            if REPEAT == END_REPEAT + 1
                warning('Unable to generate pair of points in different clusters with requested separation');
            end
        else
            break;
        end
    end
    
    y = y'; %Nx1 
    lb1 = lb1'; %Nx1
    lb2 = lb2'; %Nx1
end