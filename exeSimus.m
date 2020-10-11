% Test the semi-supervised methods

% This scripts generates a file each different
% (a) dataset, (b) clustering algorithm, (c) clustering initialization, 
% (d) constraints types.

% Each file contains a 1xC cell vector where NC is the different number of
% constrints used.

% Each cell containts a IxF cell array where I is the number of iterations
% (how many times this setting was executed changing the picked constraints
% each time) and F is the number of folds used for cross validation.

% Each sub-cell contains a 1x1 of 1xP structure that hold the clustering
% and evaluation for this fold. In case the PCSK-Means is used then each P
% specifies a different value for the sparcity parameter 's' from 1 to 
% sqrt(p) with step 'sstep' (p = the dimensionality of the dataset).

% To compute the k-fold cross validation the mean of an evaluation metric 
% over each row of the IxF cell array should be taken.


DETERM = 0;  %start with a fixed random seed
JMPCKM_OVERLOAD = 0; %0: non overloaded java MPCK-Means, 1: overloaded
CONSTR_PERC = 1; %0: flat number of constraints, 1: (%) of constraints

LOG = 2; % 0: no log file and no display, 1: log file only, 
         % 2: display only, else: both

datasets = {'40_3_3_7_0_0_0_0','40_3_5_5_0_0_0_0','real','MWM',...
    'digits048','digits389'};

% Algorithms options
constraints_type = {[1,1],[0,1],[1,0],[-1,-1]}; %Activate [ML,CL] constraints
if ~CONSTR_PERC
    constraints_number = [10:10:90,100:100:1000]; %Number of constraints to use
else
    constraints_number = 1:0.5:10; %Percentage of constraints to use
end
citer = 10; %iterations per number of constraints
sstep = 0.2; %sparsity parameter values to be tested [1.1 : sstep : sqrt(numFeatures)]
maxIter = 25; %iterations for algorithm to reach convergence
method_centers = {'Density K-Means++','MPCK-Means','ROBIN(D)'};
method_clustering = {'Lloyd','SK-Means','MPCK-Means','PCK-Means','PCSK-Means'};
method_centers_str = {'DKMPP','MPCK','ROBIN'};
method_clustering_str = {'LKM','SKM','MPCKM','PCKM','PCSKM'};

% Cross validation options
kfolds = 10; %Number of folds

% Output
startFolder = fullfile(pwd,'GenRes');
resFolder = 'results';

% Add to path
addpath(genpath(fullfile(pwd,'clustering')))            %algorithms+inits
addpath(genpath(fullfile(pwd,'cross_validation')))      %cv functions
addpath(genpath(fullfile(pwd,'utilities')))             %misc
%Pick which java version of MPCK-Means to run
if JMPCKM_OVERLOAD
    addpath(fullfile(pwd,'tests_hackedWeka'));
    addpath(genpath(fullfile(pwd,'tests_hackedWeka','weka_overload'))) 
    weka_init;
else
    addpath(fullfile(pwd,'tests_hackedWeka'));
    addpath(fullfile(pwd,'tests_hackedWeka','weka'));
    weka_init;
end



%% EXE
for I = 1:length(datasets)
    fdname = datasets{I};
    fullPath = (fullfile(startFolder,fdname,resFolder));
    if ~exist(fullPath,'dir')
        mkdir(fullPath);
    end

    if LOG == 0
        DISPLAY = 0;
        fileID = 0;
    elseif LOG == 1
        DISPLAY = 0;
        fileID = fopen(fullfile(fullPath,'log.txt'),'w');
    elseif LOG == 2
        DISPLAY = 1;
        fileID = 0;
    else
        DISPLAY = 1;
        fileID = fopen(fullfile(fullPath,'log.txt'),'w');   
    end

    psets = dir(fullfile(pwd,'datasets',fdname,'*.mat'));

    for i = 1:length(psets)
        % Load the dataset + constraints
        load(fullfile(psets(i).folder,psets(i).name));
        un = unique(x_labs);
        for L = 1:length(un)
            x_labs(x_labs==un(L)) = L;
        end
        if isequal(fdname,'real')
            tmp = [psets(i).name(1:end-4),'_constr.mat'];
            load(fullfile(pwd,'datasets_extras','constraints',tmp));
        end
        file_name = strsplit(psets(i).name,'.mat');
        file_name = file_name{1};
        
        % Generate cross validation folds
        indices = crossvalind('Kfold',x_labs,kfolds);

        for j = 1:length(method_clustering)
            %For every algorithm
            for jj = 1:length(method_centers)
                %For every init method
                for z = 1:length(constraints_type)
                    %For every constraints type
                    ML = constraints_type{z}(1);
                    CL = constraints_type{z}(2);
                    res_folds_nc = cell(1,length(constraints_number));
                    str = sprintf('%s-%s_%s_%s_ML%d_CL%d.mat',...
                        fdname,file_name,method_centers_str{jj},method_clustering_str{j},...
                        ML,CL);
                    if DISPLAY
                        str_ = strcat('Creating: ',str);
                        disp(str_);
                    end
                    if fileID
                        str_ = strcat('Creating: ',str);
                        fprintf(fileID,'%s\n',str_);                    
                    end
                    if ~exist(fullfile(fullPath,str),'file')
                        for zz = 1:length(constraints_number)
                            %For every constraints number
                            %Run the cross validation
                            res_folds = semisupervised_test_exe(x,x_labs,constr,ML,CL,...
                                method_clustering{j},method_centers{jj},indices,constraints_number(zz),...
                                citer,sstep,maxIter, DETERM,DISPLAY,fileID,CONSTR_PERC);

                            res_folds_nc{zz} = res_folds;

                            if isequal(method_clustering{j},'Lloyd') || isequal(method_clustering{j},'SK-Means')
                                if ~isequal(method_centers{jj},'MPCK-Means')
                                    break;
                                end
                            end
                        end
                        save(fullfile(fullPath,str),'res_folds_nc','-v7.3');
                    end
                end
            end
        end
    end
    fclose('all');
end