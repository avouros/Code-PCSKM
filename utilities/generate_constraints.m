function [constraints,constr] = generate_constraints(data,labels,varargin)
%GENERATE_CONSTRAINTS generates a structure containing constraints

%INPUT:
% - data:    n-by-p array where n is the number of observations and p the
%            number of atributes/features.
%
% - labels:  
%            (a) A two-column array specifying the constraint type between
%            elements.
%              If col1 < col2 => MUST-LINK constraint
%              If col1 > col2 => CANNOT-LINK constraint
%              If col1 = col2 => ignore
%              E.g. 
%                   1 | 2
%                   9 | 7
%              will create a MUST-LINK constraint between the datapoints 1 
%              and 2 and a CANNOT-LINK constraint between the datapoints 9 
%              and 7.
%
%            (b) A one-column array specifying the label of each element.
%            If an element does not have label use the special class -1
%              E.g. (1)
%                   data = [1,2; 2,2; 2,1; 4,3], labels = [1;1;2;2]
%              will create MUST-LINK constraints between the pairs of 
%              elements 1-2 and 3-4 and MUST-LINK constraints between
%              1-3, 1-4, 2-3, 2-4. 
%              E.g. (2)
%                   data = [1,2; 2,2; 2,1; 4,3], labels = [1;-1;2;2]
%              will create MUST-LINK constraints between the pair of 
%              elements 3-4 and MUST-LINK constraints between 1-3, 1-4.
%
% - varargin: name-value pair arguments. 
%            (a) 'thresholdML': MUST-LINK constraints between pairs of 
%                elements will be created only if the elements of the pairs
%                have less or equal of the specified distance. If a vector
%                of two numerics [n1,n2] is specfied then constraints will 
%                be created between pairs with distance less or equal to n2
%                and more or equal to n1. Default: +inf. Ignore: +inf
%            (b) 'thresholdCL': Same as thresholdML but in regards to the 
%                CANNOT-LINK constraints. Default: +inf. Ignore: +inf
%            (c) 'metricc': metric to be considered for calculating the
%                distance between a pair. 
%                Metrics: 'euclidean' (default), 'squaredeuclidean'

%OUTPUT:
%constraints: A structure constaining constraints:
%             x1, x2: Indeces of elements between which a constrain is
%                     available.
%             link:   Type of the constrain, 
%                     MUST-LINK = 29; CANNOT-LINK = 31;



    %% Default and user-specified options
    thresholdML = +inf;
    thresholdCL = +inf;
    metric = 'euclidean';
    
    if isempty(labels)
        constraints = [];
        return
    end

    for i = 1:length(varargin)
        if isequal(varargin,'metricc')
            metric = varargin{i+1};
            if ~isequal(metric,'squaredeuclidean') || ~isequal(metric,'euclidean')
                error('Euclidean or Euclidean^2 is allowed as metric')
            end
        elseif isequal(varargin{i},'thresholdML')
            thresholdML = varargin{i+1};
        elseif isequal(varargin{i},'thresholdCL')
            thresholdCL = varargin{i+1};
        end
    end
    
    %% Check user input
    [nn,mm] = size(labels);
    [n,~] = size(data);
    option  = 0;
    
    if mm == 2
        %Pair of indexes
        t = max(mm);
        if t(1) > n || t(2) > 2
        	error('Index of labelled datapoint exceeds the total number of points in the dataset.');
        end
        t = min(mm);
        if t(1) < n || t(2) < 2
            error('Index of labelled datapoint exceeds the total number of points in the dataset.');
        end   
        option = 1;
    elseif mm == 1
        %Labels
        if nn ~= n
            error('Labels should be a Nx1 array where N equals the number of the datapoints.');
        end
        option = 2;
    end
    
    nML = length(thresholdML);
    nCL = length(thresholdCL);
    if nML > 2 || nCL > 2
        error('Distance limits for the MUST-LINK and the CANNOT-LINK constraints can be either a scalar or a 1x2 vector.');
    end
    thresholdML = sort(thresholdML);
    thresholdCL = sort(thresholdCL);
    thresholdML(thresholdML<=0) = [];
    thresholdCL(thresholdCL<=0) = [];
    if isempty(thresholdML)
        thresholdML = +inf;
    end
    if isempty(thresholdCL)
        thresholdML = +inf;
    end    
    

    %% Execute option
    switch option
        case 0 %Wrong input
            error('Labels can be either a Mx2 array where the first and the second row spesify pair of datapoints or a Nx1 array where N equals the number of the datapoints.');
        case 1 %(a) two-column array
            constraints = struct('x1',[],'x2',[],'link',[],'distance',[]);
            constraints = repmat(constraints,nn,1);  
            for i = 1:nn
                if labels(i,1) < labels(i,2)
                    constraints(i).x1 = labels(i,1);
                    constraints(i).x2 = labels(i,2);
                    constraints(i).link = 29;                    
                elseif labels(i,1) > labels(i,2)
                    constraints(i).x1 = labels(i,2);
                    constraints(i).x2 = labels(i,1);
                    constraints(i).link = 31;                    
                else
                    %skip...
                end 
                %Compute the distance
                constraints(i).distance = pdist2(data(labels(i,1),:),data(labels(i,2),:),metric);
            end
            %Check the distance
            %ML
            if any(thresholdML==+inf)==0 %go inside only if there is no +inf
                if nML > 1
                    t = find([constraints.link]==29 & ...
                             [constraints.distance]<=thresholdML(2) & ...
                             [constraints.distance]>=thresholdML(1));
                else
                    t = find([constraints.link]==29 & ...
                             [constraints.distance]<=thresholdML);                
                end
                constraints = constraints(t);
            end
            %CL
            if any(thresholdML==+inf)==0 %go inside only if there is no +inf
                if nCL > 1
                    t = find([constraints.link]==31 & ...
                             [constraints.distance]<=thresholdCL(2) & ...
                             [constraints.distance]>=thresholdCL(1));
                else
                    t = find([constraints.link]==31 & ...
                             [constraints.distance]<=thresholdCL);                
                end
                constraints = constraints(t);
            end
            constr = labels;
        case 2 %(b) one-column array
            constraints = [];
            vec = (1:n);
            t = find(labels==-1);
            %Labels is now an Nx2 array, [index, element]
            labels = [vec',labels]; 
            %Keep only the labelled elements
            labels(t,:) = [];
            % Element wise operation...
            nl = size(labels,1);
            for i = 1:nl
                for j = i+1:nl
                    c = struct('x1',[],'x2',[],'link',[],'distance',[]);
                    d = pdist2(data(labels(i,1),:),data(labels(j,1),:),metric);
                    if labels(i,2)==labels(j,2)
                        %ML
                        if any(thresholdML==+inf)==0 %go inside only if there is no +inf
                            if nML > 1
                                if (d <= thresholdML(2) && d >= thresholdML(1)) 
                                    c.x1 = labels(i,1);
                                    c.x2 = labels(j,1);
                                    c.link = 29;  
                                    c.distance = d;
                                end
                            else
                                if d <= thresholdML
                                    c.x1 = labels(i,1);
                                    c.x2 = labels(j,1);
                                    c.link = 29;  
                                    c.distance = d;
                                end                            
                            end     
                        else
                            c.x1 = labels(i,1);
                            c.x2 = labels(j,1);
                            c.link = 29;  
                            c.distance = d;                           
                        end
                    else
                        %CL
                        if any(thresholdML==+inf)==0 %go inside only if there is no +inf
                            if nCL > 1
                                if d <= thresholdCL(2) && d >= thresholdCL(1)
                                    c.x1 = labels(i,1);
                                    c.x2 = labels(j,1);
                                    c.link = 31;  
                                    c.distance = d;
                                end
                            else
                                if d <= thresholdCL
                                    c.x1 = labels(i,1);
                                    c.x2 = labels(j,1);
                                    c.link = 31;  
                                    c.distance = d;
                                end                            
                            end       
                        else
                            c.x1 = labels(i,1);
                            c.x2 = labels(j,1);
                            c.link = 31;  
                            c.distance = d;                            
                        end
                    end
                    %assert(~isempty(c.x1),'generate_constraints, bug!')
                    constraints = [constraints;c];
                end
            end
            % Also make the two-column array format
            nc = length(constraints);
            constr = nan(nc,2);
            for i = 1:nc
                c1 = constraints(i).x1;
                c2 = constraints(i).x2;
                l = constraints(i).link;
                if l == 29
                    %ML: High to Low
                    if c1 > c2
                        constr(i,1) = c1;
                        constr(i,2) = c2;
                    else
                        constr(i,1) = c2;
                        constr(i,2) = c1;            
                    end
                elseif l == 31
                    %CL: Low to High
                    if c1 < c2
                        constr(i,1) = c1;
                        constr(i,2) = c2;
                    else
                        constr(i,1) = c2;
                        constr(i,2) = c1;            
                    end        
                end
            end                
    end
           
            
end