function x = normalizations(x,option,varargin)
%NORMALIZATIONS, various normalizations for vectors / matrices

% List of available normalizations:
% - Center data around mean
% - Sum to 1
% - Euclidean or general vector norm, see MATLAB doc https://bit.ly/2AkEfX4
% - z-score, see https://bit.ly/2l5QFXV
% - Scale

%INPUT:
% x: vector or matrix: Each element of the vector or each column will be normalized.
% option: 
%  - 'mean':     center data around mean.
%  - 'one':      sum to 1.
%  - 'n-norm':   general vector norm, default: Euclidean see MATLAB doc
%                https://bit.ly/2AkEfX4 for more info
%  - 'z-score':  standarization, see https://bit.ly/2l5QFXV for more info
%  - 'z-scorem': modified z-score, https://ibm.co/2VV6IJZ
%  - 'scale':    scalling between an interval, default [0,1]

%EXTRA INPUT (varargin): 
% (1) For 'n-norm' option it specifies the power n.
%     If omitted then n = 2 (Euclidean). 
%     E.g. normalizations(x,'n-norm',5)
% (2) For 'scale' option it specifies the upper and lower bounds.
%     If omitted then [0 1] (min-max).
%     E.g. normalizations(x,'scale',[1,10])

%OUTPUT:
% x: the normalized input vector or matrix.

%NOTE:
%    (1) In case of matrix, the operation will be performed per column.  
%    (2) In case an input vector 1xN is provided then the input x will be
%        transposed then normalized and then transoded back to 1xN.


% Author:
% Avgoustinos Vouros
% avouros1@sheffield.ac.uk

%%
    
    [n,m] = size(x);
    if n==1
        x = x';
    end

    switch option
        case 'mean'
            % center data around mean (common in PCA)
            x = x - repmat(mean(x,1),n,1);
            
        case 'one'
            % sum to 1 (common in neural networks)
            x = x ./ repmat(sum(x,1),n,1);

        case 'n-norm'
            % sum(abs(x).^n)^(1/n) for 1 <= n < inf
            % if n = 2    => Euclidean
            % if n = inf  => max(abs(x))
            if isempty(varargin)
                a = 2;
            else
                a = varargin{1};
            end
            x = x ./ repmat(vecnorm(x,a,1),n,1);
            
        case 'z-score' 
            % zero mean and unit variance
            x = (x - repmat(mean(x,1),n,1)) ./ repmat(std(x,0,1),n,1);
            
        case 'z-scorem'   
            % modified z-score using median and mad
            mad = median(abs(x - median(x,1)),1);
            if any(mad==0)
                a = 1.2533;
                m = mean(x,1);
                if any(m==0)
                    %if the mean is also 0 then set use a very small number
                    x = (x - median(x,1)) ./ repmat((a.*realmin),n,1);
                else
                    x = (x - median(x,1)) ./ repmat((a.*m),n,1);
                end
            else
                a = 1.4826;
                x = (x - median(x,1)) ./ repmat((a*mad),n,1);
            end            

        case 'scale' 
            % First do min-max ([0 1], default)
            x = (x - repmat(min(x,[],1),n,1)) ./ repmat((max(x,[],1) - min(x,[],1)),n,1);
            % Then scale
            if ~isempty(varargin)
                if length(varargin{1}) == 2
                    a = sort(varargin{1},'ascend');
                    range = a(2) - a(1);
                    x = (x.*range) + a(1);
                end
            end
    end
    
    if n==1
        x = x';
    end
end

