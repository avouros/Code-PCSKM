function delta = BinarySearchDelta(obj,s,ITER)
%% BinarySearchDelta uses the binary search algorithm to find a Delta so 
% that: Delta=0 if LASSO(w)<s; otherwise, Delta>0 so that LASSO(w)=s [1].
% Original R code can be found here: https://bit.ly/2PUFKim

% Reference:
% [1] Witten, Daniela M., and Robert Tibshirani. "A framework for feature 
%     selection in clustering." Journal of the American Statistical 
%     Association 105.490 (2010): 713-726.


% Input:
% - obj  : objective that we try to maximize.
% - s    : sparsity parameter.
% - ITER : (optional) number of iterations until Delta is detected.

% Output:
% - delta : a value for the delta parameter.


    % Default number of iterations
    if nargin < 3
        ITER = 15;
    end

    if norm(obj,2)==0 || sum(abs(obj/norm(obj,2)))<=s
        % If LASSO(w)<s then Delta=0
        delta = 0;
    else
        % Else find a Delta>0 so that LASSO(w)=s 
        lam1 = 0;
        lam2 = max(abs(obj))-1e-5;
        iter = 1;       
        while (iter<=ITER) && ((lam2-lam1)>(1e-4))
            % Soft-thresholding operator, S(x,c) = sign(x)(|x|-c)+
            su = sign(obj) .* max(0,abs(obj)-((lam1+lam2)/2)); 
            % Binary Search
            if sum(abs(su/norm(su,2))) < s
                lam2 = (lam1+lam2)/2;
            else
                lam1 = (lam1+lam2)/2;
            end 
            iter = iter + 1;
        end 
        delta = (lam1+lam2)/2;
    end
    
end

