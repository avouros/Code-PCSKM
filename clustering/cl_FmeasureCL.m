function [f_score,accuracy,recall,specificity,precision,TP,TN,FP,FN] = cl_FmeasureCL(true_values, predicted_values)
%cl_FmeasureCL. Computes various measurements from the information theory.

% INPUT:
%  true_values: data labels
%  predicted_values: clustering results

% OUTPUT:
%  f_score
%  accuracy
%  recall
%  specificity
%  precision
%  TP = true positive
%  TN = true negative
%  FP = false positive
%  FN = false negative

% Formulas obtained from: https://goo.gl/M9TXfM

    if size(true_values) ~= size(predicted_values)
        error('cl_entropy error: vectors of true and predicted values needs to have the same size');
    end

    classes = unique(true_values);
    nclasses = length(classes);
    total = length(true_values); %number of datapoints
    
    total_pairs = total*(total-1)/2;
    
    TPFP = 0;
    for i = 1:nclasses
        n = length(find(predicted_values == classes(i)));
        TPFP = TPFP + (n*(n-1)/2);
    end    
    
    TP = 0;
    for i = 1:nclasses
        for j = 1:nclasses
            tc = find(true_values == classes(i));
            pc = find(predicted_values == classes(j));
            n = length(intersect(tc,pc));
            TP = TP + (n*(n-1)/2);
        end
    end
    
    FP = TPFP - TP;
    
    
    FNTN = total_pairs - TPFP;
    FN = 0;
    for i = 1:nclasses
        for j = 1:nclasses
            tc = find(true_values == classes(i));
            pc = find(predicted_values == classes(j));    
            sn = length(intersect(tc,pc));
            n = 0;
            for k = j+1:nclasses
                pc = find(predicted_values == classes(k)); 
                n = n + length(intersect(tc,pc));
            end
            FN = FN + (sn * n);
        end
    end
    
    TN = FNTN - FN;
      

    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    f_score = (2*precision*recall) / (precision+recall);
    specificity = TN/(FP+TN);
    accuracy = TP/(TP+FP);
    
    % Special cases (trick to avoid division by 0)
    %(reference: https://goo.gl/hyv6e8)
    if TP==0 && FP==0 && FN==0
        precision = 1;
        recall = 1;
    end
    if TP==0
        if FP==0 || FN==0
            precision = 0;
            recall = 0;
        end
    end
    if precision == 0 && recall == 0
        f_score = 0;
    elseif precision == 1 && recall == 1
        f_score = 1;   
    end    
    
    % Assert final value
    if max(f_score) > 1 && min(f_score) < 0
        error('cl_FmeasureCL error: Bug found!');
    end
end
