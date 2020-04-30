clear
close all
clc

%% Percentage of used constraints

datastr = {'fisheriris';'ionospere';'TT-SC-ST';'Brod40-3'};

Ns = [150,351,424,120];
  
cv = 10;
nconstr = 1:0.5:10;

n = length(Ns);       %number of datasets
nc = length(nconstr); %number of constraints
total_labels = Ns;    %total number of labels
total_cv_labels = Ns.*(100-cv)./100;       %total number of cross validation labels
total_constraints = ((Ns.^2)-Ns)./2; %total number of constraints
total_cv_constraints = ((total_cv_labels.^2)-total_cv_labels)/2;
used_cv_constraints = zeros(n,nc);       %row=dataset, column=constraints
used_cv_labels = zeros(n,nc);            %row=dataset, column=constraints
used_cv_labels_per = zeros(n,nc);        %row=dataset, column=constraints
used_cv_labels_per_total = zeros(n,nc);  %row=dataset, column=constraints

for i = 1:n
    for j = 1:nc
        used_cv_constraints(i,j) = (total_cv_constraints(i)*nconstr(j))/100;
        tmp1 = ( -1 + sqrt( (1^2)-4*1*(-2)*used_cv_constraints(i,j) ) )/2;
        tmp2 = ( -1 - sqrt( (1^2)-4*1*(-2)*used_cv_constraints(i,j) ) )/2;
        used_cv_labels(i,j) = max([tmp1,tmp2]);
        used_cv_labels_per(i,j) = 100*used_cv_labels(i,j)/total_cv_labels(i);
        used_cv_labels_per_total(i,j) = 100*used_cv_labels(i,j)/Ns(i);
    end
end


f1 = figure;
ax = axes(f1);
bar(used_cv_labels_per_total(1,:));
title('Percentage of used labels per fold over dataset size');
%set(ax,'XTick',1:length(datastr),'XTickLabel',datastr,'XTickLabelRotation',45);
set(ax,'XTick',1:length(nconstr),'XTickLabel',nconstr,'XTickLabelRotation',0);
%set(get(ax,'XLabel'), 'String', 'dataset name');
set(get(ax,'XLabel'), 'String', 'constraints (%)');
set(get(ax,'YLabel'), 'String', 'labels (%)');

f2 = figure;
ax = axes(f2);
bar(used_cv_labels_per(1,:));
title('Percentage of used labels per fold over fold size');
%set(ax,'XTick',1:length(datastr),'XTickLabel',datastr,'XTickLabelRotation',45);
set(ax,'XTick',1:length(nconstr),'XTickLabel',nconstr,'XTickLabelRotation',0);
%set(get(ax,'XLabel'), 'String', 'dataset name');
set(get(ax,'XLabel'), 'String', 'constraints (%)');
set(get(ax,'YLabel'), 'String', 'labels (%)');

