% compute metrics

function [dm1,PD1,HD1]=eval_metrics(auto_seg,manualPoints,para)
   
    % convert mask to contour
    autoPoints1=contourc(double(auto_seg), [0 0]); autoPoints1=autoPoints1(:,2:end)';
    
if ~isempty(autoPoints1) 
% Dice Metric
dm1 = calc_dm(autoPoints1,manualPoints,para);

% Perpendicular Distance
%PD1 = calc_dist(autoPoints1,manualPoints,para);
PD1=inf;

% Hausdorff distance
HD1 = hausdorff( manualPoints, autoPoints1); 
else
dm1=0;
dm2=0;
PD1=100;
PD2=100;
HD1=100;
HD2=100;
end

end