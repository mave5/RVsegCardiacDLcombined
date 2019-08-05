function [auto_surf, manual_surf, DM, HD] = compare_rv_contours ( autocont, manualcont, N, M )

%  [auto_surf, manual_surf, DM, HD] = compare_rv_contours ( autocont, manualcont, N, M )
%
% This function compares two contours (autocont and manualcont), given the contours points and the image size (M-by-N).
% The function returns both contour areas, Dice metric and Hausdorff distance.
%
% INPUTS:
% -------
% autocont: automatic contour points
% manualcont: manual contour points
% N: number of columns of the image
% M: number of rows of the image
%
% OUTPUTS:
% -------
% autosurf: area of automatic contour
% manualsurf: area of manual contour
% DM: Dice metric between automatic contour and manual contour
% HD: Hausdorff Distance between automatic contour and manual contour

%   Copyright: LITIS EA 4108, Université de Rouen, France
%   Author: Caroline Petitjean (caroline.petitjean@univ-rouen.fr)
%   Revision: 1.0 - Date: April 20th, 2012

%Dice metric
auto_mask = mypoly2mask(autocont,M,N);
manual_mask = mypoly2mask(manualcont,M,N);

auto_surf = sum(auto_mask(:)>0);
manual_surf = sum(manual_mask(:)>0);

intersect_area = sum((auto_mask(:) + manual_mask(:))==2);
DM = 2*intersect_area/(auto_surf + manual_surf);

%Distance matrix
for i=1:length(autocont)
    for j=1:length(manualcont)
       distmat(i,j) = (autocont(i,1)-manualcont(j,1))*(autocont(i,1)-manualcont(j,1))+(autocont(i,2)-manualcont(j,2))*(autocont(i,2)-manualcont(j,2));
    end
end

distmat = sqrt(distmat);
d_m = min(distmat);
d_a = min(distmat,[],2);

%Hausdorff distance
HD = max([max(d_a) max(d_m)]);


