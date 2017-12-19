
function [subI_n,cnt_n]=edit_subI(I,subI,cnt,Mroi,reg_nn)

[ymax,xmax]=size(subI);

% otsu
o1=otsu(subI)>1;
subplot(2,1,1);
imshow(o1);
cc = bwconncomp(o1,4);
stats = regionprops(cc, 'Area','Centroid');
cnt1 =sort( [stats.Area]);
BW1 = ismember(labelmatrix(cc), reg_nn);
BW1 = imfill(BW1,'holes');

BW2= remap_mask_cnt(BW1,I,cnt);
c1=contourc(double(BW2), [0 0]); c1=c1(:,2:end)';
cnt_n=region_center(c1)';
subI_n=center2subI(I,cnt_n,Mroi);
subplot(2,1,2);
showCurveAndPhi(I,BW2);
end