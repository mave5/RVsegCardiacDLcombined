% RV segmentation
clc 
close all
clear all
addpath('functions');

%imset='TrainingSet';
imset='Test2Set';

% active contour parameters
max_its=250;
intEweight=0.5;
DLweight=.15;

% patient number
patient=48;

% slice number
psn=7

save_ena=1;
edit_ena=0;


%% convert patient number to string: XX 
pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% read contours and images
[contours,I,cnt_xy,diff_xym] =read_cont_imgs(patient,imset);
contours_endo=contours.endo;
endo_cont_names=contours.names;

% ROI size
Mroi=176;

disp(['processing patient:',num2str(patient),', slice number:',num2str(psn)]); 
I1=I(:,:,psn);    

% read contours
path2contours=['../python/results/',imset,'/pateint',pnstr]
%path2contours=['./Results/semantic/',imset,'/pateint',pnstr]
contours_cnn =read_contours(path2contours);
C1xy=contours_cnn{psn};
[m,n,z]=size(I);
cnt_xy=[m,n]/2;
init_mask=poly2mask(C1xy(:,1),C1xy(:,2),m,n);
showCurveAndPhi(I(:,:,psn),contours_cnn{psn});

% perform deformable segmentation
[RV_seg1,phi] = region_seg_subPhi(I1,cnt_xy,init_mask,max_its,intEweight,DLweight,1);
RV_seg_auto=clean_segs(RV_seg1);


subplot(1,2,1)
%showCurveAndPhi(I(:,:,psn),contours_endo{psn},init_mask)
showCurveAndPhi(I(:,:,psn),init_mask)
subplot(1,2,2)
%showCurveAndPhi(I(:,:,psn),contours_endo{psn},RV_seg_auto);
showCurveAndPhi(I(:,:,psn),RV_seg_auto);

%% save results
fname1=['Results/semantic/',imset,'/pateint',pnstr,'/'];
if exist(fname1,'dir')==0
    mkdir (fname1);
end

% save auto contour
endo_cn=char(endo_cont_names(psn));
endo_cnm = strrep(endo_cn, 'manual', 'auto');
name2=[fname1,endo_cnm];
if save_ena==1
    save_contours(double(RV_seg_auto),name2);
end


