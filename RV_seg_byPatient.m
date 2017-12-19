% RV segmentation
clc 
close all
clear all
addpath('functions');

%imset='TrainingSet';
imset='Test1Set';


% active contour parameters
max_its=100;
intEweight=0.5;
DLweight=.35;

% patient number
patient=28;

% slice number
psn=15

save_ena=0;
edit_ena=0;


%% convert patient number to string: XX 
pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% read contours and images
[contours,I,cnt_xy,diff_xym] =read_cont_imgs(patient,imset);
contours_endo=contours.endo;
endo_cont_names=contours.names;

% get sub images
Mroi=diff_xym;
thw=51;
Mroi(diff_xym>thw)=171;
Mroi(diff_xym<=thw)=91;
subI=center2subI(I,cnt_xy,Mroi);
%disImgs(subI)

%for k=1:1%length(subI)
disp(['processing patient:',num2str(patient),', slice number:',num2str(psn)]); 

I1=I(:,:,psn);    
subI1=subI{psn};

%load temp1
%showCurveAndPhi(I1,prev_cont,contours_endo{psn});
%prev_cont2=center2subI(prev_cont,cnt_xy,Mroi);
%showCurveAndPhi(subI1,prev_cont2);


% check if it large contour or small contour
if size(subI1,1)>91
    large_contour=1;
else
    large_contour=0;
end

% load DL parameters
if large_contour==1
    filename='matFiles/DLconfigure/LargeContours/Rot1Cont_V_57_RVseg_H1_300_H2_300_rho1_10div100_rho2_10div100_lambda_100';
    %filename='matFiles/DLconfigure/LargeContours/Rot1Cont_V_114_RVseg_H1_50_H2_50_rho1_10div100_rho2_10div100_lambda_100';
else
    filename='matFiles/DLconfigure/SmallContours/Rot1SmallCont_V_91_RVseg_H1_50_H2_50_rho1_10div100_rho2_10div100_lambda_100';
    %filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_10div100_rho2_10div100_lambda_100';
    %filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_30div100_rho2_30div100_lambda_100';
end
load (filename,'stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig','meanPatch','stdPatch');
patchsize=sqrt(inputSize);

% normalize image
nsubI1=normalize_data(subI1,patchsize,meanPatch,stdPatch);

% perform deep learning segmentation
init_mask1{psn}=DLN(subI1,nsubI1,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
init_mask2{psn}=clean_segs(init_mask1{psn});
%showCurveAndPhi(subI1,init_mask2{psn});

if edit_ena==1
temp0=init_mask2{psn};
temp1=edit_prior_shape(subI1,init_mask2{psn},0);
init_mask2{psn}=clean_segs(temp1);
%showCurveAndPhi(subI1,temp0,init_mask2{psn});
end

% retrun to I size
init_mask(:,:,psn)= remap_mask_cnt(init_mask2{psn},I1,cnt_xy(:,psn));

% perform deformable segmentation
[RV_seg1{psn},phi{psn}] = region_seg_subPhi(subI1,cnt_xy(:,psn),init_mask2{psn},max_its,intEweight,DLweight,0);
RV_seg2{psn}=clean_segs(RV_seg1{psn});
RV_seg_auto(:,:,psn) = remap_mask_cnt(RV_seg2{psn},I(:,:,psn),cnt_xy(:,psn));

%end

% dispaly image and segmentation
%l1=size(RV_seg_auto,3);
%nf_row=ceil(l1/4);nf_col=ceil(l1/nf_row);
%subplot1(nf_row, nf_col, 'Gap', [-.07 -.02]);
%for k=1:size(RV_seg_auto,3)
%    subplot1(k);
    %showCurveAndPhi(I(:,:,k),contours_endo{k},init_mask(:,:,k),RV_seg_auto(:,:,k))
%    showCurveAndPhi(I(:,:,psn),contours_endo{psn},RV_seg_auto(:,:,psn));
%end
subplot(1,2,1)
showCurveAndPhi(I(:,:,psn),contours_endo{psn},init_mask(:,:,psn))
subplot(1,2,2)
showCurveAndPhi(I(:,:,psn),contours_endo{psn},RV_seg_auto(:,:,psn));
h=figure
showCurveAndPhi(I(:,:,psn),RV_seg_auto(:,:,psn));
%%
%% save results
fname1=['Results/',imset,'/patient',pnstr,'/'];
if exist(fname1,'dir')==0
mkdir (fname1);
end

% save auto contour
endo_cn=char(endo_cont_names(psn));
endo_cnm = strrep(endo_cn, 'manual', 'auto');
name2=[fname1,endo_cnm];
if save_ena==1
save_contours(RV_seg_auto(:,:,psn),name2);
end

% save shape contour
endo_cnm2 = strrep(endo_cn, 'manual', 'auto');
fname2=[fname1,'shape/'];
if exist(fname2,'dir')==0
    mkdir(fname2);
end
name3=[fname2,endo_cnm2];
if save_ena==1
save_contours(init_mask(:,:,psn),name3);
end

% save mat file
fname3=[fname1,'figs/'];
if exist(fname3,'dir')==0
    mkdir(fname3);
end
endo_fignm = strrep(endo_cnm, '.txt', '');
name4=[fname3,endo_fignm];
saveas(h, name4, 'fig');

%%
para.width=size(I,1);para.height=size(I,2);
temp=(contourc(RV_seg_auto(:,:,psn),[0 0]))';
auto_contour=temp(2:end,:);
dm = calc_dm(auto_contour,contours_endo{psn},para)
