% RV segmentation
clc 
close all
clear all
addpath('functions');

%imset='Test1Set';
imset='Test2Set';


% active contour parameters
max_its=150;
intEweight=0.5;
DLweight=.1/20;

% patient number
patient=42;

% slice number
psn=20

save_ena=0;
edit_ena=0;

%% convert patient number to string: XX 
pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% read images
[contours,I,~,~] =read_cont_imgs(patient,imset);
%contours_endo=contours.endo;
endo_cont_names=contours.names;

% crop squared images from the original images 
I1=I(:,:,psn);
%I1=imrotate(I1,90);
[Isq,rd1]=rnd_crop_images(I1);

% contour size
large_contour=0;
if large_contour==1
    Mroi=171;
else
    Mroi=91;
end

%for k=1:1%length(subI)
disp(['processing patient:',num2str(patient),', slice number:',num2str(psn)]); 

% get ROI using CNN
disp('localization')
%cnnfn='matFiles/DLconfigure/LargeContours/cnnRVsims_CVTest1_BP_filterDim7_poolDim2imageDim_32';
%cnnfn='matFiles/DLconfigure/LargeContours/cnnRVsims_CVTest1_BP_filterDim11_poolDim6imageDim_64';
cnnfn='matFiles/DLconfigure/LargeContours/cnnRV_BP_I216_filterDim10_poolDim5imageDim_54';
load (cnnfn,'OptTheta','outputSize','filterDim','numFilters','poolDim','meanItr','stdItr','imr1');
rI_sq = resize_normI(Isq,imr1,1,meanItr,stdItr);
[~, ~, cnt_pred] = cnnCost(OptTheta,rI_sq, randi(100,size(rI_sq,3),2),outputSize,filterDim, numFilters, ...
                                     poolDim,1);
ROI_auto=center2mask(Isq,cnt_pred,Mroi);
Iroi_pred=center2subI(I1,cnt_pred,Mroi);
%showCurveAndPhi(I1,ROI_auto);
%figure
%showCurveAndPhi(Iroi_pred);
subI1=Iroi_pred{1};
cnt_xy=cnt_pred;

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
init_mask1=DLN(subI1,nsubI1,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
init_mask2=clean_segs(init_mask1);
%showCurveAndPhi(subI1,init_mask2{psn});

if edit_ena==1

    temp1=edit_prior_shape(subI1,init_mask2,0);
    init_mask2=clean_segs(temp1);
    %showCurveAndPhi(subI1,temp0,init_mask2);
end

% retrun to I size
init_mask= remap_mask_cnt(init_mask2,I1,cnt_xy);

% perform deformable segmentation
[RV_seg1,phi] = region_seg_subPhi(subI1,cnt_xy,init_mask2,max_its,intEweight,DLweight,0);
RV_seg2=clean_segs(RV_seg1);
RV_seg_auto = remap_mask_cnt(RV_seg2,I1,cnt_xy);

subplot(1,2,1)
showCurveAndPhi(I1,init_mask)
title('initial contour');
subplot(1,2,2)
showCurveAndPhi(I1,RV_seg_auto);
title('final contour');
h=figure;
showCurveAndPhi(I1,RV_seg_auto);
title('final contour');
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
save_contours(RV_seg_auto,name2);
end

% save shape contour
endo_cnm2 = strrep(endo_cn, 'manual', 'auto');
fname2=[fname1,'shape/'];
if exist(fname2,'dir')==0
    mkdir(fname2);
end
name3=[fname2,endo_cnm2];
if save_ena==1
save_contours(init_mask,name3);
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
%para.width=size(I,1);para.height=size(I,2);
%temp=(contourc(RV_seg_auto(:,:,psn),[0 0]))';
%auto_contour=temp(2:end,:);
%dm = calc_dm(auto_contour,contours_endo{psn},para)
