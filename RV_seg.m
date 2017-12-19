% two-dimensional segmenter combined with Deep Learning
clear all   
close all
clc
addpath('functions')

%% load the image 
disp('Load MRI images');
%imset='TrainingSet';
imset='Test1Set';
%fname1=['matFiles/',imset,'/',imset,'_all'];
fname1=['matFiles/',imset,'/Test1SetLargeCVrot'];
load (fname1)

% patient number
if strcmp(imset,'TrainingSet')
    last_patient=0;
elseif strcmp(imset,'Test1Set')
    last_patient=16;
end
patient_m=2;
patient=last_patient+patient_m;

% set slince number, for one patient
psn=1;

% set large or small contour
large_contour=1;

edit_ena=0;
save_ena=0;

% max and min slice number 
min_sn=sum(slice_per_patient(1:patient-1))+1;
max_sn=sum(slice_per_patient(1:patient-1))+slice_per_patient(patient);
slice_num=min_sn:max_sn;

% maximum iteration of active contour
max_its = 150;
intEweight=0.5;
DLweight=.1/2;

% enable plots
dis_ena=1;

% region of interest size
if large_contour==1
    Mroi=171;
else
    Mroi=91;
end
Mroi2=floor(Mroi/2);
%% get dicom info

% convert patient number to string: XX 
pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% patient folder
patient_folder=['dcom/',imset,'/patient',pnstr];

% contours folder
%contours_folder=[patient_folder,'/','P',pnstr,'contours-manual'];

% dcom folder
dicom_folder=[patient_folder,'/','P',pnstr,'dicom'];
para=get_dicominfo(dicom_folder);

display(['processing slice# ',num2str(psn)])

% original image
I=t_I{slice_num(psn)};
Isq=t_Isq(:,:,slice_num(psn));

% get center of ROI based on manual contour
C_endo=t_contours_endo{slice_num(psn)};
cnt_xy=region_center(C_endo);

%% get ROI using CNN
if large_contour==1
%cnnfn='matFiles/DLconfigure/LargeContours/cnnRVsims_CVTest1_BP_filterDim7_poolDim2imageDim_32';
cnnfn='matFiles/DLconfigure/LargeContours/cnnRVsims_CVTest1_BP_filterDim11_poolDim6imageDim_64';
load (cnnfn,'OptTheta','outputSize','filterDim','numFilters','poolDim','meanItr','stdItr','imr1');
rI_sq = resize_normI(Isq,imr1,1,meanItr,stdItr);
[~, ~, cnt_pred] = cnnCost(OptTheta,rI_sq, randi(100,1,2),outputSize,filterDim, numFilters, ...
                                     poolDim,1);
ROI_auto=center2mask(Isq,cnt_pred,Mroi);
ROI_man=center2mask(Isq,cnt_xy',Mroi);
%showCurveAndPhi(I,ROI_man,ROI_auto);
subI_pred=center2subI(I,cnt_pred,Mroi);
%showCurveAndPhi(subI_pred);
else
    
end

%% sub image, ROI
subI=t_Iroi(:,:,slice_num(psn));
%subI=subI_pred;
%cnt_xy=round(cnt_pred);
%[subI,cnt_xy]=edit_subI(I,subI,cnt_xy,Mroi,13);

% changing the size of the ROI if we want, useful for apex slices
subI2=extract_subI(subI,Mroi);

% load DL-RV parameters
if large_contour==1
filename='matFiles/DLconfigure/Normalized/NormSuffle_V_57_RVseg_H1_500_H2_500_rho1_10div100_rho2_10div100';
%filename='matFiles/DLconfigure/LargeContours/NormPCA_V_114_RVseg_H1_50_H2_50_rho1_30div100_rho2_30div100';
%filename='matFiles/DLconfigure/LargeContours/Rot1Cont_V_57_RVseg_H1_300_H2_300_rho1_10div100_rho2_10div100_lambda_100';
else
filename='matFiles/DLconfigure/SmallContours/Rot1SmallCont_V_91_RVseg_H1_50_H2_50_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1SmallCont_V_91_RVseg_H1_100_H2_100_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_91_RVseg_H1_100_H2_100_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_100_H2_100_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_30div100_rho2_30div100_lambda_100';
end
load (filename,'stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig','meanPatch','stdPatch');

% get initial mask using DL
patchsize=sqrt(inputSize);

if large_contour==1
   
    % normalize data
    nsubI=normalize_data(subI,patchsize,meanPatch,stdPatch);
    init_mask1=DLN(subI,nsubI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
    
    % clean and remove samll islands
    init_mask2=clean_segs(init_mask1);
    
    %init_mask22=edit_prior_shape(subI,init_mask2,1);
       
    % make shape same size of subI2
    init_mask3=extract_subI(init_mask2,Mroi);

    subplot(1,2,1); showCurveAndPhi(subI, [],init_mask1);title(['prior shape, I=',num2str(size(subI))]);
    contour(init_mask2,[0 0],'b','LineWidth',2);
else
    nsubI2=normalize_data(subI2,patchsize,meanPatch,stdPatch);
    init_mask1=DLN(subI2,nsubI2,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);

    % clean and remove samll islands
    init_mask2=clean_segs(init_mask1);
    %init_mask2=edit_prior_shape(subI2,init_mask2,1);
    
    
    % here, the shape is same size of subI2 
    init_mask3=init_mask2;

    subplot(1,2,1); showCurveAndPhi(subI2, [],init_mask1);title(['prior shape, I=',num2str(size(subI2))]);
    contour(init_mask2,[0 0],'b','LineWidth',2);
end

% retrun to I size
init_mask4= remap_mask_cnt(init_mask2,I,cnt_xy);

% dispaly initial mask
subplot(1,2,2); showCurveAndPhi(subI2, [],init_mask3); title(['prior shape, I=',num2str(size(subI2))]);

% run segmentation
figure;subplot(1,3,1);
[RV_seg1,phi,cnt_xy] = region_seg_subPhi(subI2,cnt_xy,init_mask3,max_its,intEweight,DLweight,1);
RV_seg2=clean_segs(RV_seg1);title(['auto contour, I=',num2str(size(subI2))]);
if edit_ena==1
RV_seg2=edit_prior_shape(subI2,RV_seg2,0);
RV_seg2=clean_segs(RV_seg2);
end
% resize mask to size subI
out_mask = remap_mask(RV_seg2,subI);
RV_seg3=out_mask;
subplot(1,3,2);showCurveAndPhi(subI, [],RV_seg3);title(['auto contour, I=',num2str(size(subI))]);
 
% compute metrics
manualPoints=t_contours_endo{slice_num(psn)};
% resize seg to I size
RV_seg4 = remap_mask_cnt(RV_seg2,I,cnt_xy);
subplot(1,3,3);showCurveAndPhi(I, manualPoints,RV_seg4);title(['auto and manual, I=',num2str(size(I))]);

% performance metrics
[dm_prior,PD_prior,HD_prior]=eval_metrics(init_mask4,manualPoints,para); 
[dm1,PD1,HD1]=eval_metrics(RV_seg4,manualPoints,para); 

%% evaluate segmentations 

% good contours
GP1=find(PD1<5);

PerDist=[PD_prior ,PD1]

DiceMetric=[dm_prior dm1]

% Hausdorff distance
HausD=[HD_prior HD1]

figure;
Ct=scaleContour(C_endo,cnt_xy,size(subI,1));
showCurveAndPhi(subI,Ct ,RV_seg3)

%% save results
fname1=['Results/',imset,'/patient',pnstr,'/'];
if exist(fname1,'dir')==0
mkdir (fname1);
end

% save auto contour
endo_cn=char(t_endo_cont_names(slice_num(psn)));
endo_cnm = strrep(endo_cn, 'manual', 'auto');
name2=[fname1,endo_cnm];
if save_ena==1
save_contours(RV_seg4,name2);
end

% save shape contour
endo_cnm2 = strrep(endo_cn, 'manual', 'shape');
fname2=[fname1,'shape/'];
if exist(fname2,'dir')==0
    mkdir(fname2);
end
name3=[fname2,endo_cnm2];
if save_ena==1
save_contours(init_mask4,name3);
end

% save mat file
modifiedStr1 = strrep(endo_cn, 'manual', 'auto');
endo_cnm = strrep(modifiedStr1, 'txt', 'mat');
name4=[fname1,endo_cnm];
%save (name4)
