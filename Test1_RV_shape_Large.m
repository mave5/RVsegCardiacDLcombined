% RV shape Test
clear all   
close all
clc
addpath('functions')

%% load the image 
disp('Load MRI images');
imset='Test1Set';
fname1=['matFiles/',imset,'/Test1SetLargeCVrot'];
%imset='TrainingSet';
%fname1=['matFiles/TrainingSet/LargeContours/TrainingSetLarge'];
load (fname1,'t_I','t_contours_endo')
[t_Isq,rd1]=rnd_crop_images(t_I);

%endo_masks=contour2mask(t_contours_endo,t_I);
%t_endo_sq=rnd_crop_images(endo_masks,rd1);
%k=randi(size(t_Isq,3));
%showCurveAndPhi(t_Isq(:,:,k),t_endo_sq(:,:,k));

% patient number
if strcmp(imset,'TrainingSet')
    last_patient=0;
elseif strcmp(imset,'Test1Set')
    last_patient=16;
end
patient_m=1;
patient=last_patient+patient_m;

% set large or small contour
large_contour=1;

save_ena=0;
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

% dcom folder
dicom_folder=[patient_folder,'/','P',pnstr,'dicom'];
para=get_dicominfo(dicom_folder);


% get center of ROI based on manual contour
for k=1:length(t_contours_endo)
    C_endo=t_contours_endo{k};
    cnt_xy(:,k)=region_center(C_endo);
end

%% get ROI using CNN
disp('localization')
if large_contour==1
%cnnfn='matFiles/DLconfigure/LargeContours/cnnRVsims_CVTest1_BP_filterDim7_poolDim2imageDim_32';
%cnnfn='matFiles/DLconfigure/LargeContours/cnnRVsims_CVTest1_BP_filterDim11_poolDim6imageDim_64';
cnnfn='matFiles/DLconfigure/LargeContours/cnnRV_BP_I216_filterDim10_poolDim5imageDim_54';
load (cnnfn,'OptTheta','outputSize','filterDim','numFilters','poolDim','meanItr','stdItr','imr1');
rI_sq = resize_normI(t_Isq,imr1,1,meanItr,stdItr);
[~, ~, cnt_pred] = cnnCost(OptTheta,rI_sq, randi(100,size(rI_sq,3),2),outputSize,filterDim, numFilters, ...
                                     poolDim,1);
ROI_auto=center2mask(t_Isq,cnt_pred,Mroi);
ROI_man=center2mask(t_Isq,cnt_xy,Mroi);
%showCurveAndPhi(I,ROI_man,ROI_auto);
%showCurveAndPhi(subI_pred);
mean((abs(cnt_xy-cnt_pred)./abs(cnt_xy))')*100;
end

for k=1:length(t_I)
t_Iroi_pred(:,:,k)=center2subI(t_I{k},cnt_pred(:,k),Mroi);
end


%% sub image, ROI

% load DL-RV parameters
if large_contour==1
%filename='matFiles/DLconfigure/Normalized/NormSuffle_V_57_RVseg_H1_500_H2_500_rho1_10div100_rho2_10div100';
%filename='matFiles/DLconfigure/LargeContours/NormPCA_V_114_RVseg_H1_50_H2_50_rho1_30div100_rho2_30div100';
filename='matFiles/DLconfigure/LargeContours/Rot1Cont_V_57_RVseg_H1_300_H2_300_rho1_10div100_rho2_10div100_lambda_100';
else
filename='matFiles/DLconfigure/SmallContours/Rot1SmallCont_V_91_RVseg_H1_50_H2_50_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1SmallCont_V_91_RVseg_H1_100_H2_100_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_91_RVseg_H1_100_H2_100_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_100_H2_100_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_10div100_rho2_10div100_lambda_100';
%filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_30div100_rho2_30div100_lambda_100';
end
load (filename,'stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig','meanPatch','stdPatch');

for k=1:size(t_Iroi_pred,3)
disp('processing images');k
I=t_I{k};
subI=t_Iroi_pred(:,:,k);
cnt_xy1=cnt_pred;
%[subI_n,cnt_n]=edit_subI(I,subI,cnt_xy1(:,k),Mroi);

%subI=t_Iroi(:,:,k);
%cnt_xy1=cnt_xy;



% changing the size of the ROI if we want, useful for apex slices
subI2=extract_subI(subI,Mroi);

% get initial mask using DL
patchsize=sqrt(inputSize);

if large_contour==1
   
    % normalize data
    nsubI=normalize_data(subI,patchsize,meanPatch,stdPatch);
    init_mask1=DLN(subI,nsubI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
    
    % clean and remove samll islands
    init_mask2=clean_segs(init_mask1);
         
    % make shape same size of subI2
    init_mask3=extract_subI(init_mask2,Mroi);

    %subplot(1,2,1); showCurveAndPhi(subI, [],init_mask1);title(['prior shape, I=',num2str(size(subI))]);
    %contour(init_mask2,[0 0],'b','LineWidth',2);
else
    nsubI2=normalize_data(subI2,patchsize,meanPatch,stdPatch);
    init_mask1=DLN(subI2,nsubI2,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);

    % clean and remove samll islands
    init_mask2=clean_segs(init_mask1);
    %init_mask2=edit_prior_shape(subI2,init_mask2,1);
    
    
    % here, the shape is same size of subI2 
    init_mask3=init_mask2;

    %subplot(1,2,1); showCurveAndPhi(subI2, [],init_mask1);title(['prior shape, I=',num2str(size(subI2))]);
    %contour(init_mask2,[0 0],'b','LineWidth',2);
end


% retrun to I size
init_mask4{k}= remap_mask_cnt(init_mask2,I,cnt_xy1(:,k));

[dm_prior(k),PD_prior(k),HD_prior(k)]=eval_metrics(init_mask4{k},t_contours_endo{k},para); 
end

for k=1:length(t_contours_endo)
    auto_seg=init_mask4{k};
    autoPoints1=contourc(double(auto_seg), [0 0]); autoPoints1=autoPoints1(:,2:end)';
    cnt_pred2(:,k)=region_center(autoPoints1);
end
%%
%% dispaly segmentation
k1=0;
rndk1=randi(length(init_mask4),1,30);
for k=1:30
    I1=t_I{rndk1(k)};
    k1=k1+1;
    subplot(5,6,k1)
    showCurveAndPhi(I1,t_contours_endo{rndk1(k)},init_mask4{rndk1(k)});      
end

mean(dm_prior)
mean(HD_prior)

k1=0;
figure
%rndk1=randi(length(init_mask4),1,30);
for k=1:30
    I1=t_I{rndk1(k)};
    k1=k1+1;
    subplot(5,6,k1)
    showCurveAndPhi(I1,ROI_man(:,:,rndk1(k)),ROI_auto(:,:,rndk1(k)));      
end
