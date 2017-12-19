% two-dimensional segmenter combined with Deep Learning
clear all   
close all
clc
addpath('functions')


%% load images 
disp('Load MRI images');
%imset='TrainingSet';
imset='Test1Set';
fname1=['matFiles/',imset,'/',imset,'_all.mat'];
%fname1=['matFiles/',imset,'/cnn_roi_test'];
load (fname1)

% patient
if strcmp(imset,'TrainingSet')
    last_patient=0;
elseif strcmp(imset,'Test1Set')
    last_patient=16;
end
patient_m=6;
patient=last_patient+patient_m;

disp(['Number of slices for patient',num2str(patient), ' => ', num2str(slice_per_patient(patient_m)) ])

% slice number for one patient
psn=1;
large_contour=1;
edit_ena=0;
save_ena=0;
clean_ena=1;

% max and min slice number 
min_sn=sum(slice_per_patient(1:patient_m-1))+1;
max_sn=sum(slice_per_patient(1:patient_m-1))+slice_per_patient(patient-last_patient);
slice_num=min_sn:max_sn;

% maximum iteration of active contour
max_its = 150;
intEweight=0.5;
DLweight=.1;

% enable plots
dis_ena=1;

% region of interest size
if large_contour==1
    Mroi=171;
    Mroi2=floor(Mroi/2);
else
    Mroi=91;
    Mroi2=floor(Mroi/2);
end
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
I=t_I{patient_m}(:,:,psn);

showCurveAndPhi(I);
h = imfreehand(gca,'closed',false); 
prior_shape_manual=createMask(h);
pos1 = getPosition(h);
[Iroi,~,prior_sm,~,contour_cnt,~,~]=imcont2sub(I,pos1,[],Mroi);
%cnt_xy =  round(ginput(1)); 
cnt_xy =contour_cnt;
%yroi=cnt_xy(2)-Mroi2:cnt_xy(2)+Mroi2;
%xroi=cnt_xy(1)-Mroi2:cnt_xy(1)+Mroi2;

% sub image, ROI
%subI=I(yroi,xroi);
subI=Iroi;
%subI=t_Iroi(:,:,slice_num(psn));
figure
showCurveAndPhi(subI,prior_sm);

% changing the size of the ROI if we want, useful for apex slices
subI2=extract_subI(subI,Mroi);

% get center of ROI based on manual contour
%cnt_xy=cnt(psn,:);

% load DL-RV parameters
if large_contour==1
    filename='matFiles/DLconfigure/LargeContours/Rot1Cont_V_57_RVseg_H1_300_H2_300_rho1_10div100_rho2_10div100_lambda_100';
    %filename='matFiles/DLconfigure/LargeContours/Rot1Cont_V_114_RVseg_H1_50_H2_50_rho1_10div100_rho2_10div100_lambda_100';
else
    filename='matFiles/DLconfigure/SmallContours/Rot1SmallCont_V_91_RVseg_H1_50_H2_50_rho1_10div100_rho2_10div100_lambda_100';
    %filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_10div100_rho2_10div100_lambda_100';
    %filename='matFiles/DLconfigure/SmallContours/Rot1Cont_V_45_RVseg_H1_200_H2_200_rho1_30div100_rho2_30div100_lambda_100';
end
load (filename,'stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig','meanPatch','stdPatch','Iroi_cv');

% get initial mask using DL
patchsize=sqrt(inputSize);

if large_contour==1
    %subI=Iroi_cv(:,:,randi(size(Iroi_cv,3)));
    % normalize data
    nsubI=normalize_data(subI,patchsize,meanPatch,stdPatch);
    %nsubI=normalize_data(subI,patchsize,mean(pmean),pstd);
    init_mask1=DLN(subI,nsubI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
    
    % clean and remove samll islands
    init_mask2=clean_segs(init_mask1);
       
    % make shape same size of subI2
    init_mask3=extract_subI(init_mask2,Mroi);

    subplot(1,2,1); showCurveAndPhi(subI, [],init_mask1);title(['prior shape, I=',num2str(size(subI))]);
    contour(init_mask2,[0 0],'b','LineWidth',2);
else
    nsubI2=normalize_data(subI2,patchsize,meanPatch,stdPatch);
    %nsubI2=sampleIMAGES(subI2,patchsize,1);
    %nsubI2=normalize_data(subI2,patchsize,mean(pmean),pstd);
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

%init_mask3=prior_sm;
% dispaly initial mask
subplot(1,2,2); showCurveAndPhi(subI2, [],init_mask3); title(['prior shape, I=',num2str(size(subI2))]);

% run segmentation
figure;subplot(1,3,1);
[RV_seg1,phi,cnt_xy] = region_seg_subPhi(subI2,cnt_xy,init_mask3,max_its,intEweight,DLweight,1);
if clean_ena==1
    RV_seg2=clean_segs(RV_seg1);
else 
    RV_seg2=RV_seg1;
end
if edit_ena==1
RV_seg2=edit_prior_shape(subI2,RV_seg2,0);
RV_seg2=clean_segs(RV_seg2);
end
% resize mask to size subI
out_mask = remap_mask(RV_seg2,subI);
RV_seg3=out_mask;
subplot(1,3,2);showCurveAndPhi(subI, [],RV_seg3);title(['auto contour, I=',num2str(size(subI))]);
 
% compute metrics
% resize seg to I size
RV_seg4 = remap_mask_cnt(RV_seg2,I,cnt_xy);
subplot(1,3,3);showCurveAndPhi(I, [],RV_seg4);title(['auto and manual, I=',num2str(size(I))]);


%% evaluate segmentations 
figure;
showCurveAndPhi(subI,[] ,RV_seg3)

%% save results
fname1=['Results/',imset,'/patient',pnstr,'/'];
if exist(fname1,'dir')==0
mkdir (fname1);
end

fname2=['matFiles/',imset,'/Test1Set_all'];
load (fname2, 't_listoffiles');
% save auto contour
temp1=t_listoffiles{patient-last_patient};
cname_i=char(temp1(2*(psn-1)+1));
endo_cn = [cname_i(32:end)];
endo_cnm = strrep(endo_cn, 'manual', 'auto');
name2=[fname1,endo_cnm];
if save_ena==1
save_contours(RV_seg4,name2);
end

%%

% convert patient number to string: XX 
%pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% patient folder
%patient_folder=['dcom/',imset,'/patient',pnstr];

% contours folder
contours_folder=[patient_folder,'/','P',pnstr,'contours-manual'];

% dcom folder
%dicom_folder=[patient_folder,'/','P',pnstr,'dicom'];

% get current folder
curr_dir=pwd;

% read contours folder
cd(contours_folder);
contours_list = dir('*.txt'); 
numFiles = size(contours_list,1);


for k=1:2:numFiles
    % k1=k/2
    k1=floor(k/2)+1;

    % collect contour file names
    cname_i=contours_list(k).name;
    endo_contours_names{k1}=cname_i;
   
    cname_o=contours_list(k+1).name;
    epi_contours_names{k1}=cname_o;
        
    % read endo contours files
    contours_endo{k1}=load (cname_i);
    
    % read epi contours files
    contours_epi{k1}=load (cname_o);
    
    % extract the number part of the name: PXX-YYYY
    dcom_names(k1,:) = [cname_i(1:8),'.dcm'];
end

% return to current folder
cd(curr_dir);

manC=contours_endo{psn};
Ct=scaleContour(manC,cnt_xy,Mroi);
plot(Ct(:,1),Ct(:,2));
