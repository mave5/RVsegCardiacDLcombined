% two-dimensional segmenter combined with Deep Learning
clear all   
close all
clc
addpath('functions')
%% load the image 
disp('Load MRI images');
imset='TrainingSet';
imsubset='Small';
fname1=['matFiles/',imset,'/',imsubset,'Contours/TrainingSet',imsubset,'.mat'];
load (fname1)

% patient
last_patient=0;
patient=last_patient+1;

% slince number for one patient
psn=1;
if strcmp(imsubset,'Large')
    large_contour=1;
else
    large_contour=0;
end
edit_ena=0;
save_ena=0;

% max and min slice number 
min_sn=sum(slice_per_patient(1:patient-1-last_patient))+1;
max_sn=sum(slice_per_patient(1:patient-1-last_patient))+slice_per_patient(patient-last_patient);
slice_num=min_sn:max_sn;

% maximum iteration of active contour
max_its = 150;
intEweight=0.5;
DLweight=.2;

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

% dcom folder
dicom_folder=[patient_folder,'/','P',pnstr,'dicom'];
para=get_dicominfo(dicom_folder);

display(['processing slice# ',num2str(psn)])

% original image
I=t_I{patient}(:,:,psn);

% get center of ROI based on manual contour
C_endo=t_contours_endo{slice_num(psn)};
cnt_xy=region_center(C_endo);

% sub image, ROI
ddy=cnt_xy(2)+randi([-10 10]);ddx=cnt_xy(1)+randi([-10 10]);
%ddy=cnt_xy(2)-7;ddx=cnt_xy(1)+10;
subI1=I(ddy-Mroi2:ddy+Mroi2,ddx-Mroi2:ddx+Mroi2);
subI2=I(cnt_xy(2)-Mroi2:cnt_xy(2)+Mroi2,cnt_xy(1)-Mroi2:cnt_xy(1)+Mroi2);


% load DL-RV parameters
%filename='matFiles/DLconfigure/Normalized/NormSuffle_V_57_RVseg_H1_500_H2_500_rho1_20div100_rho2_20div100';
%filename='matFiles/DLconfigure/Normalized/NormSuffle_V_57_RVseg_H1_100_H2_100_rho1_10div100_rho2_10div100';
filename='matFiles/DLconfigure/SmallContours/PCA3_AugSmallCont_V_91_RVseg_H1_50_H2_50_rho1_30div100_rho2_30div100_lambda_300.mat';

%load alltraindata_mean_std;
load (filename,'stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig','meanPatch','stdPatch');

% get initial mask using DL
patchsize=sqrt(inputSize);

% normalize data
nsubI=normalize_data(subI1,patchsize,meanPatch,stdPatch);
init_mask1=DLN(subI1,nsubI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
    
% clean and remove samll islands
init_mask2=clean_segs(init_mask1);
       
subplot(2,2,1); showCurveAndPhi(subI1, [],init_mask1);title(['trained movement, cnt=',num2str([ddx ddy])]);
contour(init_mask2,[0 0],'b','LineWidth',2);



% normalize data
nsubI=normalize_data(subI2,patchsize,meanPatch,stdPatch);
init_mask1=DLN(subI2,nsubI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
    
% clean and remove samll islands
init_mask2=clean_segs(init_mask1);
       
subplot(2,2,2); showCurveAndPhi(subI2, [],init_mask1);title(['trained movement, cnt=',num2str([cnt_xy(1) cnt_xy(2)])]);
contour(init_mask2,[0 0],'b','LineWidth',2);


%% using DL configure that did not train movement
% filename='matFiles/DLconfigure/InputSize_57_RVseg_H1_750_H2_750_rho1_10div100_rho2_10div100';
% load (filename,'stackedAEOptTheta','inputSize','hiddenSizeL1','hiddenSizeL2','outputSize','netconfig');
% load alltraindata_mean_std;
% 
% nsubI=normalize_data(subI1,patchsize,mean(pmean),pstd);
% init_mask1=DLN(subI1,nsubI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
%     
% clean and remove samll islands
% init_mask2=clean_segs(init_mask1);
%        
% 
% subplot(2,2,3); showCurveAndPhi(subI1, [],init_mask1);title(['did not train movement, cnt=',num2str([ddx ddy])]);
% contour(init_mask2,[0 0],'b','LineWidth',2);
% 
% %
% nsubI=normalize_data(subI2,patchsize,mean(pmean),pstd);
% init_mask1=DLN(subI2,nsubI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
%     
% clean and remove samll islands
% init_mask2=clean_segs(init_mask1);
%        
% subplot(2,2,4); showCurveAndPhi(subI2, [],init_mask1);title(['did not train movement, cnt=',num2str([cnt_xy(1) cnt_xy(2)])]);
% contour(init_mask2,[0 0],'b','LineWidth',2);
% 
% 
