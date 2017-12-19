% show manual segmentations
clear all   
close all
clc
addpath('functions')
%% load the image 
disp('Load MRI images');
imset='TrainingSet';
%imset='Test1Set';
fname1=['matFiles/',imset,'/',imset,'_all'];
%fname1=['matFiles/',imset,'/cnn_roi_test'];
load (fname1)

% patient
last_patient=0;
patient=last_patient+1;

% slince number for one patient
psn=1;

% max and min slice number 
min_sn=sum(slice_per_patient(1:patient-1-last_patient))+1;
max_sn=sum(slice_per_patient(1:patient-1-last_patient))+slice_per_patient(patient-last_patient);
slice_num=min_sn:max_sn;

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



subplot1(1,5,'Gap',[-0.1 -0.03])
k=0;
for psn=1:2:length(slice_num)
display(['processing slice# ',num2str(psn)])
% original image
I=t_I{patient}(:,:,psn);

% get center of ROI based on manual contour
C_endo=t_contours_endo{slice_num(psn)};

k=k+1;
subplot1(k)
showCurveAndPhi(I,C_endo);

end
%%
k=0;
subplot1(3,1,'Gap',[-0.1 -0.03])
for psn=1:1:3
display(['processing slice# ',num2str(psn)])
% original image
I=t_I{patient}(:,:,psn);

% get center of ROI based on manual contour
C_endo=t_contours_endo{slice_num(psn)};

k=k+1;
subplot1(k)
roi1=t_yROI{1}{1};
roi2=roi1(:,:,psn);
showCurveAndPhi(I,roi2);
Irot=imrotate(I,180);
showCurveAndPhi(Irot);

end

psn=1;
roi1=t_yROI{1}{1};
roi2=roi1(:,:,psn);
showCurveAndPhi(I,roi2);
showCurveAndPhi(t_Iroi(:,:,1),t_y_endo(:,:,1));