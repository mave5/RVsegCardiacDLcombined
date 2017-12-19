% create empty text files for epi
clear all   
close all
clc
addpath('functions')
%% load the image 
disp('Load MRI images');
imset='TrainingSet';
fname1=['matFiles/',imset,'/all_patients'];
load (fname1)

% patient
patient=1;


for psn=1:slice_per_patient(patient)
% max and min slice number 
min_sn=sum(slice_per_patient(1:patient-1))+1;
max_sn=sum(slice_per_patient(1:patient-1))+slice_per_patient(patient);
slice_num=min_sn:max_sn;


% get dicom info
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


% save results
fname1=['Results/',imset,'/patient',pnstr,'/'];
if exist(fname1,'dir')==0
mkdir (fname1);
end

% save auto contour
epi_cn=char(t_epi_cont_names(slice_num(psn)));
epi_cnm = strrep(epi_cn, 'manual', 'auto');
name2=[fname1,epi_cnm];

save_contours(t_y_epi(:,:,slice_num(psn)),name2);

end
