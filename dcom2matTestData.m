% read dcom images and contours and save them into mat files
clc 
close all
clear all
addpath('functions');

%imset='TrainingSet';
imset='Test1Set';

save_ena=0;

%% patient number
patient=32;

% convert patient number to string: XX 
pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% patient folder
patient_folder=['dcom/',imset,'/patient',pnstr];

% dcom folder
dicom_folder=[patient_folder,'/','P',pnstr,'dicom'];

% get current folder
curr_dir=pwd;

% read list file
cd(patient_folder);
contours_list = dir('*.txt'); 
cname_i=contours_list.name;

% read list file
listoffiles = importdata(cname_i);
numFiles = length(listoffiles);

cd(curr_dir);

for k=1:2:numFiles
    % k1=k/2
    k1=floor(k/2)+1;

    % collect contour file names
    cname_i=listoffiles{k};  
   
    % extract the number part of the name: PXX-YYYY
    dcom_names(k1,:) = [cname_i(32:39),'.dcm'];
end

% get current directory
curr_dir=pwd;

% colect dicom images based on the list
[IMAGES, xthickness, ythickness, zthickness] = gatherImages(dicom_folder,dcom_names);

% make images squred size
[y_max,x_max,z_max]=size(IMAGES);
imgDim=max([x_max,y_max]);
temp=zeros(imgDim,imgDim,z_max);
temp(1:y_max,1:x_max,:)=IMAGES;
IMAGES_sq=temp;

for k=1:z_max
    % show images
    subplot(3,ceil(z_max/3),k)
    I1=IMAGES_sq(:,:,k);
    imagesc(I1); colormap(gray);
    title(dcom_names(k,:))
    hold on
   
end
%% store images and region of interest on disk

I=IMAGES;
I_sq=IMAGES_sq;
filename=['matFiles/',imset,'/patient',pnstr];

if save_ena==1
    save (filename, 'I','I_sq' ,'patient','xthickness', 'ythickness' ,'zthickness','listoffiles');
    disp('results saved')
end

disp(['patient:',num2str(patient),'    list of files:',num2str(size(I,3))])
