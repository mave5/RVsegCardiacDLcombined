% read dcom images and contours and save them into .mat files
clc 
close all
clear all
addpath('functions');

save_ena=0;

imset='TrainingSet';
%imset='Test1Set';

%% patient number
patient=1;

% region of interest
Mroi_x=171;
Mroi_y=171;

% convert patient number to string: XX 
pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% patient folder
patient_folder=['dcom/',imset,'/patient',pnstr];

% contours folder
contours_folder=[patient_folder,'/','P',pnstr,'contours-manual'];

% dcom folder
dicom_folder=[patient_folder,'/','P',pnstr,'dicom'];

% get current folder
curr_dir=pwd;

% read contours folder
cd(contours_folder);
contours_list = dir('*.txt'); 
numFiles = size(contours_list,1);

nLarge2x=20;
%nLarge2x=numFiles;
for k=1:2:nLarge2x%numFiles
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

% get current directory
curr_dir=pwd;

% colect dicom images based on the list
[IMAGES, xthickness, ythickness, zthickness] = gatherImages(dicom_folder,dcom_names);

% make image square
[y_max,x_max,z_max]=size(IMAGES);
%imgDim=min([x_max,y_max]);
imgDim=max([x_max,y_max]);
%IMAGES=IMAGES(1:imgDim,1:imgDim,:);
%IMAGES=imresize(IMAGES,[256 256]);
temp=zeros(imgDim,imgDim,z_max);
temp(1:y_max,1:x_max,:)=IMAGES;
IMAGES=temp;

nLarge=z_max;
yROI=zeros(imgDim,imgDim,nLarge);
Iroi=zeros(Mroi_x,Mroi_x,nLarge);

for k=1:z_max
    % show images
    figure (1)
    subplot(3,ceil(z_max/3),k)
    I1=IMAGES(:,:,k);
    imagesc(I1); colormap(gray);
    title(dcom_names(k,:))
    hold on
   
    % read contour 
    C_endo=contours_endo{k};
    C_epi=contours_epi{k};
    Cx_endo=C_endo(:,1);Cy_endo=C_endo(:,2);
    Cx_epi=C_epi(:,1);Cy_epi=C_epi(:,2);
    
    % plot it
    [junk,xcnt,ycnt]=polycenter(Cx_endo,Cy_endo);
    plot(Cx_endo,Cy_endo,'r','LineWidth',2);   
    plot(Cx_epi,Cy_epi,'g','LineWidth',2);   
    plot(xcnt,ycnt,'r*','markerSize',12);
        
    % create segmentation mask and display it
    figure(4)
    subplot(3,ceil(z_max/3),k)
    endo_seg(:,:,k)= roipoly(I1,Cx_endo, Cy_endo);    
    epi_seg(:,:,k)= roipoly(I1,Cx_epi, Cy_epi);    
    imagesc(endo_seg(:,:,k));
    colormap(gray)
    title(['mask',num2str(k)])
    hold on
    plot(Cx_endo,Cy_endo,'b','LineWidth',2)
    plot(xcnt,ycnt,'b*','markerSize',12); 
   
    x_cntr=round(xcnt);
    y_cntr=round(ycnt);
    contour_center{k}=[x_cntr,y_cntr];
    
    % define a rectangle centered at contour
    cx_minmax=[round(min(C_endo(:,1))),round(max(C_endo(:,1)))];
    cx_cnt2=round(mean(cx_minmax));
    cy_minmax=[round(min(C_endo(:,2))),round(max(C_endo(:,2)))];
    cy_cnt2=round(mean(cy_minmax));
    
    cnt_wx=floor(Mroi_x/2);
    cnt_wy=floor(Mroi_y/2);
    wx_corners=[cx_cnt2-cnt_wx,cx_cnt2+cnt_wx,cx_cnt2+cnt_wx,cx_cnt2-cnt_wx,cx_cnt2-cnt_wx];
    if cx_cnt2>cnt_wx
        wx_roi=cx_cnt2-cnt_wx:cx_cnt2+cnt_wx;
    else
        wx_roi=1:Mroi_x;
    end
    wy_corners=[cy_cnt2-cnt_wy,cy_cnt2-cnt_wy,cy_cnt2+cnt_wy,cy_cnt2+cnt_wy,cy_cnt2-cnt_wy];
    if cy_cnt2>cnt_wy
        wy_roi=cy_cnt2-cnt_wy:cy_cnt2+cnt_wy;
    else
        wy_roi=1:Mroi_y;
    end
    % create ROI mask: this will be used for training of Automatic
    % detection network
    figure(2)
    subplot(3,ceil(z_max/3),k)
    yROI(:,:,k)=poly2mask(wx_corners,wy_corners,imgDim,imgDim);
    imshow(yROI(:,:,k));hold on
    contour(yROI(:,:,k),[0 0],'r');hold on
    plot(xcnt,ycnt,'b*','markerSize',12);
    
    % find ROI in the mask: this will be used as training data for shape
    % infering
    figure(5)
    subplot(3,ceil(z_max/3),k)
    y_endo(:,:,k)=endo_seg(wy_roi,wx_roi,k);
    y_epi(:,:,k)=epi_seg(wy_roi,wx_roi,k);
    imshow(y_endo(:,:,k));
    hold on
    plot(Mroi_x/2,Mroi_x/2,'r*','markerSize',12);

    % ROI in the image: this will be used as input training data for shape infering 
    figure(3)
    subplot(3,ceil(z_max/3),k)
    Iroi(:,:,k)=I1(wy_roi,wx_roi);
    imagesc(Iroi(:,:,k));
    colormap(gray);
    hold on
    plot(Mroi_x/2,Mroi_x/2,'r*','markerSize',12);
    contour(y_endo(:,:,k),[0 0],'r','LineWidth',2)
    contour(y_epi(:,:,k),[0 0],'g','LineWidth',2)
    %plot(Cx_LV-x_cnt+Mroi/2+1,Cy_LV-y_cnt+Mroi/2+1,'r')      
end
%% store images and region of interest on disk

I=IMAGES(:,:,1:nLarge);
filename=['matFiles/',imset,'/LargeContours','/patient',pnstr];

if save_ena==1
save (filename, 'I', 'yROI','Iroi','y_epi','y_endo','contours_endo','contours_epi',...
 'patient','xthickness', 'ythickness' ,'zthickness','contour_center','endo_contours_names','epi_contours_names');
 disp('results saved')
end
 
% make sure sizes are matched
sizes=[size(I,3),size(Iroi,3),size(y_endo,3),size(yROI,3)]
