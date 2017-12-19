
% read contours and images of a patient
function [contours,I,cnt_xy,diff_xym] =read_cont_imgs(patient,imset)

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

% flag to check if the contour-folder exists or not
check_cont_fold=0;

% read contours folder
if exist(contours_folder,'dir')
    check_cont_fold=1;
    cd(contours_folder);
    contours_list = dir('*.txt'); 
    numFiles = size(contours_list,1);
else
% get current folder
curr_dir=pwd;

% read list file
cd(patient_folder);

% there is one file here that lists the contour names
contours_list_file = dir('*.txt'); 
cl_name=contours_list_file.name;

% read list file
listoffiles = importdata(cl_name);
numFiles = length(listoffiles);

end

contours=struct; 
for k=1:2:numFiles
    % k1=k/2
    k1=floor(k/2)+1;

    % collect contour file names
    if check_cont_fold==1
        cname_i=contours_list(k).name;
        
        % read endo contours files
        contours_endo{k1}=load (cname_i);
        
        % get the center of contours
        [cnt_xy(:,k1),diff_xym(k1)]=region_center(contours_endo{k1});
        
    else
        cname_kt=listoffiles{k};  
        % extract the number part of the name: PXX-YYYY
        cname_i = [cname_kt(32:end)];
        
        contours_endo{k1}=[];
        cnt_xy(:,k1)=[0 ;0];
        diff_xym(k1)=0;
    end
    endo_contours_names{k1}=cname_i;

    % extract the number part of the name: PXX-YYYY
    dcom_names(k1,:) = [cname_i(1:8),'.dcm'];
end
contours.names=endo_contours_names;
contours.endo=contours_endo;

% return to current folder
cd(curr_dir);

% colect dicom images based on the list
[I, xthickness, ythickness, zthickness] = gatherImages(dicom_folder,dcom_names);

end

%%

% find the center of a region based on Endo contour    
function [cnt_xy,diff_xym]=region_center(C_endo)
    
    % define a rectangle centered at contour
    cx_minmax=[round(min(C_endo(:,1))),round(max(C_endo(:,1)))];   
    cx_diff=round(max(C_endo(:,1)))-round(min(C_endo(:,1)));   
    cx_cnt2=round(mean(cx_minmax));
    
    
    cy_minmax=[round(min(C_endo(:,2))),round(max(C_endo(:,2)))];
    cy_diff=round(max(C_endo(:,2)))-round(min(C_endo(:,2)));
    cy_cnt2=round(mean(cy_minmax));

    cnt_xy=[cx_cnt2,cy_cnt2]; 
    temp=max([cx_diff,cy_diff]); 
    diff_xym = (ceil((ceil(temp)/2)+0.5)*2)-1;
    
end

