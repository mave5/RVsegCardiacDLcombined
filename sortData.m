% read mat files and store them as a 3D matrix for further processing
clc 
close all
clear all
%%
imset='TrainingSet';

folder=['matFiles/',imset];
D=dir([folder,'\*.mat']);
numpatients = length(D(not([D.isdir])));

t_I=[];
t_yROI=[];
t_Iroi=[];
t_y_endo=[];
t_y_epi=[];
t_contours_endo=[];
t_contours_epi=[];
t_centers=[];
t_endo_cont_names=[];
t_epi_cont_names=[];

for k=1:numpatients
    patient=k
    % convert patient number to string: XX 
    pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];
    
    imgname=['matFiles/',imset,'/patient',pnstr];
    load (imgname);

    t_I{k}=I;
    t_yROI{k}={yROI};
    t_Iroi=cat(3,t_Iroi,Iroi);
    
    t_y_endo=cat(3,t_y_endo,y_endo);
    t_y_epi=cat(3,t_y_epi,y_epi);
    
    t_contours_endo=[t_contours_endo,contours_endo(1:size(I,3))];
    t_contours_epi=[t_contours_epi,contours_epi(1:size(I,3))];
    
    t_centers=[t_centers,contour_center(1:size(I,3))];
    slice_per_patient(k)=size(I,3);
    t_endo_cont_names=[t_endo_cont_names,endo_contours_names];
    t_epi_cont_names=[t_epi_cont_names,epi_contours_names];
end

filename=['matFiles/',imset,'/all_patients'];

save (filename,'t_I','t_yROI','t_Iroi','t_y_endo','t_y_epi','t_contours_endo','t_contours_epi','t_centers','numpatients',...
'slice_per_patient','t_endo_cont_names','t_epi_cont_names');
    
display('file has been saved')    