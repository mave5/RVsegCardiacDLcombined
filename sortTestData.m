% read mat files and store them as a 3D matrix for further processing
clc 
close all
clear all
%%
imset='Test1Set';

save_ena=1;

folder=['matFiles/',imset];
D=dir([folder,'\*.mat']);
numpatients = length(D(not([D.isdir])));

t_I=[];
t_Isq=[];
t_endo_cont_names=[];

for k=1:numpatients
    patient=k+17-1
    % convert patient number to string: XX 
    pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];
    
    imgname=['matFiles/',imset,'/patient',pnstr];
    load (imgname);

    t_I{k}=I;
    slice_per_patient(k)=size(I,3);
    t_Isq=cat(3,t_Isq,I_sq);   
    
    t_listoffiles{k}=listoffiles;
    
    xyz_thickness(k,:)=[xthickness,ythickness,zthickness];
    
end

filename=['matFiles/',imset,'/',imset,'_all.mat'];

if save_ena==1
save (filename,'t_I','t_Isq','numpatients',...
'slice_per_patient','t_listoffiles','xyz_thickness');
display('file has been saved')    
end    

