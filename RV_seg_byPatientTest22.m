% RV segmentation
clc 
close all
clear all
addpath('functions');

%imset='TrainingSet';
%imset='Test1Set';
imset='Test2Set';


% active contour parameters
max_its=150;
intEweight=0.5;
DLweight=.2;

% patient number
patient=42;

% slice number
psn=20

save_ena=0;
edit_ena=0;

%% convert patient number to string: XX 
pnstr=[num2str(floor(patient/10)),num2str(rem(patient,10))];

% read images
[contours,I,~,~] =read_cont_imgs(patient,imset);
%contours_endo=contours.endo;
endo_cont_names=contours.names;

% crop squared images from the original images 
I1=I(:,:,psn);
[Isq,rd1]=rnd_crop_images(I1);


% contour size
large_contour=1;
if large_contour==1
    Mroi=171;
else
    Mroi=91;
end

%for k=1:1%length(subI)
disp(['processing patient:',num2str(patient),', slice number:',num2str(psn)]); 


showCurveAndPhi(I1);
h = imfreehand(gca,'closed',false); 
prior_shape_manual=createMask(h);
pos1 = getPosition(h);
[subI1,~,init_mask2,~,contour_cnt,~,~]=imcont2sub(I,pos1,[],Mroi);
cnt_xy =contour_cnt;


% retrun to I size
init_mask= remap_mask_cnt(init_mask2,I1,cnt_xy);

% perform deformable segmentation
[RV_seg1,phi] = region_seg_subPhi(subI1,cnt_xy,init_mask2,max_its,intEweight,DLweight,0);
RV_seg2=clean_segs(RV_seg1);
RV_seg_auto = remap_mask_cnt(RV_seg2,I1,cnt_xy);

subplot(1,2,1)
showCurveAndPhi(I1,init_mask)
title('initial contour');
subplot(1,2,2)
showCurveAndPhi(I1,RV_seg_auto);
title('final contour');
h=figure;
showCurveAndPhi(I1,RV_seg_auto);
title('final contour');
%%
%% save results
fname1=['Results/',imset,'/patient',pnstr,'/'];
if exist(fname1,'dir')==0
mkdir (fname1);
end

% save auto contour
endo_cn=char(endo_cont_names(psn));
endo_cnm = strrep(endo_cn, 'manual', 'auto');
name2=[fname1,endo_cnm];
if save_ena==1
save_contours(RV_seg_auto,name2);
end

% save shape contour
endo_cnm2 = strrep(endo_cn, 'manual', 'auto');
fname2=[fname1,'shape/'];
if exist(fname2,'dir')==0
    mkdir(fname2);
end
name3=[fname2,endo_cnm2];
if save_ena==1
save_contours(init_mask,name3);
end

% save mat file
fname3=[fname1,'figs/'];
if exist(fname3,'dir')==0
    mkdir(fname3);
end
endo_fignm = strrep(endo_cnm, '.txt', '');
name4=[fname3,endo_fignm];
saveas(h, name4, 'fig');

%%
%para.width=size(I,1);para.height=size(I,2);
%temp=(contourc(RV_seg_auto(:,:,psn),[0 0]))';
%auto_contour=temp(2:end,:);
%dm = calc_dm(auto_contour,contours_endo{psn},para)
