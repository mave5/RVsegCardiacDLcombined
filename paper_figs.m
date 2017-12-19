clc
clear all
close all
addpath functions
%%
fn1='matFiles/TrainingSet/TrainingSet_all.mat';
load (fn1,'t_Iroi','t_y_endo','t_I')

nf_row=2;nf_col=4;
subplot1(nf_row, nf_col, 'Gap', [-.07 -.02]);
kk=[5,45,95,125];

for k=1:length(kk)
    subplot1(k)
    showCurveAndPhi(t_Iroi(:,:,kk(k)));
    subplot1(k+4)
    imshow(t_y_endo(:,:,kk(k)));
end
%%
close all
Ip1=t_I{10};
I1=Ip1(:,:,1);
r1=randi(40)
I2=I1(r1:r1+216-1,:);
showCurveAndPhi(I2);
showCurveAndPhi(I1);

% sub Image
figure
showCurveAndPhi(subI);

% sub image with initial contour
figure
showCurveAndPhi(subI, [],init_mask2)

% original image with final contour
figure
showCurveAndPhi(I, [],RV_seg4)


% resized image 
Isq=zeros(256,256);
Isq(:,1:216)=I;
Isqr=imresize(I,[64 64]);
showCurveAndPhi(Isq)


contour_cntx=contour_cnt(1);
contour_cnty=contour_cnt(2);
X=contour_cntx-Mroi2:contour_cntx+Mroi2;
Y=contour_cnty-Mroi2:contour_cnty+Mroi2;
ROI=zeros(256,256);
ROI(Y,X)=1;
showCurveAndPhi(I,ROI)
plot(contour_cntx,contour_cnty,'+')
