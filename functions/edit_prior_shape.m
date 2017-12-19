% prior shape might be inaccurate for bottom slices
% to get better prior shape, we use intersection between otsu and prior
% shape

function BW=edit_prior_shape(I,prior,disp_ena)

prior=logical((prior));

I1=I.*prior;

% otsu thresholding
ot1=otsu(I1)>1;

% intersection of otsu and prior shape
com1=logical(ot1).*prior;

% find the center contour
BW2=imfill(com1,'holes');
cc = bwconncomp(BW2,4);
stats = regionprops(cc, 'Area','Centroid');
[~,st_idx]=sort([stats.Area],'descend');
for k=1:min([2,length(stats)])%length(stats)
    cntk=stats(st_idx(k)).Centroid;
    diff_cnt(k)=norm(cntk-size(BW2,1));
end
[diffmin,idx]=min(diff_cnt);

com2= ismember(labelmatrix(cc), st_idx(idx));

if sum(com2(:))<5
    com2=prior;
end

BW=com2;


 if disp_ena==1
     imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
     contour(prior, [0 0], 'g','LineWidth',2);
     contour(BW, [0 0], 'r','LineWidth',2);
     contour(com1, [0 0], 'b','LineWidth',2);
     legend('original','edited');
 end

end

