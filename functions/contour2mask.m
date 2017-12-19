% convert contours to binary masks

function masks=contour2mask(conts,tI)


for k=1:length(conts)
    I=tI{k};
    [ymax,xmax]=size(I);    
    C1=conts{k};
    C1x=C1(:,1);C1y=C1(:,2);
    masks{k}=poly2mask(C1x,C1y,ymax,xmax);
end


end