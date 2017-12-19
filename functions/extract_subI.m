% extract a sub image from an original image
function subI=extract_subI(I,Mroi)
    for k=1:size(I,3)
        roi_x=round((size(I,1)-Mroi)/2)+1:round((size(I,1)+Mroi)/2);
        roi_y=round((size(I,2)-Mroi)/2)+1:round((size(I,2)+Mroi)/2);
        subI(:,:,k)=I(roi_y,roi_x,k);
    end
end
