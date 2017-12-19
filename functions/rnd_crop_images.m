% randomly crop square images make images square to the min dimension
function [Im_sq,rd1]=rnd_crop_images(inImage,rd1)


if iscell(inImage)
    num_imgs=length(inImage);
else
    num_imgs=size(inImage,3);
end

for k=1:num_imgs
    
    if iscell(inImage)
        Im1=inImage{k};
    else
        Im1=inImage(:,:,k);
    end
    [y_max,x_max]=size(Im1);
    if y_max>x_max
        d1=y_max-x_max;
        if nargin==1
            rd1(k)=randi(d1/2);
        end
        temp=Im1(rd1:x_max+rd1-1,:);
    else
        d1=x_max-y_max;
        if nargin==1
            rd1(k)=randi(d1/2);
        end
        temp=Im1(:,rd1:y_max+rd1-1);
    end
    Im_sq(:,:,k)=temp;
end

end