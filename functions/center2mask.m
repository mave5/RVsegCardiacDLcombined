% get centers and extract a mask

function mask=center2mask(I,cnt,Mroi)

cnt_round=round(cnt);
cnt_y=cnt_round(2,:);
cnt_x=cnt_round(1,:);

[y_max,x_max,z_max]=size(I);
mask=zeros(size(I));
M2=floor(Mroi/2);

for k=1:size(I,3)
    if cnt_y(k)>M2 && cnt_y(k)+M2<=y_max
        yroi=cnt_y(k)-M2:cnt_y(k)+M2;
    elseif cnt_y(k)<=M2
        yroi=1:Mroi;
    else
        yroi=y_max-Mroi+1:y_max;
    end
    
    if cnt_x(k)>M2 && cnt_x(k)+M2<=x_max
        xroi=cnt_x(k)-M2:cnt_x(k)+M2;
    elseif cnt_x(k)<=M2
        xroi=1:Mroi;
    else
        xroi=x_max-Mroi+1:x_max;
    end

    mask(yroi,xroi,k)=ones(Mroi,Mroi);
end

end
