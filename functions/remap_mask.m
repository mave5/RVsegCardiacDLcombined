%-- re-map mask to I

function out_mask = remap_mask(in_mask,I)

    [x_max, y_max]=size(I);
    m_cnt=ceil(x_max/2);
    
    M=size(in_mask,1);
    M2=floor(M/2);
    

    
    % top left corner
    x1=m_cnt-M2;
    y1=m_cnt-M2;

    % bottom right corner
    x4=m_cnt+M2;
    y4=m_cnt+M2;

    out_mask=zeros(x_max,y_max);
    out_mask(y1:y4,x1:x4)=in_mask;

    
end  