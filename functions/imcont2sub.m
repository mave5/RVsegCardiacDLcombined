function [Iroi,yROI,y_endo,y_epi,contour_cnt,endo_seg,epi_seg]=imcont2sub(Iorig,C_endo,C_epi,Mroi)
% takes the original images and corresponding contours of endo and epi and
% extracts a sub image and corresponding binary masks


    [y_max,x_max]=size(Iorig);
    Mroi2=floor(Mroi/2);

    
    % read contour 
    Cx_endo=C_endo(:,1);Cy_endo=C_endo(:,2);

    % create segmentation mask from contoru
    endo_seg= roipoly(Iorig,Cx_endo, Cy_endo);    

    if ~isempty(C_epi)
        Cx_epi=C_epi(:,1);Cy_epi=C_epi(:,2);
        epi_seg= roipoly(Iorig,Cx_epi, Cy_epi);    
    else
            epi_seg=[];
    end
    


    % define a rectangle centered at contour
    cx_minmax=[round(min(C_endo(:,1))),round(max(C_endo(:,1)))];
    cx_cnt=round(mean(cx_minmax));
    cy_minmax=[round(min(C_endo(:,2))),round(max(C_endo(:,2)))];
    cy_cnt=round(mean(cy_minmax));
    contour_cnt=[cx_cnt,cy_cnt];
    
    % corners
    
    wx_corners=[cx_cnt-Mroi2,cx_cnt+Mroi2,cx_cnt+Mroi2,cx_cnt-Mroi2,cx_cnt-Mroi2];  
    if cx_cnt>Mroi2 && cx_cnt+Mroi2<=size(Iorig,2)
        wx_roi=cx_cnt-Mroi2:cx_cnt+Mroi2;
    elseif cx_cnt<=Mroi2
        wx_roi=1:Mroi;
    elseif cx_cnt+Mroi2> size(Iorig,2)
        wx_roi=size(Iorig,2)-Mroi+1:size(Iorig,2);
    end

    % corners
    wy_corners=[cy_cnt-Mroi2,cy_cnt-Mroi2,cy_cnt+Mroi2,cy_cnt+Mroi2,cy_cnt-Mroi2];
    if cy_cnt>Mroi2 &&  cy_cnt+Mroi2<=size(Iorig,1)
        wy_roi=cy_cnt-Mroi2:cy_cnt+Mroi2;
    elseif cy_cnt<=Mroi2
        wy_roi=1:Mroi;
    elseif cy_cnt+Mroi2> size(Iorig,1)     
        wy_roi=size(Iorig,1)-Mroi+1:size(Iorig,1);
    end
    
    % binary mask of ROI
    yROI=poly2mask(wx_corners,wy_corners,y_max,x_max);

    % binary masks of segmentation
    y_endo=endo_seg(wy_roi,wx_roi);
    if ~isempty(C_epi)
    y_epi=epi_seg(wy_roi,wx_roi);
    else
        y_epi=[];
    end
    % extract Iroi from original image    
    Iroi=Iorig(wy_roi,wx_roi);
    
end