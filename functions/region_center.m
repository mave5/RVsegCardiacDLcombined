% find the center of a region based on Endo contour    
function cnt_xy=region_center(C_endo)

    % define a rectangle centered at contour
    cx_minmax=[round(min(C_endo(:,1))),round(max(C_endo(:,1)))];   
    cx_cnt2=round(mean(cx_minmax));
    
    cy_minmax=[round(min(C_endo(:,2))),round(max(C_endo(:,2)))];
    cy_cnt2=round(mean(cy_minmax));

    cnt_xy=[cx_cnt2,cy_cnt2];

%    Mroi_x=Mroi;
%    Mroi_y=Mroi;

%    cnt_wx=floor(Mroi_x/2);
%    cnt_wy=floor(Mroi_y/2);
    
%     wx_corners=[cx_cnt2-cnt_wx,cx_cnt2+cnt_wx,cx_cnt2+cnt_wx,cx_cnt2-cnt_wx,cx_cnt2-cnt_wx];
%     if cx_cnt2>cnt_wx
%         wx_roi=cx_cnt2-cnt_wx:cx_cnt2+cnt_wx;
%     else
%         wx_roi=1:Mroi_x;
%     end
%     wy_corners=[cy_cnt2-cnt_wy,cy_cnt2-cnt_wy,cy_cnt2+cnt_wy,cy_cnt2+cnt_wy,cy_cnt2-cnt_wy];
%     if cy_cnt2>cnt_wy
%         wy_roi=cy_cnt2-cnt_wy:cy_cnt2+cnt_wy;
%     else
%         wy_roi=1:Mroi_y;
%     end
   

end