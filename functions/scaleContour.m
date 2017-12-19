
% scale a given contour based on a center and widow size
function   Ct=scaleContour(C,cnt,M)
% C : input contour
% cnt : center
% M   : window size

m1=ceil(M/2);

Ctx=C(:,1)-cnt(1)+m1;
Cty=C(:,2)-cnt(2)+m1;
Ct=[Ctx,Cty];

end
