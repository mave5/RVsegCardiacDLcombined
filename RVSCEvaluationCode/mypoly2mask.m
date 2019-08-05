function BW = mypoly2mask(pts,M,N)
%MYPOLY2MASK Convert region polygon to region mask.
%   BW = MYPOLY2MASK(pts,M,N) computes a binary region-of-interest mask,
%   BW, from a region-of-interest polygon, represented by the matrice pts.
%   The size of BW is M-by-N.  Pixels in BW that are inside
%   the polygon are 1; pixels outside the polygon are 0.  The class
%   of BW is logical.
%
%   MYPOLY2MASK closes the polygon automatically if it isn't already
%   closed.
%
%   MYPOLY2MASK is an adaptation from POLY2MASK, a MathWorks program.

%  It is assumed that number of points is > 2
if size(pts,1) == 2
    x = pts(1,:);
    y = pts(2,:);
else
    x = pts(:,1);
    y = pts(:,2);
end

if (x(end) ~= x(1)) | (y(end) ~= y(1))
    x(end+1) = x(1);
    y(end+1) = y(1);
end

[xe,ye] = poly2edgelist(x,y);
BW = edgelist2mask(M,N,xe,ye);
