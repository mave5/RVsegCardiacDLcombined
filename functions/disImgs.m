% display gray images
function disImgs(varargin)
% I : input images as a 2D/3D matrix

I=varargin{1};
if length(varargin)==2
    C=varargin{2};
else
    C=I;
end

if iscell(I)
     z=length(I);
else
    [y,x,z]=size(I);
end

npx=ceil(4);
npy=ceil(z/4);

figure
for k=1:z
    if iscell(I)
        I1=I{k};
    else
        I1=I(:,:,k);
    end

    if iscell(I)
        C1=C{k};
    else
        C1=C(:,:,k);
    end

    subplot(npx,npy,k);
    showCurveAndPhi(I1,C1);
end


end