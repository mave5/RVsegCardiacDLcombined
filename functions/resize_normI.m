function [Ir,meanIr,stdIr] = resize_normI(Im,scale,norm_ena,meanIr,stdIr)

% get size and number of images
[yn, xn, zn]=size(Im);

xn_s=xn*scale;
yn_s=yn*scale;
visibleSize = xn_s*yn_s;   % number of input units

imgs = imresize(Im, scale);
Ir=(reshape(imgs,visibleSize,zn));

%% ---------------------------------------------------------------
% For the autoencoder to work well we need to normalize the data
% Specifically, since the output of the network is bounded between [0,1]
% (due to the sigmoid activation function), we have to make sure 
% the range of pixel values is also bounded between [0,1]
if norm_ena==1 
    if nargin==3
        [Ir,meanIr,stdIr] = normalizeData(Ir);
    elseif nargin==5
        [Ir,meanIr,stdIr] = normalizeData(Ir,meanIr,stdIr);
    end
end

Ir=reshape(Ir,yn_s,xn_s,[]);
end

%% ---------------------------------------------------------------
function [patches,meanPatches,pstd] = normalizeData(patches,meanPatches,pstd)

% Squash data to [0.1, 0.9] since we use sigmoid as the activation
% function in the output layer

if nargin==1
% Remove DC (mean of images). 
meanPatches=mean(patches,2);
end
patches = bsxfun(@minus, patches, meanPatches);

% Truncate to +/-3 standard deviations and scale to -1 to 1
if nargin==1
pstd = 3 * std(patches(:));
end
patches = max(min(patches, pstd), -pstd) / pstd;

% Rescale from [-1,1] to [0.1,0.9]
patches = (patches + 1) * 0.4 + 0.1;

end

