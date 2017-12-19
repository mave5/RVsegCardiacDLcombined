function patches = normalize_data(I,patchsize,pmean,pstd)

visibleSize = patchsize*patchsize;   % number of input units

% get size and number of images
[xn, yn, zn]=size(I);

scale=patchsize/xn;

imgs = imresize(I, scale);
patches=(reshape(imgs,visibleSize,zn));

%% ---------------------------------------------------------------
% For the autoencoder to work well we need to normalize the data
% Specifically, since the output of the network is bounded between [0,1]
% (due to the sigmoid activation function), we have to make sure 
% the range of pixel values is also bounded between [0,1]

patches = normalizeData(patches,pmean,pstd);

end

%% ---------------------------------------------------------------
function patches = normalizeData(patches,pmean,pstd)

% Squash data to [0.1, 0.9] since we use sigmoid as the activation
% function in the output layer

% Remove DC (mean of images). 
patches = bsxfun(@minus, patches, pmean);

% Truncate to +/-3 standard deviations and scale to -1 to 1
patches = max(min(patches, pstd), -pstd) / pstd;

% Rescale from [-1,1] to [0.1,0.9]
patches = (patches + 1) * 0.4 + 0.1;

end
