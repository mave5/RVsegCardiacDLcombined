% this function takes the image, and optimized parameters of the deep
% learning network and outputs a mask

function y=DLN(Ima,normImg,parameters,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig)
% Img: orginal image
%normImg: normalized image
% parameters: the learned/optimized parameters 
% inputSize : the visible size
% hiddenSizeL1 and L2: the size for layer 1 and 2
% outputSize
% netconfig

patchsize=sqrt(inputSize);
stackedAEOptTheta=parameters;

[pred_y] = stackedAEPredict(stackedAEOptTheta, inputSize, hiddenSizeL2, ...
                          outputSize, netconfig, normImg);

% mask output
y=reshape(pred_y,patchsize,patchsize,[]);

%showCurveAndPhi(reshape(normImg,patchsize,patchsize,[]),y);

scale=size(Ima,1)/patchsize;
y=imresize(y,scale);

end