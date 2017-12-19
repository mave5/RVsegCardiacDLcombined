%% shape prior using stacked autoencoder
clc;
clear all;
close all;
addpath('functions');
addpath('minFunc');
%% STEP 0: 
% parameters
patchsize = 91;
visibleSize = patchsize*patchsize;   % number of input units 
hiddenSizeL1 = 700;     % number of hidden units 
hiddenSizeL2=700;
sparsityParam1 = 0.1;   % desired average activation of the hidden units.
sparsityParam2=0.1;
lambda = 3e-3;       % weight decay parameter       
beta = 3;            % weight of sparsity penalty term       
outputSize=visibleSize; % number of output units

% STEP 1: laod training inputs and labels from mat files
%load matFiles/TrainingSet/all_patientsLarge; 
load matFiles/TrainingSet/all_patientsSmall;

% random permutation
p = randperm(size(t_Iroi_cat,3));
t_Iroi_cat_s=t_Iroi_cat(:,:,p);
t_y_endo_cat_s=t_y_endo_cat(:,:,p);

train_input=sampleIMAGES(t_Iroi_cat_s,patchsize);
train_labels=sampleIMAGES(t_y_endo_cat_s,patchsize);

%% train sparse Auto Encoder 1
%  Randomly initialize the parameters
saeTheta1 = initializeParameters(hiddenSizeL1, visibleSize);

%  Use minFunc to minimize the function
addpath minFunc/
options.Method = 'lbfgs'; % Here, we use L-BFGS to optimize our cost
                          % function. Generally, for minFunc to work, you
                          % need a function pointer with two outputs: the
                          % function value and the gradient. In our problem,
                          % sparseAutoencoderCost.m satisfies this.
options.maxIter = 400;	  % Maximum number of iterations of L-BFGS to run 
options.display = 'on';

[sae1OptTheta, cost] = minFunc( @(p) sparseAutoencoderCost(p, ...
                                   visibleSize, hiddenSizeL1, ...
                                   lambda, sparsityParam1, ...
                                   beta, train_input), ...
                                   saeTheta1, options);
%% STEP 5: Visualization of AE1
W1 = reshape(sae1OptTheta(1:hiddenSizeL1*visibleSize), hiddenSizeL1, visibleSize);
display_network(W1', 12); 
%% compute activations from layer 1
[sae1Features] = feedForwardAutoencoder(sae1OptTheta, hiddenSizeL1, ...
                                        visibleSize, train_input);

%% train sparse Auto Encoder 2                                   
%  Randomly initialize the parameters
sae2Theta = initializeParameters(hiddenSizeL2, hiddenSizeL1);

[sae2OptTheta, costL2] = minFunc( @(p) sparseAutoencoderCost(p, ...
                                  hiddenSizeL1, hiddenSizeL2, ...
                                  lambda, sparsityParam2, ...
                                  beta, sae1Features), ...
                                  sae2Theta, options);
W2 = reshape(sae2OptTheta(1:hiddenSizeL2*hiddenSizeL1), hiddenSizeL2, hiddenSizeL1);
%display_network(W2', 12); 
%% compute activation from layer 2
[sae2Features] = feedForwardAutoencoder(sae2OptTheta, hiddenSizeL2, ...
                                        hiddenSizeL1, sae1Features);

%% train multi outputs logstic regression                                    
lambda_mr=1e-4;
options_mr.maxIter = 100;
trainLabels=train_labels;
mrModel = mrTrain(hiddenSizeL2, outputSize, lambda_mr, ...
                            sae2Features, trainLabels, options_mr);

saeMultRegOptTheta = mrModel.optTheta(:);

%% fine tuning

% Initialize the stack using the parameters learned
stack = cell(2,1);
inputSize=visibleSize;

stack{1}.w = reshape(sae1OptTheta(1:hiddenSizeL1*inputSize), ...
                     hiddenSizeL1, inputSize);
stack{1}.b = sae1OptTheta(2*hiddenSizeL1*inputSize+1:2*hiddenSizeL1*inputSize+hiddenSizeL1);

stack{2}.w = reshape(sae2OptTheta(1:hiddenSizeL2*hiddenSizeL1), ...
                     hiddenSizeL2, hiddenSizeL1);
stack{2}.b = sae2OptTheta(2*hiddenSizeL2*hiddenSizeL1+1:2*hiddenSizeL2*hiddenSizeL1+hiddenSizeL2);

% Initialize the parameters for the deep model
[stackparams, netconfig] = stack2params(stack);
stackedAETheta = [ saeMultRegOptTheta ; stackparams ];


[stackedAEOptTheta, loss] = minFunc( @(x) stackedAECost(x, ...
      inputSize, hiddenSizeL2, outputSize, netconfig, ...
      lambda, train_input, train_labels), ...
      stackedAETheta, options);

%% test 
% load test data
load matFiles/TrainingSet/all_patientsSmall; 

testI=t_Iroi;
test_input=sampleIMAGES(testI,patchsize);

% predict
[pred_y_endo] = stackedAEPredict(stackedAEOptTheta, inputSize, hiddenSizeL2, ...
                          outputSize, netconfig, test_input);

% the final output is a mask of RV segmentation                      
yRVhr=reshape(pred_y_endo,patchsize,patchsize,[]);

% scale to image size
scale=size(t_Iroi,1)/patchsize;
for k=1:size(yRVhr,3)
    y1=yRVhr(:,:,k);
    yRV_h(:,:,k)=imresize(y1,scale);
    yRV_hc(:,:,k)=clean_segs(yRV_h(:,:,k));
end


%% dispaly segmentation
k1=0;
for k=1:30
    I1=testI(:,:,k);
    figure(1);
    k1=k1+1;
    subplot(5,6,k1)
    imagesc(I1);
    colormap(gray);hold on
    contour(yRV_h(:,:,k),[0 0],'r','LineWidth',2); 
    contour(t_y_endo(:,:,k),[0 0],'g','LineWidth',2);    
    contour(yRV_hc(:,:,k),[0 0],'b','LineWidth',2);    
end

% one title for all subplots
set(gcf,'NextPlot','add');
axes;
h = title(['HiddenSize=',num2str(hiddenSizeL1),' sparsity=',num2str(sparsityParam1)]);
set(gca,'Visible','off');
set(h,'Visible','on');

%% save results
filename=['matFiles/SimResults/SuffleSmallCont_V_',num2str(patchsize),'_RVseg_H1_',num2str(hiddenSizeL1),'_H2_',num2str(hiddenSizeL1),'_rho1_',num2str(sparsityParam1*100),'div100','_rho2_',num2str(sparsityParam2*100),'div100'];
save (filename);

