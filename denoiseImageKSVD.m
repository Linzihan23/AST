function [IOut,output] = denoiseImageKSVD(Image,sigma,K,varargin)
%==========================================================================
%compressed sensing overcomplete chunked adaptive sparse dictionary             
% ===================================================================
% input parameters : Image - noisy image (grayscale)
%                    sigma - Gaussian random noise window
%                    K - indicates the number of atoms in the dictionary 256.
%                    varargin - list of variable length input parameters.
% Adjustable parameters : blockSize - the size of the block on which the algorithm works.
%                         errorFactor - The allowed error factor to be represented. 
%                         maxBlocksToConsider - the maximum number of blocks that can be processed.
%                         slidingFactor - The sliding distance between processed blocks. 
%                         waitBarOn - Shows the progress of the algorithm.
% output parameters :IOut - Reconstructs the image.
%                    output - String structure:
%                    D - dictionary for denoising
% =========================================================================
%% parameter setting
reduceDC = 1;
[NN1,NN2] = size(Image);
waitBarOn = 1;
C = 1.15;
maxBlocksToConsider = 512*512;
slidingDis = 1;
bb = 8;
maxNumBlocksToTrainOn = 65000
displayFlag = 1;
%% Number of dictionary updates
numIterOfKsvd = 10;
% % if (sigma > 5)
% %     numIterOfKsvd = 10;
% % else
% %     numIterOfKsvd = 5;
% % end
%% Preprocessing
for argI = 1:2:length(varargin)
    if (strcmp(varargin{argI}, 'slidingFactor'))
        slidingDis = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'errorFactor'))
        C = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'maxBlocksToConsider'))
        maxBlocksToConsider = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'numKSVDIters'))
        numIterOfKsvd = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'blockSize'))
        bb = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'maxNumBlocksToTrainOn'))
        maxNumBlocksToTrainOn = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'displayFlag'))
        displayFlag = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'waitBarOn'))
        waitBarOn = varargin{argI+1};
    end
end

tic

%% Training dictionary data preparation
if(prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn)
    randPermutation =  randperm(prod([NN1,NN2]-bb+1));
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);

    blkMatrix = zeros(bb^2,maxNumBlocksToTrainOn);
    for i = 1:maxNumBlocksToTrainOn
        [row,col] = ind2sub(size(Image)-bb+1,selectedBlocks(i));
        currBlock = Image(row:row+bb-1,col:col+bb-1);
        blkMatrix(:,i) = currBlock(:);
    end
else
    blkMatrix = im2col(Image,[bb,bb],'sliding');
end

param.K = K;
param.numIteration = numIterOfKsvd ;
param.errorFlag = 1;
param.errorGoal = sigma*C;
param.preserveDCAtom = 0;

%% DCT Framework Initial Dictionary
Pn=ceil(sqrt(K));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

param.initialDictionary = DCT(:,1:param.K );
param.InitializationMethod =  'GivenMatrix';

if (reduceDC)
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end

if (waitBarOn)
    counterForWaitBar = param.numIteration+1;
    h = waitbar(0,'Denoising In Process ...');
    param.waitBarHandle = h;
    param.counterForWaitBar = counterForWaitBar;
end

param.displayProgress = displayFlag;
toc
fprintf('Database pre-processing time is %s s.\n',toc);

%% The most central function
tic
[Dictionary,output] = KSVD(blkMatrix,param);
toc
fprintf('Dictionary learning time is %s s.\n',toc);
% % load('Dictionary_Noiseless_new1');
output.D = Dictionary;
if (displayFlag)
    disp('finished Trainning dictionary');
end

%% Image Processing
tic
errT = sigma*C;
IMout=zeros(NN1,NN2);
Weight=zeros(NN1,NN2);
%blocks = im2col(Image,[NN1,NN2],[bb,bb],'sliding');
while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end

[blocks,idx] = my_im2col(Image,[bb,bb],slidingDis);

if (waitBarOn)
    newCounterForWaitBar = (param.numIteration+1)*size(blocks,2);
end

for jj = 1:30000:size(blocks,2)
    if (waitBarOn)
        waitbar(((param.numIteration*size(blocks,2))+jj)/newCounterForWaitBar);
    end
    jumpSize = min(jj+30000-1,size(blocks,2));
    if (reduceDC)
        vecOfMeans = mean(blocks(:,jj:jumpSize));
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
    end
    
    %Coefs = mexOMPerrIterative(blocks(:,jj:jumpSize),Dictionary,errT);
    Coefs = OMPerr(Dictionary,blocks(:,jj:jumpSize),errT);
    if (reduceDC)
        blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    else
        blocks(:,jj:jumpSize)= Dictionary*Coefs ;
    end
end
toc
fprintf('Image reconstruction time is %s s.\n',toc);

count = 1;
Weight = zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(Image)-bb+1,idx);
for i  = 1:length(cols)
    col = cols(i); 
    row = rows(i);        
    block =reshape(blocks(:,count),[bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
    count = count+1;
end;

if (waitBarOn)
    close(h);
end

%% Output
IOut = (Image+0.034*sigma*IMout)./(1+0.034*sigma*Weight);


