% =========================================================================
% Test code for digital holographic microscopy noise reduction method
%Implementation of image noise reduction using sparse dictionaries
%Environment: Win10, Matlab2020a
%Time: 2023-7-20
% =========================================================================

%% Clear residual data
clear all;
close all;
clc;

%% Data reading and pre-processing
load('date1.mat')
bb=8; 
RR=4; 
K=RR*bb^2; 
sigma = 25; 

% % newrealphase=zeros(512);
% % for i=1:1:512
% %     x=realphase(:,i);
% %     t=1:1:512;    
% %     snr=5;
% %     px_dBW=1;
% %     y1=awgn(x,snr,px_dBW);
% %     newrealphase(:,i)=y1;
% % end
load('date2.mat')

PSNRIn = 20*log10(255/sqrt(mean((newrealphase(:)-realphase(:)).^2)));

%% Gaussian random noise
% % newrealphase=realphase+sigma*randn(size(realphase));
% % PSNRIn = 20*log10(255/sqrt(mean((newrealphase(:)-realphase(:)).^2)));
% % figure(),imshow(newrealphase,[]);
% % figure,surf(newrealphase);
% % shading interp;
%% AST Denoising Test
tic
[IoutAdaptive,output] = denoiseImageKSVD(newrealphase,sigma,K);

% fprintf('Total time is %s s.\n',toc);

PSNROut = 20*log10(255/sqrt(mean((IoutAdaptive(:)-realphase(:)).^2)));
figure;
subplot(1,2,1); imshow(newrealphase,[]); title(strcat(['Noisy image, ',num2str(PSNRIn),'dB']));
subplot(1,2,2); imshow(IoutAdaptive,[]); title(strcat(['Clean Image by AST dictionary, ',num2str(PSNROut),'dB']));


figure;
I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
title('The AST dictionary');

%% Detail Retention
figure
plot(newrealphase(:,370),'--')
hold on
scatter(1:512,IoutAdaptive(:,370),'.','green')
axis tight
xlabel('Pixel'), ylabel('Phase');
legend( 'Original','Noise reduction');
%% Evaluation Indicators
% RMES
RMES = sqrt(mean((IoutAdaptive(:)-realphase(:)).^2))
% PSNR
PSNR = 20*log10(255/sqrt(mean((IoutAdaptive(:)-realphase(:)).^2)))
% ENL
enl = ENL(IoutAdaptive)
% SSI
ssl = SSI(newrealphase,IoutAdaptive)
% NC
nc = NC(realphase,IoutAdaptive)
% SSIM
sslm = SSIM(IoutAdaptive,realphase)
% FSIM
FSIM = FeatureSIM(realphase,IoutAdaptive)