% =========================================================================
% Test code for digital holographic microscopy noise reduction method
% =========================================================================

%% Clear residual data
close all;
clear;
clc;

%% Read holographic interference images
holography = double(imread('USAF1951.bmp'));
holography = holography(:,:,1);
[M,N] = size(holography);
figure(),imshow(holography,[]);title('Interference Map');

%% Fourier transform
FFTholography0 = fftshift(fft2(holography));
FFTholography = abs(FFTholography0);
figure(),imshow(FFTholography,[0,max(max(FFTholography))/50]);
title('log_Fourier spectrum');

%% Digital holographic off-axis reconstruction
[M,N] = size(holography)
FFTholography1 = zeros(M,N);
gaussfilter = zeros(M,N);

   for i = 1:M
       for k = 1:N
           distance = sqrt((k-1434).^2+(i-952).^2);
           gaussfilter(i,k) = exp(-(distance.*distance)./(2*20*50));
   
       end
end

FFTholography1=FFTholography0.*gaussfilter;
FFTholography2=ifft2(fftshift(FFTholography1));

[y,x] = size(FFTholography2);
deta = 3.45*10^-6;
d = 0*10^-3;
lx0 = x*deta;
ly0 = y*deta;
lamda = 532*10^-9;
k = 2*pi/lamda;
umax = x/(2*lx0);
vmax = y/(2*ly0);
u = linspace(-umax,umax,x);
v = linspace(-vmax,vmax,y);
[u,v] = meshgrid(u,v);
H = exp(j*k*d*(1-lamda*lamda*((u).^2+(v).^2)/2));
fftFFTholography2 = fftshift(fft2(FFTholography2));
fftFFTholography2H = fftFFTholography2.*H;
U = ifft2(fftshift(fftFFTholography2H));
I = U.*conj(U);
I = abs(U);
figure(),imshow(I,[]);
title('Reconstructed images')

%% Screenshot image
phase = atan2(imag(U),real(U));
phase = phase(1349:1860,969:1480);
figure(),imshow(phase,[]);
title('Screenshot image')

%% Phase distortion compensation
expphase = exp(j*phase);
fftphase = fftshift(fft2(expphase));
figure(),imshow(abs(fftphase),[]);
title('Center point')

[MM,NN] = size(phase);
xx = -NN/2:1:NN/2-1;
yy = -MM/2:1:MM/2-1;
[xx,yy] = meshgrid(xx,yy);
 
p10 = 0.00001;
p01 = 0.00001;
nhqx = p10/2.*xx.*xx+p01/2.*yy.*yy;
nhxwhs = exp(j*nhqx);
nhxw = atan2(imag(nhxwhs),real(nhxwhs));

phase3 = exp(j*phase+j*2*pi*((NN/2-299)/NN).*xx+j*2*pi*((MM/2-237.1)/MM).*yy-j*nhxw);
realphase = atan2(imag(phase3),real(phase3));

% % figure(),imshow(realphase,[]);
% % figure(),surf(realphase);
% % shading interp;

%% Adjustment of singularities
  [m,n]=size(realphase);
  for i=1:m
      for j=1:n
          if  realphase(i,j)>2
              realphase(i,j)= realphase(i,j)-2*pi;
          end
      end
  end

  realphase = realphase+2*pi;
  


%% show results
  figure(),imshow(realphase,[]);
  figure(),surf(realphase);
  shading interp;
  
%   save date1 realphase;
  
% % load('date.mat')
% % figure,plot(realphase(:,round(NN/2)+1))
 