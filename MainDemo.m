%====================================================
% demo
%==================================================
% If you find our code helpful in your resarch or work, please cite our paper.
% Shaobo Li , Jianhu Zhao ,Hongmei Zhang, Zijun Bi1, SiHeng Qu , "A Non-local Low Rank Algorithm for Subbottom Profile Sonar Image Denoising "

%%
baseDirectory = fileparts(mfilename('fullpath'));
addpath(genpath(baseDirectory));
%%
clear all;
%%
parpool('local',6);
%%
% load test image
x=double((imread('simulated4_2_noise.bmp')));
%%
%Filter
Id=x;
Ivessel = Id;
Id = double((Id));
Id = 255 - double(Id);
Id = medfilt2(Id, [3,3]); 
h = fspecial('gaussian',[3,3], 2);
Id = imfilter(Id, h);
sigmas = [7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14,14.5,15, 20,30];
[Ivessel0, Direction, whatScale]= FrangiFilter2D(Id,sigmas);%, Direction
Ivessel = Direction;
%%
Id = 255 - double(Id);
Id = medfilt2(Id, [3,3]); 
h = fspecial('gaussian',[3,3], 2);
Id = imfilter(Id, h);
sigmas = [6.5, 10, 30,40,50];
[Ivessel1, Direction2, whatScale2]= FrangiFilter2D(Id, sigmas);%, Direction
Ivessel2 = Direction2;
Ivessel(find(Ivessel0<Ivessel1)) = Ivessel2(find(Ivessel0<Ivessel1));
%%
%--------parameters-------------
l = 29 ;
s = 0.8;
f=3; 
t = 36 ;
sigma=25;
%--------------------------------
[Output]=FNLM(x,f,t,sigma, Ivessel, l, s);
figure, imshow(Output,[])
figure, imshow(x,[])
%=======================================
Clean=double((imread('simulated4_2_clean.bmp')));
snr=psnr(Output,Clean);
[mssim, ssim_map] = ssim(Output,Clean);
Mssim = mssim;
%%
%close the pool
delete(gcp);
