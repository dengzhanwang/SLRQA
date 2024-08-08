clear;
clc;
format long;
% GROPTBase=pwd;
% fprintf('Adding GROPT paths from %s...\n',GROPTBase);
addpath(genpath('..\'));
%图片
fpath = fullfile('..\imagedata\','*.png');

im_dir  = dir(fpath);
im_num = length(im_dir);

% parameters for denoising
nSig = 50;    % STD of the noise image
params.strname = 'Schatte';
% params.strname = 'Laplace';
params = SLRQAdnRGB_params(nSig,params.strname);

% params.strname = 'gamma';
params.phi=@Schatten_gamma;
% params.phi = @Laplace;

% record all the results in each iteration
params.PSNR = zeros(im_num,1);
params.SSIM = zeros(im_num,1);
% params.err = zeros(im_num,100);
% 
% params.n1 = zeros(im_num,100);
% params.n2 = zeros(im_num,100);
% params.n3 = zeros(im_num,100);
Output_dir = 'Results\SLRQArgb\';
mkdir(Output_dir); 
filename  = [params.strname,'_nsig_',num2str(params.nSig)];
fn_txt = strcat( filename, '_SLRQA_Results.txt' ); 

fid = fopen( fullfile(Output_dir, fn_txt), 'a+');
dat = datestr(datetime);
fprintf(fid,strcat('\r\n\r\n',dat,'\r\n'));
fprintf(fid,['Method: ' params.strname, '|gamma =' num2str(params.gamma) '\n']);
% for i = 1:1
for i=1:im_num
    
    params.image = i;
    % params.imgname = im_dir(i).name;
    params.Output_dir = Output_dir;
    params.I = double( imread(fullfile('../imagedata/', im_dir(i).name)) );
    % figure; image(params.I./255);title('原图');
    [h, w, ch] = size(params.I);
    params.h = h;
    params.w = w;
    params.ch = ch;
    randn('seed',0);
    params.nim = params.I+params.nSig*randn(size(params.I));
    % figure; image(params.nim./255);title('噪声图像');
    fprintf('%s :\n',im_dir(i).name);
    fprintf(fid,'%s :\n',im_dir(i).name);
   
    PSNR = csnr( params.nim, params.I, 0, 0 );
    SSIM = cal_ssim( params.nim, params.I, 0, 0 );
    savepath = [Output_dir,'nSig' num2str(params.nSig) '_' im_dir(i).name];
    imwrite(params.nim/255,savepath);
    fprintf('The initial value of PSNR = %2.2f, SSIM = %2.2f \n', PSNR,SSIM);
    fprintf(fid,'The initial value of PSNR = %2.2f, SSIM = %2.2f \n', PSNR,SSIM);
 
    rI = zeros(size(params.nim));
    %SLRQA不加非局部自相似性
    for k = 1:3
    [T,params,info] = ADMMRGB(params.nim(:,:,k),params);
    rI(:,:,k) = T;
    end
    rI(rI>255) = 255;
    rI(rI<0) = 0;
    PSNR = csnr( rI, params.I, 0, 0 );
    SSIM = cal_ssim( rI, params.I, 0, 0 );
    fprintf('PSNR = %2.2f, SSIM = %2.2f \n', PSNR,SSIM);
    savepath = [Output_dir, 'nSig' num2str(params.nSig) '_' params.strname '_' im_dir(i).name];
    imwrite(rI./255,savepath);
    fprintf(fid,'PSNR = %2.2f, SSIM = %2.2f \n', PSNR, SSIM);
    params.PSNR(i,1) = PSNR;
    params.SSIM(i,1) = SSIM;
    % params.err(i,1:info.iter) = info.nrm;
    % params.n1(i,1:info.iter) = info.n1;
    % params.n2(i,1:info.iter) = info.n2;
    % params.n3(i,1:info.iter) = info.n3;

    
end
apsnr = params.PSNR;
bssim = params.SSIM;
% err = params.err;
% n1 = params.n1;
% n2 = params.n2;
% n3 = params.n3;
fname = [Output_dir, 'nSig' num2str(params.nSig) '_' params.strname '_Psnr_Ssim_err.mat'];
save(fname,"bssim","apsnr");
% save(fname,'apsnr','bssim','err','n1','n2','n3');
fclose all;

function [y,grad]=Schatten_gamma(x,gamma)
    y=x.^gamma;
    grad=gamma*x.^(gamma-1);
end

function [y,grad] = Laplace(x,gamma)
y = 1-1./(exp(x./gamma));
grad = (1-y)./gamma;
end

