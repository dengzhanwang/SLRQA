clear;
clc;
format long;
GROPTBase=pwd;
fprintf('Adding GROPT paths from %s...\n',GROPTBase);

addpath(genpath('../'));

%图片
fpath = fullfile('../imagedata/','*.png');

im_dir  = dir(fpath);
im_num = length(im_dir);

% parameters for denoising
nSig = 10;    % STD of the noise image
% params.strname = 'Schatte';
% params.strname = 'Laplace';
params.strname = 'wSchatt';
params = SLRQAexperiment2_NSS_params(nSig,params.strname);
params = NSS_params(nSig,params);
% params = structAssign(params, temp);
params.phi = @Schatten;
% params.phi = @Laplace;
% params.phi = @wSchatten;

% record all the results in each iteration
params.PSNR = zeros(params.Iter+1, im_num);
params.SSIM = zeros(params.Iter+1, im_num);
Output_dir = 'Results\SLRQAexperiment2_NSS\';
params.Output_dir = Output_dir;
mkdir(Output_dir);

filename  = [params.strname,'_nsig_',num2str(params.nSig)];
fn_txt = strcat( filename, '_SLRQAexperiment2_NSS_Results.txt' ); 

fid = fopen( fullfile(Output_dir, fn_txt), 'a+');
dat = datestr(datetime);
fprintf(fid,strcat('\r\n\r\n',dat,'\r\n'));
fprintf(fid,['Method: ' params.strname, '|gamma =' num2str(params.gamma) '\n']);
fprintf(fid,'window = %2.1f, patchsize = %2.1f, step = %2.1f , \n', params.win, params.ps, params.step);
fprintf(fid,'lambda = %5.4f, mu = %4.2f, rho = %4.2f ,\n',params.lambda,params.mu,params.rho );
params.fileid = fid;

% for i=1:im_num
for i=2
    params.image = i;
    params.imgname = im_dir(i).name;
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
    % PSNR = psnr( params.nim, params.I);
    % SSIM = ssim( params.nim, params.I);
    params.PSNR(1,i) = PSNR;
    params.SSIM(1,i) = SSIM;

    savepath = [params.Output_dir, 'nSig' num2str(params.nSig) '_' im_dir(i).name ];
    imwrite(params.nim/255,savepath);
    fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
    fprintf(params.fileid,'The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
  
    [im_out,params] = SLRQA(params.nim,params.I,params);
   
    
    fprintf(params.fileid , 'The max PSNR is %2.4f , which in %2.1f th iteration.', max(params.PSNR(:,i)), find(params.PSNR(:,i)==max(params.PSNR(:,i))));
    fprintf(params.fileid , 'The max SSIM is %2.4f , which in %2.1f th iteration.', max(params.SSIM(:,i)), find(params.SSIM(:,i)==max(params.SSIM(:,i))));
    fprintf(params.fileid,'\n');
end
apsnr = params.PSNR;
bssim = params.SSIM;
fname = [params.Output_dir,'nSig' num2str(params.nSig) 'Psnr_Ssim.mat'];
save(fname,'apsnr','bssim');
fclose all;

function [y,grad]=Schatten(x,gamma)
    y=x.^gamma;
    grad=gamma*x.^(gamma-1);
end

function [y,grad] = Laplace(x,gamma)
y = 1-1./(exp(x./gamma));
grad = (1-y)./gamma;
end

function [y,grad] = wSchatten(x,gamma,w)
y = w.*(x.^gamma);
grad = w.*(gamma*x.^(gamma-1));
end


