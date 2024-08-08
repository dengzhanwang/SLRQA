function [rI,params] = SLRQA(nI, I ,params)
A = 255;
rI = nI;
   
params = SearchNeighborIndex( params );
for iter = 1 : params.Iter
    params.iter = iter;
    % iterative regularization
    rI = rI + params.delta * (nI - rI);
    rI = rI./A;
    % image to patch
    CurPat = Image2Patch( rI, params );

    if (mod(iter-1, params.Innerloop) == 0)
        params.nlsp = params.nlsp - 10;  % Lower Noise level, less NL patches
        NL_mat  =  Block_Matching(CurPat, params);% Caculate Non-local similar patches for each
    end
    % NL_mat  =  Block_Matching(CurPat, params);
    
    Y_hat = zeros(size(CurPat));
    W_hat = zeros(size(CurPat));
    for i = 1 : params.lenrc   % For each keypatch group
        Y = CurPat(:, NL_mat(1:params.nlsp,i)); % Non-local similar patches to the keypatch
   
        mY = repmat(mean( Y, 2 ),1,params.nlsp);
        Y = Y-mY;
        %转换为四元数矩阵
        Y = quaternion(zeros(params.ps2,params.nlsp),Y(1:params.ps2,:),Y(params.ps2+1:params.ps2*2,:), Y(params.ps2*2+1:params.ps2*3,:));
        [params.m,params.n] = size(Y);
        
        X = ADMM( Y, params); 
        % X = wADMM(Y,params);
        
        X = [X.x;X.y;X.z]+mY;
        
        Y_hat(:,NL_mat(1:params.nlsp,i)) = Y_hat(:,NL_mat(1:params.nlsp,i))+X;
        W_hat(:,NL_mat(1:params.nlsp,i)) = W_hat(:,NL_mat(1:params.nlsp,i))+ones(params.ps2ch, params.nlsp);
    end
    
    rI = PGs2Image(Y_hat, W_hat, params);
    rI = rI.*A;
    rI(rI>255)=255;
    rI(rI<0)=0;
    figure; image(rI./255);
    savepath = [params.Output_dir, 'nSig' num2str(params.nSig) '_' params.strname '_gamma=' num2str(params.gamma) '_' num2str(iter) '_' params.imgname];
    imwrite(rI./255,savepath);
    PSNR = csnr( I, rI, 0, 0 );
    SSIM = cal_ssim( I, rI, 0, 0 );

    fprintf(params.fileid,'Iter = %6.2f, PSNR = %2.3f, SSIM = %2.3f \n', iter, PSNR, SSIM);
    fprintf( 'Iter = %6.2f, PSNR = %2.3f, SSIM = %2.3f \n', iter, PSNR, SSIM );
    params.PSNR(iter, params.image) = PSNR;
    params.SSIM(iter, params.image) = SSIM;
   
end

end