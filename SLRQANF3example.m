clear;
clc;

fprintf('Adding GROPT paths from %s...\n',pwd);
addpath(genpath(pwd));
%图片
address=['./imagedata/'];
ssimq=zeros(8,8);
ssimsp=zeros(8,8);
psnrq=zeros(8,8);
psnrsp=zeros(8,8);
infomat=cell(8,8);
rng('default')

% ii=1;

%% 纯四元化(初始化)
% for ii = [4 6 7 8]
% dataadress = strcat('./imagedata/image',num2str(ii),'.png');
% A=cast(imreadq(dataadress), 'double') ./ 255;
% params.m=size(A,1);
% params.n=size(A,2);
% 
% A_q=quaternion(zeros(params.m,params.n),x(A),y(A),z(A));%转化为四元组向量
% 
% [U,S,V] = qsvd(A);
% rk=rank(S);
% 
% compress_ratio =0.06;
% params.rk_opt=floor(rk*compress_ratio);
% params.U_r=U(:,1:params.rk_opt);
% params.V_r=V(:,1:params.rk_opt);
% S_r=zeros(params.rk_opt);
% for jj=1:params.rk_opt
%    S_r(jj,jj)=S(jj,jj); 
% end
% params.S_r=S_r;
% A0=params.U_r*S_r*params.V_r';
% params.x0=A0;
% 
% params.gap=5;
% params.ratio = 1/4;
% params.FirstIterInitStepSize=1;
% params.lsmaxiter=100000;
% % params.linesearch = @LinesearchBacktracking;
% params.initstepsize = @InitStepSizeBB;
% params.lsmaxiter=10000;
% params.KeepNum=10;
% params.OutputGap=10;
% params.step=1;
% params.minstepsize=eps;%Wolf步长范围
% params.maxstepsize=1e7;
% params.w1=0.001;
% params.w2=1000;%BB步长范围
% params.alpha=0.0001;%c1
% params.beta=0.999;
% params.tau1=0.02;%Wolfe算法参数tau1
% params.tau2=0.98;%Wolfe算法参数tau2
% params.verbose=2;
% params.compress_ratio=0.06;
% [A info]=lowrank_pure(A,params);
% A = qtoc(A);
% B = zeros(512,512,4);
% B(:,:,2) = A(:,:,1);
% B(:,:,3) = A(:,:,2);
% B(:,:,4) = A(:,:,3);
% A = B;
% adress2=strcat('data',num2str(ii),'lr30');
% save(strcat(address,adress2),'A');
% 
% 
% end
%%
close all
draw_option = 0;
%load A_t2.mat;
lowrank = 0
if lowrank == 0
load SLRQA3info.mat
else
load SLRQA3_lowrankinfo60.mat
end

for jj= [3 5 7 8]
    corrupt=jj/10;
    seed=2099;
    rng(seed);
    for ii= 1:8
        if lowrank == 0
            dataadress = strcat('./imagedata/image',num2str(ii),'.png');
            A=cast(imreadq(dataadress), 'double') ./ 255;
            A=quaternion(zeros(512,512),x(A),y(A),z(A));
        else
            % dataaddress=strcat(address,'image',num2str(ii),'.png');
            % A=cast(imreadq('E:\毕业论文\代码\imagedata\Lenna.png'), 'double') ./ 255;
            dataaddress=strcat(address,'data',num2str(ii),'lr30');
            A=load(dataaddress);
            A=A.A;
            A=quaternion(zeros(512,512),A(:,:,2),A(:,:,3),A(:,:,4));
        end

        params.m=size(A,1);
        params.n=size(A,2);
        % set(gcf,'unit','centimeters','position',[1,2,10,10]);
        % dataaddress=strcat(address,'image',num2str(ii),'.png');
        % A=cast(imreadq(dataaddress), 'double') ./ 255;
        % [U S V]=qsvd(A);
        % S2=zeros(size(S));
        % for i=1:30
        %  S2(i,i)=S(i,i);
        % end
        % A=U*S2*V';%低秩化原矩阵
        %A的中间0.4-0.6比例将被挖掉
        % m_low=ceil(params.m*0.45);
        % m_high=ceil(params.m*0.55);
        % n_low=ceil(params.n*0.45);
        % n_high=ceil(params.n*0.55);
        % A(m_low:m_high,n_low:n_high)=0;
        %  real_support=ones(size(A));
        % real_support=A~=0;
        params.Eeps=0.1;
        % imageq(A);
        % image(A);
        [U S V]=qsvd(A);
        params.Anorm=S(1,1);
        A=A/params.Anorm;%单位化矩阵A
        % set(gcf,'unit','centimeters','position',[1,2,20,20]);
        B=rand(params.m,params.n);
        omega_real=B>corrupt;
        sum(sum(omega_real))/size(omega_real,1)/size(omega_real,2);
        % A2=quaternion(zeros(params.m,params.n),xA,xA,xA);
        % image(A2);
        % omega_real=real_support.*omega_real;
        dbstop if error
        %  input_image=A+(1-omega_real).*randq(size(A))/10;%图像加噪点
        input_image=A.*double(omega_real);
        % A1=input_image(:,:,2);
        % A2=input_image(:,:,3);
        % A3=input_image(:,:,4);
        A1=x(input_image);
        A2=y(input_image);
        A3=z(input_image);

        temp=input_image*params.Anorm;
        % imageq(temp)
        % qdraw(temp,[200 300],[200 300]);
        % imagename=strcat('image',num2str(ii),'cor.png');
        % axis off
        % exportgraphics(gcf,imagename,'Resolution',1000);
        % image(input_image)
        params.x0=input_image;
        params.omega=omega_real;
        params.loss=1-omega_real;
        params.eps=1e-4;
        params.epsilon = 1e-4;
        params.lambda=0.1;
        params.alpha=0.85;
        params.eta1=3;
        params.rho=1.1;
        params.Lk1 = 0.01;
        params.Lk2 = 0.01;
        params.gamma=0.7;
        params.mu=5;
        params.verbose=2;%意味着有输出
        params.eta2=3;
        params.phi=@WSchatten_gamma;
        params.lambda_y1=quaternion(zeros(params.m,params.n),zeros(params.m,params.n),zeros(params.m,params.n),zeros(params.m,params.n));
        params.lambda_y2=quaternion(zeros(params.m,params.n),zeros(params.m,params.n),zeros(params.m,params.n),zeros(params.m,params.n));
        % params.lambda_y1=zeros(params.m,params.n,4);
        % params.lambda_y2=zeros(params.m,params.n,4);
        params.B1=creatematrix(params.m);%已转置
        params.B2=creatematrix(params.n);
        % temp1=qmultiplication(params.B1',params.x0);
        % params.w0=qmultiplication(temp1,params.B2);
        params.w0=params.B1'*params.x0*params.B2;
        params.E=zeros(params.m,params.n);
        % params.E=1-real_support;
        % D=abs(params.w0)>0.2;
        % sum(sum(D))/size(D,1)/size(D,2)
        % C=params.w0.*D;
        % sparseA=params.B1'*C*params.B2;
        % figure(2)
        % image(sparseA)
        params.type=1;
        Ac = qtoc(A*params.Anorm);
        params.isq = 1;

%         tic
%         [Xopt,info]=SLRQA_slover(params);
%         toc
%         infomat{ii,jj}=info;
%         Xoptc=qtoc(Xopt*params.Anorm);
%         psnrq(ii,jj) = psnr(Xoptc,Ac);
%         ssimq(ii,jj)=ssim(Ac,Xoptc);

%         if draw_option == 1
%             close all
%             figure(2)
%             set(gcf,'Position',[300 300 500 500]);
%             set(gca,'Position',[0 0 1 1]);
%             imagesc(Xoptc)
%             qdraw(Xopt*params.Anorm,[200 300],[200 300]);
%             axis off
%             adress = strcat('./result/inpaint/SLRQA3/imageq',num2str(ii),'_',num2str(corrupt*100),'.png');
%             saveas(gcf, adress)
%             close all
%         end
% 
%         imageq(Xopt*params.Anorm)
%         tempc = qtoc(temp);
%         qdraw(Xopt*params.Anorm,[200 300],[200 300]);
%         imagename=strcat('image',num2str(ii),'with',num2str(corrupt*100),'SLRQA3q.png');
%         exportgraphics(gcf,imagename,'Resolution',1000);
%         axis off
        params.type=2;
        params.lambda_y2=zeros(params.m,params.n);
        params.lambda_y1=zeros(params.m,params.n);
        params.x0=A1;
        params.w0=A1;
        params.isq = 0;
        tic
        [Xopt1,info]=SLRQA_slover(params);
        params.x0=A2;
        params.w0=A2;

        [Xopt2,info]=SLRQA_slover(params);
        params.x0=A3;
        params.w0=A3;
        [Xopt3,info]=SLRQA_slover(params);
        toc
        %
        Xopt_sp =quaternion(zeros(params.m,params.n),Xopt1, Xopt2, Xopt3);
        Xopt_sp=Xopt_sp*params.Anorm;
        Xopt_spc=qtoc(Xopt_sp);
        % imageq(Xopt_sp);

%         if draw_option == 1
            % figure(3)
            % set(gcf,'Position',[300 300 500 500]);
            % set(gca,'Position',[0 0 1 1]);
            % qdraw(Xopt_sp,[200 300],[200 300]);
            % imagesc(Xopt_sp)
            % axis off
            % adress = strcat('./result/inpaint/SLRQA/imagesp',num2str(ii),'_',num2str(corrupt*100),'.png');
            % saveas(gcf, adress)

            % close all
%         end

        psnrsp(ii,jj) = psnr(Xopt_spc,Ac);
        ssimsp(ii,jj) = ssim(Ac,Xopt_spc);
        % imagename=strcat('image',num2str(ii),'with',num2str(corrupt*100),'LADMsp.png');
        % exportgraphics(gcf,imagename,'Resolution',1000);
        axis off
%         [U S V]=qsvddiy(Xopt);
%         rank(S)
%         diff_maxXopt=max(max(abs(A-Xopt)));
%         % diff_maxXoptsp=max(max(abs(A-Xopt_sp)))
%         temp=Xopt*params.Anorm;
%         Output_image=quaternion(zeros(params.m,params.n),x(temp),y(temp),z(temp));
%         % imageq(Output_image)
%         norm(Xopt-A ,'fro')/norm(A,'fro');
        % norm(Xopt_sp-A ,'fro')/norm(A,'fro')
        % psnrsp(ii)=psnr(A,Xopt_sp);

        % Xopt_sp=qtoc(Xopt_sp*params.Anorm);

        
        if lowrank == 0
            save(strcat('./result/inpaint/SLRQA3/SLRQA3info.mat'),'psnrq','ssimq','psnrsp','ssimsp',"infomat")
        else
            save(strcat('./result/inpaint/SLRQA3/SLRQA3_lowrankinfo60.mat'),'psnrq','ssimq','psnrsp','ssimsp',"infomat")
        end
    end
end
1;
if lowrank == 0
    save(strcat('./result/inpaint/SLRQA3/SLRQA3info.mat'),'psnrq','ssimq','psnrsp','ssimsp',"infomat")
else
    save(strcat('./result/inpaint/SLRQA3/SLRQA3_lowrankinfo60.mat'),'psnrq','ssimq','psnrsp','ssimsp',"infomat")
end

index = [3 5 7];
close all
figure(1)
plot(psnrq(:,3),'b-','linewidth',2)
hold on
plot(psnrq(:,5),'b-.','linewidth',2)
plot(psnrq(:,7),'b:','linewidth',2)
plot(psnrsp(:,3),'r-','linewidth',2)
plot(psnrsp(:,5),'r-.','linewidth',2)
plot(psnrsp(:,7),'r:','linewidth',2)
axis([1 8 15 45])
legend('SLRQA(Q)-30','SLRQA(Q)-50','SLRQA(Q)-70','SLRQA(RGB)-30','SLRQA(RGB)-50','SLRQA(RGB)-70','Location','southeast')
xlabel('Image Index','FontSize',15)
ylabel('PSNR','FontSize',15)
% save(strcat('./result/inpaint/SLRQA/SLRQAinfo.mat'),'psnrq','ssimq','psnrsp','ssimsp')
saveas(gcf,'./result/inpaint/SLRQA3/comppsnr.png')
figure(2)
plot(ssimq(:,3),'b-','linewidth',2)
hold on
plot(ssimq(:,5),'b-.','linewidth',2)
plot(ssimq(:,7),'b:','linewidth',2)
plot(ssimsp(:,3),'r-','linewidth',2)
plot(ssimsp(:,5),'r-.','linewidth',2)
plot(ssimsp(:,7),'r:','linewidth',2)
axis([1 8 0.75 1])
legend('SLRQA(Q)-30','SLRQA(Q)-50','SLRQA(Q)-70','SLRQA(RGB)-30','SLRQA(RGB)-50','SLRQA(RGB)-70','Location','southeast')
xlabel('Image Index','FontSize',15)
ylabel('SSIM','FontSize',15)
saveas(gcf,'./result/inpaint/SLRQA3/compssim.png')
% aa=load("LADM50.mat");
% psnrq=aa.psnrq;
% psnrsp=aa.psnrsp;
% ssimq=aa.ssimq;
% ssimsp=aa.ssimsp;
% bb=load("SLRQA50.mat");
% bb.ssimq-aa.ssimq
% bb.psnrq-aa.psnrq
% function [y,grad]=Schatten_gamma(x,gamma)
% y=x.^gamma;
% grad=gamma*x.^(gamma-1);
% end
function [y,grad]=WSchatten_gamma(x,gamma,sigmay)
    n = length(x);
    w = sqrt(n)/100/sqrt(2)./(sigmay+0.1);
%     w = ones(size(n));
    y=w.*x.^gamma;
    grad=w.*gamma.*x.^(gamma-1);
end

function [y,grad]=Schatten_gamma(x,gamma)
y=x.^gamma;
grad=gamma*x.^(gamma-1);
end

function [y,grad]=Laplace(x,gamma)
    y=1-exp(-x./gamma);
    grad=1./gamma*exp(-x./gamma);
end

