function [xopt, info] = RSD(fns, params)
% Solver of the steepest descent method
%
% INPUT:
% fns : a struct that contains required function handles
%     fns.f(x) : return objective function value at x.
%     fns.Grad(x) : return the gradient of objection function at x.
%
% params : a struct that contains parameters that are used in the solver.
%     params.x0 : initial approximation of minimizer.
%     params.maxiter : the maximum number of iterations
%     params.verbose             : '0' means silence, '1' means output information of initial and final iterate.
%                                  '2' means output information of every iterate.
%
%     params.ratio:rho
% OUTPUT:
% xopt : the last iterate
% info : informtion generated during the algorithms
%
%
% By Wen Huang
tic;
x1 = params.x0;%取初始值
fvs = fns.f(x1);%此时函数值
nabla_1 = fns.Grad(x1);%下降方向
[U, S, V]=svd(x1);
r=params.rk_opt;
S_r=zeros(r);
for jj=1:r
    S_r(jj,jj)=S(jj,jj);
end
params.U_r=U(:,1:r);
params.V_r=V(:,1:r);
params.S_r=S_r;
[params,gradf1]=proj(x1,nabla_1,params);% Riemannian gradient
xopt=[];
eta1=-gradf1;
% alpha=params.alpha;%c1
initslope = -prod(gradf1,gradf1);%h’（0）
initstepsize = 1;%alpha
itermain = 0;%当前迭代次数
% initslope2 = initslope; % Stored for true initslope
fnumber = 0; % The number of times for computing function value
gradnumber = 0; % The number of times for computing gradient value
info.f = [];
info.gradf = [];
f1 = fns.f(x1); % function value
fvs = f1; ngfs = norm(gradf1,'fro'); lf = 1; lgf = 1;% information
times = toc; Xs{1} = x1; GFs{1} = gradf1; Etas{1} = eta1;sx1{1}=norm(s(x1),'fro');
status = 0;
iter=1;
OutputGap=params.gap;
flag1=1;
stop=0;
btiter=0;
flag2=0;
while(~stop && status == 0&& flag1==1)
    %initstepsize
    if flag2==0
      [initstepsize info]=params.initstepsize(iter, initslope, fvs, Xs, GFs, Etas, fns, params);
    else
        initstepsize=1;
        flag2=0;
    end
    
    % Line search algorithm
% [eta2, stepsize, x2, f2, gradf2, slopex2, LSinfo, status] = params.linesearch(eta1, x1, fvs, gradf1, initslope, initstepsize, iter, fns, params);


%  [eta2, stepsize, x2, f2, gradf2, slopex2, LSinfo, status] = params.linesearch(eta1, x1, fvs, gradf1, initslope, initstepsize, iter, fns, params);
%      [Q_U_r, S_r ,Q_V_r , x2]=retraction(stepsize,params);%inter
  [Q_U_r, S_r ,Q_V_r,x2]=retraction(initstepsize,params);%sp

    f2=fns.f(x2);%sp
    lf=lf+1;%sp

    stepsize1=initstepsize;
    while(f2>f1+1e-4*stepsize1*initslope && btiter < params.lsmaxiter)%sp
        btiter=btiter+1;
        stepsize1=stepsize1*params.ratio;
        [Q_U_r, S_r ,Q_V_r , x2]=retraction(stepsize1,params);
        f2=fns.f(x2);
        lf=lf+1;
    end
    params.U_r=Q_U_r;
    params.V_r=Q_V_r;
    params.S_r=S_r;
    if(btiter >= params.lsmaxiter) % line search fails
        fprintf('warning: line search fails at iter:%d!\n', iter);
        status = 1;
    end

%     if(norm(x2 - x1,'fro') <= 1e-7)
%         status = 2;
%     end
    LSinfo.lf = btiter;
    LSinfo.lgf = 1;
    stepsize=stepsize1;
    
    % finish line search
    %         xopt=[xopt norm(x2.main-params.acc)];
    if(norm(x1-x2,'fro')<params.eps)
        flag1=0;
    end
%     if(norm(gradf1,'fro')<params.eps)
%         flag1=0;
%     end
    x1 = x2; f1 = f2; 
     if (norm(gradf1,'fro')<params.etanorm)
        if  (norm(s(x2))<params.eps_stop)
            flag1=0;
        end
        params.lambda=params.lambda+params.mu*s(x2);
        params.mu=params.rho*params.mu;
        fns = ProbEucQuadratic2(params.A_q,params.lambda,params.mu);
        f1=fns.f(x1);
        flag2=1;
    end

    
    if(length(Xs) >= params.KeepNum)
        Xs(1) = []; GFs(1) = []; Etas(1) = [];
    end
    %     stop = params.IsStopped( iter,Xs , fvs , ngfs , ngf0 , times , fns ,params );
    % Get ready for the next iteration
    
    
     nabla_1 =fns.Grad(x1);%sp
%     nabla_1=gradf2;%inter

    [params,gradf1]=proj(x1,nabla_1,params);% Riemannian gradient
    
      slopex2 =norm(gradf1*gradf1','fro');%sp
     
    ngf1 = norm(gradf1, 'fro');
    initslope=-slopex2;
    initstepsize=1;
    eta1=-gradf1;
    gradf2=gradf1;
    btiter=0;
    %output
    ngf2 = norm(gradf2, 'fro');
    fvs(end + 1) = f1;
    ngfs(end+1) = ngf2;
    times(end + 1) = toc;
    Xs{end + 1} = x1;
    GFs{end + 1} = gradf2;
    Etas{end + 1} = eta1;
    sx1{end+1}=norm(s(x1),'fro');
    iter = iter + 1;
    
    lgf=lgf+1;
    if(params.verbose >= 2) && (mod(iter, OutputGap) == 0)
        fprintf('iter:%d,f:%.8e,|gf|:%.3e,normsx:%.1e,snew:%.1e,t0:%.1e,tnew:%.1e,lf:%d,lgf:%d\n', iter, f1, ngf1, norm(s(x1),'fro'), slopex2, initstepsize, stepsize, sum(lf), sum(lgf));
    end
end
%     if(params.verbose >= 1)
%         fprintf('iter:%d,f:%.3e,|gf|:%.3e,|gf|/|gf0|:%.3e,lf:%d,lgf:%d,status:%d\n', iter, f1, ngf1, ngf1 / ngf0, sum(lf), sum(lgf), status);
%     end


info.iter = iter;
info.ngfs = ngfs;
info.fvs = fvs;
info.lf = lf;
info.status=status;
info.lgf = lgf;
info.Xs=Xs;
info.GFs=GFs;
info.sx1=sx1;
xopt=Xs{end};
end

function [params,grad]=proj(x,Z,params)
% [U, S, V]=svd(x);
% r=params.rk_opt;
% S_r=zeros(r);
% for jj=1:r
%     S_r(jj,jj)=S(jj,jj);
% end
% U_r=U(:,1:r);
% V_r=V(:,1:r);
U_r=params.U_r;
V_r=params.V_r;
S_r=params.S_r;
M=U_r'*Z*V_r;
U_p=Z*V_r-U_r*M;
V_p=Z'*U_r-V_r*M';
grad=U_r*M*V_r'+U_p*V_r'+U_r*V_p';
params.U_r=U_r;
params.V_r=V_r;
params.M=M;
params.U_p=U_p;
params.V_p=V_p;
params.S_r=S_r;
[Q_U,R_U]=qr([U_r U_p]);
[Q_V,R_V]=qr([V_r V_p]);

Q_U=Q_U(:,1:2*params.rk_opt);
R_U=R_U(1:2*params.rk_opt,:);
Q_V=Q_V(:,1:2*params.rk_opt);
R_V=R_V(1:2*params.rk_opt,:);
params.Q_U=Q_U;
params.Q_V=Q_V;
params.R_U=R_U;
params.R_V=R_V;
end

%%
function [Q_U_r, S_r ,Q_V_r ,Re]=retraction(alpha,params)%X+H
U_r=params.U_r;
V_r=params.V_r;
S_r=params.S_r;
U_p=params.U_p;
V_p=params.V_p;
M=params.M;
Q_U=params.Q_U;
Q_V=params.Q_V;
R_U=params.R_U;
R_V=params.R_V;
z=quaternion(zeros(size(S_r,1)),zeros(size(S_r,1)),zeros(size(S_r,1)),zeros(size(S_r,1)));
X=R_U*[S_r-alpha*M ,-alpha*eyeq(size(S_r,1));-alpha*eyeq(size(S_r,1)), z]*R_V';
[U_n,S_n,V_n]=svd(X);
r=params.rk_opt;
S_r=zeros(r);
for jj=1:r
    S_r(jj,jj)=S_n(jj,jj);
end
U_r=U_n(:,1:r);
V_r=V_n(:,1:r);

Q_U_r=Q_U*U_r;
Q_V_r=Q_V*V_r;

Re=Q_U_r*S_r*Q_V_r';
end

function a=prod(x,y)
% a=s(sum(sum(x.*y)));
a=sum(sum(s(x'*y)));
end
% a=Q_U*U_r;
% b=V_r'*Q_V';
% b=b';
