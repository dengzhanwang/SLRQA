function [rI,params,info] = ADMMRGB(nI  ,params)
% Initializing optimization variables

l1 = params.l1;
l2 = params.l2;
tau = params.nSig;
rho = params.rho;   
mu = params.mu;   %惩罚因子
lambda=params.lambda;  %稀疏惩罚因子，1范数的系数

A = 255;
nI = nI./A;

X = nI;  %初始值
Y = X;
W = mirt_dctn(Y);
lambda_y=zeros(params.h,params.w);%拉格朗日乘子

iter = 0;

while(1)
    iter=iter+1;
    %更新X
    T = mirt_idctn(W-lambda_y/mu);
    % T = quaternion(zeros(params.h,params.w),mirt_idctn(T.x),mirt_idctn(T.y),mirt_idctn(T.z));
    d = tau^2/(tau^2.*(mu+l1)+1);
    temp = d * (Y./(tau^2)+mu*T+l1.*X);
    Xnew = DC(temp,d,params);
    n1(iter) = norm(X-Xnew,'fro');
    X = Xnew;
    %更新W
    t = mirt_dctn(X);
    % t = quaternion(zeros(params.h,params.w),mirt_dctn(X.x),mirt_dctn(X.y),mirt_dctn(X.z));
    temp2 = (mu.*(t)+lambda_y+l2*W)./(mu+l2);
    % Wnew = shrinkage(temp2,lambda/(mu+l2));
    if params.epsilon > 1e-3
        params.epsilon = params.epsilon./2;
    end

    gW = gradientW(W,params);
    Wnew = shrinkageW(temp2,lambda/(mu+l2),gW,params);
    
   
    n2(iter) = norm(W-Wnew,'fro');
    W = Wnew;
    
    t1=t-W;
    lambda_y=lambda_y+mu*t1;
    mu = mu*rho;
    n3(iter) = norm(t1,'fro');

    nrm(iter) = n1(iter)+n2(iter)+n3(iter);

    if (nrm(iter)<params.eps)%||iter>30
        break;
    end

end
info.iter = iter;
info.nrm = nrm;
info.n1 = n1;
info.n2 = n2;
info.n3 = n3;

rI = X.*A;

end

function X=DC(Y,lambda,params)
[U,S,V]=qsvddiy(Y);
rk=rank(S);
% sigmac=zeros(rk,1);
S2=zeros(size(S));
sigmay=diag(S(1:rk,1:rk));
sigmac=ones(rk,1);
while(1)
    [~,grad]=params.phi(sigmac,params.gamma);
    sigma_new=max(sigmay-lambda*grad,0);
          if norm(sigma_new-sigmac)<1e-4
               break;
          end
     sigmac=sigma_new;
end
S2(1:rk,1:rk)=diag(sigma_new);
X=U*S2*V';
end
