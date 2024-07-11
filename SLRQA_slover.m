function [Xopt,info]=SLRQA_slover(params)
eta1=params.eta1;
eta2=params.eta2;
omega=params.omega;
loss=params.loss;
A=params.w0;
E=params.E;
D=params.x0;
rho=params.rho;
mu=params.mu;
lambda_y1=params.lambda_y1;
lambda_y2=params.lambda_y2;
W=params.w0;
iter=0;
B1=params.B1;
B2=params.B2;
lambda=params.lambda;
alpha=params.alpha;
%% 计算f0
% [~ S ~]=svd(A1);
% f0=trace(S);
f0_value=[];
fw_value=[];
ny1=[];
nA_W=[];
ny2=[];
% f0=value_A(A,mu,W,lambda_y1);
% f0_value=[f0];
% if params.type==1
%     params.eps=1e-3;
% else
%     params.eps=1e-3;
% end
while(1)
    %A1
    iter=iter+1;
    temp=1/(mu+params.Lk1)*((W-1/mu*lambda_y1)*mu+params.Lk1*A);
    %     [U S V]=qsvd(temp);
    %      S= shrinkage(S,1/mu);
    %      k=rank(S);
    %      U=U(:,1:k);
    %      V=V(:,1:k);
    %      sd=diag(S);
    %      sd=sd(1:k);
    %      S=diag(sd);
    %      A=U*S*V';
    A=DC(temp,1/(mu+params.Lk1),params);
    %      fw=value_W(A,mu,W,lambda_y1,lambda_y2,B1,B2,D,E,lambda,omega);
    %      temp1=qmultiplication(B1,W);
    %      temp2_w=qmultiplication(temp1,B2')+E-D;
    temp2_w=B1*W*B2'+E-D;
    %      temp2=proj(temp2_w+lambda_y2/mu,omega);
    %      temp2=qmultiplication(B1',temp2);
    %      temp2=W-(qmultiplication(temp2,B2)+W-A-lambda_y1/mu)/eta1;
    temp2=W-(B1'*proj(temp2_w+lambda_y2/mu,omega)*B2+W-A-lambda_y1/mu)/eta1;
    %     W=shrinkage(temp2,lambda/mu/eta1,params);
    gW = gradientW(W,params);
    W=shrinkageW(temp2,lambda/mu/eta1,gW,params);
    %     W=shrinkage(temp2,lambda/mu/eta1,params);



    nA_W=[nA_W norm(A-W,'fro')];
    %      nA_W=[nA_W qfro(A-W)];
    %      fw=value_W(A,mu,W,lambda_y1,lambda_y2,B1,B2,D,E,lambda,omega);
    %      fw_value=[fw_value fw];
    %
    %      temp2_e=E+B1*W*B2'-D;
    %      temp3=E-(proj(temp2_e,omega)+lambda_y2/mu)/eta2;
    %      E=shrinkage(temp3,alpha/mu/eta2);
    %      max(max(abs(E)));
    %
    t1=(A-W);
    lambda_y1=lambda_y1+mu*t1;
    t2=proj(B1*W*B2'+E-D,omega);
    lambda_y2=lambda_y2+mu*t2;
    %      lambda_y2=lambda_y2+q3multiplication(B1,W,B2');
    %      ny1=[ny1 qfro(lambda_y1)];
    %      ny2=[ny2 qfro(lambda_y2)];
    ny1=[ny1 norm(lambda_y1,'fro')];
    ny2=[ny2 norm(lambda_y2,'fro')];

    mu=mu*rho;
    nrm1(iter)=norm(norm(t1,'fro'))
    nrm2(iter)=norm(norm(t2,'fro'))
    max(norm(t1,'fro'),norm(t2,'fro'));

    if nrm1(iter)+nrm2(iter) < 1 &&  nrm1(iter)+nrm2(iter) >  5e-1
        params.epsilon = 1e-2;
    elseif nrm1(iter)+nrm2(iter) < 5e-1 && nrm1(iter)+nrm2(iter) > 1e-1
        params.epsilon = 1e-3;
    elseif nrm1(iter)+nrm2(iter) < 1e-1
        params.epsilon = 1e-4;
    end

    if (nrm1(iter)+nrm2(iter)<1e-2)%||iter>30
        %      if (qfro(A-W)<params.eps)%||iter>30
        break;
    end
    %      if (norm(A-W,'fro')<params.eps)%||iter>30
    %         break;
    %      end
    %omega
    %      temp4=abs(E)<params.Eeps;
    %      omega_new1=temp4.*omega;
    %      omega_new2=(1-temp4).*(1-omega);
    %      omega=omega_new1+omega_new2;
end
% A=rand(5,4);
% epsilon=0.1;
% output=shrinkage(A,epsilon)
info.iter=iter;
% Xopt=q3multiplication(B1,A,B2');
Xopt=B1*A*B2';
% tempX=(x(Xopt)-min(min(x(Xopt))))/(max(max(x(Xopt)))-min(min(x(Xopt))));
% Xopt=quaternion(zeros(params.m,params.n),tempX,tempX,tempX);
info.nrm1=nrm1;
info.nrm2=nrm2;
info.omega=omega;
info.E=E;
info.W=W;
info.A1=A;
info.ny1=ny1;
info.ny2=ny2;
info.lambda_y2=lambda_y2;
info.fw_value=fw_value;
info.f0_value=f0_value;
info.nA_W=nA_W;
info.omega=omega;
end


function output=shrinkage(A,epsilon,params)
if params.isq ==1
    As=s(A);
    signAs= sign(As);
    As=signAs.*max(abs(As)-epsilon,0);
    Ax=x(A);
    signAx= sign(Ax);
    Ax=signAx.*max(abs(Ax)-epsilon,0);
    Ay=y(A);
    signAy= sign(Ay);
    Ay=signAy.*max(abs(Ay)-epsilon,0);
    Az=z(A);
    signAz= sign(Az);
    Az=signAz.*max(abs(Az)-epsilon,0);
    output=quaternion(As,Ax,Ay,Az);
else
    output=sign(A).*max(abs(A)-epsilon,0);
end
%         Ax=x(A);
%         Ay=y(A);
%         Az=z(A);
%     sign_A=sign(A);
%     output=sign_A.*max(abs(A)-epsilon,0);
end

function output=shrinkageW(A,epsilon,W,params)
if params.isq ==1
    As=s(A);
    Ws=s(W);
    signAs= sign(As);
    As=signAs.*max(abs(As)-epsilon*Ws,0);
    Ax=x(A);
    Wx=x(W);
    signAx= sign(Ax);
    Ax=signAx.*max(abs(Ax)-epsilon*Wx,0);
    Ay=y(A);
    Wy=y(W);
    signAy= sign(Ay);
    Ay=signAy.*max(abs(Ay)-epsilon*Wy,0);
    Az=z(A);
    Wz=z(W);
    signAz= sign(Az);
    Az=signAz.*max(abs(Az)-epsilon*Wz,0);
    output=quaternion(As,Ax,Ay,Az);
else
    output=sign(A).*max(abs(A)-epsilon*W,0);
end
end

function output = gradientW(W,params)
if params.isq ==1
    Ws=s(W);
    Ws1 = abs(Ws) < params.epsilon;
    % gradientWs1 = params.epsilon*2 * abs(Ws);
    gradientWs1 =  1/params.epsilon*5 * abs(Ws);
    gradientWs2 = 1;
    % gradientWs1 = params.epsilon^(3/2)/4 * Ws;
    % gradientWs2 = 0.5 * (abs(Ws)+1e-6).^-0.5;
    gradientWs = Ws1.* gradientWs1 + (1 - Ws1).* gradientWs2;

    Wx=x(W);
    Wx1 = abs(Wx) < params.epsilon;
    % gradientWx1 = params.epsilon*2 * abs(Wx);
    gradientWx1 =  1/params.epsilon*5* abs(Wx);
    gradientWx2 = 1;
    % gradientWx1 = params.epsilon^(3/2)/4 * Wx;
    % gradientWx2 =0.5 * (abs(Wx)+1e-6).^-0.5;
    gradientWx = Wx1.* gradientWx1 + (1 - Wx1).* gradientWx2;

    Wy=y(W);
    Wy1 = abs(Wy) < params.epsilon;
    % gradientWy1 = params.epsilon*2 * abs(Wy);
    gradientWy1 =  1/params.epsilon*5 * abs(Wy);
    gradientWy2 = 1;
    % gradientWy1 = params.epsilon^(3/2)/4 * Wy;
    % gradientWy2 = 0.5 * (abs(Wy)+1e-6).^-0.5;
    gradientWy = Wy1.* gradientWy1 + (1 - Wy1).* gradientWy2;

    Wz=z(W);
    Wz1 = abs(Wz) < params.epsilon;
    % gradientWz1 = params.epsilon*2 * abs(Wz);
    gradientWz1 = 1/params.epsilon*5 * abs(Wz);
    gradientWz2 = 1;
    % gradientWz1 = params.epsilon^(3/2)/4 * Wz;
    % gradientWz2 = 0.5 * (abs(Wz)+1e-6).^-0.5;
    gradientWz = Wz1.* gradientWz1 + (1 - Wz1).* gradientWz2;

    output=quaternion(gradientWs,gradientWx,gradientWy,gradientWz);
else
    output = 0;
end
end

function y=proj(A,omega)
y=A.*double(omega);
end



% function f=value_A(A,mu,W,Y1)
% [~,S,~]=qsvd(A);
% f=sum(sum(S))+mu/2*norm(A-W+Y1/mu,'fro');
% end

% function f=value_W(A,mu,W,Y1,Y2,B1,B2,D,E,lambda,omega)
% f=lambda*sum(sum(abs(W)))+mu/2*norm(W-A-Y1/mu,'fro')+mu/2*norm(proj(B1*W*B2'+E-D,omega)+Y2/mu,'fro');
% end
%
function X=DC(Y,lambda,params)
simga=0;
if params.isq == 1
    [U,S,V]=qsvddiy(Y);
else
    [U,S,V]=svd(Y);
end
rk=rank(S);
sigmac=zeros(rk,1);
S2=zeros(size(S));
sigmay=diag(S(1:rk,1:rk));
sigmac=ones(rk,1);
while(1)
    [~,grad]=params.phi(sigmac,params.gamma,sigmay);
    sigma_new=max(sigmay-lambda*grad,0);
    if norm(sigma_new-sigmac)<1e-4
        break;
    end
    sigmac=sigma_new;
end
S2(1:rk,1:rk)=diag(sigma_new);
% end
% X=q3multiplication(U,S2,qtranspose(V));
X=U*S2*V';
end

