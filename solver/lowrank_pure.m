function [Xopt info]=lowrank_pure(A,params)
compress_ratio=params.compress_ratio;
A_q=quaternion(zeros(params.m,params.n),x(A),y(A),z(A));%转化为四元组向量
[U,S,V] = svd(A);
rk=rank(S);
params.rk_opt=floor(rk*compress_ratio);
params.U_r=U(:,1:params.rk_opt);
params.V_r=V(:,1:params.rk_opt);
S_r=zeros(params.rk_opt);
for jj=1:params.rk_opt
   S_r(jj,jj)=S(jj,jj); 
end
params.S_r=S_r;
A0=params.U_r*S_r*params.V_r';
params.x0=A0;
%确保秩为整数
params.mu=20000;
params.lambda=quaternion(zeros(params.m,params.n),zeros(params.m,params.n),zeros(params.m,params.n),zeros(params.m,params.n));
params.rho=1.1;
params.eps_stop=1e-4;
params.etanorm=1;
params.eps=1e-4;
params.maxlambda=20000;
params.A_q=A_q;
fns = ProbEucQuadratic2(params.A_q,params.lambda,params.mu);
%  fns = ProbEucQuadratic(params.A_q,params.mu);
% fns.f(0)
% fns.Grad(0)
% tic;
[Xopt, info] = RSD(fns, params);
info.rank=params.rk_opt;
end