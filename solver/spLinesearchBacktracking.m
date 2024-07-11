function [eta2, stepsize, x2, f2, gradf2, slopex2, LSinfo, status] = splinesearch(eta1, x1, fvs, gradf1, initslope, initstepsize, iter, fns, params)
% This function apply the simple backtracking algorithm to find an appropriate step size.
%
% INPUT:
% eta1 : search direction
% x1   : current iterate
% fvs  :  previous function values, the last one is the function value at the current iterate
% gradf1 : gradient at the current iterate
% initslope : the initial slope
% initstepsize : the initial step size
% fns : a struct that contains required function handles
%     fns.f(x) : return objective function value at x.
%     fns.Grad(x) : return the gradient of objection function at x.
% params: : a struct that contains parameters that are used
%     params.x0 : initial approximation of minimizer.
%     params.maxiter : the maximum number of iterations
%     params.verbose             : '0' means silence, '1' means output information of initial and final iterate. 
%                                  '2' means output information of every iterate.
% 
%     params.alpha:c1 

% OUTPUT:
% eta2 : step size * search direction
% stepsize : desired step size
% x2  : next iterate x + eta2
% f2  : function value at the next iterate
% gradf2 : gradient at the next iterate
% slopex2 : slope at the accepted step size
% LSinfo : Debug information in line search algorithm
%      LSinfo.lf : the number of function evaluations
%      LSinfo.lgf: the number of gradient evaluations
% status : 0 means success
%          1 means line search fails with inner iterations reaches its max
%          number of iterations达到最大迭代次数（失败）
%          2  意味着达到最大精度
% 
% By Wen Huang
status = 0;%
LSinfo.lf = 0;%函数计算次数
LSinfo.lgf = 0;%梯度计算次数
f1 = fvs(end); % function value at current iterate
f1= fns.f(x1);
stepsize1 = initstepsize;%1
rho=0.25;
x2.main=x1.main+stepsize1*eta1;
[f2,x2]=fns.f(x2);
btiter=1;%函数已使用一次
    status = 0; % line search success by default
    while(f2>f1+1e-4*stepsize1*initslope&&btiter < params.lsmaxiter)
        stepsize1=stepsize1*rho;
        x2.main=x1.main+stepsize1*eta1;
        [f2,x2]=fns.f(x2);
        btiter = btiter + 1;%每计算一次函数加1
    end
    % 更新点
   
        
        
    
    
    if(btiter >= params.lsmaxiter) % line search fails
        fprintf('warning: line search fails at iter:%d!\n', iter);
        status = 1;
    end
    
    eta2 = -fns.Grad(x2);
    if(norm(x2.main - x1.main) == 0)
        eta2 = -fns.Grad(x2);
    end
    
    if(norm(x2.main - x1.main) <= 1e-8)
        status = 2;
    end
    
    LSinfo.lf = btiter;
    LSinfo.lgf = 1;
    stepsize=stepsize1;
    
    gradf2=fns.Grad(x2);
    slopex2 = eta1' * gradf2;
 
end

