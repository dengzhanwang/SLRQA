function [eta2, stepsize, x2, f2, gradf2, slopex2, LSinfo, status] = LinesearchNonMonotonic(eta1, x1, fvs, gradf1, initslope, initstepsize, iter, fns, params)
% This function apply the backtracking algorithm to find an appropriate step size.
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
%     alpha : the coefficient in the Armijo condition
%     rho : the shrinking parameter
%     m   : the number of previous functions for nonmonotonic line search
%     lsmaxiter : the max numbers of iteration in the backtracking algorithm
%
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
%          1 means line search fails with inner iterations reaches its max number of iterations
%          2 means the numerical errors dominates the computation
% 
% By Wen Huang

    btiter=0;

        m=params.m;
        rho=params.rho;
        stepsize=initstepsize;           
        max=fvs(end);
        if(iter<10)
        for i=1:iter-1
            if(fvs(end-i)>max)
                max=fvs(end-i);
            end
        end
        else
             for i=1:9
            if(fvs(end-i)>max)
                max=fvs(end-i);
            end
                end
        end
        x2.main=x1.main+stepsize*eta1;
        while(fns.f(x2)>max+params.c1*stepsize*initslope)
        stepsize=rho*stepsize;
        btiter=btiter+1;
        x2.main=x1.main+stepsize*eta1;
        end
%         x2.main=x1+stepsize*eta1;
        [f2 x2]=fns.f(x2);
        eta2=-fns.Grad(x2);
        gradf2=-eta2;
        slopex2=gradf1'*eta2;
        status=0;
            if(btiter >= params.lsmaxiter)
        fprintf('warning: line search fails at iter:%d!\n', iter);
        status = 1;
    end
    if(norm((x1.main + eta2) - x1.main) == 0)
        status = 2;
    end
    LSinfo.lf = 1 + btiter; LSinfo.lgf = 1;
    
    
%     status = 1;
% LSinfo.lf = 0;
% LSinfo.lgf = 0;
% m = params.m;
% rho = params.rho;
% alpha = params.alpha;
% lsmaxiter = params.lsmaxiter;
% L = length(fvs);
% stepsize = initstepsize;
% eta2 = eta1 * stepsize;
% x2.main = x1.main + eta2;
% [f2 x2]= fns.f(x2);
% gradf2 = fns.Grad(x2);
% LSinfo.lf = LSinfo.lf + 1;
% LSinfo.lgf = LSinfo.lgf + 1;
% slopex2 = -gradf2' * gradf2;
% iter = 0;
% while 1
%     Max = -inf;
%     for i = 1:m
%         if L - i + 1 >= 1 && fvs(L - i + 1) >= Max
%             Max = fvs(L - i + 1);
%         end
%     end
%         if f2 <= Max + alpha * stepsize * initslope
%             gradf2 = fns.Grad(x2);
%             LSinfo.lgf = LSinfo.lgf + 1;
%             slopex2 = -gradf2' * gradf2;
%             status = 0;
%             return;
%         end
%     stepsize = stepsize * rho;
%     eta2 = eta1 * stepsize;
%     x2.main = x1.main + eta2;
%     [f2 x2] = fns.f(x2);
%     LSinfo.lf = LSinfo.lf + 1;
%     iter = iter + 1;
%     if iter > lsmaxiter
%         status = 1;
%         break;
%     end
% end
end

