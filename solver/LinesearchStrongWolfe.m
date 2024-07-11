function [nabla_2, stepsize, x2, f2, gradf2, slopex2, LSinfo, status] = LinesearchStrongWolfe(eta1, x1, fvs, gradf1, initslope, initstepsize, iter, fns, params)
% This function apply the backtracking algorithm to find an appropriate step size.
%
% INPUT:
% eta1 : search direction
% x1   : current iterate
% fvs  : previous function values, the last one is the value at the current iterate
% gradf1 : gradient at the current iterate
% initslope : the initial slope
% initstepsize : the initial step size
% fns : a struct that contains required function handles
%     fns.f(x) : return objective function value at x.
%     fns.Grad(x) : return the gradient of objection function at x.
% params: : a struct that contains parameters that are used
%     alpha: the coefficient in the Armijo condition c1
%     beta : the coefficient in the curvature condition c2
%     lsmaxiter : the max numbers of inner iteration
%     minstepsize : minimum step size tau1
%     maxstepsize : maximum step size tau2
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
% status    : 4 means some other serious errors happened,
%             3 means step size excesses the maximum step size
%             2 means step size reaches the minimum step size
%             1 means unable to find a point that satisfies the curvature condition
%             0 means a desired step size is found
%
% By Wen Huang

% phase one algorithm
status = 0;%

f1 = fvs(end); % function value at current iterate
stepsize1 = initstepsize;%1
c1=params.alpha;
c2=params.beta;
[Q_U_r,S_r,Q_V_r, x2]=retraction(stepsize1,params);
% x2=x1+stepsize1*eta1;
f2=fns.f(x2);
nabla_2 = -fns.Grad(x2);%下降方向
eta2=-nabla_2;
LSinfo.lf = 1;%函数计算次数
LSinfo.lgf = 1;%梯度计算次数
slopex1=initslope;
slopex2=prod(eta1,nabla_2);%h'(alpha)
btiter=1;%函数已使用一次
alo=0.0001;
i=0;
falo=f1;
while(1)
    if(f2>f1+c1*initslope*stepsize1||(i>=1&&f2>=f1))
        aup=stepsize1;
        fv=f1;
        falo=f1;
        faup=f2;
        aloslope=initslope;
        status=0;
        LSinfo.lf = btiter;
        [nabla_2, stepsize, x2, f2, gradf2, status, slopex2, LSinfo] = Zoom(fns, params, initslope, x1, eta1, fv, alo, falo, aloslope, aup, faup, LSinfo);
        break;
    elseif(abs(slopex2)<=-c2*initslope)
        status=0;
        
        LSinfo.lf = btiter;
        stepsize=stepsize1;
         if stepsize<=params.minstepsize
            status=2;
            stepsize=params.minstepsize;
        elseif stepsize>=params.maxstepsize
             status=3;
            stepsize=params.maxstepsize;
        else
            status=0;
        end
        gradf2=-nabla_2;
        break;
    elseif(slopex2>=0)
        aup=stepsize1;
        fv=f1;
        falo=f1;
        faup=f2;
        aloslope=initslope;
        LSinfo.lf = btiter;
        [nabla_2, stepsize, x2, f2, gradf2, status, slopex2, LSinfo] = Zoom(fns, params, initslope, x1, eta1, fv, alo, falo, aloslope, aup, faup, LSinfo);
    end
    %更新初始条件，为下一次迭代做准备
    alo=stepsize1;
    falo=f2;
    stepsize1=-initslope*stepsize1^2/(2*(f2-f1-initslope*stepsize1));
%     eta1=eta2;
    x2.main=x1.main+stepsize1*eta1;%更新x2
    f2=fns.f(x2);
%     gradf2=fns.Grad(x2);
%     eta2 = -gradf2*stepsize1;
%     falo=f2;
    
%     slopex2=-eta1'*eta2;%h'(alpha)
%     LSinfo.lf =LSinfo.lf+ 1;%函数计算次数
% LSinfo.lgf =LSinfo.lgf+ 1;%梯度计算次数
    i=i+1;
    stepsize=stepsize1;
  if iter > params.lsmaxiter
        status = 1;
        stepsize = stepsize1;
        nabla_2 = stepsize * eta1;
        x2 = x1 + nabla_2;
        f2 = fns.f(x2);
        gradf2 = fns.Grad(x2);
        LSinfo.lf = LSinfo.lf + 1;
        LSinfo.lgf = LSinfo.lgf + 1;
        slopex2 = -prod(gradf2,gradf2);
        return;
  end
    
    
    
 end
end




function [eta2, stepsize, x2, f2, gradf2, status, slopex2, LSinfo] = Zoom(fns, params, initslope, x1, eta1, fv, alo, falo, aloslope, aup, faup, LSinfo)
% INPUT:
% initslope : the initial slope
% eta : search direction
% x   : current iterate
% fv  : the function value at the current iterate
% alo : alpha_low
% falo : function value when step size is alpla_low
% aloslope : h'(alpha_low)
% aup : alpha_up
% faup : function value when step size is alpha_up
% fns : a struct that contains required function handles
%     fns.f(x) : return objective function value at x.
%     fns.Grad(x) : return the gradient of objection function at x.
% params: : a struct that contains parameters that are used
%     alpha: the coefficient in the Armijo condition
%     beta : the coefficient in the curvature condition
%     lsmaxiter : the max numbers of inner iteration
%     minstepsize : minimum step size tau1
%     maxstepsize : maximum step size tau2
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
% status    : 4 means some other serious errors happened,
%             3 means step size excesses the maximum step size
%             2 means step size reaches the minimum step size
%             1 means unable to find a point that satisfies the curvature condition
%             0 means a desired step size is found
%
% By Wen Huang

% phase two algorithm
c1=params.alpha;
c2=params.beta;
i=0;
while(1)
    
    aplus=alo+(-aloslope)*(aup-alo)^2/(2*(faup-falo- aloslope*(aup-alo))) ;
    aplus=min(max(aplus,alo+params.tau1*(aup-alo)),min(aplus,alo+params.tau2*(aup-alo)));

    xplus.main=x1.main+aplus*eta1;
    [faplus,xplus]=fns.f(xplus);
    dfaplus=fns.Grad(xplus);
        LSinfo.lf=LSinfo.lf+1;
    LSinfo.lgf=LSinfo.lgf+1;
    if (aup==aplus)<0.00001
        stepsize=aplus;
         if stepsize<=params.minstepsize
            status=2;
            stepsize=params.minstepsize;
        elseif stepsize>=params.maxstepsize
            status=3;
            stepsize=params.maxstepsize;
        else
            status=0;
        end
             status=0;% w
        eta2=stepsize*eta1;
        x2.main=x1.main+eta2;
        [f2 x2]=fns.f(x2);
        LSinfo.lf=LSinfo.lf+1;
        gradf2=fns.Grad(x2);
        LSinfo.lgf=LSinfo.lgf+1;
        slopex2=-eta1'*gradf2;
        return; 
    end
    if(faplus>fv+c1*aplus*initslope||faplus>=falo)
        aup=aplus;
        continue;
    elseif(abs(dfaplus)<=-c2*initslope)
        stepsize=aplus;
        if stepsize<=params.minstepsize
            status=2;
            stepsize=params.minstepsize;
        elseif stepsize>=params.maxstepsize
            status=3;
            stepsize=params.maxstepsize;
        else
            status=0;
        end
         status=0;
        eta2=stepsize*eta1;
        x2.main=x1.main+eta2;
        [f2 x2]=fns.f(x2);
        LSinfo.lf=LSinfo.lf+1;
        gradf2=fns.Grad(x2);
        LSinfo.lgf=LSinfo.lgf+1;
        slopex2=-eta1'*gradf2;
        return;
    elseif(dfaplus>0)
        aup=aplus;
    else
        alo=aplus;
    end
        i=i+1;
    if i==params.lsmaxiter
        stepsize=aplus;
         if stepsize<=params.minstepsize
            status=2;
            stepsize=params.minstepsize;
        elseif stepsize>=params.maxstepsize
            status=3;
            stepsize=params.maxstepsize;
        else
            status=0;
        end
        status=1;
        eta2=stepsize*eta;
        x2=x1+eta2;
        f2=fns.f(x2);
        LSinfo.lf=LSinfo.lf+1;
        gradf2=fns.Grad(x2);
        LSinfo.lgf=LSinfo.lgf+1;
        slopex2=-eta1'*gradf2;
        return;

    end
    
end



end

function a=prod(x,y)
a=real(sum(sum(x.*y)));
end
