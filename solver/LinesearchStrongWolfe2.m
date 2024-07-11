function [eta2, stepsize, x2, f2, gradf2, slopex2, LSinfo, status] = LinesearchStrongWolfe2(eta1, x1, fvs, gradf1, initslope, initstepsize, iter, fns, params)
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
f1= fns.f(x1);%h(0)
stepsize1 = initstepsize;%1
c1=params.alpha;
c2=params.beta;
x2.main=x1.main+stepsize1*eta1;
[f2,x2]=fns.f(x2);
gradf2 = fns.Grad(x2);%下降方向

LSinfo.lf = 1;%函数计算次数
LSinfo.lgf = 1;%梯度计算次数
slopex1=initslope;
slopex2=eta1'*gradf2;%h'(alpha)
btiter=1;%函数已使用一次
alo=0;
i=0;
falo=f1;
aloslope=slopex1;
eta2=-gradf2;
while(1)
    if(f2>f1+c1*initslope*stepsize1||(i>=1&&f2>=falo))
        aup=stepsize1;
        xaup.main=x1.main+aup*eta1;
        [faup xaup]=fns.f(xaup);
        [eta2, stepsize, x2, f2, gradf2, status, slopex2, LSinfo] = Zoom(fns, params, initslope, x1, eta1, f1, alo, falo, aloslope, aup, faup, LSinfo);
        break;
    elseif(abs(slopex2)<=-c2*initslope)
        status=0;
        stepsize=stepsize1;
        if stepsize<=params.minstepsize
            status=2;
            stepsize=params.minstepsize;
        elseif stepsize>=params.maxstepsize
            %              status=3;
            stepsize=params.maxstepsize;
        else
            status=0;
        end
        break;
    elseif(slopex2>=0)
        aup=stepsize1;
        xaup.main=x1.main+aup*eta1;
        [faup xaup]=fns.f(xaup);
        [eta2, stepsize, x2, f2, gradf2, status, slopex2, LSinfo] = Zoom(fns, params, initslope, x1, eta1, f1, alo, falo, aloslope, aup, faup, LSinfo);
    end
    %更新初始条件，为下一次迭代做准备
    alo=stepsize1;
    xalo.main=x1.main+eta1*alo;
    [falo,xalo]=fns.f(xalo);
    dfalo=fns.Grad(xalo);
    aloslope=dfalo'*eta1;
    stepsize2=-initslope*stepsize1^2/(2*(falo-f1-initslope*alo));
            if numel(find(isnan(stepsize2)))
                stepsize=stepsize1;
        status=1;
        return;
            end
    stepsize1=max([stepsize2 1.1*alo]);

    
    x2.main=x1.main+stepsize1*eta1;%更新x2
    [f2,x2]=fns.f(x2);
    gradf2=fns.Grad(x2);
    slopex2=gradf2'*eta1;
    eta2=-gradf2;
    status=0;
    i=i+1;
    
    %   if iter > params.lsmaxiter
    %         status = 1;
    %         stepsize = stepsize1;
    %         eta2 = stepsize * eta1;
    %         x2.main = x1.main +stepsize*eta2;
    %         f2 = fns.f(x2);
    %         gradf2 = fns.Grad(x2);
    %         LSinfo.lf = LSinfo.lf + 1;
    %         LSinfo.lgf = LSinfo.lgf + 1;
    %         slopex2 = -gradf2' * gradf2;
    %         return;
    %   end
    
    
    
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
    aplus=min(max(aplus,alo+params.tau1*(aup-alo)),alo+params.tau2*(aup-alo));
    
    if numel(find(isnan(aplus)))
        status=1;
        eta2=stepsize*eta1;
        x2.main=x1.main+eta2;
        [f2 x2]=fns.f(x2);
        LSinfo.lf=LSinfo.lf+1;
        gradf2=fns.Grad(x2);
        LSinfo.lgf=LSinfo.lgf+1;
        slopex2=-eta1'*gradf2;
        return;
    end
    xplus.main=x1.main+aplus*eta1;
    [faplus,xplus]=fns.f(xplus);
    dfaplus=fns.Grad(xplus);
    aplusslope=eta1'*dfaplus;
    LSinfo.lf=LSinfo.lf+1;
    LSinfo.lgf=LSinfo.lgf+1;
    if (aup-aplus)<0.00001
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
        
        x2.main=x1.main+stepsize*eta1;
        [f2 x2]=fns.f(x2);
        LSinfo.lf=LSinfo.lf+1;
        gradf2=fns.Grad(x2);
        eta2=-gradf2;
        LSinfo.lgf=LSinfo.lgf+1;
        slopex2=-eta1'*gradf2;minstepsize
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
        
        x2.main=x1.main+stepsize*eta1;
        [f2 x2]=fns.f(x2);
        
        LSinfo.lf=LSinfo.lf+1;
        gradf2=fns.Grad(x2);
        eta2=-gradf2;
        LSinfo.lgf=LSinfo.lgf+1;
        slopex2=eta1'*gradf2;
        return;
    elseif(aplusslope>0)
        aup=aplus;
    else
        alo=aplus;
    end
    i=i+1;
    %     if i==params.lsmaxiter
    %         stepsize=aplus;
    %          if stepsize<=params.minstepsize
    %             status=2;
    %             stepsize=params.minstepsize;
    %         elseif stepsize>=params.maxstepsize
    %             status=3;
    %             stepsize=params.maxstepsize;
    %         else
    %             status=0;
    %         end
    %         status=1;
    %         eta2=stepsize*eta;
    %         x2=x1+eta2;
    %         f2=fns.f(x2);
    %         LSinfo.lf=LSinfo.lf+1;
    %         gradf2=fns.Grad(x2);
    %         LSinfo.lgf=LSinfo.lgf+1;
    %         slopex2=-eta1'*gradf2;
    %         return;
    %
    %     end
    
end



end
