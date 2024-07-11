function [initstepsize, info] = InitStepSizeBB(iter, initslope, fvs, Xs, GFs, Etas, fns, params)
% Find initial step size
%
% INPUT:
% iter : the current iteration
% initslope : the initial slope
% fvs  : previous function values, the last one is the value at the current iterate
% Xs   : previous iterates
% GFs  : previous gradients
% Etas : previous search directions
% fns : a struct that contains required function handles
% params: : a struct that contains parameters that are used
%     FirstIterInitStepSize: the initial step size at the very first iteration
%
% OUTPUT:
% initstepsize : initial step size
% info: a struct that contains debug information
% 
% By Wen Huang

    if(iter == 1)
        initstepsize = params.FirstIterInitStepSize; % initial step size at the very first iteration
    else
        % BB stepsize with safeguard. Note that the safeguard for
        % monotonic line search and non-monotonic linesearch is different!
        if(params.step==1)%µ¥µ÷ËÑË÷

      ak=2*(fvs(end)-fvs(end-1))/-sum(sum(s(Etas{end}'*Etas{end})));
        else
            ak = params.FirstIterInitStepSize;
        end
   sk=Xs{end}-Xs{end-1};
yk=GFs{end}-GFs{end-1};
initstepsize=abs(sum(sum(s(sk'*yk)))/sum(sum(s(yk'*yk))));
initstepsize=max(min(initstepsize,params.w2*ak),params.w1*ak);
    end
    info.lf = 0; info.lgf = 0;

end
