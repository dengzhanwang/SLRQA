function stop = IsStopped(iter, Xs, fvs, ngfs, ngf0, times, fns, params)
% Check stopping criterion
%
% INPUT:
% iter : the current iteration
% Xs   : previous iterates
% fvs  : previous function values, the last one is the value at the current iterate
% ngfs : previous norms of gradients, the last one is the norm at the current iterate
% times : computational times until every iteration
% fns : a struct that contains required function handles
% params: : a struct that contains parameters that are used
%     params.Tolerance          : tolerance of stopping criterion.
%     params.Max_Iteration      : the maximum number of iterations.
%
% OUTPUT:
% stop : stop the algorithm or not
%     0 : not stop
%     1 : stop
% 
% By Wen Huang
stop=0;
% if sum(times)>10
%     stop =1 ;
% end
if(abs(ngfs(end)/ngfs(1))<eps)
    stop =1;
end
if(iter>params.lsmaxiter)
    stop=1;
end
end
