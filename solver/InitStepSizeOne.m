function [initstepsize, info] = InitStepSizeOne(iter, initslope, fvs, Xs, GFs, Etas, fns, params)
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
%
% OUTPUT:
% initstepsize : initial step size
% info: a struct that contains debug information
% 
% By Wen Huang

    initstepsize = 1;
    info.lf = 0; info.lgf = 0;
end
