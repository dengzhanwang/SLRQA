function  par  =  SearchNeighborIndex(par)
% This Function Precompute the all the patch indexes in the Searching window
% -NeighborIndex is the array of neighbor patch indexes for each keypatch
% -NumIndex is array of the effective neighbor patch numbers for each keypatch
% -SelfIndex is the index of keypatches in the total patch index array
par.maxr = par.h - par.ps + 1;
par.maxc = par.w - par.ps + 1;
r          =  1:par.step:par.maxr;  %keypatch 行指标 
par.r          =  [r r(end) + 1:par.maxr];  % 不能整除时，余数取重叠大patch
c          =  1:par.step:par.maxc;   %keypatch 的列指标
par.c          =  [c c(end) + 1:par.maxc];
par.lenr = length(par.r);
par.lenc = length(par.c);
par.ps2 = par.ps^2;
par.ps2ch = par.ps2 * par.ch;
% Total number of patches in the test image
par.maxrc = par.maxr * par.maxc;
% Total number of seed patches being processed
par.lenrc = par.lenr * par.lenc;   %keypatch的个数
% index of each patch in image
par.Index     =   (1:par.maxrc);
par.Index    =   reshape(par.Index, par.maxr, par.maxc);
% preset variables for all the patch indexes in the Searching window
par.NeighborIndex    =   int32(zeros(4 * par.win^2, par.lenrc));
par.NumIndex        =   int32(zeros(1, par.lenrc));
par.SelfIndex   =   int32(zeros(1, par.lenrc));

for  i  =  1 : par.lenr
    for  j  =  1 : par.lenc   %第i行第j列的patch
        row = par.r(i);   %当前keypatch所在的行列
        col = par.c(j);
        off = (col-1) * par.maxr + row;   %按列数，keypatch中第一个元素在maxr*maxc中的指标
        off1 = (j-1) * par.lenr + i;  %按列数keypatch的指标，即是第几个keypatch
        
        % the range indexes of the window for searching the similar patches
        rmin    =   max( row - par.win, 1 );  %keypatch的指标变化范围
        rmax    =   min( row + par.win, par.maxr );
        cmin    =   max( col - par.win, 1 );
        cmax    =   min( col + par.win, par.maxc );
        
        idx     =   par.Index(rmin:rmax, cmin:cmax);
        idx     =   idx(:);  %拉成一列
        
        par.NumIndex(off1)  =  length(idx);   %每个keypatch搜索临近像素的个数
        par.NeighborIndex(1:par.NumIndex(off1), off1)  =  idx;   %每个keypatch搜索临近像素的指标
        par.SelfIndex(off1) = off;
    end
end
