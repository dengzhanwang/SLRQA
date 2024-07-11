clear;
clc;
qlena = cast(imreadq('D:\毕业论文\代码\Lenna.png'), 'double') ./ 255;
% image(qlena);
[U,S,V]=svd(qlena);
S2=zeros(size(S));
for jj=1:100
   S2(jj,jj)=S(jj,jj); 
end
q2=U*S2*V';
image(q2);
a=randq(10);
[U,S,V]=svd(a);
% clear;
% clc;
% a=randq(128,128);
% tic;
% [U,S,V]=svd(a);
% toc