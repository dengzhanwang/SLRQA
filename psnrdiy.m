function p=psnr(rim,im)

%The function calculate the peak signal-to-noise ration.
% p=10*(log10(size(im,1)*size(im,2))-log10(norm(rim(:)-im(:))^2));

[m,n] = size(im);
MSE = norm(rim-im,'fro')^2/(m*n); 
p = -10*log10(MSE); 