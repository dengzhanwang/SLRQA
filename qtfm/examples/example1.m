% a=imread('D:\BaiduNetdiskDownload\A01建模视频及播放器（下载不了或速度慢可以用群下载的方式）\上课用的课件和代码(下载后记得解压，所有视频配套的都在里面) (1)\第13讲.奇异值分解SVD和图形处理\代码和例题数据\千与千寻.jpg');
% a=a;
% image(a);
clear all;clc;close all;  
img = imread('Lenna.png');  
% my picture is named lenna.bmp while yours may be not  
I = rgb2gray(img);  
% Attention: we use the axis off to get rid of the axis.  
figure(1),image(I); %equals to imagesc(I,[1 64]);you can try it.  
colorbar,title('show by image in figure1');axis off;  
figure(2),imagesc(I);  
%equals to imagesc(I,[min(I(:)) max(I(:))]);you can try it.  
colorbar,title('show by imagesc in figure2');axis off;  
%colormap(gray) %use this statement you can get a gray image.  
figure(3),imshow(I),colorbar,title('show by imshow in figure3'); 
imwrite(uint8(I), 'E:\毕业论文\代码\qtfm\examples\img_gray.jpg');
set(gcf,'unit','centimeters','position',[1,2,200,140]);
imshow(I)