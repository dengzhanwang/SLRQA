function imageql(A,x1,y1)
%This function aims to image the quaternion matrix A
[m,n]=size(A);
B=zeros(m,n,3);
B(:,:,1)=x(A);
B(:,:,2)=y(A);
B(:,:,3)=z(A);
image(x1,y1,B);
end