function imageq(A)
%This function aims to image the tensor matrix A(m*n*4)
[m,n]=size(A);
B=zeros(m,n,3);
B(:,:,1)=x(A);
B(:,:,2)=y(A);
B(:,:,3)=z(A);
image(B);
end