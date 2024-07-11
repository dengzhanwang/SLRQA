function [U S V]=qqsvd(A)
%输入的是以三维矩阵表示的四元数组
if size(A,3)==4
    [m,n,~]=size(A);
    A1=A(:,:,1);
    if strcmp(class(A1),'uint8')==1
        A1=zeros(m,n);
    end
    A2=A(:,:,2);
    A3=A(:,:,3);
    A4=A(:,:,4);
    A_2=[A1+A2*i,A3+A4*i;-A3+A4*i,A1-A2*i];
    [U S V]=svd(A_2);
    U1_r=real(U(1:m,2:2:end));
    U1_i=imag(U(1:m,2:2:end));
    U2_r=-real(U(m+1:end,2:2:end));
    U2_i=imag(U(m+1:end,2:2:end));
    V1_r=real(V(1:n,2:2:end));
    V1_i=imag(V(1:n,2:2:end));
    V2_r=-real(V(n+1:end,2:2:end));
    V2_i=imag(V(n+1:end,2:2:end));
    U=zeros(size(A));
    V=zeros(size(A));
    U(:,:,1)=U1_r;
    U(:,:,2)=U1_i;
    U(:,:,3)=U2_r;
    U(:,:,4)=U2_i;
    V(:,:,1)=V1_r;
    V(:,:,2)=V1_i;
    V(:,:,3)=V2_r;
    V(:,:,4)=V2_i;
%     U=quaternion(U1_r,U1_i,U2_r,U2_i);
% V=quaternion(V1_r,V1_i,V2_r,V2_i);
S=S(2:2:end,2:2:end);
else
    [U S V]=svd(A);
end