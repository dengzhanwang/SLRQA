function [X]=creatematrix(m) %%生成小波变换矩阵
    k=1;
    r2=sqrt(2);
    a=1/r2;
    while(2*k<=m)
        k=k*2;
    end
    for j=k+1:m
       X(j,j)=1;
    end
    block=1;
    while(1<k)
        k2=k/2;
        for j=1:k2
            X(1+(j-1)*2*block:1+(j-1)*2*block+block-1,k2+j)=a;
            X(1+(j-1)*2*block+block:1+(j-1)*2*block+block+block-1,k2+j)=-a;
        end
        k=k/2;
        a=a/r2;
        block=block*2;
    end
    X(1:block,1)=a*r2;
    Y=0;
end