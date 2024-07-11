function out=color2complex(A)
%This function aims to  transform a color image data into its complex adjoint form
if size(A,3) ==3
[m n]=size(A);
A1=A(:,:,1);
A2=A(:,:,2);
A3=A(:,:,3);
out=[A1 ,A2+i*A3;A3*i-A2 , A1];
end

end