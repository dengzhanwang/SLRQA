function output=qtranspose(A)
output=zeros(size(A));
output(:,:,1)=A(:,:,1)';
output(:,:,2)=-A(:,:,2)';
output(:,:,3)=-A(:,:,3)';
output(:,:,4)=-A(:,:,4)';
end