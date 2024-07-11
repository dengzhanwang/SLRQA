function output=qtoc(quat)
% this function aims to convert a pure quaternion matrix to an 3d matrix.
if strcmp(class(quat),'quaternion')==1
[m,n]=size(quat);
xx=x(quat);
xy=y(quat);
xz=z(quat);
output=zeros(m,n,3);
output(:,:,1)=xx;
output(:,:,2)=xy;
output(:,:,3)=xz;
end

end