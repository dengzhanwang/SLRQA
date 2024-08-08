function output=shrinkageW(A,epsilon,W,params)
if params.isq ==1
As=s(A);
Ws=s(W);
signAs= sign(As);
As=signAs.*max(abs(As)-epsilon*Ws,0);
Ax=x(A);
Wx=x(W);
signAx= sign(Ax);
Ax=signAx.*max(abs(Ax)-epsilon*Wx,0);
Ay=y(A);
Wy=y(W);
signAy= sign(Ay);
Ay=signAy.*max(abs(Ay)-epsilon*Wy,0);
Az=z(A);
Wz=z(W);
signAz= sign(Az);
Az=signAz.*max(abs(Az)-epsilon*Wz,0);
output=quaternion(As,Ax,Ay,Az);
else
output=sign(A).*max(abs(A)-epsilon*W,0);
end
end
