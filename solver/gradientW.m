function output = gradientW(W,params)
if params.isq ==1
Ws = W.s;
% Ws=s(W);
Ws1 = abs(Ws) < params.epsilon;
% gradientWs1 = params.epsilon*2 * abs(Ws);
gradientWs1 = 1/params.epsilon*5 * abs(Ws);
gradientWs2 = 1;
% gradientWs1 = params.epsilon^(3/2)/4 * Ws;
% gradientWs2 = 0.5 * (abs(Ws)+1e-6).^-0.5;
gradientWs = Ws1.* gradientWs1 + (1 - Ws1).* gradientWs2;
Wx=x(W);
Wx1 = abs(Wx) < params.epsilon;
% gradientWx1 = params.epsilon*2 * abs(Wx);
gradientWx1 = 1/params.epsilon*5* abs(Wx);
gradientWx2 = 1;
% gradientWx1 = params.epsilon^(3/2)/4 * Wx;
% gradientWx2 =0.5 * (abs(Wx)+1e-6).^-0.5;
gradientWx = Wx1.* gradientWx1 + (1 - Wx1).* gradientWx2;
Wy=y(W);
Wy1 = abs(Wy) < params.epsilon;
% gradientWy1 = params.epsilon*2 * abs(Wy);
gradientWy1 = 1/params.epsilon*5 * abs(Wy);
gradientWy2 = 1;
% gradientWy1 = params.epsilon^(3/2)/4 * Wy;
% gradientWy2 = 0.5 * (abs(Wy)+1e-6).^-0.5;
gradientWy = Wy1.* gradientWy1 + (1 - Wy1).* gradientWy2;
Wz=z(W);
Wz1 = abs(Wz) < params.epsilon;
% gradientWz1 = params.epsilon*2 * abs(Wz);
gradientWz1 = 1/params.epsilon*5 * abs(Wz);
gradientWz2 = 1;
% gradientWz1 = params.epsilon^(3/2)/4 * Wz;
% gradientWz2 = 0.5 * (abs(Wz)+1e-6).^-0.5;
gradientWz = Wz1.* gradientWz1 + (1 - Wz1).* gradientWz2;
output=quaternion(gradientWs,gradientWx,gradientWy,gradientWz);
else
 
W1 = abs(W) < params.epsilon;
% gradientWs1 = params.epsilon*2 * abs(Ws);
gradientW1 = 1/params.epsilon*5 * abs(W);
gradientW2 = 1;
% gradientWs1 = params.epsilon^(3/2)/4 * Ws;
% gradientWs2 = 0.5 * (abs(Ws)+1e-6).^-0.5;
gradientW = W1.* gradientW1 + (1 - W1).* gradientW2;
output = gradientW;
end
end
