function  params =  SLRQAdnRGB_params( nSig,str)
params.strname = str;
params.nSig = nSig;
params.isq = 2;
params.epsilon = 1;
if str == 'Schatte'
if nSig <= 10
    params.gamma = 0.5;   
    params.lambda = 0.01;   
    params.eps  =  1e-4;   
    params.rho = 2;
    params.mu = 2;
    params.l1 = 0.01;
    params.l2 = 0.01;
    
elseif nSig <= 30
    params.gamma = 0.5;  
    params.lambda = 0.3;   
    params.eps  =  1e-3;   
    params.rho = 2;
    params.mu = 1;
    params.l1 = 0.1;
    params.l2 = 0.1;   
elseif nSig <= 50
    params.gamma = 0.5;
    params.lambda = 0.5;   
    params.eps  =  1e-3;   
    params.rho = 2;
    params.mu = 2;
    params.l1 = 0.01;
    params.l2 = 0.01;   
end
else
    if nSig <= 10
    params.gamma = 0.1;
    params.lambda = 0.1;   
    params.eps  =  1e-4;   
    params.rho = 2;
    params.mu = 2;
    params.l1 = 0.1;
    params.l2 = 0.1;
    end
end

end

