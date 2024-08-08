function  [params]= NSS_params(nSig,params)

params.nSig = nSig;
params.Innerloop = 2 ;
params.delta = 0.01; %迭代正则项
if nSig <= 10
    params.win = 30;  % Non-local patch searching window
   
    params.ps = 10;    % Patch size, larger values would get better performance, but will be slower
    params.Iter = 10;  % total iter numbers
    params.nlsp = 70;   % Initial Non-local Patch number

elseif nSig <= 30
    
    params.win = 30;
    params.ps = 12;   
    params.Iter = 12;
    params.nlsp = 80;
   
elseif nSig <= 50
    
    params.win = 30;
    params.ps = 14;   
    params.Iter = 14;
    params.nlsp = 90;
     
end

params.step  = params.ps-1; % The step between neighbor patches, smaller values would get better performance, but will be slower

end

