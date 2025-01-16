%% Applies the nonlinear operator F to given input waves and PWV and returns the result as a vector
function [y] = operatorF(x1,x2,u,param)
y=zeros(param.Nwaves*param.m,1);
for it=1:param.Nwaves
    tau_pos=param.distance(it)/param.effectivedT/u;
    tau_neg=(param.distance(param.Nwaves)-param.distance(it))/param.effectivedT/u;
    W_pos = transpose(exp(-1i * 2 * pi * tau_pos * [0:floor(param.m/2)-1 floor(-param.m/2):-1] / param.m)); 
    W_neg = transpose(exp(-1i * 2 * pi * tau_neg * [0:floor(param.m/2)-1 floor(-param.m/2):-1] / param.m)); 
    y(((it-1)*param.m+1):(it*param.m))=x1.*W_pos+x2.*W_neg;
end