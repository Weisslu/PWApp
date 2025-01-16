%% Applies the nonlinear operator F to given input waves and PWV and returns the result as a vector
function [y] = operatorF_k(x1,x2,u,k,param)
y=zeros(param.m,1);
tau_pos=param.distance(k)/param.effectivedT/u;
tau_neg=(param.distance(param.Nwaves)-param.distance(k))/param.effectivedT/u;
W_pos = transpose(exp(-1i * 2 * pi * tau_pos * [0:floor(param.m/2)-1 floor(-param.m/2):-1] / param.m)); 
W_neg = transpose(exp(-1i * 2 * pi * tau_neg * [0:floor(param.m/2)-1 floor(-param.m/2):-1] / param.m)); 
y=x1.*W_pos+x2.*W_neg;
end