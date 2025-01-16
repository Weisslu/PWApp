%% Splitting operator for N waves as matrix. 
%% input: 2 "split" waves
%% Output: N- summed waves
function [A] = operatorA(param)
A=zeros(param.Nwaves*param.m,2*param.m);
for it=1:param.Nwaves
    W_pos = transpose(exp(-1i * 2 * pi * param.tau(it) * [0:floor(param.m/2)-1 floor(-param.m/2):-1] / param.m)); 
    W_neg = transpose(exp(-1i * 2 * pi * (param.tau(param.Nwaves)-param.tau(it)) * [0:floor(param.m/2)-1 floor(-param.m/2):-1] / param.m)); 
    A(((it-1)*param.m+1):(it*param.m),1:param.m)=diag(W_pos);
    A(((it-1)*param.m+1):(it*param.m),param.m+1:2*param.m)=diag(W_neg);      
end
% figure
% surf(real(A));
% view([0 0 90])
% shading('interp');
end