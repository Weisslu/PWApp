function [rho_rec,data_error,pwv_rec,alpha] = solveLinearReal_tikh_min(rho,PWVs,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
steps=length(PWVs);
data_error=zeros(steps,1);
alpha=zeros(steps,1);
for run=1:steps
    param.tau=param.distance/param.effectivedT/PWVs(run);
    [rho_rec,alpha(run)]=solveLinear_tikh(rho,param);
    data_error(run)=dataError(rho,rho_rec,param.Nwaves);
end
j=data_error==min(data_error);
pwv_rec=PWVs(j);
param.tau=param.distance/param.effectivedT/pwv_rec;
[rho_rec,alpha_min] = solveLinear_tikh(rho,param);
end