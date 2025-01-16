function [rho_rec,fitting_error,data_error,pwv_rec,alpha] = solveLinear_tikh_min(rho,PWVs,param,tikh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
steps=length(PWVs);
fitting_error=zeros(steps,1);
data_error=zeros(steps,1);
alpha=zeros(steps,1);
for run=1:steps
    param.tau=param.distance/param.effectivedT/PWVs(run);
    if tikh==1
        [rho_rec,alpha(run)]=solveLinear_tikh(rho,param);
    elseif tikh==2
        [rho_rec,alpha(run)]=solveLinear_FISTA(rho,param);
    else
        [rho_rec,alpha(run)]=solveLinearSim_voss(rho,param);
    end
    if strcmp(str1,mode)
        fitting_error(run)=fittingError(rho.for{1},rho.back{param.Nwaves},rho_rec.for{1},rho_rec.back{param.Nwaves});
    end
    data_error(run)=dataError(rho,rho_rec,param.Nwaves);
end
j=data_error==min(data_error);
pwv_rec=PWVs(j);
param.tau=param.distance/param.effectivedT/pwv_rec;
if tikh==1
    [rho_rec,alpha_min] = solveLinear_tikh(rho,param);
elseif tikh==2
    [rho_rec,alpha_min] = solveLinear_FISTA(rho,param);
else
    [rho_rec,alpha_min]=solveLinearSim_voss(rho,param);
end

end