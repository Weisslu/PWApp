function [rho_rec,fitting_error,data_error,pwv_rec,alpha] = solveLinear_tikh_min(rho,PWVs,param,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
steps=length(PWVs);
fitting_error=zeros(steps,1);
data_error=zeros(steps,1);
alpha=zeros(steps,1);
str1='sim';
str2='tikh';
str5='tikh + bw cond';
str3='FISTA';
str6='FISTA + bw cond';
str4='direct';
if param.showwaitbar~=0
    h=waitbar(0,"Estimating PWV");
end
for run=1:steps
    if param.showwaitbar~=0
        waitbar(run/steps,h,"Estimating PWV");
    end
    param.tau=param.distance/param.effectivedT/PWVs(run);
    if strcmp(method,str2)
        [rho_rec,alpha(run)]=solveLinear_tikh(rho,param);
    elseif strcmp(method,str3)
        [rho_rec,alpha(run)]=solveLinear_FISTA(rho,param);
    elseif strcmp(method,str6)
        [rho_rec,alpha(run)]=solveLinear_FISTAbw(rho,param);
    elseif strcmp(method,str4)
        [rho_rec,alpha(run)]=solveLinearSim_voss(rho,param);
    elseif strcmp(method,str5)
        [rho_rec,alpha(run)]=solveLinear_tikh_GD(rho,param);
    end
    if strcmp(str1,param.mode)
        fitting_error(run)=fittingError(rho.for{1},rho.back{param.Nwaves},rho_rec.for{1},rho_rec.back{param.Nwaves});
    end
    data_error(run)=dataError(rho,rho_rec,param.Nwaves);
end
j=find(data_error==min(data_error),1);
pwv_rec=PWVs(j);
param.tau=param.distance/param.effectivedT/pwv_rec;
if strcmp(method,str2)
    [rho_rec,alpha_min] = solveLinear_tikh(rho,param);
elseif strcmp(method,str3)
    [rho_rec,alpha_min] = solveLinear_FISTA(rho,param);
elseif strcmp(method,str5)
    [rho_rec,alpha_min] = solveLinear_tikh_GD(rho,param);
elseif strcmp(method,str6)
    [rho_rec,alpha_min] = solveLinear_FISTAbw(rho,param);
elseif strcmp(method,str4)
    [rho_rec,alpha_min]=solveLinearSim_voss(rho,param);
end
if param.showwaitbar~=0
    close(h);
end
end