%% Author: LWeissinger
%% Solve linear Pulse wave splitting problem for a set of predefined PWVs.
%% Best fitting solution is the reconstructed PWV.
function [rho_rec,fitting_error,data_error,pwv_rec,success,alpha] = solveLinear_tikh_min(rho,PWVs,param,method)
    success = 0;
    steps=length(PWVs);
    fitting_error=zeros(steps,1);
    data_error=zeros(steps,1);
    alpha=zeros(steps,1);
    str1='sim';
    str2='tikh';
    str3='FISTA';
    str4='FISTA + bw cond';
    str5='tikh + bw cond';
    
    if param.showwaitbar~=0
        h=waitbar(0,"Estimating PWV", 'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
        setappdata(h, 'canceling', 0);
    end
    for run=1:steps
        if param.showwaitbar~=0
            if getappdata(h, 'canceling')
                disp('User cancelled the operation.');
                delete(h)
                break;
            end
            waitbar(run/steps,h,"Estimating PWV");
        end
        param.tau=param.distance/param.effectivedT/PWVs(run);
        if strcmp(method,str2)
            [rho_rec,alpha(run)]=solveLinear_tikh(rho,param);
        elseif strcmp(method,str3)
            [rho_rec,alpha(run)]=solveLinear_FISTA(rho,param);
        elseif strcmp(method,str4)
            [rho_rec,alpha(run)]=solveLinear_FISTAbw(rho,param);
        elseif strcmp(method,str5)
            [rho_rec,alpha(run)]=solveLinear_tikh_GD(rho,param);
        end
        if strcmp(str1,param.mode)
            fitting_error(run)=fittingError(rho.for{1},rho.back{param.Nwaves},rho_rec.for{1},rho_rec.back{param.Nwaves});
        end
        data_error(run)=dataError(rho,rho_rec,param.Nwaves);
        if run == steps
            success = 1;
            if param.showwaitbar~=0
                delete(h);
            end 
        end
    end
    j=find(data_error==min(data_error),1);
    pwv_rec=PWVs(j);
    param.tau=param.distance/param.effectivedT/pwv_rec;
    if strcmp(method,str2)
        [rho_rec,~] = solveLinear_tikh(rho,param);
    elseif strcmp(method,str3)
        [rho_rec,~] = solveLinear_FISTA(rho,param);
    elseif strcmp(method,str4)
        [rho_rec,~] = solveLinear_tikh_GD(rho,param);
    elseif strcmp(method,str5)
        [rho_rec,~] = solveLinear_FISTAbw(rho,param);
    end
end