%% Author: LWeissinger
%% Description: All physical and discretization parameters are specified and
%% stored in the struct param
%% Output: Struct containing all involved physical parameters of the problem
function param = paramSettingsAPP(tRes,nframes,Nwaves,distances)
    param.mode='real';
    param.parameter_choice_rule='apriori'; %fixed here - other options only available in simulation setting
    % max iterations in Gradient Descent/FISTA
    param.max_it=5000;
    param.Nwaves=Nwaves;
    param.m=nframes;
    param.effectivedT=tRes;
    param.distance=distances;
    param.taxis=transpose(0:param.m-1)*param.effectivedT;
end

