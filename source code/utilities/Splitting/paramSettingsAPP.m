function param = paramSettingsAPP(tRes,nframes,Nwaves,distances)
%Author: LWeissinger
%Description: All physical and discretization parameters are specified and
%stored in the struct param
%Input: Optional: Discretization level, actual PWV
%Output: Struct containing all involved physical parameters of the problem
param.Nwaves=Nwaves;
param.mode='real';

%smoothing-parameter multiplicand (Recommended: Do not change)
param.beta=1; 
%Parameter in velocity part in norm of X (Recommended: Do not change - value is not used in current methods)
param.c=4000; 
% Parameters for (linear) Regularization
% either 'apriori'(specify below in this case), 'dp' or 'minerr'(minerr only available in simulation)
param.parameter_choice_rule='apriori'; 
 
% 
param.reg_man=1;
%Parameters for (nonlinear) Landweber part
%activate nesterov acceleration? (1/0)
param.nesterov=1; 
% max iterations in landweber (inner) iteration
param.max_it=5000; 

%max iterations in split (outer) iteration
param.stop=500;
%Set ~=0 to let it run until max iterations (fixed stopping rule)
param.justmax=1;

%% We recommend no changes below this line. Parameters here mostly belong to physical quantities to create realistic simulation waves
%% ----------------------------------------------------

param.m=nframes;
param.effectivedT=tRes;
%fprintf("Calculating length of artery segment...\n")
%[all_distances,length_of_total_artery_segment]=curvedistance(pixdim.*branch,3,8);
%fprintf("Length of artery-segment via Bezier-curves: %.4f m \n",length_of_total_artery_segment/1000)
param.distance=distances;
param.taxis=transpose(0:param.m-1)*param.effectivedT;
end

