function [Frechet] = FreDerAdjoint(x1,x2,u,param)
%Adjoint of Frechet derivative in a point (x1,x2,u) as matrix form
%Without adjoint embedding!
%WL 29.11.2022
c=100; %constant in velocity component of norm
%Is alternatively applied in solveNonlinearLandweber with the adjoint
%embedding - if the Operator ' is used to compute the adjoint (for steepest
%descent recommended.
Frechet=zeros(2*param.m+1,param.Nwaves*param.m);
omega=[0:floor(param.m/2)-1 floor(-param.m/2):-1]';
for it=1:param.Nwaves
    tau_pos=param.distance(it)/param.effectivedT/u;
    tau_neg=(param.distance(param.Nwaves)-param.distance(it))/param.effectivedT/u;
    W_pos = exp(1i * 2 * pi * tau_pos * omega / param.m); 
    W_neg = exp(1i * 2 * pi * tau_neg * omega / param.m); 
    a_pos=2*pi*1i*omega*param.distance(it)/param.effectivedT/param.m;
    a_neg=2*pi*1i*omega*(param.distance(param.Nwaves)-param.distance(it))/param.effectivedT/param.m;
    Frechet(1:param.m,((it-1)*param.m+1):(it*param.m))=diag(W_pos);
    Frechet(param.m+1:2*param.m,((it-1)*param.m+1):(it*param.m))=diag(W_neg);
    Frechet(2*param.m+1,((it-1)*param.m+1):(it*param.m))=((x1.*a_pos/(u^2)).*W_pos+(x2.*a_neg/(u^2)).*W_neg)/c;
end
end