function [rho_rec,pwv] = bjornfotsmethod(rho,param)
%Bjornfots method assumes that there is no backwards wave. In our
%simulation data, there is the assumtion of a backwardswave, so we want to
%implement bjornfots methods as described in their paper to check its
%accuracy in presence of a backwards wave
% WL 19.04.2023, Ithaca
x_0=zeros(param.m+1,1);
for j=1:param.Nwaves
    x_0(1:param.m)=x_0(1:param.m)+rho.sum{j}/param.Nwaves;
end
x_0(param.m+1)=10;
%Function-handle here!
%TODO
fun=@(x)MLE(x,rho,param);
x=fminunc(fun,x_0);
rho_rec=x(1:param.m);
pwv=x(param.m+1);
end

function f=MLE(x,rho,param)
data=zeros(param.Nwaves,param.m);
V=zeros(param.Nwaves,param.m);
tau=param.distance/param.effectivedT/x(param.m+1);
for j=1:param.Nwaves
    data(j,:)=rho.sum{j};
    V(j,:)=FourierShift(x(1:param.m),tau(j));
end
f=sum((V-data).^2,'all');
end