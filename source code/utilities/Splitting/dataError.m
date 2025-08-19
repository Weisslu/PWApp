%Author: LWeissinger, 31.07.2025
function [rel_error] = dataError(rho,rho_rec,Nwaves)
% data error, but not in fourier space. Should be equal due to parseval,
% plancherel?
abs=0;
ref=0;
for it=1:Nwaves
    abs=abs+norm(rho.sum{it}-rho_rec.sum{it})^2;
    ref=ref+norm(rho.sum{it})^2;
end
rel_error=sqrt(abs/ref);
end
