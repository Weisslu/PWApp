function [rel_error] = dataError_man(rho1,rho2,rho3,rho1_rec,rho2_rec,rho3_rec)
% fitting error, but not in fourier space. Only possible for simulation
% data, when original split waves are known.
abs=norm(rho1-rho1_rec)^2+norm(rho2-rho2_rec)^2+norm(rho3-rho3_rec)^2;
ref=norm(rho1,2)^2+norm(rho2)^2+norm(rho3)^2;
rel_error=sqrt(abs/ref);
end