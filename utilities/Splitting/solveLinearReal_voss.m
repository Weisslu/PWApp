function [rho1f_rec,rho2b_rec,rho1b_rec,rho2f_rec,rho1_rec,rho2_rec,alpha] = solveLinearReal_voss(rho1,rho2,param)
%OUTDATED
%solve Pulsewave splitting problem with direct regularization approach
%according to [Voss]
rho1=rho1(:);
rho2=rho2(:);
Rho1=fft(rho1);
Rho2=fft(rho2);
alpha=param.alpha;
%Important difference to solveLinear_tikh: Method works only for rounded shifts. Reason still unclear 
tau=floor(param.tau);


W = exp(-1i * 2 * pi * tau * transpose([0:floor(param.m/2)-1 floor(-param.m/2):-1]) / param.m); 

denominator=directRegDenominator(W,param.alpha);
Rho1f=(Rho1-Rho2.*W)./denominator;
Rho2b=(Rho2-Rho1.*W)./denominator;
rho1f_rec=ifft(Rho1f);
rho2b_rec=ifft(Rho2b);
rho2f_rec=FourierShift(rho1f_rec,param.tau);
rho1b_rec=FourierShift(rho2b_rec,param.tau);
rho1_rec=rho1f_rec+rho1b_rec;
rho2_rec=rho2f_rec+rho2b_rec;

end