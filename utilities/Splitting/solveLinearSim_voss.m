function [rho_rec,alpha] = solveLinearSim_voss(rho,param)
%solve Pulsewave splitting problem with direct regularization approach
%according to [Voss]
%% --------------------------------------------
% CARE: not updated for the general N-wave case
%% --------------------------------------------
Rho1=rho.data{1};
Rho2=rho.data{2};
alpha=param.alpha;
%Important difference to solveLinear_tikh: Method works only for rounded shifts. Reason still unclear 
tau=floor(param.tau(2));


W = exp(-1i * 2 * pi * tau * transpose([0:floor(param.m/2)-1 floor(-param.m/2):-1]) / param.m); 

denominator=directRegDenominator(W,alpha);
Rho1f=(Rho1-Rho2.*W)./denominator;
Rho2b=(Rho2-Rho1.*W)./denominator;
rho_rec.for{1}=ifft(Rho1f);
rho_rec.back{2}=ifft(Rho2b);

%% Compute best possible parameter
endit=0;
    if param.delta~=0 
        a_start=0.1*param.delta;
    else 
        a_start=0.1;
    end
    while ~endit
        if param.reg_man 
            break; 
        end
        a1=0;
        a3=a_start;
        a2=a3/2;
        x1=zeros(param.m,2);
        x2=zeros(param.m,2);
        x3=zeros(param.m,2);
        %Compute solutions here!TODO
        denominator=directRegDenominator(W,a1);
        X1_1f=(Rho1-Rho2.*W)./denominator;
        X1_2b=(Rho2-Rho1.*W)./denominator;
        x1(:,1)=ifft(X1_1f);
        x1(:,2)=ifft(X1_2b);
        denominator=directRegDenominator(W,a2);
        X2_1f=(Rho1-Rho2.*W)./denominator;
        X2_2b=(Rho2-Rho1.*W)./denominator;
        x2(:,1)=ifft(X2_1f);
        x2(:,2)=ifft(X2_2b);
        denominator=directRegDenominator(W,a3);
        X3_1f=(Rho1-Rho2.*W)./denominator;
        X3_2b=(Rho2-Rho1.*W)./denominator;
        x3(:,1)=ifft(X3_1f);
        x3(:,2)=ifft(X3_2b);
        err1=relativeError(rho1f,rho2b,x1(:,1),x1(:,2));
        err2=relativeError(rho1f,rho2b,x2(:,1),x2(:,2));
        err3=relativeError(rho1f,rho2b,x3(:,1),x3(:,2));
        for k=1:15
            errv=[err1 err2 err3];
            if err1==min(errv)
                a3=a2;
                err3=err2;
                x3=x2;
                a2=(a3+a1)/2;
                x2=zeros(param.m,2);
                denominator=directRegDenominator(W,a2);
                X2_1f=(Rho1-Rho2.*W)./denominator;
                X2_2b=(Rho2-Rho1.*W)./denominator;
                x2(:,1)=ifft(X2_1f);
                x2(:,2)=ifft(X2_2b);
                err2=relativeError(rho1f,rho2b,x2(:,1),x2(:,2));
            elseif err3==min(errv)
                a1=a2;
                err1=err2;
                x1=x2;
                a2=(a3+a1)/2;
                x2=zeros(param.m,2);
                denominator=directRegDenominator(W,a2);
                X2_1f=(Rho1-Rho2.*W)./denominator;
                X2_2b=(Rho2-Rho1.*W)./denominator;
                x2(:,1)=ifft(X2_1f);
                x2(:,2)=ifft(X2_2b);
                err2=relativeError(rho1f,rho2b,x2(:,1),x2(:,2));
            elseif err2==min(errv)
                a4=(a2+a1)/2;
                a5=(a3+a2)/2;
                x4=zeros(param.m,2);
                x5=zeros(param.m,2);
                denominator=directRegDenominator(W,a4);
                X4_1f=(Rho1-Rho2.*W)./denominator;
                X4_2b=(Rho2-Rho1.*W)./denominator;
                x4(:,1)=ifft(X4_1f);
                x4(:,2)=ifft(X4_2b);
                err4=relativeError(rho1f,rho2b,x4(:,1),x4(:,2));
                denominator=directRegDenominator(W,a5);
                X5_1f=(Rho1-Rho2.*W)./denominator;
                X5_2b=(Rho2-Rho1.*W)./denominator;
                x5(:,1)=ifft(X5_1f);
                x5(:,2)=ifft(X5_2b);    
                err5=relativeError(rho1f,rho2b,x5(:,1),x5(:,2));
                errv=[ err4 err2 err5];
                if err4==min(errv)
                    a3=a2;
                    err3=err2;
                    x3=x2;
                    a2=a4;
                    err2=err4;
                    x2=x4;
                elseif err5==min(errv)
                    a1=a2;
                    err1=err2;
                    x1=x2;
                    a2=a5;
                    err2=err5;
                    x2=x5;
                elseif err2==min(errv)
                    a1=a4;
                    err1=err4;
                    x1=x4;
                    a3=a5;
                    err3=err5;
                    x3=x5;
                end
            end
        end
        errv=[err1 err2 err3];
        if err1==min(errv)
            alpha=a1;
            rho_rec.for{1}=x1(:,1);
            rho_rec.back{2}=x1(:,2);
        elseif err2==min(errv)
            alpha=a2;
            rho_rec.for{1}=x2(:,1);
            rho_rec.back{2}=x2(:,2);
        elseif err3==min(errv)
            alpha=a3;
            rho_rec.for{1}=x3(:,1);
            rho_rec.back{2}=x3(:,2);
        end 
        if alpha==a_start
            a_start=2*alpha;
        else
            fprintf("alpha found as %.8f with min-error after %i iterations\n",alpha,k);
            endit=1;
        end
    end
rho_rec.for{2}=FourierShift(rho_rec.for{1},param.tau(2));
rho_rec.back{1}=FourierShift(rho_rec.back{2},param.tau(2));
for it=1:param.Nwaves
    rho_rec.sum{it}=rho_rec.for{it}+rho_rec.back{it};
end
end