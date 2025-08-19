%Author: LWeissinger
function [rho_rec,alpha] = solveLinear_tikh(rho,param)
%% Solve Pulsewave splitting problem with Tikhonov Regularization
%% Transform Problem into Frequency Domain
if ~(isfield(param,'tau')) 
    error("No PWV known to run linear skript")
end
% Force the input to be a column vector.
for it=1:param.Nwaves
    rho.sum{it} = rho.sum{it}(:);
end

%% Get discrete matrix A corresponding to operator A in the fourier domain

A=operatorA(param);

%% Precomputations

Y=zeros(param.Nwaves*param.m,1);
for it=1:param.Nwaves
    Y(((it-1)*param.m+1):(it*param.m))=rho.data{it};
end
%Include adjoint embedding operator for weighted spaces
d=(param.beta*abs([0:floor(param.m/2)-1 floor(-param.m/2):-1]).^2+1).^(-param.s);
Adj_embedding=diag([d d]);
Aadj=Adj_embedding*A';


%% Parameter choice rule
startbed=0;
endbed=0;    
endit=0;
gen=0;
kit=0;
alpha=param.alpha;
str1='sim';
if strcmp(param.mode,str1)
    noise=param.delta;%Kann bei bedarf verändert werden: siehe auch heuristische parameter choice rules
else
    noise=0.01;
end
switch param.parameter_choice_rule
    case 'apriori'
        X=(Aadj*A+param.alpha*eye(2*param.m))\(Aadj*Y);
        rho_rec.for{1}=real(ifft(X(1:param.m)));
        rho_rec.back{param.Nwaves}=real(ifft(X(param.m+1:2*param.m)));
        %fprintf("alpha chosen apriori as %.8f \n",param.alpha);
    case 'dp'
        alpha=0.1*param.delta;
        while ~endit
            alpha
            kit=kit+1;
            X=(Aadj*A+alpha*eye(2*param.m))\(Aadj*Y);
            Y_rec=A*X;
            dis=0;
            for it=1:param.Nwaves
                rho_rec.sum{it}=real(ifft(Y_rec(((it-1)*param.m+1):it*param.m)));
                dis=dis+norm(rho.sum{it}-rho_rec.sum{it},'fro')^2;
            end
            discrepancy=sqrt(dis);
            if startbed==0
                if discrepancy>param.tau_dp*noise
                    startbed=1;
                    alpha_start=alpha;
                else
                    alpha=2*alpha;
                end
            end
            if startbed==1
                gen=gen+1;
                if discrepancy<=param.tau_dp*noise
                    alpha=alpha+alpha_start/(2^gen);
                    endbed=1;
                else
                    alpha=alpha-alpha_start/(2^gen);
                end
            end
            if gen>=14 && endbed==1
                rho_rec.for{param.Nwaves}=real(ifft(X(1:param.m)));
                rho_rec.back{param.Nwaves}=real(ifft(X(param.m+1:2*param.m)));
                endit=1;
            end
        end
        fprintf("alpha found with discrepancy Principle as %.8f  after %i iterations\n",alpha,kit);           
    case'minerr'
        if param.delta~=0 
            a_start=0.1*param.delta;
        else 
            a_start=0.1;
        end
        while ~endit
            a1=0;
            a3=a_start;
            a2=a3/2;
            x1=zeros(param.m,2);
            x2=zeros(param.m,2);
            x3=zeros(param.m,2);
            X1=(Aadj*A+a1*eye(2*param.m))\(Aadj*Y);
            x1(:,1)=real(ifft(X1(1:param.m)));
            x1(:,2)=real(ifft(X1(param.m+1:2*param.m)));
            X2=(Aadj*A+a2*eye(2*param.m))\(Aadj*Y);
            x2(:,1)=real(ifft(X2(1:param.m)));
            x2(:,2)=real(ifft(X2(param.m+1:2*param.m)));
            X3=(Aadj*A+a3*eye(2*param.m))\(Aadj*Y);
            x3(:,1)=real(ifft(X3(1:param.m)));
            x3(:,2)=real(ifft(X3(param.m+1:2*param.m)));
            err1=fittingError(rho.for{1},rho.back{2},x1(:,1),x1(:,2));
            err2=fittingError(rho.for{1},rho.back{2},x2(:,1),x2(:,2));
            err3=fittingError(rho.for{1},rho.back{2},x3(:,1),x3(:,2));
            for it=1:15
                errv=[err1 err2 err3];
                if err1==min(errv)
                    a3=a2;
                    err3=err2;
                    x3=x2;
                    a2=(a3+a1)/2;
                    x2=zeros(param.m,2);
                    X2=(Aadj*A+a2*eye(2*param.m))\(Aadj*Y);
                    x2(:,1)=real(ifft(X2(1:param.m)));
                    x2(:,2)=real(ifft(X2(param.m+1:2*param.m)));      
                    err2=fittingError(rho.for{1},rho.back{2},x2(:,1),x2(:,2));
                elseif err3==min(errv)
                    a1=a2;
                    err1=err2;
                    x1=x2;
                    a2=(a3+a1)/2;
                    x2=zeros(param.m,2);
                    X2=(Aadj*A+a2*eye(2*param.m))\(Aadj*Y);
                    x2(:,1)=real(ifft(X2(1:param.m)));
                    x2(:,2)=real(ifft(X2(param.m+1:2*param.m)));      
                    err2=fittingError(rho.for{1},rho.back{2},x2(:,1),x2(:,2));
                elseif err2==min(errv)
                    a4=(a2+a1)/2;
                    a5=(a3+a2)/2;
                    x4=zeros(param.m,2);
                    x5=zeros(param.m,2);
                    X4=(Aadj*A+a4*eye(2*param.m))\(Aadj*Y);
                    x4(:,1)=real(ifft(X4(1:param.m)));
                    x4(:,2)=real(ifft(X4(param.m+1:2*param.m)));      
                    err4=fittingError(rho.for{1},rho.back{2},x4(:,1),x4(:,2));
                    X5=(Aadj*A+a5*eye(2*param.m))\(Aadj*Y);
                    x5(:,1)=real(ifft(X5(1:param.m)));
                    x5(:,2)=real(ifft(X5(param.m+1:2*param.m)));      
                    err5=fittingError(rho.for{1},rho.back{2},x5(:,1),x5(:,2));
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
                rho_rec.back{param.Nwaves}=x1(:,2);
            elseif err2==min(errv)
                alpha=a2;
                rho_rec.for{1}=x2(:,1);
                rho_rec.back{param.Nwaves}=x2(:,2);
            elseif err3==min(errv)
                alpha=a3;
                rho_rec.for{1}=x3(:,1);
                rho_rec.back{param.Nwaves}=x3(:,2);
            end 
            if alpha==a_start
                a_start=2*alpha;
            else
                fprintf("alpha found as %.8f with min-error after %i iterations\n",alpha,it);
                endit=1;
            end
        end
    otherwise
        error("not a valid stopping rule");
end
%Create whole set of reconstructed waves
for it=2:param.Nwaves
rho_rec.for{it}=FourierShift(rho_rec.for{1},param.tau(it));
rho_rec.back{param.Nwaves+1-it}=FourierShift(rho_rec.back{param.Nwaves},param.tau(param.Nwaves)-param.tau(param.Nwaves+1-it));
end
for it=1:param.Nwaves
    rho_rec.sum{it}=rho_rec.for{it}+rho_rec.back{it};
end
end

