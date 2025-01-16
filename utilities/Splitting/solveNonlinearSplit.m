function [rho_rec,u] = solveNonlinearSplit(rho,param,initial_pwv,save)
%% Solve Pulsewave splitting problem via Separating the problem into looking for waves and velocity alternatingly
%pwv=zeros(param.stop,1);
pwv(1)=initial_pwv; % start velocity
%% Precomputations
Y=zeros(param.Nwaves*param.m,1);
conv=zeros(10,1);
for it=1:param.Nwaves
    Y(((it-1)*param.m+1):(it*param.m))=rho.data{it};
end
d=(param.beta*abs([0:floor(param.m/2)-1 floor(-param.m/2):-1]).^2+1).^(-param.s);
Adj_embedding=diag([d d]);
k=0;
while 1
    k=k+1;
    %% Tikhonov step for waves
    param.tau=param.distance/param.effectivedT/pwv(k); %tau=L/u delay in pixel
    A=operatorA(param);
    Aadj=Adj_embedding*A';
    X=(Aadj*A+param.alpha*eye(2*param.m))\(Aadj*Y);
    x1=X(1:param.m);
    x2=X(param.m+1:2*param.m);
    %figure 
    %plot(abs(x1)-abs(x2))
    if 0 %Additional assumption that fourier coeffs of backwave are smaller than fourier coeffs of forwave
        j=abs(x1)-abs(x2) <0;
        rel_diff=((abs(x1)-abs(x2))./abs(x2));
        x2(j)=x2(j).*(1+rel_diff(j));
    end
    %figure
    %plot(abs(x1)-abs(x2))

    %% Landweber Iteration for velocity
    [pwv(k+1)] = solveNonlinearLandweberVelocityOnly(pwv(k),rho,x1,x2,param);
    conv=circshift(conv,1);
    conv(1)=abs(pwv(k+1)-pwv(k));
    if k>10 && sum(conv)<0.0001
        break
    end
end  
fprintf("Outer Iterations: %i\n",k+1)
u=pwv(k+1);
if save
    figure
    plot(1:(k+1),pwv)
    yline(param.PWV,'g')
    legend("computed pwv","real pwv")
    xlabel("Iterations")
    ylabel("pwv")
    str="ADM_u"+num2str(param.PWV)+"_d"+num2str(param.delta)+"_N"+num2str(param.Nwaves);
    savefig(str+".fig")
    saveas(gcf,str+".png")
end
param.tau=param.distance/param.effectivedT/u; %tau=L/u delay in pixel
[rho_rec,~] = solveLinear_tikh(rho,param);

    
