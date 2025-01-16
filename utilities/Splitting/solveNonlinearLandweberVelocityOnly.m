function [pwv] = solveNonlinearLandweberVelocityOnly(pwv,rho,x1,x2,param)
%% Solve Pulsewave splitting problem with Tikhonov Regularization
%% Transform Problem into Frequency Domain

% Force the input to be a column vector.
Rho=zeros(param.Nwaves*param.m,1);
for i=1:param.Nwaves
    Rho(param.m*(i-1)+1:param.m*i)=rho.data{i};
end
%% Landweber Iteration
it=0;
pwv_old=pwv;
conv=zeros(10,1);
%Include adjoint embedding operator for weighted spaces
while it<param.max_it
    it=it+1;
    if param.nesterov
        z=pwv+(it-1)/(it+2)*(pwv-pwv_old);
    else
        z=pwv;
    end
    residual=Rho - operatorF(x1,x2,z,param);
    %===Plot y_k==
%     Y=operatorF(x1,x2,z,param);
%     figure
%     hold on
%     for k=1:param.Nwaves
%         plot(1:param.m,ifft(Y(((k-1)*param.m+1):(k*param.m),1)))
%     end
%     hold off
%     figure
%     hold on
%     for k=1:param.Nwaves
%         plot(1:param.m,ifft(residual(((k-1)*param.m+1):(k*param.m),1)))
%     end
%     hold off
    %=============
    FrDer=Frechetderivative_x(x1,x2,z,param);
    s=FrDer'*residual;
    %Omega chosen via steepest descent
    omega(it)=1*(abs(s)/norm(FrDer*s))^2;
    %u(it)=pwv;
    pwv_old=pwv;
    pwv=real(z+omega(it)*s);
    if pwv<1
        pwv=1;
    elseif pwv>10
        pwv=10;
    end
    conv=circshift(conv,1);
    conv(1)=abs(pwv-pwv_old);
    if sum(conv)<0.00001 && it>10
        break
    end
    %Heuristic discrepancy principle, avoiding minimum for small iterations
    %(first 5)
%     if ~param.justmax
%     HDP_res=Rho-operatorF(x1,x2,pwv,param);
%     res_norm(it)=norm(HDP_res,2);
%     HDP(it)=sqrt(it)*norm(HDP_res,2); 
%         if it>5
%             if HDP(it)>HDP(it-1) && HDP(it-1)<HDP(it-2)
%                 break
%             end
%         end
%     end
end
% it
% if it==param.max_it
%     fprintf("Maximum iterations reached at %i iterations \n",it)
% else
%     fprintf("Landweber stopped after %i iterations via the Heuristic DP \n",it);
% end
% figure
% plot(1:it,res_norm,1:it,omega)
% legend("residual norm","step size")
%  figure
%  plot(1:it,u)
%  legend("u")
end

