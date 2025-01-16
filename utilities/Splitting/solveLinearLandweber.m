function [rho_rec,pwv] = solveLinearLandweber(rho,param)
%% Solve Pulsewave splitting problem with Tikhonov Regularization
%% Transform Problem into Frequency Domain
% Force the input to be a column vector.
Rho=zeros(param.Nwaves*param.m,1);
for i=1:param.Nwaves
    Rho(param.m*(i-1)+1:param.m*i)=fft(rho.sum{i});
end
if ~(isfield(param,'tau'))
    error("No specified pwv!")
end
%% Landweber Iteration
it=0;
stop=0;
sol=zeros(2*param.m,1);
sol_old=sol;
%Include adjoint embedding operator for weighted spaces
d=(param.beta*abs([0:floor(param.m/2)-1 floor(-param.m/2):-1]).^2+1).^(-param.s);
Adj_embedding=diag([d d]);
A=operatorA(param);
while it<param.max_it && ~stop
    it=it+1;
    if param.nesterov
        z=sol+(it-1)/(it+2)*(sol-sol_old);
    else
        z=sol;
    end
    residual=Rho -A*z;
    s=Adj_embedding*A'*residual;
    %Omega chosen via steepest descent
    omega(it)=(norm([d.^(-1);d.^(-1)].*s)/norm(A*s))^2;
    sol_old=sol;
    sol=z+omega(it)*s;
    %Heuristic discrepancy principle, avoiding minimum for small iterations
    %(first 5)
    HDP_res=Rho-A*sol;
    HDP(it)=sqrt(it)*norm(HDP_res,2); 
    if it>5
        if HDP(it)>HDP(it-1) && HDP(it-1)<HDP(it-2) 
            stop=1;
        end
    end
end
if it==param.max_it
    fprintf("Maximum iterations reached at %i iterations \n",it)
else
    fprintf("Landweber stopped after %i iterations via the Heuristic DP \n",it);
end
%figure
%plot(1:it,res_norm,1:it,HDP)
%legend("residual norm","step size")
rho_rec.for{1}=real(ifft(sol(1:param.m)));
rho_rec.back{param.Nwaves}=real(ifft(sol((param.m+1):(2*param.m))));

%Create whole set of reconstructed waves
for it=2:param.Nwaves
rho_rec.for{it}=FourierShift(rho_rec.for{1},param.tau(it));
rho_rec.back{param.Nwaves+1-it}=FourierShift(rho_rec.back{param.Nwaves},param.tau(param.Nwaves)-param.tau(param.Nwaves+1-it));
end
for it=1:param.Nwaves
    rho_rec.sum{it}=rho_rec.for{it}+rho_rec.back{it};
end
end

