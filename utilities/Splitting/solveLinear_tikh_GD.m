function [rho_rec,alpha] = solveLinear_tikh_GD(rho,param)
%% Solve Pulsewave splitting problem with Tikhonov Regularization
%% Implemented via Gradient Descent algorithm to incorporate constraint x2<x1
%% Transform Problem into Frequency Domain
alpha=param.alpha;
if param.parameter_choice_rule~='apriori'
    error("GD only with a-priori choice rule implemented")
end
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
omega=1/(2*double(param.Nwaves)+alpha);%Operator norm known for Lipschitz-constant
x=zeros(2*param.m,1);
x_old=x;
t=1;
it=0;
while it~=param.max_it
    it=it+1;
    t_old=t;
    t=(1+sqrt(1+4*t^2))/2;
    z=x+(t_old-1)/(t)*(x-x_old);
    x_old=x;
    gradg=Aadj*(A*z-Y)+alpha*z;
    x=z-omega*gradg;
    x1=x(1:param.m);
    x2=x(param.m+1:2*param.m);
    if 1 %Additional assumption that fourier coeffs of backwave are smaller than fourier coeffs of forwave
        j=abs(x1)-abs(x2) <0;
        rel_diff=((abs(x1)-abs(x2))./abs(x2));
        x2(j)=x2(j).*(1+rel_diff(j));
    end
    x(1:param.m)=x1;
    x(param.m+1:2*param.m)=x2;
end
rho_rec.for{1}=real(ifft(x(1:param.m)));
rho_rec.back{param.Nwaves}=real(ifft(x(param.m+1:2*param.m)));
        
%Create whole set of reconstructed waves
for it=2:param.Nwaves
rho_rec.for{it}=FourierShift(rho_rec.for{1},param.tau(it));
rho_rec.back{param.Nwaves+1-it}=FourierShift(rho_rec.back{param.Nwaves},param.tau(param.Nwaves)-param.tau(param.Nwaves+1-it));
end
for it=1:param.Nwaves
    rho_rec.sum{it}=rho_rec.for{it}+rho_rec.back{it};
end
end

