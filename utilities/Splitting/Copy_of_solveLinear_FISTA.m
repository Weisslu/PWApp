function [rho_rec,alpha] = solveLinear_FISTA()
%% Solve Pulsewave splitting problem with Tikhonov Regularization

% Force the input to be a column vector.
%% Get discrete matrix A corresponding to operator A in the fourier domain
theta=0:179;
[A,~,x_true]=paralleltomo(100,theta);%operatorA(param);
% Add noise to the data.
x_true=reshape(x_true,100,[]);
[Fx,Fy]=gradient(x_true);
x_true=reshape(boolean(Fx+Fy),[],1);
b_ex=A*x_true;
rng(0);
e = randn(size(b_ex));
Y = b_ex + 0.01*norm(b_ex)*e/norm(e);


%% Precomputations


%Include adjoint embedding operator for weighted spaces
%d=(param.beta*abs([0:floor(param.m/2)-1 floor(-param.m/2):-1]).^2+1).^(-param.s);
%Adj_embedding=diag([d d]);
Aadj=A';%Adj_embedding*A';
omega=1/(normest(A)^2);%Operator norm known
stop=0;
x=zeros(100^2,1);
x_old=x;
t=1;
it=0;
HDP=zeros(3,1);
nY=norm(Y);
nX=norm(double(x_true));
while stop~=1 && it~=5000
    it=it+1;
    t_old=t;
    t=(1+sqrt(1+4*t^2))/2;
    z=x+(t_old-1)/(t)*(x-x_old);
    x_old=x;
    gradg=Aadj*(A*z-Y);
    %omega(it)=(norm([d.^(-1);d.^(-1)].*gradg)/norm(A*gradg))^2;
    
    x=softThresholding(z-omega*gradg,10*omega);
%     HDP=circshift(HDP,1);
%     HDP_res=Y-A*x;
%     HDP(1)=sqrt(it)*norm(HDP_res,2); 
%     if it>10
%         if (HDP(1)>HDP(2) && HDP(2)<HDP(3)) || it==param.stop
%             stop=1;
%         end
%     end
    %e_fit(it)=norm(x-x_true)/nX;
    %e_res(it)=norm(A*x-Y)/nY;
end
figure
imagesc(reshape(x,100,[]))
%xt=(A'*A+0.01*eye(100^2))^(-1)*(A'*Y);
%figure
%imagesc(reshape(xt,100,[]))
end

