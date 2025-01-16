function val = norm_X(s,param)
%LW 01.12.2022
% Norm in (L^2_s)^2xR, taking in account a multiplication of a constant with the velocity part. 
% 
weight=(1+abs([0:floor(param.m/2)-1 floor(-param.m/2):-1]').^param.s);
val=norm([[weight;weight].*s(1:2*param.m);sqrt(param.c)*s(2*param.m+1)]);
end