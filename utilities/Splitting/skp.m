function [val] = skp(x1,x2)
% Just L2 inner product for two vectors
%intervallbreite ist param.m 
% there should be periodic continuation of the involved functions before
% applying trapz! For starters, assuming pw constant basis functions should
% suffice.
%val=trapz(x1.*conj(x2));
x1=x1(:);
x2=x2(:);
val=x1'*x2;
end