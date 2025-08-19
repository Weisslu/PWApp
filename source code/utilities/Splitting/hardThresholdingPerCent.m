%Author: LWeissinger, 31.07.2025
function H=hardThresholdingPerCent(x,t)
H=zeros(size(x));
n=numel(x);
s=sort(abs(x));
if floor(n*t)<=0
        H=x;
elseif floor(n*t)>n
else
    k=abs(x)>s(floor(n*t));
    H(k)=x(k);
end
end