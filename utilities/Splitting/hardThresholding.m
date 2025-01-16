function H=hardThresholding(x,t)
H=zeros(size(x));
k=abs(x)>t;
H(k)=x(k);
end