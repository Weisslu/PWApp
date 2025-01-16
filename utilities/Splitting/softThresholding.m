function S=softThresholding(x,t)
S=sign(x).*max(abs(x)-t,0);
end