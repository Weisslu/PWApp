function [magWeightVel, angio] = calc_angio(MAG, v, Venc)

% time resolved
Vmag = squeeze(sqrt( sum( v.^2,4)));
Vmag(Vmag > Venc) = Venc;

magWeightVel = 32000*MAG.*sin( pi/2*Vmag / Venc);

% time averaged
mm = mean(MAG,4); 
vMean = mean(v,5);
Vmag = squeeze(sqrt( sum( vMean.^2,4)));
Vmag(Vmag > Venc) = Venc;

%According to algorithm I in Anderson et. al. 2008
angio = 32000*mm.*sin( pi/2*Vmag / Venc);
return