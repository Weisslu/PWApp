%Author: LWeissinger 23.11.2022
%Converts the data set into a useable set of waves, adaptable number of
%waves!
function rho = getRealWave(waveforms,param)
    for it=1:param.Nwaves
        rho.sum{it}=waveforms(it,:)';
        rho.data{it}=fft(rho.sum{it}); 
    end
end