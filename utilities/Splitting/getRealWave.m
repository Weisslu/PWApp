%LW 23.11.2022
%Converts the data set into a useable set of waves, adaptable number of
%waves!
function rho = getRealWave(waveforms,param)
%color=["k" "g" "b" "r" "c" "m" "y"];
%figure
%hold on
for it=1:param.Nwaves
    rho.sum_raw{it}=waveforms(it,:)';
    rho.sum{it}=rho.sum_raw{it};
    rho.data{it}=fft(rho.sum_raw{it}); 
    %p1(it)=plot(param.taxis,rho.sum_raw{it},color(1+mod(it,7))+':','linewidth',1);
    %p1label(it)="\rho_{"+num2str(it)+"} raw";
    %legend(cat(1,p1label),'Location','northeast');
end
%hold off
%axis tight
%title("Input waves from test data")
end