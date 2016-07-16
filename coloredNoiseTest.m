import dsp.*;
Fs = 5e6;
H = dsp.ColoredNoise('InverseFrequencyPower',0.7,'SamplesPerFrame',Fs);
Pxx = zeros(Fs,1e3);
for i=1:1e1
    noise=step(H);
    nfft = length(noise);
    Pxx(:,i) = abs(fft(noise,nfft)).^2/length(noise)/Fs;
end
plot(mean(Pxx,2));

%%
% Create a single-sided spectrum
Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
plot(Hpsd)

%%
