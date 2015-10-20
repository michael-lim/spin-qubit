import dsp.*;

Fs = 5009999;
t=tic;
H = dsp.ColoredNoise('InverseFrequencyPower',0.7,'SamplesPerFrame',N_points);
noise=step(H);
toc(t)

rng default
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t) + randn(size(t));