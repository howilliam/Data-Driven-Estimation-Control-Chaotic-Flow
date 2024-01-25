function [Y,f]=fourier_transform(t,y)

T=(t(2)-t(1));
N=length(y);

Fs=1/T;
NFFT = 2^nextpow2(N); % Next power of 2 from length of y
Y = fft(y,NFFT)/N;
Y = 2*abs(Y(1:NFFT/2+1));
f = Fs/2*linspace(0,1,NFFT/2+1);

end